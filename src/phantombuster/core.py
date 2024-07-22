import logging
import os
import os.path
import sys
import multiprocessing as mp
import glob

import ntpath
import pandas as pd
#from phantombuster import old
from dataclasses import dataclass
import json

import phantombuster.remoter
from phantombuster import porcelain, plumbing, stores
from phantombuster.stores import deduplicator_to_pyarrow_table
from phantombuster.remoter import Scheduler
from phantombuster.io_ import PathsAndFiles, write_parquet
from phantombuster.project import Project
from phantombuster.config_files import read_barcode_hierarchy_file, read_input_files_file, read_regex_file
import click
from typing import Optional, List
import pyarrow.parquet
import pyarrow.csv

from pathlib import Path


def demultiplex(input_files_file, regex_file, barcode_hierarchy_file, project, debug=False, show_qc=False):
    project.create()

    input_groups = read_input_files_file(input_files_file)
    barcodes = read_barcode_hierarchy_file(barcode_hierarchy_file)
    barcodes_names = [bc['name'] for bc in barcodes]
    regexes = read_regex_file(regex_file)

    for bc_name in barcodes_names:
        names_with_numbers = []
        if bc_name.rstrip("0123456789") != bc_name:
            names_with_numbers.append(bc_name)
    if len(names_with_numbers) > 0:
        logging.error(f'Barcode names missconfigured. barcode names cannot have numbers at the end. Problematic barcode names: {names_with_numbers}.')
        raise Exception(f'Barcodes missconfigured: {names_with_numbers}')
           
    missconfigured_groups = []
    for input_group in input_groups:
        rs = regexes.get_regexes_for_group(input_group.group)
        r_tags = [name for r in rs.values() for name in r.groupindex.keys()]
        r_tags  = list(set([tag.rstrip('0123456789') for tag in r_tags]))

        missing_tags = [bc for bc in barcodes_names if bc not in r_tags]
        superflous_tags = [tag for tag in r_tags if tag not in barcodes_names]
        if len(missing_tags) > 0 or len(superflous_tags) > 0:
            missconfigured_groups.append((input_group, missing_tags, superflous_tags))

    if len(missconfigured_groups) > 0:
        for g, missing, superflous in missconfigured_groups:
            logging.error(f'Group {g.group} with example file {g.files[0].path} missconfigured. Missing tags: {missing}. Unused tags: {superflous}.')
        raise Exception('Groups missconfigured')


    phantombuster.remoter.load(stores.load, dir=str(project.tmp_dir))(plumbing.deduplicate_section)
    phantombuster.remoter.save(stores.save, dir=str(project.tmp_dir))(plumbing.deduplicate_section)

    phantombuster.remoter.load(stores.load, dir=str(project.tmp_dir))(plumbing.combine)
    phantombuster.remoter.save(stores.save, dir=str(project.tmp_dir))(plumbing.combine)

    scheduler = project.get_scheduler()

    with scheduler:
        maintable, (reads, remaining_reads, counters_success, counters_fail) = porcelain.demultiplex(input_files_file, project.demultiplex_output_path, project.demultiplex_stats_path, regex_file, barcode_hierarchy_file, debug, show_qc, scheduler)
    m = max([len(str(count)) for count in counters_fail.values()] + [len(str(reads))])

    logging.info(f"{str(reads).rjust(m)} Reads processed")
    logging.info(f"{str(remaining_reads).rjust(m)} Reads Remain ({100*remaining_reads/reads}%)")
    for reason, count in counters_fail.items():
        logging.info(f"{str(count).rjust(m)} Reads were unfit due to {reason}")

    r = (maintable, (reads, remaining_reads, counters_success, counters_fail))
    return r


def error_correct(project, error_threshold, barcode_hierarchy_file):
    project.create()

    phantombuster.remoter.load(stores.load, dir=str(project.tmp_dir))(porcelain.error_correct_partition)
    phantombuster.remoter.save(stores.save, dir=str(project.tmp_dir))(porcelain.error_correct_partition)

    barcode_hierarchy = read_barcode_hierarchy_file(barcode_hierarchy_file)

    table_file = project.demultiplex_output_path
    out_file = project.error_correct_output_path

    error_corrected = porcelain.error_correct(table_file, out_file, error_threshold, barcode_hierarchy, project)

    print(f"LIDs after correction: {error_corrected['lids_corrected']}")
    print(f"Reads after correction: {error_corrected['reads_corrected']}")

    with open(project.error_correct_stats_path, mode="w") as f:
        json.dump(error_corrected, f)

    logging.info('Error correction done')


def hopping_removal(project, hopping_barcodes, threshold):
    project.create()

    table_file = project.error_correct_output_path

    r, stats = porcelain.hopping_removal(table_file, threshold, hopping_barcodes)
    write_parquet(r, project.hopping_removal_output_path)
    with open(project.hopping_removal_stats_path, mode='w') as f:
        json.dump(stats, f)


def threshold(project, threshold_file):
    project.create()

    # read in threshold file
    thresholds = porcelain.read_threshold_file(threshold_file)

    porcelain.threshold(project, threshold_file)
