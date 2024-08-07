import regex
import json
import phantombuster.plumbing as plumbing
import phantombuster.config_files

from .plumbing import pairwise
from phantombuster import bamindexer
from phantombuster.io_ import write_parquet

import os
import logging
from collections import defaultdict, namedtuple, Counter
import gzip

import pyarrow as pa
import pyarrow.types
import pyarrow.csv
import pyarrow.parquet
import numpy as np

import polars as pl


class PhantomBusterException(Exception):
    pass

class PhantomBusterUnknownFileType(PhantomBusterException):
    pass


def index_group(group, every):
    if len(group.files) == 1:
        return index_file(group.files[0].path, every)
    else:
        idxs = [[0]* len(group.files)] + [[np.inf]* len(group.files)]
        return list(zip(*idxs))


def index_file(file, every):
    if file.endswith('.bam'):
        return index_bamfile(file, every)
    elif file.endswith('.fastq.gz'):
        return index_fastqgz_file(file, every)

def index_bamfile(bamfile, every):
    with open(bamfile, mode="rb") as f:
        return bamindexer.index(f, every)


def _read_complete_read(f, read_size=1024):
    buffer = []
    lines = 0
    tries = 0
    while lines < 7 and tries < 100:
        read = f.read(read_size)
        if len(read) == 0:
            raise EOFError()
        lines += read.count(b'\n')
        buffer.append(read)
        tries += 1
    if tries >= 100:
        raise ValueError("Can't find reads in fastq file")
    complete_read = b''.join(buffer)
    lines = complete_read.split(b'\n')
    for line_idx in range(5):
        if len(lines[line_idx]) > 0 and lines[line_idx][0] == b'@'[0]:
            break
        else:
            pass
            #print('no @ in', lines[line_idx][0])
    else:
        raise EOFError()

    single_read = b'\n'.join(lines[line_idx:line_idx+4])
    offset = complete_read.find(single_read)
    return single_read, offset


def index_fastqgz_file(file, every):
    with gzip.open(file) as f:
        read, offset = _read_complete_read(f)
        total_readsize = len(read)
        reads_analysed = 1

        readsize = int(total_readsize / reads_analysed)

        idx = 0
        indices = []
        while True:
            try:
                f.seek(idx)
                read, offset = _read_complete_read(f, readsize)

                indices.append(idx+offset)
                idx = idx + offset + every

                total_readsize += len(read)
                reads_analysed += 1
                readsize = int(total_readsize / reads_analysed)

            except EOFError:
                break
        indices.append(np.inf)
    return indices


def demultiplex(input_files_file, outpath, outpath_stats, regex_file, barcode_hierarchy_file, debug, show_qc, scheduler):

    logger = logging.getLogger("phantombuster")

    input_groups = phantombuster.config_files.read_input_files_file(input_files_file)

    logger.info(f"Deduplicating {len(input_groups)} many groups")
    logger.info("Reading mappings")

    # read mapping files
    barcode_hierarchy = phantombuster.config_files.read_barcode_hierarchy_file(barcode_hierarchy_file)

    regex_dictionary = phantombuster.config_files.read_regex_file(regex_file)

    logger.info("Scheduling indexing %s group(s)", len(input_groups))
    # BGZF blocks per task, each block is <64KB so 20000 should roughly be 1.28 GB
    if scheduler.query_workers() == 0:
        logger.warning("Work is getting scheduled and no workers are currently connected. Launch workers now.")

    idxs_promises = []
    for group in input_groups:
        promise = scheduler.schedule(index_group, group, 20000)
        idxs_promises.append(promise)

    for promise in idxs_promises:
        r = 0
        waiting_time = 0
        waiting_notice_time = 10
        while isinstance(r, int):
            waiting_time += r
            if waiting_time > waiting_notice_time:
                waiting_time = 0
                number_of_workers = scheduler.query_workers()
                if number_of_workers == 0:
                    logger.warning("Currently no workes are connected to the server. Server will make no progress without workers.")
                    if waiting_notice_time == 10:
                        waiting_notice_time = 60
                    elif waiting_notice_time == 60:
                        waiting_notice_time = 10*60
                else:
                    waiting_notice_time = 10*60
            r = scheduler.wait(promise, return_after_wait=True)

    idxs = [scheduler.wait(promise, unpack=True) for promise in idxs_promises]

    logger.info("Parsing file(s) piecewise, %s many subtasks", sum(len(idx) for idx in idxs))
    logger.debug("Index is %s", idxs)
    outs = []
    for group, index in zip(input_groups, idxs):
        for section in pairwise(index):
            logger.debug("Scheduling deduplicate section %s", section)
            t = scheduler.schedule(plumbing.deduplicate_section, group, section, regex_dictionary, barcode_hierarchy, show_qc)
            outs.append(t)
    scheduler.wait_all(outs, unpack=False)

    logger.info("Extracted reads from individual files, combining %s many files", len(outs))
    new_outs = []
    while len(outs) > 1:
        while len(outs) > 1:
            a = outs.pop(0)
            b = outs.pop(0)
            logger.debug("Waiting on two files")
            scheduler.wait(a, unpack=False)
            scheduler.wait(b, unpack=False)
            logger.debug("Scheduling merging of two files, %s remaining.", len(outs))
            result = scheduler.schedule(plumbing.combine, [a, b], barcode_hierarchy)
            new_outs.append(result)
        if len(outs) == 1:
            logger.debug("Cannot merge, only one file remaining")
            new_outs.append(outs[0])
        logger.debug("One round done, there were %s many files, now there are %s", len(outs), len(new_outs))
        outs = new_outs
        new_outs = []

    logger.debug("Waiting for last combination round")
    out = scheduler.wait(outs[-1], unpack=True)
    if out is None or out[0] is None:
        logger.error("Not a single read was extracted. This might indicates a problem with the references or a missconfiguration of the regexes. Aborting. ")
        logger.debug('Result of combining all files is: %s', out)
        raise Exception()

    write_parquet(out[0], outpath)

    with open(outpath_stats, mode="w") as f:
        json.dump(out[1], f)

    logger.info("Demultiplexing done")
    return out

def read_threshold_file(threshold_file):
    if threshold_file.endswith(".csv"):
        f = pa.csv.read_csv(threshold_file)
    elif threshold_file.endswith(".parquet"):
        f = pa.parquet.read_table(threshold_file)
    else:
        raise PhantomBusterUnknownFileType("Thresholdfile {threshold_file} has unknown file type. Supported are .csv and .parquet.")
    return f

def error_correct(samplefile, outfilename, threshold, barcode_hierarchy, project, remove_ambigious):

    logger = logging.getLogger('phantombuster')

    logger.debug(f"Reading table {samplefile}")
    df = pl.scan_parquet(samplefile)
    number_of_reads, number_of_lineages = df.select(pl.col('reads').sum(), pl.len()).collect()
    number_of_reads, number_of_lineages = number_of_reads[0], number_of_lineages[0]
    logger.debug(f"Table has {number_of_lineages} many rows and {number_of_reads} many reads")

    first_bc = barcode_hierarchy[0]

    if first_bc["type"] == 'reference':
        logger.debug(f"First barcode in hierarchy '{barcode_hierarchy[0]['name']}' is of 'reference' type, partitioning")
        partitions = df.select(pl.col(first_bc['name']).unique().alias('values')).collect()['values'].to_list()
        logger.debug(f"Partitioning Done, created {len(partitions)} many partitions.")
    else:
        logger.debug(f"First barcode in hierarchy '{barcode_hierarchy[0]['name']}' is of 'random' type, no partitioning performed")
        partitions = [None]

    scheduler = project.get_scheduler()
    with scheduler as scheduler:
        logger.debug('Distributing each partition to one worker')
        tasks = []
        for value in partitions:
            task = scheduler.schedule(error_correct_partition, samplefile, value, None, threshold, barcode_hierarchy, use_column=first_bc['name'], remove_ambigious=remove_ambigious)
            tasks.append(task)
        logger.debug('Distributed all partitions, waiting for workers')
        results = scheduler.gather(*tasks)
        logger.debug('Workers done, performing postprocessing')

    tables = [table for table, stats in results]

    cast_tables = []
    for table in tables:
        for column_name in table.column_names:
            if pa.types.is_string(table[column_name].type):
                table = table.set_column(table.column_names.index(column_name), column_name, table[column_name].cast("large_string"))
            if pa.types.is_binary(table[column_name].type):
                table = table.set_column(table.column_names.index(column_name), column_name, table[column_name].cast("large_binary"))
        cast_tables.append(table)
    tables = cast_tables

    corrected_table = pyarrow.concat_tables(tables)

    number_of_lineages_corrected = len(corrected_table)
    number_of_reads_corrected = pa.compute.sum(corrected_table["reads"]).as_py()

    write_parquet(corrected_table, outfilename)

    return {"sample": samplefile,
            "lids_uncorrected": number_of_lineages,
            "lids_corrected": number_of_lineages_corrected,
            "reads_uncorrected": number_of_reads,
            "reads_corrected": number_of_reads_corrected}

def error_correct_partition(samplefile, lower_or_value, length_opt, threshold, barcode_hierarchy, use_column=False, remove_ambigious=True):

    logger = logging.getLogger("phantombuster")

    if lower_or_value is None:
        logger.info(f'Read table from {samplefile} completly')
        table = pl.scan_parquet(samplefile).collect().to_arrow()
    elif use_column:
        assert isinstance(use_column, str)
        assert isinstance(lower_or_value, str)
        value = lower_or_value
        table = pl.scan_parquet(samplefile).filter(pl.col(use_column) == value).collect().to_arrow()
        logger.info(f'Read table from {samplefile} with {use_column} only {value}.')
    else:
        assert isinstance(lower_or_value, int)
        lower = lower_or_value
        assert length_opt is not None
        length = length_opt

        table = pa.parquet.read_table(samplefile, memory_map=True, use_threads=False)
        table = table.slice(lower, length)
        logger.info(f'Read table from {samplefile} from {lower} with length {length}.')

    number_of_lineages = len(table)
    number_of_reads = pa.compute.sum(table["reads"]).as_py()

    logger.info('Starting error correction')
    corrected_table = plumbing.error_correct(table, barcode_hierarchy, threshold)
    logger.info('Done error correction')

    if remove_ambigious:
        logger.info('Removing barcode sequences with ambigious characters (NRYSWKMBDHV).')
        corrected_table = plumbing.remove_ambigious(corrected_table, barcode_hierarchy)

    number_of_lineages_corrected = len(corrected_table)
    number_of_reads_corrected = pa.compute.sum(corrected_table["reads"]).as_py()

    stats = {"lids_uncorrected": number_of_lineages,
            "lids_corrected": number_of_lineages_corrected,
            "reads_uncorrected": number_of_reads,
            "reads_corrected": number_of_reads_corrected}

    return corrected_table, stats

def _read_mappingfile(mappingfile):
    table = pa.csv.read_csv(mappingfile)
    map = dict(zip(table["sample"].to_pylist(), table["group"].to_pylist()))
    return map

def hopping_removal(input_file, alpha_threshold, hopping_barcodes):
    logger = logging.getLogger("phantombuster")
    master_table = pa.parquet.read_table(input_file)

    number_of_lineages = len(master_table)
    number_of_reads = pa.compute.sum(master_table["reads"]).as_py()

    for column_name in master_table.column_names:
        if pa.types.is_string(master_table[column_name].type):
            logger.warning(f"Detected datatype 'string' in column {column_name}, casting to 'large_string'")
            master_table = master_table.set_column(master_table.column_names.index(column_name), column_name, master_table[column_name].cast("large_string"))
        if pa.types.is_binary(master_table[column_name].type):
            logger.warning(f"Detected datatype 'binary' in column {column_name}, casting to 'large_binary'")
            master_table = master_table.set_column(master_table.column_names.index(column_name), column_name, master_table[column_name].cast("large_binary"))

    logger.info("Starting Hopping Removal")
    threshold_list = []
    for barcode in hopping_barcodes:
        logger.info(f"Considering Barcode '{barcode}'")
        threshold = plumbing.calculate_hopping_threshold(master_table, barcode, alpha_threshold=alpha_threshold)
        threshold_list.append(threshold)
    master_table = plumbing.call_hopping(master_table, threshold_list)

    logger.info("Hopping Removal Done")

    number_of_lineages_corrected = len(master_table)
    number_of_reads_corrected = pa.compute.sum(master_table["reads"]).as_py()

    stats = {"lids_before": number_of_lineages,
             "lids_after": number_of_lineages_corrected,
             "reads_before": number_of_reads,
             "reads_after": number_of_reads_corrected}

    return master_table, stats


def threshold(project, threshold_file):
    thresholds = pl.DataFrame(read_threshold_file(threshold_file))
    table = pl.read_parquet(project.hopping_removal_output_path)

    lids_before = len(table)
    reads_before = table.select(pl.col('reads').sum())['reads'][0]

    if len(thresholds.columns) == 1:
        table = table.filter(pl.col('reads') >= thresholds['threshold'][0])
    else:
        columns = table.columns
        table = table.join(thresholds, on=[col for col in table.columns if col != 'threshold' and col in thresholds.columns]).filter(pl.col('reads') >= pl.col('threshold'))
        table = table.select(columns)

    lids_after = len(table)
    reads_after = table.select(pl.col('reads').sum())['reads'][0]

    stats = {"lids_before": lids_before,
             "lids_after": lids_after,
             "reads_before": reads_before,
             "reads_after": reads_after}
    return table.to_arrow(), stats
