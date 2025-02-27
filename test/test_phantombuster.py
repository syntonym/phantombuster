import sys
import hashlib
import subprocess
import os
import tempfile
import pytest
import shutil
import numpy.random
import numpy as np
import glob
from click.testing import CliRunner
import gzip
import regex
import polars as pl


from os.path import join as pjoin
import uuid
import pandas
from collections import Counter
import pyarrow
import pyarrow as pa
import datetime

from phantombuster import plumbing, analysis, core
import phantombuster.cli
import phantombuster.io_
import phantombuster.remoter
import phantombuster.plumbing
import phantombuster.statistics
from phantombuster.merge_cython import xmerge
from phantombuster.plumbing import pyarrow_table_to_counter
from phantombuster.plumbing import deduplicator_to_pyarrow_table
from phantombuster.plumbing import get_sectioned_iterator, pairwise
from phantombuster.porcelain import index_bamfile, index_fastqgz_file
from phantombuster.error_corrector import ErrorCorrector, Umi, UmiCount
from phantombuster.project import Project
import phantombuster.config_files
import phantombuster.remoter.serialization
import regex

import zmq
from phantombuster.handler import PushHandler
import logging
import threading
import time



phantombuster.remoter.Scheduler.DEFAULT_WORKERS = 2

ALPHABET = list("ACGT")

phantombuster.cli.configure_logging(None, False)


def generate_barcode(length, k=1, alphabet=ALPHABET):
    barcodes = numpy.random.choice(ALPHABET, (k, length))
    return ["".join(barcode) for barcode in barcodes[:, ]]


def gen_bc(barcodes, weights=None, blocksize=1024):
    while True:
        block = np.random.choice(barcodes, blocksize, p=weights)
        for lid in block:
            yield lid


class ReadGenerator:

    def __init__(self, samples=10, genes=20):
        self.fixed = "GGGC"

        self.b2regex = "(?P<umi>[ACGT]{9})"
        self.bcregex = "(?P<grna>[ACGT]{11})"
        self.queryregex = "GGGC(?P<sample>[ACGT]{5})(?P<lid>[ACGT]{30})"

        self.samples = generate_barcode(5, samples)
        self.sample_dict = {s: f"Sample{i}" for i, s in enumerate(self.samples)}

        self.genes = generate_barcode(11, genes)
        self.gene_dict = {s: f"Gene{i}" for i, s in enumerate(self.genes)}

    def generate_reads(self, N, lids=12000, umis=1000):

        lids = generate_barcode(30, lids)
        umis = generate_barcode(9, umis)

        lids_pregen = gen_bc(lids)
        umis_pregen = gen_bc(umis)
        genes_pregen = gen_bc(self.genes)
        sample_pregen = gen_bc(self.samples)

        reads = []

        for i in range(N):
            lid = next(lids_pregen)
            umi = next(umis_pregen)
            gene = next(genes_pregen)
            sample = next(sample_pregen)
            read = make_read(f"{self.fixed}{sample}{lid}", umi, gene)
            reads.append(read)

        return reads


def make_read(query, b2, bc):
    query_qc = "F"*len(query)
    b2_qc = "F"*len(b2)
    bc_qc = "F"*len(bc)

    read = f"D00689:314:CBYRBANXX:1:1104:4980:1967\t4\t*\t0\t0\t*\t*\t0\t0\t{query}\t{query_qc}\tB2:Z:{b2}\tQ2:Z:{b2_qc}\tBC:Z:{bc}\tQT:Z:{bc_qc}\tRG:Z:CBYRBANXX.1\n".encode("utf-8")
    return read


def copy_file_efficiently(src, dst):
    with open(src, mode="br") as f_src:
        count = os.fstat(f_src.fileno()).st_size
        with open(dst, mode="bw") as f_dst:
            shutil.copyfile(src, dst)


# -- Unit Tests -- #

def test_unit_serialization_extensions_1():
    input_file = phantombuster.config_files.InputFile('a', 'b', 'c')
    serialized = phantombuster.remoter.serialization.serialize_value(input_file)
    deserialized = phantombuster.remoter.serialization.deserialize_value(serialized, None)

    assert deserialized == input_file


def test_unit_serialization_extensions_2():
    input_file = phantombuster.config_files.RegexDictionary({('', '', 'b2'): '*]*([)<>*()'})
    serialized = phantombuster.remoter.serialization.serialize_value(input_file)
    deserialized = phantombuster.remoter.serialization.deserialize_value(serialized, None)

    assert deserialized == input_file


def test_unit_regex_validation_1():
    regex1 = regex.compile("AAACC(?P<sample1>[CGT]{10})")
    regex2 = regex.compile("AAACC(?P<sample2>[ACGT]{10})")
    regex3 = regex.compile("AAACC(?P<lid>[ACGT]{10})")

    assert analysis.regexes_are_valid([regex1, regex2, regex3])


def test_unit_regex_validation_2():
    regex1 = regex.compile("G(?P<grna>[ACGT]{6})")
    regex2 = regex.compile("AACC")
    regex3 = regex.compile("AAACC(?P<lid>[ACGT]{10})")

    assert not analysis.regexes_are_valid([regex1, regex2, regex3])


def test_unit_plumbing_extract_read_information_1():
    read = b"D00689:314:CBYRBANXX:1:1104:4980:1967\t4\t*\t0\t0\t*\t*\t0\t0\tCCTATGGGCGTGGAAAGGACGAAACACCGCATGGCCCTGAAGAATGATGGTT\tBBBBBFFFFFFFFFFFFFBFFFFFFFFFFF<<FFFFFFF</F/<FFFFFFFF\tB2:Z:CGCTATCAA\tQ2:Z:<BBB<///B\tBC:Z:TTTCCCACCCT\tQT:Z:BBBBBFFFFFF\tRG:Z:CBYRBANXX.1\n"
    with tempfile.NamedTemporaryFile() as f:
        f.write(read)
        f.flush()

        b2regex = regex.compile("")
        bcregex = regex.compile("(?P<grna>[ACGT]{11})")
        queryregex = regex.compile("GGGC(?P<sample>[ACGT]{5})(?P<lid>[ACGT]{30})")

        regex_dictionary = {'B2': b2regex, 'BC': bcregex, 'query': queryregex}

        bamfile = iter(open_sequencing_file(f.name, type='sam'))
        read = next(bamfile)
        features = phantombuster.plumbing.extract_read(regex_dictionary, read)

        print(features)

        assert features["grna"] == "TTTCCCACCCT"
        assert features["sample"] == "GTGGA"
        assert features["lid"] == "AAGGACGAAACACCGCATGGCCCTGAAGAA"


def test_unit_plumbing_extract_read_information_2():
    read = b"D00689:314:CBYRBANXX:1:1104:4980:1967\t4\t*\t0\t0\t*\t*\t0\t0\tCCTATGGGCGTGGAAAGGACGAAACACCGCATGGCCCTGAAGAATGATGGTT\tBBBBBFFFFFFFFFFFFFBFFFFFFFFFFF<<FFFFFFF</F/<FFFFFFFF\tB2:Z:CGCTATCAA\tQ2:Z:<BBB<///B\tBC:Z:TTTCCCACCCT\tQT:Z:BBBBBFFFFFF\tRG:Z:CBYRBANXX.1\n"
    with tempfile.NamedTemporaryFile() as f:
        f.write(read)
        f.flush()

        b2regex = regex.compile("(?P<lid1>[ACGT]{9})")
        bcregex = regex.compile("(?P<grna>[ACGT]{11})")
        queryregex = regex.compile("GGGC(?P<sample>[ACGT]{5})(?P<lid2>[ACGT]{30})")
        regex_dictionary = {'B2': b2regex, 'BC': bcregex, 'query': queryregex}

        bamfile = iter(open_sequencing_file(f.name, type='sam'))
        read = next(bamfile)
        features = phantombuster.plumbing.extract_read(regex_dictionary, read)

        assert features["grna"] == "TTTCCCACCCT"
        assert features["sample"] == "GTGGA"
        assert features["lid"][:30] == "CGCTATCAAAAGGACGAAACACCGCATGGC"


def test_unit_plumbing_extract_read_information_uncommon_names():
    read = b"D00689:314:CBYRBANXX:1:1104:4980:1967\t4\t*\t0\t0\t*\t*\t0\t0\tCCTATGGGCGTGGAAAGGACGAAACACCGCATGGCCCTGAAGAATGATGGTT\tBBBBBFFFFFFFFFFFFFBFFFFFFFFFFF<<FFFFFFF</F/<FFFFFFFF\tB2:Z:CGCTATCAA\tQ2:Z:<BBB<///B\tBC:Z:TTTCCCACCCT\tQT:Z:BBBBBFFFFFF\tRG:Z:CBYRBANXX.1\n"
    with tempfile.NamedTemporaryFile() as f:
        f.write(read)
        f.flush()

        b2regex = regex.compile("(?P<cellidentifier1>[ACGT]{9})")
        bcregex = regex.compile("(?P<cellline>[ACGT]{11})")
        queryregex = regex.compile("GGGC(?P<organoid>[ACGT]{5})(?P<cellidentifier2>[ACGT]{30})")
        regex_dictionary = {'B2': b2regex, 'BC': bcregex, 'query': queryregex}

        bamfile = iter(open_sequencing_file(f.name, type='sam'))
        read = next(bamfile)
        features = phantombuster.plumbing.extract_read(regex_dictionary, read)

        assert features["cellline"] == "TTTCCCACCCT"
        assert features["organoid"] == "GTGGA"
        assert features["cellidentifier"][:30] == "CGCTATCAAAAGGACGAAACACCGCATGGC"


def test_unit_plumbing_extract_read_information_as_lt47_experiment():
    read = b"D00689:314:CBYRBANXX:1:1104:13784:1992\t4\t*\t0\t0\t*\t*\t0\t0\tACGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGTT\tBB//BFBFF/FFFFFFFFFFFFFFFFFFB<FF<FFF///<FFFF</FBFFFF\tB2:Z:CGCTATCAA\tQ2:Z:<BBBB/F/F\tBC:Z:CGTATTTTAGT\tQT:Z:B/<<<FFFFFF\tRG:Z:CBYRBANXX.1"
    with tempfile.NamedTemporaryFile() as f:
        f.write(read)
        f.flush()

        b2regex = regex.compile("^[ACGTN]{3}(?P<sample>[ACGTN]{5})")
        bcregex = regex.compile("")
        queryregex = regex.compile("(?P<lid>[ACGTN]{5,6}(?P<lib>ACGT|GTAC){s<=1}[ACGTN]+)")
        regex_dictionary = {'B2': b2regex, 'BC': bcregex, 'query': queryregex}

        bamfile = iter(open_sequencing_file(f.name, type='sam'))
        read = next(bamfile)
        features = phantombuster.plumbing.extract_read(regex_dictionary, read)
        print(features)

        assert features["sample"] == "TATCA"
        assert features["lid"] == "CGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGTT"
        assert features["lib"] == "AAGT"


def test_unit_pyarrow_table_to_counter():
    table1 = pyarrow.Table.from_arrays(
            [
                pyarrow.array(["AAC", "ACG", "GTA"]), 
                pyarrow.array(["sample1", "sample1", "sample1"]),
                pyarrow.array(["lib1", "lib1", "lib2"]),
                pyarrow.array([1, 1, 5]),
            ],names=["lid", "samplename", "lib", "reads"])

    c = pyarrow_table_to_counter(table1, ["lid", "samplename", "lib"])
    print(c)
    assert c[("AAC", "sample1", "lib1")] == 1
    assert c[("ACG", "sample1", "lib1")] == 1
    assert c[("GTA", "sample1", "lib2")] == 5


def test_unit_error_corrector_1():
    ec = ErrorCorrector(4, 0, 0, 1)
    tags = {Umi(0, 0, "AAAA"): UmiCount(1, 1), Umi(0, 0, "AAAT"): UmiCount(10,1), Umi(0, 0, "AATA"): UmiCount(10, 1), Umi(0,0, "ATAA"): UmiCount(10, 1)}
    corr_lids = ec.process(tags)
    # AAAA can be corrected to AAAT, AATA or ATAA and is thus ambigious and should be removed.
    assert Umi(0, 0, "AAAA") not in corr_lids[1]
    assert corr_lids[0][Umi(0,0, "AAAT")].count == 10
    assert corr_lids[0][Umi(0,0, "AATA")].count == 10
    assert corr_lids[0][Umi(0,0, "ATAA")].count == 10
    assert len(corr_lids[0]) == 3


def test_unit_error_corrector_2():
    ec = ErrorCorrector(50, 0, 0, 1)
    tags = {Umi(0, 0, "AAAAGTACGTGTCAGACTGTGAGAGACTGACAGTCAGAGTCTGACAACGT"): UmiCount(112, 1), Umi(0, 0, "AAAAGTACGTGTCAGACTGTGAGAGACTAACAGTCAGAGTCTGACAACGT"): UmiCount(2,1), Umi(0, 0, "AAAAGTACGTGAGAGACTCTGTGTCACTCAGTCTGACAGTCTGTGAACGT"): UmiCount(6, 1), Umi(0,0, "AAAAGTACGTGAGAGACTCCGTGTCACTCAGTCTGACAGTCTGTGAACGT"): UmiCount(2, 1), Umi(0, 0, "AAAAGTACGTGTCAGACTGTGAGAGACTGACAGTCAGAGTCTGGCAACGT"): UmiCount(2, 1), Umi(0,0,"AAAAGTACGTGTCAGACTGTGAGAGACTGACAGTCAGAGTCTTACAACGT"): UmiCount(6, 1), Umi(0,0,"AAAAGTACGTGTGACTCTGTGTGTGTGACTGTGTCTCTGTCTCACAACGT"): UmiCount(40, 1)}
    corr_lids = ec.process(tags)

    for key, count in corr_lids[0].items():
        print(key)
        print(count)
    for lid, correct in corr_lids[1].items():
        print(f"{lid.tag} {corr_lids[0][lid]} -> {correct.tag} : {corr_lids[0][correct]}")
    assert corr_lids[0][Umi(start=0, end=0, tag='AAAAGTACGTGTCAGACTGTGAGAGACTGACAGTCAGAGTCTGACAACGT')].count == 122
    assert corr_lids[0][Umi(start=0, end=0, tag='AAAAGTACGTGAGAGACTCTGTGTCACTCAGTCTGACAGTCTGTGAACGT')].count == 8
    assert corr_lids[0][Umi(start=0, end=0, tag='AAAAGTACGTGTGACTCTGTGTGTGTGACTGTGTCTCTGTCTCACAACGT')].count == 40

    assert sum(tag.count for tag in corr_lids[0].values()) == 170


def test_unit_plumbing_error_correct_1():
    sample = ["S1", "S1", "S1", "S2", "S2", "S2"]
    tags = ["AAAAC", "AAACC", "AACCC", "ACCCC", "CCCCC", "GTAAC"]
    grna = [    "A",     "A",     "A",     "A",     "A",     "A"]
    reads = [     1,       2,     10 ,      1 ,     10 ,      1]
    table = pyarrow.table({"sample": sample, "tags": tags, "reads": reads, "grna": grna})
    for column_name in ['grna', 'tags']:
        table = table.set_column(table.column_names.index(column_name), column_name, table[column_name].cast("large_string"))
    barcode_hierarchy = [{"name": "sample", "type": "reference"}, {"name": "grna", "type": "reference"}, {"name": "tags", "type": "random"}]
    corrected_table = plumbing.error_correct(table, barcode_hierarchy, 1).to_pydict()
    assert corrected_table["sample"] == ["S1", "S2", "S2"]
    assert corrected_table["tags"]  == ["AACCC", "CCCCC", "GTAAC"]
    assert corrected_table["reads"] == [     13,      11,       1]
    assert corrected_table["grna"] == ["A", "A", "A"]


def test_unit_plumbing_error_correct_2():
    tags = ["AAAAC", "AAACC", "AACCC", "ATTTT", "TTTTT", "GTAAC"]
    reads = [     1,       2,     10 ,      1 ,     10 ,      1]
    table = pyarrow.table({"tags": tags, "reads": reads})
    for column_name in ['tags']:
        table = table.set_column(table.column_names.index(column_name), column_name, table[column_name].cast("large_string"))
    barcode_hierarchy = [{"name": "tags", "type": "random"}]
    corrected_table = plumbing.error_correct(table, barcode_hierarchy, 1).to_pydict()
    assert corrected_table["tags"]  == ["AACCC", "GTAAC", "TTTTT"]
    assert corrected_table["reads"] == [     13,        1,       11]


def test_unit_plumbing_error_correct_3_first_base():
    sample = ["S1", "S1", "S1"]
    tags = ["AAAAA", "AAACC", "CAAAA"]
    grna = [    "A",     "A",     "A"]
    reads = [    10,       10,     1 ]
    table = pyarrow.table({"sample": sample, "tags": tags, "reads": reads, "grna": grna})
    for column_name in ['grna', 'tags']:
        table = table.set_column(table.column_names.index(column_name), column_name, table[column_name].cast("large_string"))
    barcode_hierarchy = [{"name": "sample", "type": "reference"}, {"name": "grna", "type": "reference"}, {"name": "tags", "type": "random"}]
    corrected_table = plumbing.error_correct(table, barcode_hierarchy, 1).to_pydict()
    assert corrected_table["sample"] == ["S1", "S1"]
    assert corrected_table["tags"]  == ["AAAAA", "AAACC"]
    assert corrected_table["reads"] == [     11,      10]
    assert corrected_table["grna"] == ["A", "A"]


def test_unit_plumbing_apply_correction():
    tags = ["AAAAC", "AAACC", "AACCC", "ACCCC", "CCCCC", "GTAAC"]
    reads = [     1,       2,     10 ,      1 ,     10 ,      1]
    corrections = {"AAAAC": "AACCC", "AAACC": "AACCC", "AACCC": "AACCC", "CCCCC": "CCCCC", "GTAAC": "GTAAC"}
    table = pyarrow.table({"tags": tags, "reads": reads})
    table = plumbing.apply_correction(table, "tags", corrections)
    d = table.to_pydict()
    assert d["tags"]  == ["AACCC", "CCCCC", "GTAAC"]
    assert d["reads"] == [     13,      10,       1]


def test_unit_plumbing_sort_and_partition():
    table = pyarrow.Table.from_pydict({
        'sample': ['A', 'A', 'B', 'B'],
        'reads' : [ 1 ,  1,   2 ,  3 ],
        'lib'   : ['a', 'b', 'a', 'a'],
        })

    partitions = plumbing.sort_and_partition(table, [('sample', 'ascending'), ('lib', 'ascending')], ('sample', 'lib'))
    print(partitions)
    assert ('A', 'a') in partitions
    assert ('A', 'b') in partitions
    assert ('B', 'a') in partitions
    assert ('B', 'b') not in partitions

    assert np.sum(partitions[('A', 'a')]['reads'].to_numpy()) == 1
    assert np.sum(partitions[('A', 'b')]['reads'].to_numpy()) == 1
    assert np.sum(partitions[('B', 'a')]['reads'].to_numpy()) == 5

    assert len(partitions[('A', 'a')]['reads'].to_numpy()) == 1
    assert len(partitions[('A', 'b')]['reads'].to_numpy()) == 1
    assert len(partitions[('B', 'a')]['reads'].to_numpy()) == 2


def test_unit_extract_section():

    fixed = "GGGC"

    b2regex = "(?P<umi>[ACGT]{9})"
    bcregex = "(?P<grna>[ACGT]{11})"
    queryregex = "(?P<lib>GGGC)(?P<sample>[ACGT]{5})(?P<lid>[ACGT]{30})"
    regex_dictionary = phantombuster.config_files.RegexDictionary()
    regex_dictionary.add_regex('B2', b2regex)
    regex_dictionary.add_regex('BC', bcregex)
    regex_dictionary.add_regex('query', queryregex)

    samples = ["GTGGA", "TTTGA"]
    sample_dict = {s: f"Sample{i}" for i, s in enumerate(samples)}

    samples0_mutated = "GGGGA"

    genes = ["TTTCCCACCCT"]
    gene_dict = {s: f"Gene{i}" for i, s in enumerate(genes)}

    lib_dict = {"GGGC": "lib1"}

    lid1 = "AAGGACGAAACACCGCATGGCCCTGAAGAA"
    lid2 = "CAGGACGAAACACCGCATGGCCCTGAAGAA" # hamming distance 1 to lid1

    umi1 = "CGCTATCAA"
    umi2 = "AAGTATCAA"
    umi3 = "CGCTAGCCT"

    reads = [
            make_read(f"{fixed}{samples[0]}{lid1}", umi1, genes[0]),
            make_read(f"{fixed}{samples0_mutated}{lid1}", umi1, genes[0]), # should be collapsed with read 1
            make_read(f"{fixed}{samples[0]}{lid1}", umi1, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid1}", umi2, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid2}", umi1, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid2}", umi3, genes[0]),
            make_read(f"{fixed}{samples[1]}{lid2}", umi3, genes[0]),
            ]


    barcode_hierarchy = [{"name": "lib", "type": "reference", "reference": lib_dict, "referencefile": None, "threshold": 0, "min_length": None, "max_length": None},
                         {"name": "sample", "type": "reference", "reference": sample_dict, "threshold": 1, "min_length": None, "max_length": None},
                         {"name": "grna", "type": "reference", "reference": gene_dict, "threshold": 0, "min_length": None, "max_length": None},
                         {"name": "lid", "type": "auto", "threshold": "auto", "min_length": None, "max_length": None},
                         {"name": "umi", "type": "auto", "threshold": "auto", "min_length": None, "max_length": None},
                         ]


    with tempfile.NamedTemporaryFile(suffix='.sam') as f:
        for read in reads:
            f.write(read)
        f.flush()

        lidonly = False
        count_qc = False

        input_file = phantombuster.config_files.InputFile(f.name, '', '')

        f = phantombuster.io_.SamSequencingFile(input_file.path)
        assert len(list(f)) == len(reads)

        input_group = phantombuster.config_files.InputGroup([f], "")

        table, stats = phantombuster.plumbing.deduplicate_section(input_group, None, regex_dictionary, barcode_hierarchy, count_qc)

    cols = table.column_names
    cols.remove("reads")

    c = pyarrow_table_to_counter(table, cols)

    assert sum(c.values()) == 7

    for k, v in c.items():
        print(k, v)

    assert c[("lib1", "Sample0", "Gene0", lid1, umi1)] == 3
    assert c[("lib1", "Sample0", "Gene0", lid1, umi2)] == 1
    assert c[("lib1", "Sample0", "Gene0", lid2, umi1)] == 1
    assert c[("lib1", "Sample0", "Gene0", lid2, umi3)] == 1
    assert c[("lib1", "Sample1", "Gene0", lid2, umi3)] == 1

    assert sum(c.values()) == len(reads)


def test_unit_deduplication_persister():

    from phantombuster import stores
    from phantombuster.remoter.task import Task

    result = [{("lid1", "lib1"): 1, ("lid2", "lib1"): 2}, [3, {"sample": {"sample1": 3}}, {"sample":0} ]]

    result[0] = deduplicator_to_pyarrow_table(Counter(result[0]), [{'name': "lid", 'type': 'random'}, {'name': "lib", 'type': 'reference'}])

    with tempfile.TemporaryDirectory() as tmpdirname:
        id = stores.save(result, dir=tmpdirname)
        v = stores.load(id, dir=tmpdirname)

    r1 = pyarrow_table_to_counter(v[0], ["lid", "lib"])
    r2 = pyarrow_table_to_counter(v[0], ["lid", "lib"])

    assert r1[("lid1", "lib1")] == r2[("lid1", "lib1")]
    assert r1[("lid2", "lib1")] == r2[("lid2", "lib1")]

    assert r1 == r2


def test_unit_combine():

    fixed = "GGGC"

    b2regex = "(?P<umi>[ACGT]{9})"
    bcregex = "(?P<grna>[ACGT]{11})"
    queryregex = "(?P<lib>GGGC)(?P<sample>[ACGT]{5})(?P<lid>[ACGT]{30})"
    regex_dictionary = phantombuster.config_files.RegexDictionary()
    regex_dictionary.add_regex('B2', b2regex)
    regex_dictionary.add_regex('BC', bcregex)
    regex_dictionary.add_regex('query', queryregex)

    samples = ["GTGGA"]
    sample_dict = {s: f"Sample{i}" for i, s in enumerate(samples)}

    genes = ["TTTCCCACCCT"]
    gene_dict = {s: f"Gene{i}" for i, s in enumerate(genes)}

    lib_dict = {"GGGC": "lib1"}

    lid1 = "AAGGACGAAACACCGCATGGCCCTGAAGAA"
    lid2 = "CAGGACGAAACACCGCATGGCCCTGAAGAA"

    umi1 = "CGCTATCAA"
    umi2 = "AAGTATCAA"
    umi3 = "CGCTAGCCT"

    reads = [
            make_read(f"{fixed}{samples[0]}{lid1}", umi1, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid1}", umi1, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid1}", umi2, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid2}", umi1, genes[0]),
            make_read(f"{fixed}{samples[0]}{lid2}", umi3, genes[0]),
            ]

    #                    error is at line in read3                                         |


    barcode_hierarchy = [{"name": "lib", "type": "reference", "reference": lib_dict, "referencefile": None, "threshold": 0, "min_length": None, "max_length": None},
                         {"name": "sample", "type": "reference", "reference": sample_dict, "threshold": 0, "min_length": None, "max_length": None},
                         {"name": "grna", "type": "reference", "reference": gene_dict, "threshold": 0, "min_length": None, "max_length": None},
                         {"name": "lid", "type": "auto", "threshold": "auto", "min_length": None, "max_length": None},
                         {"name": "umi", "type": "auto", "threshold": "auto", "min_length": None, "max_length": None},
                         ]

    with tempfile.NamedTemporaryFile(suffix='.sam') as f:
        for read in reads:
            f.write(read)
        f.flush()

        lidonly = False
        count_qc = False

        input_file = phantombuster.config_files.InputFile(f.name, '', '')

        f = phantombuster.io_.SamSequencingFile(input_file.path)
        assert len(list(f)) == len(reads)

        input_group = phantombuster.config_files.InputGroup([f], "")

        out1 = phantombuster.plumbing.deduplicate_section(input_group, None, regex_dictionary, barcode_hierarchy, count_qc)
        out2 = phantombuster.plumbing.deduplicate_section(input_group, None, regex_dictionary, barcode_hierarchy, count_qc)

    table, stats = phantombuster.plumbing.combine([out1, out2], barcode_hierarchy)

    cols = table.column_names
    cols.remove("reads")

    c = pyarrow_table_to_counter(table, cols)

    for k, v in c.items():
        print(k, v)

    assert c[("lib1", "Sample0", "Gene0", lid1, umi1)] == 4
    assert c[("lib1", "Sample0", "Gene0", lid1, umi2)] == 2
    assert c[("lib1", "Sample0", "Gene0", lid2, umi1)] == 2
    assert c[("lib1", "Sample0", "Gene0", lid2, umi3)] == 2

    assert sum(c.values()) == 2*len(reads)


def test_unit_combine_plumbing():

    table1 = pyarrow.Table.from_arrays(
            [
                pyarrow.array(["AAC", "ACG", "GTA"], type=pyarrow.large_string()), 
                pyarrow.array(["sample1", "sample1", "sample1"], type=pyarrow.large_string()),
                pyarrow.array([1, 1, 5], type=pyarrow.uint64()),
            ],names=["lid", "sample", "reads"])

    table2 = pyarrow.Table.from_arrays(
            [
                pyarrow.array(["AAC", "ACG", "GTA"], type=pyarrow.large_string()), 
                pyarrow.array(["sample1", "sample1", "sample2"], type=pyarrow.large_string()),
                pyarrow.array([1, 1, 5], type=pyarrow.uint64()),
            ],names=["lid", "sample", "reads"])
    bhc = [{"name": "sample", "type": "reference"}, {"name": "lid", "type": "random"}]
    combined = phantombuster.plumbing._combine_duplicators_cython(table1, table2, bhc)


    print(combined)
    assert len(combined) == 4
    assert combined["sample"][3].as_py() == "sample2"



def test_unit_plumbing_correction():
    tags  = ["AAAA", "AAAT", "ACGT", "AATT", "AACC", "AACT"]
    reads = [     32,      8,      4,      2,    16 ,     1]

    corrections = plumbing.calculate_corrections(tags, reads, 4, 1)

    assert corrections["AAAA"] == "AAAA"
    assert corrections["AAAT"] == "AAAA"
    assert corrections["ACGT"] == "ACGT"
    assert corrections["AATT"] == "AAAA"
    assert "AACT" not in corrections
    assert corrections["ACGT"] == "ACGT"
    assert len(corrections) == 5


def test_unit_plumbing_calculate_overlap():
    table1_lids = ["AAAA", "ACGT", "GGTT"]
    table2_lids = ["AAGT", "ACGT", "GCTT", "TTTG"]
    table1_reads= [11, 10, 8]
    table2_reads= [40, 13, 3, 1]

    table1 = pa.Table.from_pydict({"lid": np.array(table1_lids), "reads": table1_reads, "samplename": np.array(["s1"]*3)})
    table2 = pa.Table.from_pydict({"lid": np.array(table2_lids), "reads": table2_reads, "samplename": np.array(["s2"]*4)})

    mastertable = pa.concat_tables([table1, table2])
    sorted_table, overlap = plumbing.calculate_overlap(mastertable, "samplename")

    print("sorted_table", sorted_table.to_pandas())
    print("overlap", overlap)

    assert len(overlap) == len(table1) + len(table2)
    assert len(sorted_table) == len(overlap)
    out_table = plumbing.hopping_removal(sorted_table, overlap, 0.5)
    partitions = plumbing.sort_and_partition(out_table)
    for p in partitions:
        print(p, partitions[p].to_pandas())

    assert len(partitions["s1"]) == 3
    assert len(partitions["s2"]) == 4
    assert all(partitions["s1"]["reads"].to_numpy() == np.array(table1_reads))
    assert all(partitions["s2"]["reads"].to_numpy() == np.array(table2_reads))
    assert all(partitions["s1"]["lid"].to_numpy() == np.array(table1_lids))
    assert all(partitions["s2"]["lid"].to_numpy() == np.array(table2_lids))

    out_table = plumbing.hopping_removal(sorted_table, overlap, 1.0)
    print(out_table.to_pandas())
    partitions = plumbing.sort_and_partition(out_table)

    assert len(partitions["s1"]) == 2
    assert len(partitions["s2"]) == 4
    assert all(partitions["s1"]["reads"].to_numpy() == np.array([11, 8]))
    assert all(partitions["s2"]["reads"].to_numpy() == np.array(table2_reads))
    assert all(partitions["s1"]["lid"].to_numpy() == np.array(["AAAA", "GGTT"]))
    assert all(partitions["s2"]["lid"].to_numpy() == np.array(table2_lids))


def test_unit_plumbing_calculate_overlap_multiple():
    table1_lids = ["AAAA", "ACGT", "GCTT"]
    table1_lib  = [   "A",   "A",     "A"]
    table2_lids = ["AAGT", "ACGT", "GCTT", "TTTG"]
    table2_lib  = [   "A",   "A",     "B",    "B"]
    table1_reads= [11, 10, 8]
    table2_reads= [40, 13, 3, 1]

    table1 = pa.Table.from_pydict({"lid": np.array(table1_lids), "reads": table1_reads, "lib": table1_lib, "samplename": np.array(["s1"]*3)})
    table2 = pa.Table.from_pydict({"lid": np.array(table2_lids), "reads": table2_reads, "lib": table2_lib, "samplename": np.array(["s2"]*4)})

    mastertable = pa.concat_tables([table1, table2])
    sorted_table, overlap = plumbing.calculate_overlap(mastertable, "samplename", ["lid", "lib"])

    print("sorted_table", sorted_table.to_pandas())
    print("overlap", overlap)

    assert len(overlap) == len(table1) + len(table2)
    assert len(sorted_table) == len(overlap)
    out_table = plumbing.hopping_removal(sorted_table, overlap, 0.5)
    partitions = plumbing.sort_and_partition(out_table)
    for p in partitions:
        print(p, partitions[p].to_pandas())

    assert len(partitions["s1"]) == 3
    assert len(partitions["s2"]) == 4
    assert all(partitions["s1"]["reads"].to_numpy() == np.array(table1_reads))
    assert all(partitions["s2"]["reads"].to_numpy() == np.array(table2_reads))
    assert all(partitions["s1"]["lid"].to_numpy() == np.array(table1_lids))
    assert all(partitions["s2"]["lid"].to_numpy() == np.array(table2_lids))

    out_table = plumbing.hopping_removal(sorted_table, overlap, 1.0)
    print(out_table.to_pandas())
    partitions = plumbing.sort_and_partition(out_table)

    assert len(partitions["s1"]) == 2
    assert len(partitions["s2"]) == 4
    assert all(partitions["s1"]["reads"].to_numpy() == np.array([11, 8]))
    assert all(partitions["s2"]["reads"].to_numpy() == np.array([40, 13, 3, 1]))
    assert all(partitions["s1"]["lid"].to_numpy() == np.array(["AAAA", "GCTT"]))
    assert all(partitions["s2"]["lid"].to_numpy() == np.array(table2_lids))


def test_unit_plumbing_calculate_jaccard():

    table = pyarrow.Table.from_arrays(
            [
                pyarrow.array(["AAC", "ACG", "AAC"]), 
                pyarrow.array(["sample1", "sample1", "sample2"]),
                pyarrow.array(["lib2", "lib1", "lib2"]),
                pyarrow.array([1, 1, 5]),
            ],names=["lid", "samplename", "lib", "reads"])

    jaccard = plumbing.calculate_jaccard(table, 'samplename')

    assert jaccard[('sample1', 'sample2')] == 1/2

    jaccard = plumbing.calculate_jaccard(table, 'lib')

    assert len(jaccard) == 0


def test_unit_calculate_correction():
    lids = ['AAACCCGGGTTT', 'GAACCCGGGTTT', 'GGACCCGGGTTT', 'GGGCCCGGGTTT', 'GGGTTTAAACCC']
    reads = [10, 3, 1, 10, 11]

    corrections = plumbing.calculate_corrections(np.array(lids), np.array(reads), len(lids[0]), 1)

    assert len(corrections) == 4
    assert corrections['AAACCCGGGTTT'] == 'AAACCCGGGTTT'
    assert corrections['GAACCCGGGTTT'] == 'AAACCCGGGTTT'
    assert corrections['GGGCCCGGGTTT'] == 'GGGCCCGGGTTT'
    assert corrections['GGGTTTAAACCC'] == 'GGGTTTAAACCC'
    assert corrections.get('GGACCCGGGTTT', None) == None


def test_unit_xmerge_plumbing():

    table1 = pyarrow.Table.from_arrays(
            [
                pyarrow.array(["AAC", "ACG", "GTA"], type=pyarrow.large_string()), 
                pyarrow.array(["sample1", "sample1", "sample1"], type=pyarrow.large_string()),
                pyarrow.array([1, 3, 5], type=pyarrow.uint64()),
            ],names=["lid", "sample", "reads"])

    table2 = pyarrow.Table.from_arrays(
            [
                pyarrow.array(["AAC", "ACG", "GTA"], type=pyarrow.large_string()), 
                pyarrow.array(["sample1", "sample1", "sample2"], type=pyarrow.large_string()),
                pyarrow.array([9, 17, 33], type=pyarrow.uint64()),
            ],names=["lid", "sample", "reads"])
    bhc = [{"name": "sample", "type": "reference"}, {"name": "lid", "type": "random"}]

    length = len(table1) + len(table2)
    idx = np.zeros(length, dtype=np.uint)
    counts = np.zeros(length, dtype=np.uint)
    column_names = table1.column_names
    column_names.remove('reads')
    cols = [table1.column_names.index(bc['name']) for bc in bhc]
    print('cols', cols)
    idx_r, counts_r, idx_i, i, j, state = xmerge(table1, table2, idx, counts, cols)

    assert state == 0
    idx = idx[:idx_i]
    counts = counts[:idx_i]


    print("index", idx)
    print('counts', counts)
    print('idx_i', idx_i)

    assert idx_i == 4


def test_unit_deduplicate_single():
    sample = ['s1', 's1', 's1', 's2', 's2']
    lids =   [ 'A',  'A',  'A',  'A',  'B']
    umis =   [ 'X',  'X',  'Y',  'X',  'X']
    reads =  [  1 ,   2 ,   1 ,   2,    5 ]


    sample, lids, umis, reads = np.array(sample), np.array(lids), np.array(umis), np.array(reads)

    table = pyarrow.Table.from_arrays([sample, lids, umis, reads], names=['sample', 'lid', 'umi', 'reads'])
    t = plumbing.deduplicate_single(table)

    c = {col: t[col].to_numpy() for col in t.column_names}

    assert np.sum(c['reads']) == np.sum(reads)
    assert len(c['reads']) == 4


def test_unit_remove_ambigious():
    sample = [ 's1',  's2']
    lids =   ['ACG', 'NNN']
    reads =  [   3 ,    2 ]

    barcode_list= [{"name": "sample", "type": "reference", "reference": None, "referencefile": None, "threshold": 0, "min_length": None, "max_length": None},
                   {"name": "lid", "type": "random", "threshold": "auto", "min_length": None, "max_length": None},
                  ]

    sample, lids, reads = np.array(sample), np.array(lids), np.array(reads)

    table = pyarrow.Table.from_arrays([sample, lids, reads], names=['sample', 'lid', 'reads'])
    table = plumbing.remove_ambigious(table, barcode_list)

    assert np.sum(table['reads']) == 3
    assert len(table) == 1


@pytest.mark.slow
def test_unit_combine_large_tables():
    bc_hierarchy = [{"name": "LID"}]
    N = 10000000
    a = pyarrow.table({"reads": np.array([1]*N), "LID": np.array(["A"*20]*N)})
    a = a.set_column(a.column_names.index("LID"), "LID", a["LID"].cast("large_string"))
    b = pyarrow.table({"reads": np.array([1]*N), "LID": np.array(["A"*20]*N)})
    b = b.set_column(b.column_names.index("LID"), "LID", b["LID"].cast("large_string"))
    d = phantombuster.plumbing._combine_duplicators_cython(a, b, bc_hierarchy)
    assert len(d) == N


def test_unit_statistics():
    df = pl.DataFrame({'grna': ['A', 'A', 'A', 'B', 'C', 'C'],
                       'category': ['control', 'control', 'control', 'test', 'test', 'test'],
                       'reads': [10, 20, 30, 25, 5, 5]})

    r = phantombuster.statistics.calculate_pvalues(df, column='reads', group_by=['grna'], N=10000)
    for row in r.iter_rows(named=True):
        if row['grna'] == 'B':
            assert row['pvalue'] - 0.66 < 0.05
            assert row['lineages'] == 1
        elif row['grna'] == 'C':
            assert row['pvalue'] == 0.0
            assert row['lineages'] == 2


# -- Integration test helpers -- #


def prepare_dir(template):
    # make temporary directory for all data
    tmpdir = pjoin("/tmp", uuid.uuid4().hex)
    os.mkdir(tmpdir)

    template_dir = pjoin(os.path.dirname(os.path.abspath(__file__)), 'test_data_files', template)
    files = glob.glob(pjoin(template_dir, '*'))

    for p in files:
        p = os.path.basename(p)
        copy_file_efficiently(pjoin(template_dir, p), pjoin(tmpdir, p))

    return tmpdir


def execute_in(tmpdir):
    old_wd = os.getcwd()
    os.chdir(tmpdir)
    yield 'input_files.csv', 'regexes.csv'
    os.chdir(old_wd)


# -- Integration Test Fixtures -- #

@pytest.fixture()
def setup_lineage_data_small():
    tmp_dir = prepare_dir('small')
    yield from execute_in(tmp_dir)


@pytest.fixture()
def setup_lineage_data_small_fastq():
    tmp_dir = prepare_dir('small_fastq')
    yield from execute_in(tmp_dir)


@pytest.fixture()
def setup_lineage_data_medium():
    tmp_dir = prepare_dir('medium')
    yield from execute_in(tmp_dir)


@pytest.fixture()
def setup_lineage_data_big():
    tmp_dir = prepare_dir('big')
    yield from execute_in(tmp_dir)


@pytest.fixture()
def setup_lineage_data_1010():
    tmp_dir = prepare_dir('1010')
    yield from execute_in(tmp_dir)



# -- Integration Tests -- #


def test_integration_read_regex_file(setup_lineage_data_small):
    regex_dict = phantombuster.config_files.read_regex_file('regexes.csv')
    regexes = regex_dict.get_regexes_for_group('*')
    assert regexes['B2'] == regex.compile('^[ACGTN]{3}(?P<sample>[ACGTN]{5})')
    assert regexes['query'] == regex.compile('(?P<lid>[ACGTN]{5,6}(?P<lib>ACGT|GTAC){s<=1}[ACGTN]+)')

    assert regexes['B2'].match('AAACGTAGC').group('sample') == 'CGTAG'
    assert regexes['query'].match('AAAAAGTACAAAAAAAAAAA').group('lib') == 'GTAC'


def test_integration_read_input_files_file(setup_lineage_data_small):
    input_groups = phantombuster.config_files.read_input_files_file('input_files.csv')
    assert phantombuster.config_files.InputFile('101.bam', '101.bam', '') in input_groups[0].files


def test_integration_read_barcode_hierarchy(setup_lineage_data_small):
    bch = phantombuster.config_files.read_barcode_hierarchy_file('barcode_hierarchy.csv')
    assert [bc["name"] for bc in bch] == ['lib', 'sample', 'lid']
    assert [bc["type"] for bc in bch] == ['reference', 'reference', 'random']


def test_integration_indexing_bam(setup_lineage_data_1010):
    bamfiles, regex_file = setup_lineage_data_1010
    idx = index_bamfile("1010.bam", 1)
    every_read_index = []
    f = open_sequencing_file('1010.bam')

    for read in f:
        every_read_index.append(f.tell())

    with_index = []
    for start, end in pairwise(idx):
        with open_sequencing_file('1010.bam') as f:
            f.seek(start)
            for read in f:
                with_index.append(f.tell())
                if f.tell() >= end:
                    break

    assert with_index == every_read_index


def test_integration_indexing_fastq(setup_lineage_data_1010):
    bamfiles, regex_file = setup_lineage_data_1010
    idx = index_fastqgz_file("1010.fastq.gz", 1)
    with gzip.open("1010.fastq.gz") as f:
        lines = f.readlines()

    every_read_index = []
    i = 0
    index = 0
    while True:
        every_read_index.append(index)
        if i+4 < len(lines):
            index += sum(len(l) for l in lines[i:i+4])
            i += 4
        else:
            break

    assert idx[:-1] == every_read_index
   

from phantombuster.io_ import open_sequencing_file

def test_integration_read_fastq(setup_lineage_data_1010):
    bamfiles, regex_file = setup_lineage_data_1010


    with open_sequencing_file('1010.fastq.gz') as f:
        fi = iter(f)
        read = next(fi)
        assert read['name'] == '@D00689:314:CBYRBANXX:1:1104:13784:1992'
        assert read['seq'] == 'ACGGGTGAAGTGGAAAGGACGACACACCGGTACAAGATCAACAACAAAGGTT'


def test_integration_read_bam(setup_lineage_data_1010):
    bamfiles, regex_file = setup_lineage_data_1010

    with open_sequencing_file('1010.bam') as f:
        fi = iter(f)
        read = next(fi)
        assert read['query'] == 'ACGGGTGAAGTGGAAAGGACGACACACCGGTACAAGATCAACAACAAAGGTT'
        assert read['B2'] == 'CGCTATCAA'


@pytest.mark.slow
def test_integration_fullrun_only(setup_lineage_data_small):
    input_files_file, regex_file = setup_lineage_data_small
    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)

    result = pa.parquet.read_table(pjoin("out", 'data', "demultiplex.parquet")).to_pandas()
    r1 = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGT"]
    r2 = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGACACACCGGTACAAGATCAACAACAAAGGT"]

    assert r1["lib"].iloc[0] == "lib1"
    assert r1["reads"].iloc[0] == 100

    assert r2["lib"].iloc[0] == "lib1"
    assert r2["reads"].iloc[0] == 1


@pytest.mark.slow
def test_integration_fullrun_fastq_only(setup_lineage_data_small_fastq):
    input_files_file, regex_file = setup_lineage_data_small_fastq

    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)

    result = pa.parquet.read_table(pjoin("out", 'data', "demultiplex.parquet")).to_pandas()
    r1 = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGT"]
    r2 = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGACACACCGGTACAAGATCAACAACAAAGGT"]

    assert r1["lib"].iloc[0] == "lib1"
    assert r1["reads"].iloc[0] == 100

    assert r2["lib"].iloc[0] == "lib1"
    assert r2["reads"].iloc[0] == 1


@pytest.mark.slow
def test_integration_fullrun_many_only(setup_lineage_data_medium):
    input_files_file, regex_file = setup_lineage_data_medium
    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)

    result = pa.parquet.read_table(pjoin("out", 'data', "demultiplex.parquet")).to_pydict()
    print(result)
    idx = result["lid"] == "GGGCGCCGAAAGGGACATTCCAGCAGCTCT"
    assert result["lib"][idx] == "lib1"
    assert result["reads"][idx] == 1025

    assert sum(result['reads']) == 10000


@pytest.mark.slow
def test_integration_fullrun_no_prefix(setup_lineage_data_small):
    input_files_file, regex_file = setup_lineage_data_small

    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)

    result = pa.parquet.read_table(pjoin("out", 'data', "demultiplex.parquet")).to_pandas()
    r1 = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGT"]
    r2 = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGACACACCGGTACAAGATCAACAACAAAGGT"]
    assert r1["lib"][0] == "lib1"
    assert r1["reads"][0] == 100

    assert r2["lib"].iloc[0] == "lib1"
    assert r2["reads"].iloc[0] == 1


@pytest.mark.slow
def test_integration_fullrun_multiple_files(setup_lineage_data_big):
    input_files_file, regex_file = setup_lineage_data_big
    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)

    result = pa.parquet.read_table(pjoin("out", 'data', "demultiplex.parquet")).to_pandas()
    for idx, data in result.iterrows():
        print(*[data[col] for col in result.columns])
    r = result[result["lid"] == "CGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGT"]
    assert r["lib"][0] == "lib1"
    assert r["reads"][0] == 500

@pytest.mark.slow
def test_integration_pvalue(setup_lineage_data_big):
    input_files_file, regex_file = setup_lineage_data_big
    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)
    core.error_correct(project, error_threshold=1, barcode_hierarchy_file="barcode_hierarchy.csv", remove_ambigious=True)
    core.hopping_removal(project, ['lid'], 0.05)
    core.threshold(project, 'thresholds.csv')
    core.pvalue(project, project.threshold_output_path, 'categories.csv', group_by=['sample'])


@pytest.mark.slow
def test_integration_fullrun_with_ec(setup_lineage_data_big):
    input_files_file, regex_file = setup_lineage_data_big
    project = Project("out")

    core.demultiplex(input_files_file, regex_file=regex_file, barcode_hierarchy_file="barcode_hierarchy.csv", project=project, debug=True, show_qc=False)

    print("Finished phantombuster demultiplex")

    print("E16: All files in temporary directory")
    for root, dirs, files in os.walk("."):
        print(root)
        for f in files:
            print(pjoin(root, f))

    core.error_correct(project=project, error_threshold=1, barcode_hierarchy_file='barcode_hierarchy.csv', remove_ambigious=True)

    print("E15: All files in temporary directory")
    for root, dirs, files in os.walk("."):
        print(root)
        for f in files:
            print(pjoin(root, f))

    before_ec = pa.parquet.read_table(pjoin("out", 'data', "demultiplex.parquet"))

    for lib, sample, lid, reads in zip(before_ec["lib"].to_pylist(), before_ec["sample"].to_pylist(), before_ec["lid"].to_pylist(), before_ec["reads"].to_pylist()):
        print(lib, sample, lid, reads)

    result_orig = pa.parquet.read_table(pjoin("out", 'data', "error_correct.parquet"))
    df = pl.from_arrow(result_orig)

                                               # ERROR ->@<-# 
                                   #'CGGGTGAAGTGGAAAGGACGACACACCGGTACAAGATCAACAACAAAGGT'
    r = df.filter(pl.col('lid') == "CGGGTGAAGTGGAAAGGACGAAACACCGGTACAAGATCAACAACAAAGGT")

    assert r["lib"][0] == "lib1"
    assert r["reads"][0] == 508

    assert len(df['sample'].unique()) == 2

    for column_name in result_orig.column_names:
        print(column_name, result_orig[column_name].type)

    #subprocess.run(["phantombuster", "--verbose", "hopping_removal", "--outdir=out"], timeout=20, check=True)


    print("E17 All files in temporary directory")
    for root, dirs, files in os.walk("."):
        print(root)
        for f in files:
            print(pjoin(root, f))

    #result = pa.parquet.read_table(pjoin("out", "data", "hopping_removal", "25-H9-day00-01.parquet")).to_pandas()

@pytest.mark.slow
def test_integration_logger_handler():

    zmq_log_handler = None

    try:
        logger = logging.getLogger()
        #configure_logging(outputlog="worker" + name, verbose=True)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(process)s - %(levelname)s - %(message)s")
        logger.setLevel(logging.DEBUG)

        context = zmq.Context()

        zmq_log_handler = PushHandler("inproc://test_logger", context=context)
        zmq_log_handler.setLevel(logging.DEBUG)
        logger.addHandler(zmq_log_handler)

        thread = threading.Thread(target=logging.debug, args=["hello"])
        thread.start()

        socket = context.socket(zmq.PULL)
        socket.bind("inproc://test_logger")
        msg = None
        count = 0
        while msg is None and count <= 10:
            time.sleep(1)
            try:
                msg = socket.recv(flags=zmq.NOBLOCK)
            except zmq.error.Again as e:
                count += 1
        if msg is None:
            msg = socket.recv(flags=zmq.NOBLOCK)

        assert b"hello" in msg
    finally:
        if zmq_log_handler:
            logger.removeHandler(zmq_log_handler)
