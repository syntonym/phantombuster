from typing import Optional
import os.path
import logging
import pysam
import gzip
from dataclasses import dataclass
import pyarrow.parquet


STAGES = ["deduplication", "error_corrected", "hopping_removal", "thresholded"]

@dataclass
class PathsAndFiles:
    outdir: str
    indexdir: str
    summarydir: str
    datadir: str
    tmpdir: str
    prefix: Optional[str]
    bamfiles: Optional[list[str]]
    tsv_filename: Optional[str]
    sort_filename: Optional[str]
    rest_filename: Optional[str]
    grna_filename: Optional[str]
    indexfile: Optional[str]
    task_db_filename: Optional[str]

    def __init__(self, outdir: str, prefix: Optional[str], bamfiles: Optional[list[str]]) -> None:

        assert prefix not in ["index", "summary", "data", "plot", "tmp"]

        self.outdir = outdir
        self.indexdir = os.path.join(self.outdir, "index")
        self.summarydir = os.path.join(self.outdir, "summary")
        self.datadir = os.path.join(self.outdir, "data")
        self.plotdir = os.path.join(self.outdir, "plot")
        self.tmpdir = os.path.join(self.outdir, "tmp")

        self.prefix = prefix
        self.bamfiles = bamfiles

        if self.prefix:
            self.tsv_filename = self.prefixs + ".tsv"
        else:
            if self.bamfiles:
                self.tsv_filename = self.bamfiles[0].replace(".bam", ".tsv")
            else:
                self.tsv_filename = None

        if self.tsv_filename is not None:
            self.indexfile = self.tsv_filename.replace(".tsv", ".index")
            self.sort_filename = self.tsv_filename.replace(".tsv", ".sort.tsv.gz")
            self.rest_filename = self.tsv_filename.replace(".tsv", ".rest.tsv.gz")
            self.grna_filename = self.tsv_filename.replace(".tsv", ".grna.tsv.gz")
        else:
            self.indexfile = None
            self.sort_filename = None
            self.rest_filename = None
            self.grna_filename = None

        self.task_db_filename = os.path.join(self.tmpdir, ".tasks.sqlite3")

    @property
    def prefixs(self):
        if self.prefix is None:
            return ""
        else:
            return self.prefix

    def in_dir(self, filename):
        path = [self.outdir]
        if self.prefix:
            path.append(self.prefix)
        path.apend(filename)
        return os.path.join(*path)

    def in_data(self, filename):
        path = [self.datadir]
        if self.prefix:
            path.append(self.prefix)
        path.apend(filename)
        return os.path.join(*path)

    def stage_path(self, stage):
        if stage not in STAGES:
            raise Exception(f"Invalid stage {stage}, valid stages are {STAGES}")
        path = [self.datadir]
        if self.prefix:
            path.append(self.prefix)
        path.append(stage)
        return os.path.join(*path)

    def create(self) -> None:
        for d in [self.indexdir, self.summarydir, self.datadir, self.plotdir, self.tmpdir]:
            os.makedirs(d, exist_ok=True)
        for rel_path in STAGES:
            f = [self.datadir]
            if self.prefix:
                f.append(self.prefix)
            f.append(rel_path)
            d = os.path.join(*f)
            os.makedirs(d, exist_ok=True)


def generate_directory_structure(dirs):
    # generating directory structure
    logging.info("Creating output directories")
    for direct in [
        dirs.outdir,
        dirs.datadir,
        dirs.indexdir,
        dirs.summarydir,
        dirs.plotdir,
        dirs.tmpdir,
    ]:
        logging.info(f"Creating directory {direct}.")
        # creating output directory if not already created
        if not os.path.exists(direct):
            os.mkdir(direct)
        else:
            logging.info("output directory already exists, no need to create it")


def open_sequencing_file(path, type=None):
    if type is None:
        if path.endswith('.bam'):
            type = 'bam'
        elif path.endswith('.sam'):
            type = 'sam'
        elif path.endswith('.fastq.gz'):
            type = 'fastq'
        else:
            raise ValueError(f'Cannot detect type for file {path}')

    if type == 'bam':
        return BamSequencingFile(path)
    elif type == 'sam':
        return SamSequencingFile(path)
    elif type == 'fastq':
        return FastqSequencingFile(path)
    else:
        raise ValueError(f'Unknown file type {type}')


class SequencingFile:

    def __init__(self, path):
        pass

    def seek(self, idx):
        pass

    def tell(self):
        pass

    def __iter__(self):
        pass


class SamSequencingFile(SequencingFile):

    def __init__(self, path):
        self.path = path
        self.file = pysam.Samfile(path, "r", check_sq=False)
        self._last_tell = None

    def seek(self, idx):
        self.file.seek(idx)
        self._last_tell = None

    def tell(self):
        return self.file.tell()

    def __iter__(self):
        while True:
            try:
                read = next(self.file)
            except StopIteration:
                break

            try:
                b2 = read.get_tag('B2')
                b2_q = read.get_tag('Q2')
                bc = read.get_tag('BC')
                bc_q = read.get_tag('QT')
                query = read.query_sequence
                query_qc = read.query_qualities
                yield {'query': query, 'query_qc': query_qc, 'b2': b2, 'b2_qc': b2_q, 'bc': bc, 'bc_qc': bc_q}
            except KeyError:
                pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()


class BamSequencingFile(SamSequencingFile):

    def __init__(self, path):
        self.path = path
        self.file = pysam.Samfile(path, "rb", check_sq=False)
        self._last_tell = None


class FastqSequencingFile(SequencingFile):

    def __init__(self, path):
        self.path = path
        self._last_tell = None
        if path.endswith('.gz'):
            self.file = gzip.open(path)
        else:
            self.file = open(path)

    def seek(self, idx):
        self.file.seek(idx)

    def tell(self):
        return self.file.tell()

    def __iter__(self):
        while True:
            name = self.file.readline().strip()
            if name == b'':
                break
            seq = self.file.readline().strip()
            _ = self.file.readline()
            qual = self.file.readline().strip()
            read = {'name': name.decode('utf8'), 'seq': seq.decode('utf8'), 'qual': qual.decode('utf8')}
            yield read

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()


def write_parquet(table, path):
    pyarrow.parquet.write_table(table, path)
