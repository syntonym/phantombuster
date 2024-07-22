#!/usr/bin/env python

import logging
import logging.config
import os
import os.path
import sys
import multiprocessing as mp
import glob
import collections

import ntpath
from dataclasses import dataclass
import json

from phantombuster import porcelain, plumbing, stores
from phantombuster import core
from phantombuster.stores import deduplicator_to_pyarrow_table
from phantombuster.remoter import Worker
from phantombuster.io_ import write_parquet
from phantombuster.project import Project
import click
import click.core
from typing import Optional, List, Mapping
import pyarrow.parquet
import pyarrow.csv

from pathlib import Path

from gettext import gettext


class OrderedGroup(click.Group):
    def __init__(self, name: Optional[str] = None, commands: Optional[Mapping[str, click.Command]] = None, **kwargs):
        super(OrderedGroup, self).__init__(name, commands, **kwargs)
        self.commands = commands or collections.OrderedDict()
        self.main_commands = collections.OrderedDict()
        self.secondary_commands = collections.OrderedDict()

    def list_commands(self, ctx: click.Context) -> Mapping[str, click.Command]:
        return self.commands

    def list_main_commands(self, ctx: click.Context) -> Mapping[str, click.Command]:
        return self.main_commands

    def list_secondary_commands(self, ctx: click.Context) -> Mapping[str, click.Command]:
        return self.secondary_commands

    def add_main_command(self, cmd, name=None):
        """Registers another :class:`Command` with this group.  If the name
        is not provided, the name of the command is used.
        """
        self.add_command(cmd, name)
        name = name or cmd.name
        if name is None:
            raise TypeError("Command has no name.")
        click.core._check_multicommand(self, name, cmd, register=True)
        self.main_commands[name] = cmd

    def add_secondary_command(self, cmd, name=None):
        """Registers another :class:`Command` with this group.  If the name
        is not provided, the name of the command is used.
        """
        self.add_command(cmd, name)
        name = name or cmd.name
        if name is None:
            raise TypeError("Command has no name.")
        click.core._check_multicommand(self, name, cmd, register=True)
        self.secondary_commands[name] = cmd

    def main_command(self, *args, **kwargs):
        """A shortcut decorator for declaring and attaching a command to
        the group. This takes the same arguments as :func:`command` and
        immediately registers the created command with this group by
        calling :meth:`add_command`.

        To customize the command class used, set the
        :attr:`command_class` attribute.

        .. versionchanged:: 8.0
            Added the :attr:`command_class` attribute.
        """
        from click.decorators import command

        if self.command_class is not None and "cls" not in kwargs:
            kwargs["cls"] = self.command_class

        def decorator(f):
            cmd = command(*args, **kwargs)(f)
            self.add_main_command(cmd)
            return cmd

        return decorator

    def secondary_command(self, *args, **kwargs):
        """A shortcut decorator for declaring and attaching a command to
        the group. This takes the same arguments as :func:`command` and
        immediately registers the created command with this group by
        calling :meth:`add_command`.

        To customize the command class used, set the
        :attr:`command_class` attribute.

        .. versionchanged:: 8.0
            Added the :attr:`command_class` attribute.
        """
        from click.decorators import command

        if self.command_class is not None and "cls" not in kwargs:
            kwargs["cls"] = self.command_class

        def decorator(f):
            cmd = command(*args, **kwargs)(f)
            self.add_secondary_command(cmd)
            return cmd

        return decorator

    def format_commands(self, ctx, formatter):
        """Extra format methods for multi methods that adds all the commands
        after the options.
        """
        main_commands = []
        for subcommand in self.list_main_commands(ctx):
            cmd = self.get_command(ctx, subcommand)
            # What is this, the tool lied about a command.  Ignore it
            if cmd is None:
                continue
            if cmd.hidden:
                continue

            main_commands.append((subcommand, cmd))

        secondary_commands = []
        for subcommand in self.list_secondary_commands(ctx):
            cmd = self.get_command(ctx, subcommand)
            # What is this, the tool lied about a command.  Ignore it
            if cmd is None:
                continue
            if cmd.hidden:
                continue

            secondary_commands.append((subcommand, cmd))



        # allow for 3 times the default spacing
        if len(main_commands) or len(secondary_commands):
            limit = formatter.width - 6 - max(len(cmd[0]) for cmd in main_commands)

            rows = []
            for subcommand, cmd in main_commands:
                help = cmd.get_short_help_str(limit)
                rows.append((subcommand, help))

            if rows:
                with formatter.section(gettext("Main Commands")):
                    formatter.write_dl(rows)

            limit = formatter.width - 6 - max(len(cmd[0]) for cmd in secondary_commands)

            rows = []
            for subcommand, cmd in secondary_commands:
                help = cmd.get_short_help_str(limit)
                rows.append((subcommand, help))

            if rows:
                with formatter.section(gettext("Secondary Commands")):
                    formatter.write_dl(rows)

def configure_logging(outputlog, verbose):
    logging_config = {
                      'version': 1,
                      'formatters':{'default': {'format': "%(asctime)s %(levelname)-8s %(name)-15s %(message)s",
                                                'datefmt': "%Y-%m-%d %H:%M:%S"}},
                      'handlers': {'console': {'class': 'logging.StreamHandler', 'formatter': 'default', 'stream': 'ext://sys.stdout'}},
                      'loggers': {'remoter': {'level': 'WARNING'}},
                      'root': {'handlers': ['console'], 'level': 'INFO'}
                      }
    if outputlog:
        logging_config['handlers']['file'] = {'class': 'logging.FileHandler', 'formatter': 'default', 'filename': outputlog}
        logging_config['root']['handlers'].append('file')
    if verbose:
        logging_config['root']['level'] = 'DEBUG'

    logging.config.dictConfig(logging_config)
    logging.info('Logging configured')


def log_call(function, **kwargs):
    """Log a CLI call with all arguments"""
    logging.info(f"PhantomBuster {function} was called with the following arguments: {kwargs}")


@click.group(cls=OrderedGroup)
@click.version_option(package_name='phantombuster')
@click.option("--verbose/--silent", default=False, help="Enable verbose debugging")
@click.option("-o", "--outputlog", type=click.Path(), help="Output file for logs")
@click.option("--save-results/--no-save-results", type=bool, default=True, help="DEBUG OPTION set no-save-results to not save which stages were already done")
def phantombuster(verbose: bool, outputlog: str, save_results: bool) -> None:
    """Bioinformatical tool to remove phantoms from barcode based NGS sequencing data.

    The tool consists of four stages, represented by the four main commands.
    """
    configure_logging(outputlog, verbose)


# -- Main Commands -- #


@phantombuster.main_command()
@click.argument("input", type=click.Path(exists=True))
@click.option("--outdir", required=True)
@click.option("--regex-file", required=True)
@click.option("--barcode-hierarchy-file", type=click.Path(exists=True), required=True)
@click.option("--debug/--production", default=False)
@click.option("--show-qc/--no-qc", default=False)
@click.option("--force/--no-force", default=False)
def demultiplex(input, regex_file, debug, outdir, show_qc, force, barcode_hierarchy_file):
    """
    Demultiplex BAM/FASTA files into parquet files according to the barcode hierarchy
    """
    log_call("demultiplex", input=input, regex_file=regex_file, debug=debug,
             outdir=outdir, show_qc=show_qc, barcode_hierarchy_file=barcode_hierarchy_file)
    project = Project(outdir)

    try:
        core.demultiplex(input, regex_file, barcode_hierarchy_file, project, debug=debug, show_qc=show_qc)
    except Exception as e:
        logging.exception("Pipeline encountered an error. Aborting.")
        raise click.Abort()
    return


@phantombuster.main_command()
@click.option("--outdir", required=True)
@click.option("--error-threshold", default=1)
@click.option("--barcode-hierarchy-file", required=True)
def error_correct(outdir, error_threshold, barcode_hierarchy_file):
    """
    Correct sequencing errors originating from single-nucleotide errors
    """
    log_call("error-correct", outdir=outdir, error_threshold=error_threshold)
    project = Project(outdir)
    core.error_correct(project, error_threshold, barcode_hierarchy_file)


@phantombuster.main_command()
@click.argument('hopping-barcodes', nargs=-1)
@click.option("--outdir", required=True)
@click.option("--threshold", default=0.05, type=float)
def hopping_removal(outdir, threshold, hopping_barcodes):
    """
    Remove phantom combinations originating from index hopping
    """
    log_call("hopping-removal", outdir=outdir, threshold=threshold, hopping_barcodes=hopping_barcodes)
    project = Project(outdir)
    hopping_barcodes = [bc.split(',') for bc in hopping_barcodes]
    core.hopping_removal(project, hopping_barcodes, threshold)

@phantombuster.main_command()
@click.option("--outdir", required=True)
@click.option("--threshold-file", required=True)
def threshold(outdir, threshold_file):
    """
    Remove combinations with a read count below a threshold
    """
    log_call("threshold", outdir=outdir, threshold_file=threshold_file)
    project = Project(outdir)
    core.threshold(project, threshold_file)

# -- Helper Commands -- #

@phantombuster.secondary_command()
@click.option("--outdir", default=None, required=True)
@click.option("--name", default=None)
def worker(outdir, name):
    """
    Start a worker process
    """
    import phantombuster as phantombuster
    import phantombuster.plumbing

    project = Project(outdir)
    project.create()
    path = project._get_server_path()

    print(f"Connecting to {path}")
    worker = Worker(path, name=name)
    worker.start_async()

@phantombuster.secondary_command()
@click.argument("parquetfile")
@click.argument("outfile", default=None, required=False)
def to_csv(parquetfile, outfile):
    """
    Convert a parquet file to a CSV file.
    """
    log_call("to_csv", sample=parquetfile, outdir=outfile)
    table = pyarrow.parquet.read_table(parquetfile)
    if outfile is None:
        outfile = parquetfile.replace(".parquet", ".csv")
    pyarrow.csv.write_csv(table, outfile)


@phantombuster.secondary_command()
@click.argument("csvfile")
@click.argument("outfile", default=None, required=False)
def to_parquet(csvfile, outfile):
    """
    Convert a CSV file to a parquet file.
    """
    log_call("to_parquet", csvfile=csvfile, outfile=outfile)
    table = pyarrow.csv.read_csv(csvfile)
    if outfile is None:
        outfile = csvfile.replace(".csv", ".parquet")
    write_parquet(table, outfile)


@phantombuster.secondary_command()
@click.argument("prefixes", nargs=-1)
@click.option("--outdir", required=True)
@click.option("--prefix")
@click.option("--barcode-hierarchy-file", type=click.Path(exists=True), required=True)
def merge(prefixes, outdir, prefix, barcode_hierarchy_file):
    """
    Merge multiple prefixes under one prefix
    """
    # Log the call of this function with all parameters to the logfile
    log_call("merge", prefixes=prefixes, outdir=outdir, prefix=prefix, barcode_hierarchy_file=barcode_hierarchy_file)

    master_paths = PathsAndFiles(outdir, prefix, None)
    master_paths.create()

    try:
        barcode_hierarchy = plumbing.read_barcode_hierarchy_file(barcode_hierarchy_file)
    except Exception:
        raise Exception("Could not read barcode hierarchy file correctly")

    to_merge = [PathsAndFiles(outdir, prefix, None) for prefix in prefixes]

    results = [stores.load(('deduplication', True), paths.stage_path('deduplication')) for paths in to_merge]
    out = plumbing.combine(results, barcode_hierarchy)

    stores.save(out, master_paths.stage_path('deduplication'), id='deduplication')

if __name__ == "__main__":
    phantombuster()
