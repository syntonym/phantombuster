"""
Algorithmic steps, functions take tables, scalars, etc. as input and also output tables, scalars etc.
Does not deal with IO. Functions need to be stitched together to do something useful.
"""

import pandas as pd
import time
import logging
import pysam
from collections import Counter, defaultdict, deque
import functools
import regex
from typing import Optional, Dict
import pyarrow
import pyarrow.compute
from phantombuster.merge_cython import merge
from phantombuster.error_corrector import ErrorCorrector, Umi, UmiCount
from phantombuster.io_ import open_sequencing_file
import numpy as np
from itertools import tee
import scipy.stats
from dataclasses import dataclass
from scipy.stats._discrete_distns import binom

import pyarrow as pa

from phantombuster import bamindexer


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

_marker = object()

try:
    StringCounter = Counter[str] # type: ignore
except:
    StringCounter = Counter # type: ignore

featuretype = tuple

# copied from more-itertools https://github.com/more-itertools/more-itertools 
class peekable:
    """Wrap an iterator to allow lookahead and prepending elements.

    Call :meth:`peek` on the result to get the value that will be returned
    by :func:`next`. This won't advance the iterator:

        >>> p = peekable(['a', 'b'])
        >>> p.peek()
        'a'
        >>> next(p)
        'a'

    Pass :meth:`peek` a default value to return that instead of raising
    ``StopIteration`` when the iterator is exhausted.

        >>> p = peekable([])
        >>> p.peek('hi')
        'hi'

    peekables also offer a :meth:`prepend` method, which "inserts" items
    at the head of the iterable:

        >>> p = peekable([1, 2, 3])
        >>> p.prepend(10, 11, 12)
        >>> next(p)
        10
        >>> p.peek()
        11
        >>> list(p)
        [11, 12, 1, 2, 3]

    peekables can be indexed. Index 0 is the item that will be returned by
    :func:`next`, index 1 is the item after that, and so on:
    The values up to the given index will be cached.

        >>> p = peekable(['a', 'b', 'c', 'd'])
        >>> p[0]
        'a'
        >>> p[1]
        'b'
        >>> next(p)
        'a'

    Negative indexes are supported, but be aware that they will cache the
    remaining items in the source iterator, which may require significant
    storage.

    To check whether a peekable is exhausted, check its truth value:

        >>> p = peekable(['a', 'b'])
        >>> if p:  # peekable has items
        ...     list(p)
        ['a', 'b']
        >>> if not p:  # peekable is exhausted
        ...     list(p)
        []

    """

    def __init__(self, iterable):
        self._it = iter(iterable)
        self._cache = deque()

    def __iter__(self):
        return self

    def __bool__(self):
        try:
            self.peek()
        except StopIteration:
            return False
        return True

    def peek(self, default=_marker):
        """Return the item that will be next returned from ``next()``.

        Return ``default`` if there are no items left. If ``default`` is not
        provided, raise ``StopIteration``.

        """
        if not self._cache:
            try:
                self._cache.append(next(self._it))
            except StopIteration:
                if default is _marker:
                    raise
                return default
        return self._cache[0]

    def prepend(self, *items):
        """Stack up items to be the next ones returned from ``next()`` or
        ``self.peek()``. The items will be returned in
        first in, first out order::

            >>> p = peekable([1, 2, 3])
            >>> p.prepend(10, 11, 12)
            >>> next(p)
            10
            >>> list(p)
            [11, 12, 1, 2, 3]

        It is possible, by prepending items, to "resurrect" a peekable that
        previously raised ``StopIteration``.

            >>> p = peekable([])
            >>> next(p)
            Traceback (most recent call last):
              ...
            StopIteration
            >>> p.prepend(1)
            >>> next(p)
            1
            >>> next(p)
            Traceback (most recent call last):
              ...
            StopIteration

        """
        self._cache.extendleft(reversed(items))

    def __next__(self):
        if self._cache:
            return self._cache.popleft()

        return next(self._it)

    def _get_slice(self, index):
        # Normalize the slice's arguments
        step = 1 if (index.step is None) else index.step
        if step > 0:
            start = 0 if (index.start is None) else index.start
            stop = maxsize if (index.stop is None) else index.stop
        elif step < 0:
            start = -1 if (index.start is None) else index.start
            stop = (-maxsize - 1) if (index.stop is None) else index.stop
        else:
            raise ValueError('slice step cannot be zero')

        # If either the start or stop index is negative, we'll need to cache
        # the rest of the iterable in order to slice from the right side.
        if (start < 0) or (stop < 0):
            self._cache.extend(self._it)
        # Otherwise we'll need to find the rightmost index and cache to that
        # point.
        else:
            n = min(max(start, stop) + 1, maxsize)
            cache_len = len(self._cache)
            if n >= cache_len:
                self._cache.extend(islice(self._it, n - cache_len))

        return list(self._cache)[index]

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._get_slice(index)

        cache_len = len(self._cache)
        if index < 0:
            self._cache.extend(self._it)
        elif index >= cache_len:
            self._cache.extend(islice(self._it, index + 1 - cache_len))

        return self._cache[index]


def _get_type(c, item):
    """
    Get the type of a column depending on its name and one example item.
    Returns None for undetermined type.
    """
    name, type = c['name'], c['type']
    if type == 'reference':
        return pyarrow.large_string()
        return pyarrow.dictionary(pyarrow.int64(), pyarrow.string())
    if name == 'reads':
        return pyarrow.uint64()
    return pyarrow.large_string()


def deduplicator_to_pyarrow_table(deduplicator, barcode_hierarchy):
    """
    Transform a collections.Counter to a pyarrow table.
    The counter has a tuple as keys and the readcount as values,
    the order of the columns must be the same order as the tuple-key + "reads"
    """
    barcode_hierarchy = barcode_hierarchy+ [{'name': "reads", 'type': 'count'}]
    column_names = [c['name'] for c in barcode_hierarchy]
    rows = ((*items, readcount) for items, readcount in deduplicator.items())
    rows = peekable(rows)
    try:
        first_row = rows.peek()
    except StopIteration:
        return None

    arrays = zip(*rows)
    types = [_get_type(c, value) for c, value in zip(barcode_hierarchy, first_row)]

    assert len(first_row) == len(barcode_hierarchy)

    try:
        pyarrow_arrays = [pyarrow.array(x, type=type) for x, type in zip(arrays, types)]
    except Exception as e:
        logging.error('First row %s', first_row)
        logging.error('Column names %s', column_names)
        logging.error(f'First entry of each array: {[col[0] for col in arrays]}, types: {types}')
        raise e
    t = pyarrow.Table.from_arrays(pyarrow_arrays, column_names)
    return t

def pyarrow_table_to_counter(table, columns):
    """
    Transform a pyarrow.table to a collections.counter.
    `columns` determines the order of presence of columsn in the counter.
    """
    d = table.to_pydict()
    counter = Counter(dict(zip(zip(*(d[col] for col in columns)), d["reads"])))
    return counter


def get_sectioned_iterator(inputgroup, start, end):
    if len(inputgroup.files) == 1:
        inputfile = inputgroup.files[0]
        with open_sequencing_file(inputfile.path) as f:
            if start != 0:
                f.seek(start)
            for read in f:
                yield read
                if f.tell() >= end:
                    break
    else:
        input_files = [(iter(open_sequencing_file(inputfile.path)), inputfile.prefix) for inputfile in inputgroup.files]

        running = True

        while running:
            read = {}
            for f, prefix in input_files:
                try:
                    partial_read = next(f)
                except StopIteration:
                    running = False
                    break
                for tag in partial_read:
                    read[prefix+tag] = partial_read[tag]
            if running:
                yield read


def deduplicate_section(inputgroup, section,
        regex_dictionary,
        barcode_hierarchy,
        count_qc: bool):
    logger = logging.getLogger()

    if section is None:
        start, end = 0, np.inf
    else:
        start, end = section

    regexes = regex_dictionary.get_regexes_for_group(inputgroup.group)

    feature_deduplicator : Counter = Counter()

    # aggregate some information while iterating through all reads
    processed_reads = 0                   # count processed reads
    remaining_reads = 0                   # count surviving reads
    counters_success = {bc["name"]: StringCounter() for bc in barcode_hierarchy if bc["type"] == "reference" }
    counters_fail = StringCounter()

    logging_buffer = 5000000

    logger.debug(f"Processing group {inputgroup.group}")

    opened_input_file = get_sectioned_iterator(inputgroup, start, end)

    for read in opened_input_file:
        insufficient = False
        processed_reads += 1
        if processed_reads % logging_buffer == 0:
            logger.info(f"Processed {processed_reads} reads from group {inputgroup.group}")

        features = extract_read(regexes, read)

        for bc in barcode_hierarchy:
            tag = bc["name"]
            if tag not in features or len(features[tag]) == 0:
                counters_fail[tag] += 1
                insufficient = True

        if not insufficient:
            for bc in barcode_hierarchy:
                col = bc["name"]
                if bc["min_length"] is not None:
                    if len(features[col]) < bc["min_length"]:
                        insufficient = True
                        counters_fail[col+"_length"] += 1
                if bc ["max_length"] is not None:
                    features[col] = features[col][:bc["max_length"]]

        if not insufficient:
            for bc in barcode_hierarchy:
                if bc["type"] == "reference":
                    col = bc["name"]
                    correct = deterministic_correct(features[col], bc["reference"], None, bc["threshold"])
                    features[col] = correct

                    if correct is None:
                        counters_fail.update([col+"_correction_failed"])
                        insufficient = True

        if not insufficient:
            for bc in barcode_hierarchy:
                if bc["type"] == "reference":
                    col = bc["name"]
                    counters_success[col].update([features[col]])

            data = tuple(features.get(bc["name"], "") for bc in barcode_hierarchy)
            feature_deduplicator.update([data])
            remaining_reads += 1

    logger.info(f"Finished parsing {inputgroup.group}.")

    logger.info("Sorting each sample")

    table = deduplicator_to_pyarrow_table(feature_deduplicator, barcode_hierarchy)
    if table is not None:
        index = pyarrow.compute.sort_indices(table, sort_keys=[(bc["name"], "ascending") for bc in barcode_hierarchy])
        table = table.take(index)

    logger.info("demultiplexing thread finished")

    return table, (processed_reads, remaining_reads, counters_success, counters_fail) 


def extract_read(regex_dictionary, read):
    """Extract features from read

    Reads are made out of three regions of interest: The query, the b2tag and the bctag.
    Features can be scattered across these regions, e.g. part1 barcode is in the query region
    and part2 is in the b2tag.

    The function applies the user defined regex for each region and then aggregates
    the features. 
    """
    groups = {}

    for tag in read:
        re = regex_dictionary.get(tag, None)
        if re:
            match_ = re.search(read[tag])
        else:
            continue

        if not match_:
            continue
        groups.update(match_.groupdict())

    tag_stems = [tag.rstrip('0123456789') for tag in groups.keys()]

    # normalize ["sample", "sample1"] to ["sample0", "sample1"] so that sort handles this correctly
    for tag0 in tag_stems:
        if tag0 in groups:
            groups[tag0 + "0"] = groups[tag0]
            del groups[tag0]

    features = {}
    for tag in tag_stems:
        combined_tag = functools.reduce(
            lambda x, y: x + groups[y],
            sorted(name for name in groups.keys() if name.startswith(tag)),
            "",
        )
        features[tag] = combined_tag

    return features


def _combine_duplicators_cython(a,b, barcode_hierarchy):
    if a is None:
        return b
    if b is None:
        return a
    else:
        a_reads_as_uint64 = a["reads"].cast(pyarrow.uint64())
        b_reads_as_uint64 = b["reads"].cast(pyarrow.uint64())

        a = a.remove_column(a.column_names.index("reads"))
        b = b.remove_column(b.column_names.index("reads"))
        a = a.append_column("reads", a_reads_as_uint64)
        b = b.append_column("reads", b_reads_as_uint64)
        c = merge(a,b, barcode_hierarchy)
        if c == -1:
            raise Exception("Unsupported type in merge")

        sort_keys = [(bc["name"], "ascending") for bc in barcode_hierarchy]
        index = pyarrow.compute.sort_indices(c, sort_keys=sort_keys)
        c = c.take(index)

        return c


def combine(deduplications, barcode_hierarchy):
    main_table, processed_reads, remaining_reads, counters_success, counters_fail = None, 0, 0, defaultdict(Counter), Counter()

    for i, deduplication in enumerate(deduplications):
        partial_table, (_processed_reads, _remaining_reads,  _counters_success, _counters_fail) = deduplication

        main_table = _combine_duplicators_cython(main_table, partial_table, barcode_hierarchy)

        processed_reads += _processed_reads
        remaining_reads += _remaining_reads
        counters_success = {name: counters_success[name] + Counter(_counters_success[name]) for name in set(counters_success.keys()) | set(_counters_success.keys())}
        counters_fail += _counters_fail

        del deduplication
        del partial_table

    return main_table, (processed_reads, remaining_reads, counters_success, counters_fail)

def error_correct(table, barcode_hierarchy, threshold):
    partitions = [table]
    tag_names = [tag["name"] for tag in barcode_hierarchy[1:]]
    for tag_idx, tag in enumerate(barcode_hierarchy):
        new_partitions = []
        for partition in partitions:
            if tag["type"] == "random":
                partition = error_correct_column(partition, tag['name'], threshold)

            sort_necessary = tag['type'] == 'random'
            ps = sort_and_partition(partition, sort_columns=[(tag['name'], 'ascending')], partition_columns=[tag['name']], unwrap_if_single=False, sort=sort_necessary)
            for partition in ps.values():
                new_partitions.append(partition)
        partitions = new_partitions
    if len(partitions) > 0:
        table = pyarrow.concat_tables(partitions)
    return table

def error_correct_column(partition, tag, threshold):
    column = partition[tag].to_numpy()
    reads = partition["reads"].to_numpy()

    tag_length = len(column[0])

    tag_and_reads = partition.drop([column for column in partition.column_names if not column in [tag, 'reads']])
    tag_and_reads = deduplicate_single(tag_and_reads)

    correction = calculate_corrections(tag_and_reads[tag].to_numpy(), tag_and_reads['reads'].to_numpy(), tag_length,  threshold)
    partition = apply_correction(partition, tag, correction)
    return partition

def deduplicate_single(table):
    idx = pyarrow.compute.sort_indices(table, sort_keys=[(bc, "ascending") for bc in table.column_names])
    table = table.take(idx)

    cache = {col: table[col].to_numpy() for col in table.column_names}
    cache['reads'] = cache['reads'].copy()

    mask = np.ones(len(table), dtype=np.bool_)

    def compare(table, cache, i, j):
        for col in table.column_names:
            if col == "reads":
                continue
            if cache[col][i] != cache[col][j]:
                return False
        return True

    j = 0
    for i in range(1, len(table)):
        if compare(table, cache, i, j):
            cache['reads'][j] += cache['reads'][i]
            mask[i] = 0
        else:
            j = i

    table = table.set_column(table.column_names.index('reads'), 'reads', pyarrow.array(cache['reads']))
    table = table.filter(mask)
    return table

def calculate_corrections(tags, reads, tag_length, threshold):

    ec = ErrorCorrector(tag_length, 0, 0, threshold)
    tags = [(Umi(0, 0, tag), UmiCount(r, 1)) for tag, r in zip(tags, reads)]
    tags_dict = dict(tags)
    result = ec.process(tags_dict)

    corrections = result[1]
    corrected   = result[0]

    correction = {a.tag: b.tag for a, b in corrections.items()}

    return correction

def apply_correction(table, correction_column, correction):
    col = table[correction_column]
    mask = []
    corrected_tags = []
    for tag in col.to_numpy():
        corrected_tag = correction.get(tag, None)
        if corrected_tag:
            corrected_tags.append(corrected_tag)
            mask.append(True)
        else:
            mask.append(False)
    table = table.filter(mask)
    corrected_column = pyarrow.array(corrected_tags)
    table = table.set_column(table.column_names.index(correction_column), correction_column, corrected_column)
    table = deduplicate_single(table)
    return table


def calculate_overlap(table, barcode, lid_barcodes=('lid',), map=None):

    if map:
        mapping = lambda x: map.get(x, x)
    else:
        mapping = lambda x: x

    idx_nreads = pyarrow.compute.sort_indices(table, sort_keys=[*[(bc, "ascending") for bc in lid_barcodes], ("reads", "descending")])
    sorted_table = table.take(idx_nreads)

    sorted_lids = [sorted_table[bc] for bc in lid_barcodes]
    sorted_reads = sorted_table["reads"]
    sorted_sample = sorted_table[barcode]

    values = []

    current_lid = None
    current_n = None
    current_sample = None


    for lid, n, sample in zip(zip(*sorted_lids), sorted_reads, sorted_sample):
        lid = tuple([bc.as_py() for bc in lid])
        n = n.as_py()
        if lid == current_lid:
            assert current_n >= n
            if mapping(sample) != mapping(current_sample):
                d = n / current_n
            else:
                d = 1.1 # LIDs in the same sample group should never be removed.
            values.append(d)
        else:
            values.append(1.0)
            current_lid = lid
            current_n = n
            current_sample = sample

    values = np.array(values)

    return sorted_table, values


def calculate_hopping_threshold(table, hopping_barcodes):

    anchor_barcodes = [col for col in table.column_names if col not in hopping_barcodes]
    if 'reads' in anchor_barcodes:
        anchor_barcodes.remove('reads')

    idx_nreads = pyarrow.compute.sort_indices(table, sort_keys=[*[(bc, "ascending") for bc in anchor_barcodes], ("reads", "descending")])
    reverse_idx = pyarrow.compute.sort_indices(idx_nreads)
    sorted_table = table.take(idx_nreads)

    anchor_columns = [sorted_table[bc] for bc in anchor_barcodes]
    reads = sorted_table['reads']

    values = []

    last_anchor = None
    total_reads = 0
    seen_readcounts = []
    categories = 0
    sig_categories = 0

    for anchor_values, read_count in zip(zip(*anchor_columns), reads):
        anchor = tuple([bc.as_py() for bc in anchor_values])
        read_count = read_count.as_py()

        if last_anchor != anchor:
            if len(seen_readcounts) == 1:
                values.append(seen_readcounts[0])
            else:
                if sig_categories > 1:
                    p = 1/sig_categories if sig_categories > 0 else 1
                    threshold = binom.isf(0.05, total_reads, p)
                    for r in seen_readcounts:
                        values.append(threshold)
                        #result = scipy.stats.binomtest(read_count, total_reads, p, alternative="greater").pvalue
                        #values.append(result)
                else:
                    if len(seen_readcounts) > 0:
                        threshold = max(seen_readcounts)
                        for r in seen_readcounts:
                            values.append(threshold)
            last_anchor = anchor
            total_reads = read_count
            sig_categories = 1 if read_count > 1 else 0
            categories = 1
            seen_readcounts = [read_count]
        else:
            total_reads += read_count
            categories += 1
            if read_count > 1:
                sig_categories += 1
            seen_readcounts.append(read_count)

    # Calculate the last group
    if len(seen_readcounts) == 1:
        values.append(seen_readcounts[0])
    else:
        p = 1/sig_categories if sig_categories > 0 else 1
        threshold = binom.isf(0.05, total_reads, p)
        for r in seen_readcounts:
            values.append(threshold)
    values = np.array(values)
    values = pyarrow.compute.take(values, reverse_idx)
    return values


def call_hopping(table, threshold_list):
    reads = table["reads"].to_numpy()
    condition = reads >= threshold_list[0]
    for thresholds in threshold_list[1:]:
        condition = condition & (reads>= thresholds)
    condition = pyarrow.array(condition)
    out_table = table.filter(condition)
    return out_table


def hopping_removal(table, overlap, threshold):
    condition = pyarrow.array(overlap >= threshold)
    out_table = table.filter(condition)
    return out_table

def sort_and_partition(table, sort_columns=(("samplename", "ascending"), ("reads", "descending")), partition_columns=("samplename",), unwrap_if_single=True, return_bounds=False, sort=True):
    partition_keys, partitions = _sort_and_partition(table, sort_columns, partition_columns, sort=sort)

    if not return_bounds:
        partitions = [p for p, start, end in partitions]
    d = dict(zip([tuple(key) for key in partition_keys], partitions))

    # in case of only one partition key, make keys the values instead of tuple with length 1
    if unwrap_if_single:
        if len(list(d.keys())[0]) == 1:
            d = {key[0]: value for key, value in d.items()}

    return d

def _sort_and_partition(table, sort_columns=(("samplename", "ascending"), ("reads", "descending")), partition_columns=("samplename",), sort=True):
    if isinstance(partition_columns, str):
        partition_columns = [partition_columns]

    search_column_keys = [c[0] for c in sort_columns]
    p_col_idx = [search_column_keys.index(col) for col in partition_columns]
    table_columns = table.column_names

    # all partition columns must appear as sorting columns. More sorting columns are allowed.
    assert all(i != -1 for i in p_col_idx)
    # partition columns must be in same order as sort columns
    assert sorted(p_col_idx) == p_col_idx
    assert all(col in table_columns for col in search_column_keys)

    if sort:
        index = pyarrow.compute.sort_indices(table, sort_keys=sort_columns)
        table = table.take(index)

    partitions     = [(table, 0, len(table))]
    partition_keys = [[]]

    for col in partition_columns:
        new_partitions = []
        new_partition_keys = []

        for part, key in zip(partitions, partition_keys):
            p, orig_start, orig_length = part

            new_keys  = pyarrow.compute.unique(p[col])

            start = 0
            indexes = []
            for new_key in new_keys:
                idx = p[col].index(new_key, start=start).as_py()
                indexes.append(idx)
                start = idx
            indexes.append(len(p))

            d = {}
            for new_key, (start, end) in zip(new_keys, pairwise(indexes)):
                new_partition_key = [*key, new_key.as_py()]
                new_partition = (p.slice(start, end-start), orig_start+start, end-start)

                new_partitions.append(new_partition)
                new_partition_keys.append(new_partition_key)

        partitions = new_partitions
        partition_keys = new_partition_keys
    return partition_keys, partitions



def threshold_table(table, threshold):
    table_th = table.join(threshold, "sample")
    condition = pyarrow.compute.greater_equal(table_th["reads"], table_th["threshold"])
    return table.filter(condition)

def _metric(a, b):
    return sum([aa != bb for aa, bb in zip(a,b)])

def deterministic_correct(barcode, references, default_value, threshold):
    # TODO do this better
    r = references.get(barcode, None)
    if r is not None:
        return r
    else:
        for reference_bc, reference_name in references.items():
            if _metric(barcode, reference_bc) <= threshold:
                return reference_name
        return default_value

def calculate_threshold(value_dict):
    th = np.inf
    for v1 in value_dict.keys():
        for v2 in value_dict.keys():
            if v1 == v2:
                continue
            th = min(th, int((_metric(v1, v2)-1)/2))
    return th


def calculate_jaccard(table, barcode):

    idx_nreads = pyarrow.compute.sort_indices(table, sort_keys=[("lid", "ascending"), (barcode, 'ascending')])
    sorted_table = table.take(idx_nreads)

    sorted_lids = sorted_table["lid"]
    sorted_sample = sorted_table[barcode]

    counts = {}
    intersection_counts = {}

    current_lid = None
    current_samples = set()

    for lid, sample in zip(sorted_lids, sorted_sample):
        lid = lid.as_py()
        if lid == current_lid:
            current_samples.add(sample.as_py())
        else:
            for sample1 in current_samples:
                for sample2 in current_samples:
                    if sample1 == sample2:
                        counts[sample1] = counts.get(sample1, 0) +1
                    else:
                        key = (sample1, sample2)
                        intersection_counts[key] = intersection_counts.get(key, 0) + 1
            current_lid = lid
            current_samples = {sample.as_py()}

    # for last lid
    for sample1 in current_samples:
        for sample2 in current_samples:
            if sample1 == sample2:
                counts[sample1] = counts.get(sample1, 0) +1
            else:
                key = (sample1, sample2)
                intersection_counts[key] = intersection_counts.get(key, 0) + 1

    jaccard_index = {}
    for key in intersection_counts:
        sample1, sample2 = key
        intersection_count = intersection_counts[key]
        count1, count2 = counts[sample1], counts[sample2]
        ji = intersection_count / (count1 + count2 - intersection_count)
        jaccard_index[key] = ji

    return jaccard_index
