import pandas
import json
import hashlib
import os.path
import pyarrow
import pyarrow.parquet
from collections import Counter
from phantombuster.io_ import write_parquet
import uuid

def save(result, dir, id=None):
    table, stats  = result

    if id is None:
        id = str(uuid.uuid4())

    if table is not None:
        write_parquet(table, os.path.join(dir, id+'.parquet'))
        r = True
    else:
        r = False

    with open(os.path.join(dir, id + f"_stats.json"), mode="w") as f:
        json.dump(stats, f)

    return (id, r)

def load(id_r, dir):
    id, r = id_r

    if r is True:
        table = pyarrow.parquet.read_table(os.path.join(dir, id+'.parquet'))
    else:
        table = None

    with open(os.path.join(dir, id + f"_stats.json")) as f:
        stats = json.load(f)

    return (table, stats)


def _to_array(iterator, size, type):
    pyarrow.array(iterator, size=size, type=type)

def deduplicator_to_pyarrow_table(deduplicator, columns, lengths, types):
    size = len(deduplicator)
    rows = ((*items, readcount) for items, readcount in deduplicator.items())
    arrays = zip(*rows)
    t = pyarrow.Table.from_arrays([_to_array(iterator, size, type) for iterator, type in zip(arrays, types)], columns)
    return t

