from dataclasses import dataclass
import polars as pl
import numpy as np



def calculate_pvalue(null, value, n):
    idx = null[n-1,:].searchsorted(value)
    return idx / null.shape[1]


def null_distribution(df, column='cells', N=1_000, max_n=950):
    a = df[column].to_numpy()
    a = np.log2(a)
    r = np.random.choice(a, (max_n, N))
    np.cumsum(r, axis=0, out=r)
    ns = np.arange(max_n) + 1
    r /= ns.reshape((max_n, 1))
    r.sort(axis=1)
    assert r.shape == (max_n, N)
    return r


def calculate_pvalues(df, column="reads", group_by=['grna'], N=1_000_000, max_n=950):
    controls = df.filter(pl.col('category') == 'control')
    to_test = df

    if len(to_test) == 0:
        raise ValueError('No control lineages')

    ctrl_mean = controls[column].log(2).mean()
    to_test = to_test.group_by(*group_by).agg(pl.col(column).log(2).mean().alias('teststatistic'), pl.len().alias('lineages'))
    to_test = to_test.with_columns((pl.col('teststatistic') - ctrl_mean).alias('foldchange'))

    max_size = to_test['lineages'].max()

    if max_size < max_n:
        max_n = max_size

    null = null_distribution(controls, column=column, max_n=max_n, N=N)

    out = {'pvalue': [], 'teststatistic': [], 'foldchange': [], "lineages": []}

    for tag in group_by:
        out[tag] = []


    for row in to_test.iter_rows(named=True):
        pvalue = calculate_pvalue(null, row['teststatistic'], min(row['lineages'], max_n))
        for tag in group_by:
            out[tag].append(row[tag])
        out['teststatistic'].append(row['teststatistic'])
        out['foldchange'].append(row['foldchange'])
        out['pvalue'].append(pvalue)
        out['lineages'].append(row['lineages'])

    df = pl.DataFrame(out)
    return df


