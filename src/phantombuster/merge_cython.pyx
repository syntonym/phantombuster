# cython: language_level=3
# distutils: language = c++
# cython: language=c++

import pyarrow
cimport pyarrow
from pyarrow cimport CTable
from pyarrow.lib cimport pyarrow_wrap_array, pyarrow_unwrap_table, TableBatchReader, CRecordBatch, CStringArray, CDoubleArray, CInt64Builder, CArray, c_get_memory_pool, CDoubleBuilder, CUInt64Array, CUInt64Builder, CFixedSizeBinaryArray, pyarrow_wrap_schema, CSchema, Type, CFixedSizeBinaryType, CDataType, CLargeStringArray, CLargeBinaryArray, CBinaryArray

import numpy

cimport cython
from cython.operator cimport dereference as deref

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.cast cimport dynamic_cast, static_cast
from libc.stdint cimport uint8_t, int32_t, int64_t


ctypedef CArray* CArrayPtr
ctypedef CStringArray* CStringArrayPtr
ctypedef CLargeStringArray* CLargeStringArrayPtr
ctypedef CDoubleArray* CDoubdynamicyPtr
ctypedef CUInt64Array* CUInt64ArrayPtr
ctypedef CFixedSizeBinaryArray* CFixedSizeBinaryArrayPtr
ctypedef CFixedSizeBinaryType* CFixedSizeBinaryTypePtr
ctypedef CBinaryArray* CBinaryArrayPtr
ctypedef CLargeBinaryArray* CLargeBinaryArrayPtr

#@cython.boundscheck(False)
#@cython.wraparound(False)
#cdef int _compare(const uint8_t* a, const uint8_t* b):
    #return ((deref(a) > deref(b)) - (deref(a) < deref(b)));
    #if l1 < l2:
        #return -1
    #elif l1 > l2:
        #return 1
    #else:
        #return 0
    #return l1.compare(l2)
    #

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _compare_string(string a, string b):
    cdef int r = a.compare(b)
    if r > 0:
        return 1
    elif r < 0:
        return -1
    else:
        return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _compare_uint(const uint8_t* a, const uint8_t* b, int lena, int lenb):
    cdef unsigned int i = 0
    while i < lena:
        if a[i] < b[i]:
            return -1
        elif a[i] > b[i]:
            return 1
        i += 1
    return 0

cdef int _compare_type_aware(shared_ptr[CArray] a, shared_ptr[CArray] b, int i, int j):

    cdef Type t1 = deref(a).type_id()
    cdef Type t2 = deref(b).type_id()

    cdef int32_t length_a
    cdef int32_t length_b

    cdef int64_t length64_a
    cdef int64_t length64_b

    cdef unsigned int lena
    cdef unsigned int lenb
    cdef unsigned int offseta
    cdef unsigned int offsetb
    cdef string str_a
    cdef string str_b
    cdef const uint8_t* bin_a
    cdef const uint8_t* bin_b

    if t1 == Type._Type_STRING:
        str_a = dynamic_cast[CStringArrayPtr](a.get()).GetString(i)
        str_b = dynamic_cast[CStringArrayPtr](b.get()).GetString(j)
        return _compare_string(str_a, str_b)
    if t1 == Type._Type_LARGE_STRING:
        str_a = dynamic_cast[CLargeStringArrayPtr](a.get()).GetString(i)
        str_b = dynamic_cast[CLargeStringArrayPtr](b.get()).GetString(j)
        return _compare_string(str_a, str_b)
    if t1 == Type._Type_FIXED_SIZE_BINARY:
        bin_a = dynamic_cast[CFixedSizeBinaryArrayPtr](a.get()).GetValue(i)
        bin_b = dynamic_cast[CFixedSizeBinaryArrayPtr](b.get()).GetValue(j)
        lena = (dynamic_cast[CFixedSizeBinaryTypePtr](deref(a).type().get())).byte_width()
        lenb = (dynamic_cast[CFixedSizeBinaryTypePtr](deref(b).type().get())).byte_width()
        return _compare_uint(bin_a, bin_b, lena, lenb)
    if t1 == Type._Type_BINARY:
        bin_a = dynamic_cast[CBinaryArrayPtr](a.get()).GetValue(i, &length_a)
        bin_b = dynamic_cast[CBinaryArrayPtr](b.get()).GetValue(j, &length_b)
        return _compare_uint(bin_a, bin_b, length_a, length_b)
    if t1 == Type._Type_LARGE_BINARY:
        bin_a = dynamic_cast[CLargeBinaryArrayPtr](a.get()).GetValue(i, &length64_a)
        bin_b = dynamic_cast[CLargeBinaryArrayPtr](b.get()).GetValue(j, &length64_b)
        return _compare_uint(bin_a, bin_b, length64_a, length64_b)
    return -2

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int compare(shared_ptr[CRecordBatch] a, shared_ptr[CRecordBatch] b, unsigned int i, unsigned int j, cols, unsigned int max_i, unsigned int max_j):

    if i == max_i:
        return 1
    if j == max_j:
        return -1

    cdef int r = 0;

    cdef int ci

    for c in cols:

        ci = c

        r = _compare_type_aware(deref(a).column(ci), deref(b).column(ci), i, j)
        #r = _compare(str_a, str_b)
        if r != 0:
            return r
    return r

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef xmerge(table1, table2, unsigned long[:] idx, unsigned long[:] counts, cols):

    cdef unsigned int i = 0
    cdef unsigned int offset_i = 0
    cdef unsigned int j = 0
    cdef unsigned int offset_j = 0

    cdef unsigned int idx_i = 0
    cdef int c

    cdef int reads_column = table1.column_names.index("reads")

    cdef shared_ptr[CTable] t1 = pyarrow_unwrap_table(table1)
    cdef shared_ptr[CTable] t2 = pyarrow_unwrap_table(table2)

    cdef TableBatchReader* tbr1 = new TableBatchReader(deref(t1))
    cdef TableBatchReader* tbr2 = new TableBatchReader(deref(t2))

    cdef shared_ptr[CRecordBatch] batch1
    cdef shared_ptr[CRecordBatch] batch2

    tbr1.ReadNext(&batch1);
    tbr2.ReadNext(&batch2);

    cdef unsigned int max_i = deref(batch1).num_rows()
    cdef unsigned int max_j = deref(batch2).num_rows()
    cdef unsigned int j_offset = len(table1)

    cdef shared_ptr[CSchema] schema = deref(batch1).schema()

    cdef CUInt64ArrayPtr r1 = static_cast[CUInt64ArrayPtr](deref(batch1).column(reads_column).get())
    cdef CUInt64ArrayPtr r2 = static_cast[CUInt64ArrayPtr](deref(batch2).column(reads_column).get())

    while i < max_i or j < max_j:

        c = compare(batch1, batch2, i, j, cols, max_i, max_j)

        if c == 1:
            idx[idx_i] = j + offset_j + j_offset
            counts[idx_i] = r2.Value(j)
            j += 1
            idx_i += 1
        elif c == -1:
            idx[idx_i] = i + offset_i
            counts[idx_i] = r1.Value(i)
            i += 1
            idx_i += 1
        elif c == 0:
            idx[idx_i] = i + offset_i
            counts[idx_i] = r1.Value(i) + r2.Value(j)
            i += 1
            j += 1
            idx_i += 1
        else:
            return idx, counts, idx_i, i, j, -1

        if i != 0 and i == max_i:
            offset_i += i
            i = 0
            tbr1.ReadNext(&batch1);
            if batch1 == NULL:
                max_i = 0
            else:
                max_i = deref(batch1).num_rows()
                r1 = static_cast[CUInt64ArrayPtr](deref(batch1).column(reads_column).get())
        if j != 0 and j == max_j:
            offset_j += j
            j = 0
            tbr2.ReadNext(&batch2);
            if batch2 == NULL:
                max_j = 0
            else:
                max_j = deref(batch2).num_rows()
                r2 = static_cast[CUInt64ArrayPtr](deref(batch2).column(reads_column).get())

    return idx, counts, idx_i, i, j, 0

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef merge(table1, table2, bhc):

    length = len(table1) + len(table2)

    idx = numpy.zeros(length, dtype=numpy.uint)
    counts = numpy.zeros(length, dtype=numpy.uint)

    column_names = table1.column_names

    t1 = table1
    t2 = table2

    column_names.remove("reads")

    _, __, idx_i, i, j, r = xmerge(t1, t2, idx, counts, [table1.column_names.index(bc['name']) for bc in bhc])
    if r != 0:
        return r

    t3 = pyarrow.concat_tables([table1, table2])

    t3 = t3.remove_column(t3.column_names.index("reads"))
    idx = idx[:idx_i]
    counts = counts[:idx_i]

    a = pyarrow.array(idx)

    counts_array = pyarrow.array(counts)

    t3 = t3.take(a)

    t3 = t3.append_column("reads", counts_array)

    return t3
