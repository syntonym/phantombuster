Files
=====

The configuration of PhantomBuster for each individual experiment is controlled via three CSV files.
They are passed via command line arguments to the individual stages.

input_files.csv
---------------

CSV file that lists each input file (BAM or FASTQ) to be processed.
Each row corresponds to a single input file.
The column `file` is mandatory and must be the path to the input file.
Paths can be absolute or relative to the directory in which phantombuster is executed.
The column `group` is optional and groups different files into logical groups, in case that the sequences of a read are located in different files, as is common for FASTQ files.
In FASTQ format the read 1, read 2, index 1 and index 2 are commonly seperated into files of the schema `SampleName_S1_L001_X_001.fastq.gz` where `X` is `R1` for read 1, `R2` for read 2 (paired end sequencing), `I1` for index 1 and `I2` for index 2.
Each group must be unique and files that do not share reads/sequences must be in seperate groups.
The column `prefix` is optional and configures which regular expressions are used for the file (see also `regexes.csv`).
For the implementation that parses the file see `phantombuster.config_files.read_input_files_file`.

.. list-table:: input_files.csv columns
   :header-rows: 1

   * - Column
     - Mandatory
     - Meaning
   * - file
     - Yes
     - Path to input file
   * - group
     - No
     - Grouping of files which must be processed together
   * - prefix
     - No
     - Regexes to use for extraction

.. list-table:: input_files.csv Minimal Example
   :header-rows: 1

   * - file
   * - ``/scratch/experiment20240702/data/20240702.bam``

.. _input_complex_example:

.. list-table:: input_files.csv Complex Example
   :header-rows: 1

   * - file
     - group
     - prefix
   * - ``SampleOne_S1_L001_R1_001.fastq.gz``
     - S1
     - FASTQ_R1
   * - ``SampleOne_S1_L001_R2_001.fastq.gz``
     - S1
     - FASTQ_R2
   * - ``SampleOne_S1_L001_I1_001.fastq.gz``
     - S1
     - FASTQ_I1
   * - ``SampleOne_S1_L001_I2_001.fastq.gz``
     - S1
     - FASTQ_I2
   * - ``SampleTwo_S2_L001_R1_001.fastq.gz``
     - S2
     - FASTQ_R1
   * - ``SampleTwo_S2_L001_R2_001.fastq.gz``
     - S2
     - FASTQ_R2
   * - ``SampleTwo_S2_L001_I1_001.fastq.gz``
     - S2
     - FASTQ_I1
   * - ``SampleTwo_S2_L001_I2_001.fastq.gz``
     - S2
     - FASTQ_I2
   * - ``20240702.bam``
     - S3
     - BAM

regexes.csv
-----------

Regular expression to extract the barcodes from the sequences.
Each row specifies one regular expression to apply to a read region.
The column `regex` is mandatory and specifies the regular expression to use when extracting barcodes.
For accepted regular expression syntax see the `regex module <https://github.com/mrabarnett/mrab-regex>`_.
Named groups are extracted as barcodes (e.g. `(?P<sample>\w+)` is a regular expression that captures a barcode named `sample`).
Named groups with numbers at the end are concatenated to make a single barcode (e.g. `(?P<sample0>\w+)ACGTACGT(?P<sample1>\w+)` results in a single sample barcode).
The column `tag` is mandatory and specifies the read region the regex is applied to.
Possible values are `query`, `bc` and `b2` for BAM files, and `name` and `seq` for FASTQ files.
The column `prefix` is optional and is used to specify different regular expressions for different input files.
Each `prefix` value in the `input_files.csv` should have a corresponding entry in the `regexes.csv` file.
The column `group` is optional and can be used to specify regexes for a specific group of input files.
The use of the `group` column is deprecated, instead use the `prefix` column.
For the implementation that parses the file see `phantombuster.config_files.read_regex_file`.

.. list-table:: regexes.csv Columns
   :header-rows: 1

   * - Column
     - Mandatory
     - Meaning
   * - regex
     - Yes
     - Regular expression to extract barcodes
   * - tag
     - Yes
     - Target Read region (query, bc, b2, name, seq)
   * - prefix
     - No
     - Indicates for which input files the regex applies to
   * - group
     - No
     - Deprecated. Indicates the input file group the regex applies to

.. list-table:: regexes.csv Minimal Example
   :header-rows: 1

   * - tag
     - regex
   * - b2
     - ``"^[ACGTN]{3}(?P<sample>[ACGTN]{5})"``
   * - query
     - ``"(?P<lid>[ACGTN]{5,6}(?P<lib>ACGT|GTAC){s<=1}[ACGTN]+)"``


.. _regexes_complex_example:

.. list-table:: regexes.csv Complex Example
   :header-rows: 1

   * - prefix
     - tag
     - regex
   * - FASTQ_R1
     - seq
     - ``"^(?P<cell>[ACGTN]{12})\w*(AGGACGAAACACC){s<=1}(?P<grna>\w{20})"``
   * - FASTQ_R2
     - seq
     - ``"(?P<lid>\w+)"``
   * - FASTQ_I1
     - seq
     - ``"(?P<sample0>\w{8})"``
   * - FASTQ_I2
     - seq
     - ``"(?P<sample1>\w{8})"``
   * - BAM
     - query
     - ``"^(?P<cell>[ACGTN]{12})\w*(AGGACGAAACACC){s<=1}(?P<grna>\w{20})"``
   * - BAM
     - bc
     - ``"(?P<sample0>\w{8})"``
   * - BAM
     - b2
     - ``"(?P<sample1>\w{8})"``

The complex examples continues the example in :ref:`input_complex_example`. 
Each named group/barcode has a corresponding entry in :ref:`barcodes_complex_example`.

.. _barcodes.csv:

barcodes.csv
------------

CSV files that lists all barcodes occuring in the experiment and their type.
All columns are mandatory.
The column `barcode` is the name of the barcode, each barcode specified in `regexes.csv` must correspond to a barcode here.
The column `type` configures whether the valid sequences of the barcode are known (`reference`) or not (`random`).
The column `referencefile` is a path to a CSV file with valid sequences for `reference` barcodes and is ignored for `random` barcodes.
For the content of the reference file, see :ref:`barcode_sequence_files`.
The column `threshold` is used for the error correction of `reference` barcodes.
A barcode sequence is corrected to a reference sequence if their hamming distance is below the error threshold.
A value of `auto` will determine the largest possible error threshold that allows for unique error correction and is recommended.
For `random` barcodes the column is ignored.
The column `min_length` determines the minimal length of the barcode sequences, sequences of that barcode below the value are discarded.
The column `max_length` determines the maximal length of the barcode sequences, sequences of that barcode above the value are discarded.

The order of the barcodes is significant and used for the error correction of random barcodes.
When correcting random barcode sequences, only sequences are compared for which all barcode sequences above the barcode are equal.
Practically that means that barcodes should be specified from general to specific.
With the four barcodes `sample`, `grna`, `lineage` and `cell` two sequences of the `lineage` barcode would only be compared if they have the same `sample` and `grna` values.
Sequences with different `sample` or `grna` values can not originate from the same lineage and are thus not compared.
For two `cell` barcode squences to be compared their `lineage` sequence would also need to be the same.

The `min_length` and `max_length` columns overlap in their purpose with the length restrictions directly in the regular expression.
As the regular expression allows to configure minimal, maximal lengths and more directly, the `min_length` and `max_length` are deprecated and should be set to `-`.
Instead formulate any resctrictions on the barcodes directly in the regular expression.

For the implementation that parses the file see `phantombuster.config_files.read_barcode_hierarchy_file`.

.. list-table:: barcodes.csv Columns
   :header-rows: 1

   * - Column
     - Mandatory
     - Meaning
   * - barcode
     - Yes
     - Name of the barcode
   * - type
     - Yes
     - Barcode type, either `reference` or `random`
   * - referencefile
     - Yes
     - Path to reference file for `reference` barcodes, ignored otherwise
   * - threshold
     - Yes
     - Threshold for error correction (``auto`` or int)
   * - min_length
     - Yes
     - DEPRECATED Minimal length of barcode, set to `-` to disable
   * - max_length
     - Yes
     - DEPRECATED Maximal length of barcode, set to `-` to disable

.. list-table:: barcodes.csv Minimal Example
   :header-rows: 1

   * - barcode
     - type
     - referencefile
     - threshold
     - min_length
     - max_length
   * - sample
     - reference
     - ``/scratch/experiment20240702/sample_barcodes.csv``
     - auto
     - \-
     - \-
   * - lib
     - reference
     - ``/scratch/experiment20240702/library_barcodes.csv``
     - auto
     - \-
     - \-
   * - lid
     - random
     - \-
     - \-
     - 50
     - 50

.. _barcodes_complex_example:

.. list-table:: barcodes.csv Complex Example
   :header-rows: 1

   * - barcode
     - type
     - referencefile
     - threshold
     - min_length
     - max_length
   * - sample
     - reference
     - ``sample_barcodes.csv``
     - auto
     - \-
     - \-
   * - grna
     - random
     - \-
     - \-
     - \-
     - \-
   * - lid
     - random
     - \-
     - \-
     - 50
     - 50
   * - cell
     - random
     - \-
     - \-
     - \-
     - \-

The complex examples continues the example in :ref:`regexes_complex_example`.
The regexes extract values for `sample0` and `sample1`, which are concatenated to a single `sample` barcode.
Thus, here only a single `sample` barcode is listed.

thresholds.csv
--------------

Read threshold for the thresholding step.
Barcode combinations with a read count below the read threshold are discarded.
Only the column `threshold` is mandatory, which specifies the read threshold.
All other columns must be the name of a `reference` barcode as specified in :ref:`barcodes.csv`.
Valid values of the column are then the names in the corresponding barcode sequence file (see :ref:`barcode_sequence_files`).


.. list-table:: thresholds.csv Columns
   :header-rows: 1

   * - Column
     - Mandatory
     - Meaning
   * - threshold
     - Yes
     - Read threshold to apply
   * - `[BARCODE NAME]`
     - No
     - specifies to which barcode combinations the threshold applies
   * - ...
     - No
     - Multiple barcode names can be supplied


.. list-table:: thresholds.csv Minimal Example
   :header-rows: 1

   * - threshold
   * - 100


.. list-table:: thresholds.csv Complex Example
   :header-rows: 1

   * - sample
     - threshold
   * - Sample1
     - 80
   * - Sample2
     - 120

.. _barcode_sequence_files:

barcode sequences files
-----------------------

For each `reference` barcode in :ref:`barcodes.csv` a file with all valid barcode sequences must be provided.
These files must contain two columns.
The `name` column assigns each sequence a human readable name, e.g. in the case of sample barcodes the sample name or id.
The `barcode` column must consist of the sequence.

.. list-table:: barcode sequence files Columns
   :header-rows: 1

   * - Column
     - Mandatory
     - Meaning
   * - name
     - Yes
     - human readable name
   * - barcode
     - Yes
     - sequence (consisting of ACGT)

.. list-table:: example
   :header-rows: 1
   
   * - name
     - barcode
   * - Sample1
     - ``CGTACTAGATAGAGAG``
   * - Sample2
     - ``TCCTGAGCTCTACTCT``
