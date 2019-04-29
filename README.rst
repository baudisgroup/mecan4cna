Mecan4CNA
=========

Due to sample impurity and measurement bias, a copy number profile
usually needs to be calibrated for the position of baseline (normal copy
numbers) and DNA levels. Because each profile usually has different
signal scale, it’s also important to normalize profiles in comparitive
analysis.

Mecan4CNA (Minimum Error Calibration and Normalization for Copy Number
Analysis) uses an algebraic method to estimate the baseline and the
distance between DNA levels (refered as level distance). It can be used
for both single file analysis and multi-file normalization.

Key features:

-  Calibration of single file or normalization of multiple files
-  Only requires a segmentation file (from any platform)
-  Detailed results and plots for in depth analysis
-  Fast

How to install
--------------

The easiest way is to install through pip:

::

    pip install mecan4cna
    segment_liftover --help

How to use
----------

See the
`manual <https://github.com/baudisgroup/mecan4cna/blob/master/manual.md>`__
for details.

Quick start
~~~~~~~~~~~

::

    mecan4can -i [SEGMENT_FILE] -o [OUTPUT_PATH]

Demo mode
~~~~~~~~~

::

    mecan4can --demo

This will copy a few example files to the current directory and run with
default settings. It actually invokes the ``run_mecan`` script, which
will also be copied over and can be used as a template for customized
analysis.

General Usage
~~~~~~~~~~~~~

::

    Usage: mecan4cna [OPTIONS]

    Options:
      -i, --input_file FILENAME       The input file.
      -o, --output_path TEXT          The path for output files.
      -p, --plot                      Whether to show the signal histogram.
      -b, --bins_per_interval INTEGER RANGE
                                      The number of bins in each copy number
                                      interval.
      -v, --intervals INTEGER RANGE   The number of copy number intervals.
      -pt, --peak_thresh INTEGER RANGE
                                      The minimum probes of a peak.
      -st, --segment_thresh INTEGER RANGE
                                      The minimum probes of a segment.
      --model_steps INTEGER RANGE     The incremental step size in modeling.
      --mpd_coef FLOAT                Minimun Peak Distance coefficient in peak
                                      detection.
      --max_level_distance FLOAT      The maximum value of level distance.
      --min_level_distance FLOAT      The minimum value of level distance.
      --min_model_score INTEGER RANGE
                                      The minimum value of model score.
      --info_lost_ratio_thresh FLOAT  The threshold of information lost ratio.
      --info_lost_range_low FLOAT     The low end of information lost range.
      --info_lost_range_high FLOAT    The high end of information lost range.
      --help                          Show this message and exit.

Required options are:

-  ``-i FILENAME``
-  ``-o TEXT``

Input file format
~~~~~~~~~~~~~~~~~

The input should be a segmentation file:

-  have at least **5** columns as id, chromosome, start, end, probes and
   value (in exact order, names do not matter). Any additional columns
   will be ignored.
-  the first line of the file is assumed to be column names, and will be
   ignored. Do not put empty lines at the beginning of the file.
-  be **tab separated**, without quoted values

An example:

::

    id  chro    start   end num_probes  seg_mean
    GSM378022   1   775852  143752373   9992    0.025
    GSM378022   1   143782024   214220966   6381    0.1607
    GSM378022   2   88585000    144628991   4256    0.0131
    GSM378022   2   144635510   146290468   146 0.1432
    GSM378022   3   48603   8994748 1469    0.0544

Output files
~~~~~~~~~~~~

4 files will be created in the output path. If the mecan fails to detect
anything (not enough aberrant segments or no valid models), only the
histogram will be created:

-  baseNdistance.txt : contains the estimated baseline and level
   distance.
-  histogram.pdf : a visual illustration of signal distributions.
-  models.tsv : a tab seperated table that details all information of
   all models.
-  peaks.tsv : a tab seperated table shows the determined signal peaks
   and their relative DNA levels comparing to the baseline.

Common problems
---------------

Error of matplotlib in conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you run into this prolbem, simply re-install matplotlib with:

::

    conda install -n YOURENV matplotlib
