# Mecan4CNA Manual

mecan is a Python tool for calibrating and normalizing copy number segmentation files. See the main [README](https://github.com/baudisgroup/mecan4cna/blob/master/README.md) for introduction and installation.

## Copy number value convention

Usually, segmentation files use log2 ratios as the value of a segment. However, log values are not linear, mecan converts all values to copy number values. Ideally, a value of 2 should correspond to 2 copies of the DNA.

## Options
```
Usage: mecan4cna [OPTIONS]

Options:
  -i, --input_file FILENAME       The input file.
  -o, --output_path TEXT          The path for output files.
  -n, --normalize                 Calibrate and normalize the input file.
  -p, --plot                      Whether to save the signal histogram.
  -b, --bins_per_interval INTEGER RANGE
                                  The number of bins in each copy number
                                  interval.
  -v, --intervals INTEGER RANGE   The number of copy number intervals.
  --demo                          Copy example files and run a demo script in
                                  the current directory.
  -pt, --peak_thresh INTEGER RANGE
                                  The minimum probes of a peak.
  -st, --segment_thresh INTEGER RANGE
                                  The minimum probes of a segment.
  --model_steps INTEGER RANGE     The incremental step size in modeling.
  --mpd_coef FLOAT                Minimum Peak Distance coefficient in peak
                                  detection.
  --max_level_distance FLOAT      The maximum value of level distance.
  --min_level_distance FLOAT      The minimum value of level distance.
  --min_model_score INTEGER RANGE
                                  The minimum value of the model score.
  --info_lost_ratio_thresh FLOAT  The threshold of information lost ratio.
  --info_lost_range_low FLOAT     The low end of information lost range.
  --info_lost_range_high FLOAT    The high end of information lost range.
  --ld_scaler FLOAT               The scaler of level distance in
                                  normalization.
  --help                          Show this message and exit.
```

## Required options

### input file

```
-i, --input_file FILENAME
```
The input should be a segmentation file:

- have at least **5** columns: id, chromosome, start, end, probes and value (in exact order, names do not matter). Any additional columns will be ignored.
- the first line of the file is assumed to be column names, and will be ignored. Do not put empty lines at the beginning of the file.
- be **tab separated**, without quoted values


An example:

```
id	chro	start	end	num_probes	seg_mean
GSM378022	1	775852	143752373	9992	0.025
GSM378022	1	143782024	214220966	6381	0.1607
GSM378022	2	88585000	144628991	4256	0.0131
GSM378022	2	144635510	146290468	146	0.1432
GSM378022	3	48603	8994748	1469	0.0544
```

### output path

```
-o, --output_path TEXT
```
4 files will be created in the output path. If mecan fails to detect anything, only 1 file will be created:

- base_level.txt : contains the estimated baseline and level distance.
- histogram.pdf : a visual illustration of signal distributions.
- models.tsv : a tab separated table that details all information of all models.
- peaks.tsv : a tab separated table shows the determined signal peaks and their relative DNA levels compared to the baseline.

## Basic parameters

Users can change these parameters according to their needs. 

### calibration and normalization
```
-n, --normalize 
```
A flag option, when specified, the input file will also be normalized and saved as a separated file.

### plot
```
-p, --plot
```
A flag option, when specified, the histogram plot will show up in a separated window.

### intervals 
```
-v, --intervals INTEGER RANGE
```
The total number of DNA copy level intervals to use in modeling and plotting. The default value is 4, which means DNA copies from 0 to 4 are modeled and copies beyond 4 are ignored. Note, this is the ideal range. In practice, the measured copy number value is usually smaller than the actual values. This range is enough to include all signals in most cases, therefore, change of this value usually only has visual effects.

### bins per interval
```
-b, --bins_per_interval INTEGER RANGE
```
The number of bins between copy levels (for example, between 2 and 3 copies). It is used for both modeling and plotting. The default value is 20. Change of this value may have a dramatic impact on the estimation results. A high value preserves more details but also suffers more from noisy data; a low value is less accurate but is often effective in filtering out noisy data. User can modify this value according to the number of abnormal segments and the scale of changes. In general, a value between 10 and 25 performs well in most cases.

### peak threshold
```
-pt, --peak_thresh INTEGER RANGE
```
A threshold for the minimum number of probes a signal peak should have. The default value is 1000. This parameter can remove noisy signals and reduce computation time. Modify according to the total number of probes of the data generation platform. 

### segment threshold
```
-st, --segment_thresh INTEGER RANGE
```
A threshold for the minimum number of probes segment should have. The default value is 3. This parameter can remove noisy signals and reduce computation time. 

### demo 
```
--demo 
```
A flag option, when specified, it will copy 5 example files to the current directory and run with default settings. It invokes the ```run_mecan_example.sh``` script, which will also be copied over and can be used as a template for customized analysis.


## Advanced parameters

Modification of these parameters has a great impact on model performance. Please refer to our publication for detailed explanations.

### Minimum peak distance coefficient
```
--mpd_coef FLOAT
```
A parameter that is used as the minimum distance of two peaks. Default is 0.1.

### Maximum level distance
```
--max_level_distance FLOAT  
```
The maximum value of a level distance. When there are no better solutions, this limitation is ignored. Default is 1.3.

### Minimum level distance
```
--min_level_distance FLOAT
```
The minimum value of a level distance. When there are no better solutions, this limitation is ignored. default is 0.3.

### Minimum model score
```
--min_model_score INTEGER RANGE
```
The minimum value of a model score. Default is 9.

### Information lost range low
```
--info_lost_range_low FLOAT
```
The bottom limit of information lost range. Default is 0.2.

### Information lost range high
```
--info_lost_range_high FLOAT 
```
The upper limit of information lost range. Default is 0.8.

### Information lost ratio threshold
```
--info_lost_ratio_thresh FLOAT 
```
Models with information lost ration below this threshold will be abandoned. When there are no better solutions, this limitation is ignored.  

### Scaler of level distance in normalization
```
--ld_scaler FLOAT
```
In normalization, scaling down the level distance to have a loose condition that includes more segments, or scaling up to have a strict condition that includes fewer segments. The default value is 1. 