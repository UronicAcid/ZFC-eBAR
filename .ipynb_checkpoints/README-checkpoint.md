
# ZFC-eBAR #

ZFC-eBAR is a software tool for analyzing CRISPR-based high-throughput screening data, specifically designed to assess the enrichment of epegRNA/sgRNA from read count data. It is built on the framework of the previously developed [ZFC](https://github.com/wolfsonliu/zfc) by our lab, with the following updates:

1. eBAR (external barcode) count data may not be proportional, so zLFC (z-score of LFC) is calculated for each barcode individually.
2. The lfc-control mean modeling is done separately for each barcode, using Lowess regression for smoother fitting.
3. Although the ``--use_small_count`` option can be used to correct for low counts, we still recommend setting ``--min_ctrl_counts`` to filter out control samples (usually Day0) with average counts below a certain threshold.


## Dependency ##

ZFC-eBAR is designed with python3 and requires packages that are available
in Linux, Mac OS, and Windows.

* Python3.x
* numpy >= 1.10
* scipy >= 1.0
* pandas >= 0.16
* matplotlib >= 2.0
* sklearn >= 0.20

## Installation ##

Clone this repository, then install the software.

```{shell}
$ conda create -n zfc-ebar python==3.11.3
$ cd /path/to/ZFC-eBAR
$ python setup.py install
```

## Usage ##

The help of ZFC software:

```{shell}
usage: zfc-ebar [-h] [-i INPUT] [-o OUTPREFIX]
           [--normalization {total,median,upper_quartile,median_ratio,none}]
           [--top-n-sgrna TOP_N_SGRNA] [--top-n-gene TOP_N_GENE]
           [--null-iteration NULL_ITERATION] [--plot] [--use_small_count]
           [--min_ctrl_counts MIN_CTRL_COUNTS]
           [--min_ctrl_exp_counts MIN_CTRL_EXP_COUNTS] [--version]

Calculate fold change of screening data (zscore log2 fold change).

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Raw count table with header should include: <gene>,
                        <guide>, <barcode>, <ctrl>, <exp>.<ctrl> is the raw
                        counts of control group, <exp> is the raw counts of
                        treatment group. For screening without barcode, the
                        barcode column can be the same with guide.
  -o OUTPREFIX, --outprefix OUTPREFIX
                        Output file prefix, can be the file directory path
                        with file name prefix. The directory in the outprefix
                        should be built before analysis.
  --normalization {total,median,upper_quartile,median_ratio,none}
                        Normalization of raw count data, default is total.
                        Support method: total (Total sum normalization);
                        median (Median normalization); upper_quartile (Upper
                        quartile normalization (0.75)); median_ratio (Median
                        ratio normalization); none (No normalization).
  --top-n-sgrna TOP_N_SGRNA
                        Only consider top n barcodes for each sgRNA. Default
                        to use all the data.
  --top-n-gene TOP_N_GENE
                        Only consider top n barcodes for each gene. Default to
                        use all the data.
  --null-iteration NULL_ITERATION
                        The iteration to generate null distribution in
                        calculating the p value of genes. The larger the
                        better, but slower in calculation, default to be 100.
  --plot                Output additional QC figures.
  --use_small_count     whether to add a small count to the count matrix,
                        avoid to devide by 0, default to be False.
  --min_ctrl_counts MIN_CTRL_COUNTS
                        Min counts in ctrl, default to be 1.
  --min_ctrl_exp_counts MIN_CTRL_EXP_COUNTS
                        Min counts in ctrl + exp, default to be 1.
  --version             show program's version number and exit
```


## Example ##

Data files can be found at https://zenodo.org/.

```shell
./bin/zfc_ebar -i example/A549_D0_vs_D35.count.txt --normalization median --min_ctrl_count 100 -o A549_median_min_ctrl_counts_100 
```

If there is a bias in the fitting, you can adjust the parameters of the lowess regression in the code: ``lowess_regression`` function in ``fit_by_lowess.py``. 

The setting of ``--min_ctrl_count`` and ``--use_small_count`` parameters were different in the original HCT116 analysis. 
```shell
./bin/zfc_ebar -i example/HCT116_D0_vs_D35.count.txt --normalization median --min_ctrl_count 1 --use_small_count -o HCT116_median_use_small_count
```
Specific parameters and corresbonding quanlity control can be found in example/ZFC-eBAR_HCT116.ipynb. This Jupyter Lab can be run in the same environment and includes all the necessary functions. For HCT116, there is only a negligible difference between JupyterLab and the aforementioned command-line version, with a correlation of 1.0 between the two results.
