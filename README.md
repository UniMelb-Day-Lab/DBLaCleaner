# DBLaCleaner
A pipeline for demultiplexing and identifying VAR DBLalpha tags in Illumina paired end data. The pipeline relies on HMMER, flexbar, usearch and pear software programs. It outputs a number of report statistics as well as a fasta file of DBLalpha sequences along with a corresponding contaminant fasta file.

##Installation
The program is written in python to run on unix systems. It requires the programs HMMER, flexbar, usearch and pear to be installed seperately. The hardcoded location of these programs need to be entered into the corresponding variables at the top of the python script.

##Considerations
The pipeline uses Uchime denovo chimeric read filtering. This has been shown to reduce errors due to chimeric reads without removing too many valid recombinant reads. The Uparse approach was also investigated but found to be overly harsh for this problem. 

##Instructions
```
usage: cleanDBLalpha.py [-h] -o OUTPUTDIR -r READ1 -R READ2 -d DESC
                        [--NoFilter] [--NoChimeric] [--minSize MIN_SIZE]
                        [--barcodeThreshold BARCODE_THRESHOLD] [--perID PERID]
                        [--cpu CPU] [--verbose]

Process paired end Illumina reads to obtain DBLalpha OTUs.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        location of output directory. Will be created if it
                        doesn't exist
  -r READ1, --read1 READ1
                        location of first read file.
  -R READ2, --read2 READ2
                        location of paired read file.
  -d DESC, --desc DESC  location of desc mapping file.
  --NoFilter            turn off filtering for merged reads using usearch if
                        more than 1 expected error (default=True)
  --NoChimeric          turn off chimeric filtering using uchime denovo
                        algorithm. (default=True)
  --minSize MIN_SIZE    minimum support for a read to be kept. (default=15)
  --barcodeThreshold BARCODE_THRESHOLD
                        number of errors allowed in a barcode/primer pair.
                        (default=0)
  --perID PERID         percentage ID threshold. (default=0.96)
  --cpu CPU             number of cpus to use. (default=1)
  --verbose             print verbose output (default=False)
```
##Generating summary report
To generate a nice summary report of the cleaning run a Rmarkdown files is provided. At the moment it can be compiled using the command below assuming R is installed along with the R libraries `data.table`, `ggplot2`, `knitr` and `stringr`. Hopefully this will be made a bit easier in the future.

```
R -e 'rmarkdown::render("/path/to/summaryReport.Rmd", params = list(summary_file = "/path/to/summaryStatistics.log"), output_file="/output/path/summaryStatistics.html")'
```
