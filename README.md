# BatchAdjust() - CyTOF Batch Adjust


Harmonize all samples in a cytometry experiment using anchor samples included in each batch to compute adjustment factors for each channel in each batch.



## Installation

BatchAdjust()   runs in the R environment.

Resources for installing and getting started with R are available at the Comprehensive R Archive Network:

https://cran.r-project.org/manuals.html




### Install required R packages

Install flowCore if you haven't already.


#### flowCore
At the R command line enter:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```



## Usage

BatchAdjust() is a command line application for R. It is designed to run in a linux environment or macOS.

In an R session, navigate to where you downloaded BatchAdjust.R, and load the application by typing:

```

source("BatchAdjust.R")

```



# Arguments

```
BatchAdjust(
   basedir=".",
   outdir=".",
   channelsFile = "ChannelsToAdjust.txt",
   batchKeyword="Barcode_",
   anchorKeyword = "anchor stim",
   method="95p",
   transformation=FALSE,
   addExt=NULL,
   plotDiagnostics=TRUE)
```



###### basedir: 
directory to look for input (source) FCS files. All files to be adjusted must be in this directory.

###### outdir:  
directory to write resulting batch adjusted files. Must not be the same as basedir to avoid overwriting original data files, unless addExt is set (see below). 


###### channelsFile: 
plain text file listing channels to adjust, one per line. 
Only channels listed here will be adjusted, and only channels that should be adjusted should be listed here. E.g. open channels and barcoding channels should be omitted from this file. 
Channel names must match those in the FCS files exactly.

For an example, see ChannelsToAdjust\_example.txt.

###### batchKeyword:
"Barcode\_" (refer to File naming requirements)

###### anchorKeyword:
"anchor stim" (refer to File naming requirements)


###### method:
95p | SD | quantile

quantile: quantile normalization

SD: scaling to reference batch standard deviation

50p: scaling to reference batch 50th percentile (median)

95p: scaling to reference batch 95th percentile

Batches may be scaled to an user-defined percentile by specifying any number (1-100) followed by the letter 'p'. For example method="80p" would scale channels to the 80th percentile of the reference batch.


###### transformation:
 TRUE | FALSE

TRUE: asinh transformation is applied before batch adjustment. sinh is applied to the adjusted data before writing results.

FALSE: No transformation is applied.


###### addExt:
A character string to append to output filenames (before .fcs extension) to distinguish from input filenames.

With the default=NULL, output filenames are identical to input filenames.
If addExt is not NULL, basedir and outdir may be the same directory.


###### plotDiagnostics:
 TRUE | FALSE

Generate distribution plots for each adjusted channel for all batch anchor samples before and after adjustment.

Also generate a figure summarizing permutation test results for decreased variability in signal levels among batches for all channels.

Note: this may take longer than the adjustment process itself. It is safe to interrupt this process, and doing so will not affect the batch adjusted .fcs files.




# File naming requirements:

Filenames for FCS files to be batch adjusted must contain a reference to which batch they belong.

One sample from each batch must also contain the anchor keyword defined by the parameter "anchorKeyword" to indicate the control sample expected to be consistent across batches.

Adjustments for all batches are relative to batch 1, and samples in batch 1 are not changed. Therefore one batch must be must be labeled as batch 1. To change your reference batch, simply rename your files.




##### Non-Anchor file naming:   
xxx[batchKeyword][##]\_xxx.fcs

'xxx' is optional and may be any characters.

Note that underscore '\_' is required after batch number (to distinguish e.g. Batch10\_ from Batch1\_).



##### Anchor file naming:   
xxx[batchKeyword][##]\_[anchorKeyword]xxx.fcs

'xxx' is optional and may be any characters.

Note that no other characters are allowed between [batchKeyword][##]\_[anchorKeyword].

Note that underscore '\_' is required after batch number (to distinguish e.g. Batch10\_ from Batch1\_).

Separators such as '\_' are allowed, but must be specified in the anchor keywords.





###### File name examples 1:

011118\_Barcode\_7\_anchor stim.fcs (anchor sample)

011118\_Barcode\_7\_310A\_T0.fcs

011118\_Barcode\_7\_310A\_T6.fcs

011118\_Barcode\_7\_310A\_T6\_R.fcs



For the above file name examples, batchKeyword and anchorKeyword parameters should be set as follows:

batchKeyword = "Barcode\_"

anchorKeyword = "anchor stim"





###### File name examples 2:

Set10\_CTstim.fcs (anchor sample)

Set10\_CCP13T0.fcs

Set10\_CCP13T6LPS.fcs

Set10\_CCP13T6PBS.fcs

Set10\_CCP13T6R848.fcs



For the above file name examples, batchKeyword and anchorKeyword parameters should be set as follows:

batchKeyword = "Set"

anchorKeyword = "CTstim"










