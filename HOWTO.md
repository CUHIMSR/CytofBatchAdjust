# HOWTO


## Main command

Here is the main command to compute and apply the batch adjustment.

"." corresponds to the working directory of R. In RStudio use the files tab in the bottom right panel. Click on the "..." at the right of the tab to select the directory of FCS files. Then in the More menu, click "Set as working directory". Now, the files you see in the Files tab are those contained in the "." directory.

"batch_normalized" corresponds to a new directory in the current directory that will be created automatically (and should not exist yet). Normalized files will be located there.

Other parameters are explained in the README.

```
BatchAdjust(
  basedir = ".",
  outdir = "batch_normalized",
  channelsFile = "ChannelsToAdjust.txt",
  batchKeyword = "Batch",
  anchorKeyword = "Anchor",
  method = "95p",
  transformation = FALSE,
  addExt = NULL,
  plotDiagnostics = TRUE)
```

If the file ChannelsToAdjust.txt does not exist, it will be filled with the
channels of one of your FCS file. The program will stop to allow you to edit
the list (typically remove some of the channels that should not be adjusted).
Once done with editing, run again the same command, the process will start.


## Errors during plots

You can rerun the plot diagnostics if an error occurs. You have to run the
following command in R console. Copy the values you set in the call to
batchAdjust. Graphics will be overwritten.

```
call_plotAllAPrePost1ch(
  plotnz=TRUE,
  xlim=c(0,8), # increase xlim if it fails, ie change 8 by 12 or more
  postdir="copy_the_outdir",
  anchorKeyword="copy_the_anchorKeyword", 
  batchKeyword="copy_the_batchKeyword",
  predir="copy_the_basedir",
  colorPre="blue", colorPost="wheat", 
  addExt="copy_the_addExt", 
  channelsFile="copy_the_channelsFile")
```

## Copying FCS files with new names

The file namming convention is quite strict, and you probably will need to
rename your FCS files. There is two commands to achieve a copy and rename of
FCS files. The steps are the following. First, you will select all your FCS
files and build a list of file names. Then, you will edit this list in order
to set new file names matching the convention. Finally, you will launch the
copy and rename process.

```
build_index_file(
  base_dir = ".",
  index_file = "indexFile.csv"
)
```

At this stage, you should edit the indexFile.csv with an office spreadsheet software and fill the NEW column with file names that match the convention. Once finished, save the file in its original format and run the following.


```
copy_indexed_files(
  index_file = "indexFile.csv",
  destination_dir = "./renamed",
  dry_run = FALSE
)
```

The copy_files process only uses ORIGIN and NEW columns, all other columns
being ignored. Empty values in the NEW column are skipped, ie the original FCS
file is not copied to the new location.
