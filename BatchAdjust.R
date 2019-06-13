####         
#
# BatchAdjust.R   ##  ** modified from fcs19.R and fcs13.R **  
#
#  scale by XXp or SD:
#        batch to batch*XXp(1)/XXp(batch)
#        batch to batch*SD(1)/SD(batch)
#  if(pooled % zeros) > maxFrac0ForMedianThreshold (0.39 for hybrid) map to SD of batch 1 as reference (fcs12)
#  else map to XXp% of batch 1 
#  
#   modified from fcs17.R. Main difference: no mapping function, just straight scaling
#
#  fcs19_hybrid 

# source("https://bioconductor.org/biocLite.R")
# biocLite("flowCore")
library("flowCore") # for read.FCS

#g_asinh_b <- 1/5; # global asinh factor; transform is asinh(x*b); corresponds to a_b in cytofkit cytof_exprsExtract()
g_asinh_b <- 1; # global asinh factor; transform is asinh(x*b); corresponds to a_b in cytofkit cytof_exprsExtract()

# logToFile
# Maintain a log.
logToFile <- function(logfilename, logmessage, timestamp=FALSE, echo=TRUE, overwrite=FALSE){
   # Should we add a timestamp?
   if(timestamp){
      logmessage <- sprintf("%s  %s", Sys.time(), logmessage);
   }
   if(overwrite==TRUE){
      # Over write existing.
      logfile <- file(logfilename, open="wt");
   }else{
      # Open the file for appending, in text mode.
      logfile <- file(logfilename, open="at");
   }
   writeLines(logmessage, con=logfile);
   close(logfile);
   if(echo){
      print(logmessage, quote=FALSE);
   }
}

# Return the requested (80th) percentile of the vector vec.
#  perc <- .8; # percentile. median would be .5. Highest % zeros in any sample is 76.35%
get80thPercentile <- function(vec, perc=.8){
   npoints <- length(vec);
   vec_sorted <- sort(vec);
   perci <- ceiling(perc * npoints);
   perc_value <- vec_sorted[perci];
   return(perc_value);
}


# Find all anchor files in basedir,
# parse file names to return a vector of batch numbers.
#    xxx[batchKeyword][##][anchorKeyword]xxx.fcs
# example: 011118_Barcode_7_anchor stim.fcs:  anchorKeyword="anchor stim", batchKeyword="Barcode_"
# example: Set10_CTT0.fcs:  anchorKeyword="CTT0"; batchKeyword="Set";
listBatchesPresent <- function(basedir, batchKeyword="Barcode_", anchorKeyword="anchor stim"){
   batches_present <- c();
   basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
   grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);

   #ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
   # Exclude Barcode_5 and Barcode_6 (they weren't stimulated).
   ls_cmd <- sprintf("ls -1 %s/*%s*.fcs  | grep -v Barcode_5 | grep -v Barcode_6", basedir_escapeSpace, grepForAnchor_escapeSpace);
   anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));

   underscore_anchorKeyword <- sprintf("_%s", anchorKeyword);
   for(ananchor in anchors_list){
      first_part <- unlist(strsplit(ananchor, underscore_anchorKeyword, fixed=TRUE));
      next_part <- unlist(strsplit(first_part[1], batchKeyword, fixed=TRUE));
      this_batch_num <- as.numeric(next_part[2]); # numeric
      batches_present <- c(batches_present, this_batch_num);
   }
   return(sort(batches_present));
}


# Get the minimum number of events across anchor files.
getMinEventCount <- function(anchorKeyword="anchor stim", basedir="/Users/ron-home/projects/DATA/SLE_Malaria/Bead_Normalized_Debarcoded_Singlet_JG_unzipped"){

   #T0 <- Sys.time();
   whichlines <- NULL;
   basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
   grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
   #ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
   # Exclude Barcode_5 and Barcode_6 (they weren't stimulated).
   ls_cmd <- sprintf("ls -1 %s/*%s*.fcs | grep -v Barcode_5 | grep -v Barcode_6", basedir_escapeSpace, grepForAnchor_escapeSpace);
   anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
   anchor_counter <- 0;
   minCount <- Inf;
   for(ananchor in anchors_list){
      anchor_counter <- anchor_counter + 1;
      thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines);
      #thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines, column.pattern="Time");
      #this_data <- exprs(thisFCSobject);
      Nrows <- nrow(thisFCSobject);
      if(Nrows < minCount){
         minCount <- Nrows;
      }
      #print(sprintf("%i %i %s", anchor_counter, nrow(thisFCSobject), basename(ananchor)),q=F);
   }
   #T1 <- Sys.time();
   #print(T1-T0);
   return(minCount);
}


# List all columns in the most recently created fcs file.
# This is only used if no channelsToAdjust file is specified.
# Try to remove non-data channels:   "Time" "Event_length" "Center" "Offset" "Width" "Residual" 
# return cols_to_norm 
get_cols_to_norm <- function(basedir, anchorKeyword=c()){
   cols_to_skip <- c("Time","Event_length","Center","Offset","Width","Residual");
   if(is.null(anchorKeyword)){
      anchorKeyword <- "";
   }
   basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
   grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);

   ls_cmd <- sprintf("ls -1t %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
   anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));

   if(length(anchors_list) == 0){
      stop("Found no FCS files, no fields to adjust.");
   }
   cols_to_norm <- c();
   anchor_file <- anchors_list[1];
   anchor_object <- read.FCS(anchor_file, which.lines=10);
   anchor_data <- exprs(anchor_object);
   chi <- 0;
   chis_to_drop <- c();
   for(channelname in colnames(anchor_data)){
      chi <- chi + 1;
      for(skip in cols_to_skip){
         #if( length(grep(pattern=skip, x=channelname, ignore.case=TRUE)) > 0 ){}
         if( length(grep(pattern=skip, x=channelname, fixed=TRUE)) > 0 ){
            chis_to_drop <- c(chis_to_drop, chi);
            next;
         }
      }
   }
   cols_to_norm <- colnames(anchor_data)[-chis_to_drop];
   print(sprintf("No channels file. Using channels in most recent FCS file: %s", anchor_file), q=F);
   return(cols_to_norm);
}

getBatchNumFromFilename <- function(anchorFileName, batchKeyword="Barcode_", anchorKeyword="anchor stim"){
   underscore_anchorKeyword <- sprintf("_%s", anchorKeyword);
   first_part <- unlist(strsplit(anchorFileName, underscore_anchorKeyword, fixed=TRUE));
   next_part <- unlist(strsplit(first_part[1], batchKeyword, fixed=TRUE));
   this_batch_num <- as.numeric(next_part[2]); # numeric
   return(this_batch_num);
}

getNZ <- function(vec){
   wnz <- which(vec > 0);
   return(vec[wnz]);
}

# Drop the trailing p, return a number 1-100: "80p" -> 80 ; "50p" -> 50
parseP <- function(str="80p"){
   if(length(grep("p$", str)) == 1){
      num <- as.numeric(sub("p$", "", str));
      if(is.finite(num) && num>=1 && num <=100){
         return(num);
      }
   }
   print(sprintf("parseP: %s doesn't fit pattern", str));
   stop();
}


# Return mappingFunctionsList for quantile normalization.
# mappingFunctionsList[[batch]][[acol]]
getValueMappings <- function(anchorKeyword, batchKeyword, basedir, minCount, batches, cols_to_norm, transformation=TRUE, outputfile, nz_only=FALSE){

   mt0 <- Sys.time();
   whichlines <- minCount;
   basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
   grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
   #ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
   # Exclude Barcode_5 and Barcode_6 (they weren't stimulated).
   ls_cmd <- sprintf("ls -1 %s/*%s*.fcs | grep -v Barcode_5 | grep -v Barcode_6", basedir_escapeSpace, grepForAnchor_escapeSpace);
   anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
   anchor_counter <- 0;

   # Build a 2-level list: anchorDataListList[[col_to_norm]][[batch]]
   # Top level is indexed by column name (channel) to be adjusted: anchorDataList <- anchorDataListList[[col_to_norm]]
   # Second level is indexed by batch: anchorData_batchX <- unlist(anchorDataList[[batch]])
   # Initialize list.
   anchorDataListList <- list();
   pZeros <- list(); # percentage of zeros.
   # Read events for an anchor.
   # Load into lists

   Nbatches <- length(batches);
   for(ananchor in anchors_list){
      thisBatchNum <- getBatchNumFromFilename(basename(ananchor), batchKeyword=batchKeyword, anchorKeyword=anchorKeyword); # note this is numeric
      anchor_counter <- anchor_counter + 1;
      logToFile(outputfile, sprintf("Start loading events batch: %i (%i of %i)", thisBatchNum, anchor_counter, Nbatches), timestamp=TRUE);
      thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines);
      this_data <- exprs(thisFCSobject);
      for(acol in colnames(this_data)){
         if(!(acol %in% cols_to_norm)){
            next;
         }
         if(transformation){
            anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- asinh(this_data[,acol] * g_asinh_b);
         } else{
            anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- this_data[,acol];
         }
         pZeros[[acol]][[as.character(thisBatchNum)]] <- 100 * sum(anchorDataListList[[acol]][[as.character(thisBatchNum)]] <= 0) /
                                                           length(anchorDataListList[[acol]][[as.character(thisBatchNum)]]);

         if(nz_only){
            wgz <- which(anchorDataListList[[acol]][[as.character(thisBatchNum)]] > 0);
            anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- anchorDataListList[[acol]][[as.character(thisBatchNum)]][wgz];
         }
      }
   }


   # Create the reference quantile.
   nqpoints <- 100000;
   qtype <- 8; # Quantile algorithm. type=7 is default. type=8 recommended by Hyndman and Fan (1996)
   refq <- list();
   for(acol in cols_to_norm){
      # WAS: Data for all batches in this channel.
      #refq[[acol]] <- quantile(unlist(anchorDataListList[[acol]], use.names=FALSE), probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
      # Now, just use batch 1 as target.
      refq[[acol]] <- quantile(anchorDataListList[[acol]][["1"]], probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
   }

   # Create a mapping function for each batch for each column to norm.
   # This will map values in for the anchor to the reference values.
   mappingFunctionsList <- list();
   # mappingFunctionsList[[batch]][[acol]]
   for(abatch in batches){
      thisBatchChar <- as.character(abatch);
      thisBatchFunctionsList <- list(); # indexed by column name to adjust
      for(acol in cols_to_norm){
         ## Option to help at end?
         #maxRefValue <- max(refq[[acol]]);
         #if(max(anchorDataListList[[acol]][[thisBatchChar]]) < maxRefValue){
         #   qx <- quantile(c(anchorDataListList[[acol]][[thisBatchChar]], maxRefValue), probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
         #} else{
         #   qx <- quantile(anchorDataListList[[acol]][[thisBatchChar]], probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
         #}
         qx <- quantile(anchorDataListList[[acol]][[thisBatchChar]], probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
         spf <- splinefun(x=qx, y=refq[[acol]], method="monoH.FC", ties=min);
         #spf <- approxfun(x=qx, y=refq[[acol]], rule=2, ties=min, yleft=0); # clips
         thisBatchFunctionsList[[acol]] <- spf;
      }
      mappingFunctionsList[[thisBatchChar]] <- thisBatchFunctionsList;
      logToFile(outputfile, sprintf("Done mappingFunctionsList[[%s]]", thisBatchChar), timestamp=TRUE);
   }

   save(mappingFunctionsList, file=sprintf("%s/mappingFunctionsList.Rdata", dirname(outputfile)));

   mt1 <- Sys.time();
   logToFile(outputfile, "getValueMappings duration:", timestamp=TRUE);
   logToFile(outputfile, format(mt1-mt0), timestamp=FALSE);
   return(mappingFunctionsList);
}


# Return	scalingFactorsList 
# scalingFactorsList[[batch]][[acol]] 
# method = 80p | hybrid | SD | sd 
getScalingFactors <- function(anchorKeyword, batchKeyword, basedir, minCount, batches, cols_to_norm, transformation=TRUE, nz_only=FALSE, outputfile, method="80p"){
   
   mt0 <- Sys.time();
   #whichlines <- minCount;
   whichlines <- NULL;
   basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
   grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
   #ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
   # Exclude Barcode_5 and Barcode_6 (they weren't stimulated).
   ls_cmd <- sprintf("ls -1 %s/*%s*.fcs | grep -v Barcode_5 | grep -v Barcode_6", basedir_escapeSpace, grepForAnchor_escapeSpace);
   anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
   anchor_counter <- 0;

   # Build a 2-level list: anchorDataListList[[col_to_norm]][[batch]]
	# Top level is indexed by column name (channel) to be adjusted: anchorDataList <- anchorDataListList[[col_to_norm]]
	# Second level is indexed by batch: anchorData_batchX <- unlist(anchorDataList[[batch]])
	# Initialize list.
	anchorDataListList <- list();
	pZeros <- list(); # percentage of zeros.
   # Read events for an anchor. 
   # Load into lists
   Nbatches <- length(batches);
   for(ananchor in anchors_list){
      thisBatchNum <- getBatchNumFromFilename(basename(ananchor), batchKeyword=batchKeyword, anchorKeyword=anchorKeyword); # note this is numeric
      anchor_counter <- anchor_counter + 1;
      logToFile(outputfile, sprintf("Start loading events batch: %i (%i of %i)", thisBatchNum, anchor_counter, Nbatches), timestamp=TRUE);
      thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines);
      this_data <- exprs(thisFCSobject);
      if(!is.null(minCount)){
         if(minCount > nrow(this_data)){
            stop("minCount too high");
         }
         this_data <- this_data[sample.int(n=nrow(this_data), size=minCount),];
      }
      for(acol in colnames(this_data)){
         if(!(acol %in% cols_to_norm)){
            next;
         }
         if(transformation){
            #anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- sort(asinh(this_data[,acol] * g_asinh_b));
            anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- (asinh(this_data[,acol] * g_asinh_b));
         } else{
            #anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- sort(this_data[,acol]);
            anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- (this_data[,acol]);
         }
         #nZeros[[acol]][[as.character(thisBatchNum)]] <- sum(anchorDataListList[[acol]][[as.character(thisBatchNum)]] == 0);
         pZeros[[acol]][[as.character(thisBatchNum)]] <- 100 * sum(anchorDataListList[[acol]][[as.character(thisBatchNum)]] <= 0) / 
                                                           length(anchorDataListList[[acol]][[as.character(thisBatchNum)]]);
      }
   }
	
   # Create the reference.
   #allSD <- list();
   #allMean <- list();
   #allMedian <- list();
   maxFrac0ForMedianThreshold <- 0.39; # hybrid
   if(method == "SD" || method == "sd"){
      maxFrac0ForMedianThreshold <- -1; # do SD only
   }
   else { # method == "80p" || method == "50p" || method == "hybrid" ...
      # For % or median scaling.
      if(method == "hybrid"){
         perc <- .8; # percentile. median would be .5. Highest % zeros in any sample is 76.35%
      } else{ # method == "50p", or "80p"...
         maxFrac0ForMedianThreshold <- 1;   # do 80p only (or median)
         perc <- parseP(method)/100;
         print(sprintf("method %s using perc: %.2f", method, perc), q=F);
      }
   }
   FractionZerosPooled <- list();
	for(acol in cols_to_norm){
      zmax <- max(unlist(pZeros[[acol]]));
      zmin <- min(unlist(pZeros[[acol]]));
      zsd <- sd(unlist(pZeros[[acol]]));
      pooled_vec <- unlist(anchorDataListList[[acol]], use.names=FALSE);
      #allSD[[acol]] <- sd(pooled_vec);
      #allMean[[acol]] <- mean(pooled_vec);
      #allMedian[[acol]] <- median(pooled_vec);
      zeros_ref <- sum(pooled_vec <= 0);
      length_ref <- length(pooled_vec);
      pzeros_ref <- 100*zeros_ref/length_ref;
      fzeros_ref <- zeros_ref/length_ref;
      FractionZerosPooled[[acol]] <- fzeros_ref;
      print(sprintf("acol: %s  zsd: %.2f   zmin: %2.f%%  zmax: %.2f%%   pooled0s: %.2f%%   %i of %i", acol, zsd, zmin, zmax, pzeros_ref, zeros_ref, length_ref));
      if( FractionZerosPooled[[acol]] > maxFrac0ForMedianThreshold ){
         print(sprintf("acol: %s using SD scaling zsd: %.2f   zmin: %2.f%%  zmax: %.2f%%   pooled0s: %.2f%%   %i of %i", acol, zsd, zmin, zmax, pzeros_ref, zeros_ref, length_ref));
      } else{
         print(sprintf("acol: %s using percentile scaling zsd: %.2f   zmin: %2.f%%  zmax: %.2f%%   pooled0s: %.2f%%   %i of %i", acol, zsd, zmin, zmax, pzeros_ref, zeros_ref, length_ref));
      }
      
   }


   # Create a scaling factor for each batch for each column to norm.
   # This will scale values in for the anchor to the reference values.
   scalingFactorsList <- list();
   # scalingFactorsList[[batch]][[acol]] 
   for(abatch in batches){
      thisBatchChar <- as.character(abatch);
      thisBatchScalingFactors <- list(); # indexed by column name to adjust

      for(acol in cols_to_norm){
         thisvec <- anchorDataListList[[acol]][[thisBatchChar]];
         refvec <- anchorDataListList[[acol]][["1"]];
         if(nz_only){
            thisvec <- getNZ(thisvec);
            refvec <- getNZ(refvec);
         }
         # Make the target based on the fraction of zeros in this channel.
         if( FractionZerosPooled[[acol]] > maxFrac0ForMedianThreshold ){
            # use SD scaling
            scf <- sd(refvec) / sd(thisvec);
         } else{
            # use percentile 80 or median scaling
            refvalue <- get80thPercentile(refvec, perc);
            thisvalue <- get80thPercentile(thisvec, perc);
            if(refvalue == 0){
               stop(sprintf("zero scaling factor: batch %s channel %s", thisBatchChar, acol));
            }
            if(thisvalue == 0){
               stop(sprintf("undefined scaling factor: batch %s channel %s", thisBatchChar, acol));
            }
            scf <- refvalue / thisvalue;
         }

         thisBatchScalingFactors[[acol]] <- scf;
      }
      scalingFactorsList[[thisBatchChar]] <- thisBatchScalingFactors;
      logToFile(outputfile, sprintf("Done scalingFactorsList[[%s]]", thisBatchChar), timestamp=TRUE);
   }

   save(scalingFactorsList, file=sprintf("%s/scalingFactorsList.Rdata", dirname(outputfile)));
 
   mt1 <- Sys.time();
   logToFile(outputfile, "getScalingFactors duration:", timestamp=TRUE);
   logToFile(outputfile, format(mt1-mt0), timestamp=FALSE);
   return(scalingFactorsList); 
}

# method = 80p | hybrid | SD | sd | quantile | QN
BatchAdjust <- function(
   basedir=".",
   outdir=".",
   channelsFile = "ChannelsToAdjust.txt",
   batchKeyword="Barcode_",
   anchorKeyword = "anchor stim",
   nz_only=FALSE,
   method="80p",
   transformation=FALSE){

   whichlines <- NULL;
   timestamp <- format(Sys.time(), format="%Y.%m.%d.%H%M%S");
   outputfile <- sprintf("%s/LOG_BatchAdjust.%s.txt", outdir, timestamp);

   logToFile(logfilename=outputfile, logmessage="BatchAdjust.R", timestamp=TRUE, echo=TRUE, overwrite=FALSE);
   logToFile(outputfile, sprintf("basedir:%s",basedir));
   logToFile(outputfile, sprintf("outdir:%s",outdir));
   logToFile(outputfile, sprintf("channelsFile:%s",channelsFile));
   logToFile(outputfile, sprintf("batchKeyword:%s",batchKeyword));
   logToFile(outputfile, sprintf("anchorKeyword:%s",anchorKeyword));
   if(transformation){
      logToFile(outputfile, sprintf("transformation:%s","TRUE"));
   }else{
      logToFile(outputfile, sprintf("transformation:%s","FALSE"));
   }
   if(nz_only){
      logToFile(outputfile, sprintf("nz_only:%s","TRUE"));
   }else{
      logToFile(outputfile, sprintf("nz_only:%s","FALSE"));
   }
   logToFile(outputfile, sprintf("method:%s", method));


   t0 <- Sys.time(); # Time full duration.

   # Escape any spaces in file names.
   # Also enforce that batchKeyword doesn't contain spaces...
   basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
   outdir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=outdir, fixed=TRUE);
   anchorKeyword_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);

   # Parameter checking.
   if(basedir==outdir){
      print(basedir);
      print(outdir);
      stop("basedir=outdir will overwrite source files. Please specifiy basedir and/or outdir.");
   }
   # Enforce that batchKeyword doesn't contain spaces...
   if( length(grep(pattern=" ", x=batchKeyword, fixed=TRUE)) > 0 ){
      stop("batchKeyword (%s) must not contain spaces.");
   }
	if(is.null(channelsFile)){
      cols_to_norm <- get_cols_to_norm(basedir=basedir, anchorKeyword=anchorKeyword);
      logToFile(outputfile, "No channelsFile provided.");
   } else{
      logToFile(outputfile, sprintf("Reading channels file: %s", channelsFile));
      cols_to_norm <- unique(as.character(read.table(file=channelsFile, header=F, sep="\n", quote="", as.is=TRUE)[,1]));
   }
   N_ref_channels <- length(cols_to_norm);
   logToFile(outputfile, sprintf("Adjusting %i channels.", N_ref_channels));

   batchesToAdjust <- listBatchesPresent(basedir, batchKeyword=batchKeyword, anchorKeyword=anchorKeyword);
   logToFile(outputfile, "batchesToAdjust:");
   logToFile(outputfile, paste(batchesToAdjust, collapse=" "));


   batch_counter <- 0;
   file_counter <- 0;

   nbatches <- length(batchesToAdjust);

   minCountAnchors <-  NULL; # get all events. Doesn't need to be the same number for all.
   #minCountAnchors <-  310085; # for anchor stim files
   #if(is.null(minCountAnchors)){
   #   logToFile(outputfile, "minCount: Using all anchor events.");
   #} else{
   #   logToFile(outputfile, sprintf("minCount: Using %i events per anchor.", minCountAnchors));
   #}

   if(method == "quantile" || method == "QN" || method == "Quantile"){
      mappingFunctionsList <- getValueMappings(anchorKeyword=anchorKeyword, batchKeyword=batchKeyword, basedir=basedir, minCount=minCountAnchors, batches=batchesToAdjust, cols_to_norm=cols_to_norm, transformation=transformation, outputfile=outputfile, nz_only=nz_only);
      # mappingFunctionsList[[batch]][[colname]]
   } else{
      scalingFactorsList <- getScalingFactors(anchorKeyword=anchorKeyword, batchKeyword=batchKeyword, basedir=basedir, minCount=minCountAnchors, batches=batchesToAdjust, cols_to_norm=cols_to_norm, transformation=transformation, nz_only=nz_only, outputfile=outputfile, method=method);
      # scalingFactorsList[[batch]][[colname]] 
   }

   # For each batch...
   for(thisbatch in batchesToAdjust){
      td0 <- Sys.time(); # time this batch
      if(method == "quantile" || method == "QN"){
         #thisBatchMappingFunctions <- mappingFunctionsList[[as.character(thisbatch)]];
         thisBatchAdjustments <- mappingFunctionsList[[as.character(thisbatch)]];
      } else{
         #thisBatchScalingFactors <- scalingFactorsList[[as.character(thisbatch)]];
         thisBatchAdjustments <- scalingFactorsList[[as.character(thisbatch)]];
      }

      batch_counter <- batch_counter + 1;
      logToFile(outputfile, sprintf("Batch %i of %i", batch_counter, nbatches));
      logToFile(outputfile, as.character(thisbatch));

      # Apply to each file.
      #  Non-Anchor filenaming:   xxx[batchKeyword][##]_xxx.fcs
      ls_cmd <- sprintf("ls -1 %s/*%s%i_*.fcs", basedir_escapeSpace, batchKeyword, thisbatch);
      fcs_files <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
      for(fcsfile in fcs_files){
         file_counter <- file_counter + 1;
         logToFile(outputfile, sprintf("file %i", file_counter));

         # Make the adjustment per channel.
         tf0 <- Sys.time();
         thisFCSobject <- read.FCS(fcsfile, transformation=NULL, which.lines=whichlines);
         this_data <- exprs(thisFCSobject);
         these_parameters <- parameters(thisFCSobject);

         if(transformation){
            this_data <- asinh(this_data*g_asinh_b);
         }
         for(colname in colnames(this_data)){
            if(!(colname %in% cols_to_norm)){
               next;
            }
            if(method == "quantile" || method == "QN"){
               # Apply mapping function.
               mapFun <- thisBatchAdjustments[[colname]];
               this_data[,colname] <- mapFun(this_data[,colname]);
            } else{
               # Apply scaling factor.
               scalingFactor <- thisBatchAdjustments[[colname]];
               this_data[,colname] <- scalingFactor * (this_data[,colname]);
            }
            w0 <- which(this_data[,colname] < 0);
            if(length(w0) > 0){
               this_data[w0,colname] <- 0;
            }
         }

         if(transformation){ # undo the transform: inverse of asinh is sinh
            # asinh(ck/5) == ckt
            # 5*sinh(ckt) == ck
            this_data <- sinh(this_data) / g_asinh_b;
         }
         #exprs(thisFCSobject) <- this_data;
         # Need to set the ranges:
         new_params <- pData(these_parameters);
         new_params[,"minRange"] <- floor(apply(this_data, MARGIN=2, FUN=min));
         new_params[,"maxRange"] <- ceiling(apply(this_data, MARGIN=2, FUN=max));
         new_params[,"range"] <- new_params[,"maxRange"] - new_params[,"minRange"] + 1;
         pData(these_parameters) <- new_params;
         desc <-  thisFCSobject@description;
         for(paramrowi in 1:nrow(new_params)){
            rangeStr <- sprintf("$P%iR", paramrowi);
            desc[[rangeStr]] <- new_params[paramrowi, "maxRange"] + 1;
         }
         newFCSobject <- flowFrame(exprs=this_data, parameters=these_parameters, description=desc);
         #newFCSobject <- flowFrame(exprs=this_data, parameters=these_parameters)

         addExt <- "_BN";
         # Add an extension to the output file name to distinguish.
         #  addExt <- c(); # or don't
         if(is.null(addExt)){
            outfilename <- sprintf("%s/%s", outdir, basename(fcsfile));
         }else{
            replfcs <- sprintf("%s.fcs", addExt);
            #basenameW_BNext <- gsub(pattern="\\.fcs$", replacement="_BN.fcs", x=basename(fcsfile), fixed=F);
            basenameW_BNext <- gsub(pattern="\\.fcs$", replacement=replfcs, x=basename(fcsfile), fixed=F);
            outfilename <- sprintf("%s/%s", outdir, basenameW_BNext);
         }
         #write.FCS(thisFCSobject, filename=outfilename);
         write.FCS(newFCSobject, filename=outfilename);
         tf1 <- Sys.time();
         print("Per file read,norm,write:");
         print(tf1-tf0);
      }
      td1 <- Sys.time();
      logToFile(outputfile, "Per batch:"); # ~1-2 minutes 
      logToFile(outputfile, format(td1-td0));
   } # for each batch number

   t1 <- Sys.time();
   logToFile(outputfile, "Finished. Duration:", timestamp=TRUE);
   logToFile(outputfile, format(t1-t0)); 

} # BatchAdjust

   





