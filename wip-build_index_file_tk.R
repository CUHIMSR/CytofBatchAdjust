library(tcltk)


build_files_frame <- function(filepath, filename = NULL) {
  browser()
  files_df <- data.frame(fullpath = filepath, stringsAsFactors = FALSE)
  files_df$dirname <- dirname(files_df$fullpath)
  files_df$basename <- basename(files_df$fullpath)
  files_df$batchnum <- 1
  files_df$reference <- ""
  files_df$reference[1] <- "x"  # as an example for editing
  if (!is.null(filename)) {
    write.csv(files_df, file = filename, quote = FALSE, row.names = FALSE)
  }
}

build_index_file_tk <- function() {
  
  ##--------------------------##
  ## parameter initialization ##
  ##--------------------------##
  
  files <- data.frame(fullpath = NULL)
  cur_dir <- getwd()
  
  index_file_name <- tclVar("FCS_index.csv")
  FCS_dir <- tclVar(cur_dir)
  all_FCS <- tclVar("")

  ret_var <- tclVar("")
  
  
  ##-------------------##
  ##  button functions ##
  ##-------------------##
  
  add_FCS_dir <- function() {
    FCS_dir <- ""
    FCS_dir <- tclvalue(tkchooseDirectory(title = "Choose a directory and add all FCS to the list"))
    if (FCS_dir != "") {
      new_FCS <- dir(path = FCS_dir, pattern = "\\.fcs",
                     ignore.case = TRUE, full.names = TRUE)
      tclvalue(all_FCS) <- paste0(tclvalue(all_FCS), new_FCS, collapse = "}{")
    }
  }
  
  FCS_dir_help <- function() {
    tkmessageBox(title = "FCS_dir", message = "The directory that contains FCS files.", 
                 icon = "info", type = "ok")
  }

  index_file_name_help <- function() {
    tkmessageBox(title = "Index file name", message = "The file name to store the index of all FCS to process.", 
                 icon = "info", type = "ok")
  }
  
  
  submit <- function() {
    has_error = FALSE
    if (has_error == FALSE) {
      tclvalue(ret_var) <- "OK"
      tkdestroy(tt)
    }
  }
  
  quit <- function() {
    tkdestroy(tt)
  }
  
  
  ##----------------##
  ##  build the GUI ##
  ##--------------- ##
  
  ## head line
  
  tt <- tktoplevel(borderwidth = 20)
  tcl("wm", "attributes", tt, topmost=TRUE)
  tkwm.title(tt, "Index FCS to CSV")
  
  if(.Platform$OS.type == "windows") box_length <- 63 else box_length <- 55
  cell_width <- 3
  bt_width <- 8

  ## index_file_name
  index_file_name_label <- tklabel(tt, text = "Index file Name:")
  index_file_name_entry <- tkentry(tt, textvariable = index_file_name, width = box_length)
  index_file_name_hBut <- tkbutton(tt, text = " ? ", command = index_file_name_help)
  
  ## FCS_dir
  FCS_dir_label <- tklabel(tt, text = "FCS Directory:")
  FCS_dir_entry <- tkentry(tt, textvariable = FCS_dir, width = box_length)
  FCS_dir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, command = add_FCS_dir)
  FCS_dir_hBut <- tkbutton(tt, text = " ? ", command = FCS_dir_help)
  
  ## batch Id
  # FCS_dir_label <- tklabel(tt, text = "Raw FCS Directory:")
  # FCS_dir_entry <- tkentry(tt, textvariable = FCS_dir, width = box_length)
  
  
  ## submit / reset / quit
  submit_button <- tkbutton(tt, text = "Save", command = submit)
  quit_button <- tkbutton(tt, text = "Quit", command = quit)
  
  
  ## display GUI
  tkgrid(index_file_name_label, index_file_name_hBut, index_file_name_entry,
         padx = cell_width)
  tkgrid.configure(index_file_name_label, index_file_name_entry,
                   index_file_name_hBut, sticky = "e")
  
  tkgrid(FCS_dir_label, FCS_dir_hBut, FCS_dir_entry, FCS_dir_button, 
         padx = cell_width)
  tkgrid.configure(FCS_dir_label, FCS_dir_entry, FCS_dir_button, FCS_dir_hBut, 
                   sticky = "e")
  
  tkgrid(tklabel(tt, text = "\n"), padx = cell_width)  # leave blank line
  
  tkgrid(submit_button, tklabel(tt, text = ""), 
         quit_button, padx = cell_width)
  tkgrid.configure(quit_button, sticky = "w")
  
  tkwait.window(tt)
  
  ##-------------------##
  ## Return parameters ##
  ##-------------------##
  
  if (tclvalue(ret_var) != "OK") {
    okMessage <- "Process is cancelled."
  } else {
    FCS_files <- strsplit(tclvalue(all_FCS), "}{", fixed = TRUE)[[1]]
    build_files_frame(FCS_files, tclvalue(index_file_name))
  }
  
}
  
indexFCS_GUI()





launchShinyAPP_GUI <- function(message="cytofkit", dir = getwd(), obj = NULL){
  
  if(message == "Analysis is cancelled."){
    message("Analysis is cancelled!")
  }else{
    ifAPP <- tclVar("n")
    ss <- tktoplevel(borderwidth = 10)
    tkwm.title(ss, "cytofkit: Analysis Done")
    
    onYes <- function() {
      tclvalue(ifAPP) <- "y"
      tkdestroy(ss)
    }
    
    onNo <- function() {
      tclvalue(ifAPP) <- "n"
      tkdestroy(ss)
    }
    yesBut <- tkbutton(ss, text = " Yes ", command = onYes)
    noBut <- tkbutton(ss, text = " No ", command = onNo)
    openDirBut <- tkbutton(ss, text = "Open", command = function(){opendir(dir)})
    okBut <- tkbutton(ss, text = "OK", command = function(){tkdestroy(ss)})
    tkgrid(tklabel(ss, text = message))
    
    tkgrid(openDirBut)
    tkgrid(tklabel(ss, text = "\n"))
    tkgrid(tklabel(ss, text = "Launch Shiny APP to check your results:"))
    tkgrid(noBut, tklabel(ss, text = "    "), yesBut)
    tkgrid.configure(noBut, sticky = "e")
    tkgrid.configure(yesBut, sticky = "e")
    tkwait.window(ss)
    
    if(tclvalue(ifAPP) == "y"){
      cytofkitShinyAPP(obj)
    }
  }
}


## function for opening the results directory
opendir <- function(dir = getwd()){
  if (.Platform['OS.type'] == "windows"){
    shell.exec(dir)
  } else {
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}
