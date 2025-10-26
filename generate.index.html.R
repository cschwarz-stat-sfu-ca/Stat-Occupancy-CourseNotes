# This will generate an index file for all naked directories (i.e., missing an index.html file)
# for all directories under the current working directory.

# Based on 
# https://stackoverflow.com/questions/3785055/how-can-i-create-a-simple-index-html-file-which-lists-all-files-directories
# except I could not get the find option to work.


# get all the directories under the current working directory

dirs <- list.dirs( recursive=TRUE)

# process each directory in turn


library(plyr)
plyr::l_ply(dirs, function(directory){
  cat("\n\n\n ***Starting work on ", directory, "\n")
  
  # check if a hidden directory from git
  dir.chunks <- unlist(strsplit(directory, "[\\/]"))
  cat("dir.chunks ",unlist(dir.chunks), "\n")
  if(! any(dir.chunks == ".git")){
    cat("Not a hidden directory \n")
    
    # check if this directory already has an index.html file. If so, then skip
    files <- dir(directory)
    if(! any(files == "index.html")){
      cat("No html file - generate one \n")
      current_wd <- getwd()
      setwd(directory)
      system("tree -H '' -L 1 --noreport --houtro '' --prune --charset utf-8 --ignore-case --timefmt '%d-%b-%Y %H:%M' -I 'index.html' -T 'Downloads' -s -D -o index.html")
      setwd(current_wd)
      cat("Index file generated \n")
    } 
    if( any(files == "index.html")){
      cat("Index file already exists - left alone \n")
    }  
    
    
    
    
  } 
  if(any(dir.chunks == ".git")){
     cat("Hidden directory \n")
  }
  
})

