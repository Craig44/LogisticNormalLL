# this function in fishmod
make.filename<-function(file="",path="",add.terminal=F) {
  if(path != "") {
    plc <- substring(path, nchar(path))
    if(!(plc == "\\" | plc == "/")) path <- paste(path, "\\", sep = "")
  }
  filename <- paste(path, file, sep = "")
  if(add.terminal==T) {
    plc <- substring(filename, nchar(filename))
    if(!(plc == "\\" | plc == "/")) filename <- paste(filename, "\\", sep = "")
  }
  return(filename)
}