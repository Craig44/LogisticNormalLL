## initialisation script

# Set paths and working directories

if (Sys.info()[["user"]]=="Dell User") {
  thisPath <- "C:\\Work\\Projects\\PhD\\Logistic_normal\\"
} 

library(nzPlot)
source(paste(thisPath,"\\Rstudio\\general_Funs\\assign.directories.R",sep=""))
source(paste(thisPath,"\\Rstudio\\general_Funs\\make.filename.R",sep=""))
DIR<-assign.directories(base=thisPath)


