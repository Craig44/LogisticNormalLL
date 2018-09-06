assign.directories<-function(base="")
    {
    DIR<-list()
    
    DIR[["Base"]]<-base
    
    DIR[["R"]]<-make.filename("Rstudio",DIR[["Base"]],T)
    
    DIR[["TMB"]]<-make.filename("TMB",DIR[["Base"]],T)

    DIR[["General fun"]]<-make.filename("general_Funs",DIR[["R"]],T)

    DIR[["logistic fun"]]<-make.filename("general_Funs",DIR[["R"]],T)
    
    DIR[["csl_files"]]<-make.filename("CSL_files",DIR[["Base"]],T)

    DIR[["figures"]]<-make.filename("Figures",DIR[["Base"]],T)
  
    DIR[["C++"]]<-make.filename("C++",DIR[["Base"]],T)
    
    DIR[["C++_work"]]<-make.filename("Debug",DIR[["C++"]],T)
    
    return(DIR)
} 
