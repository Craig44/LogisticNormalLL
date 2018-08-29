#
#
#  prune.pop
#  prune.est
#  prune.out
#
#  prune.pars
#
#  remove.SA.from.est
#
#  make.moab.script



prune.pop = function(pop, current.year){
  
#-------------------------------------------------------------
#
# This functions takes a CASAL population file and prunes it 
# back to a different fishing year from what it is set up for
#  
# Arguments
#
# est:    CASAL population file read in by extract.csl.file
#
# current.year:  year to prune the estimation file back to 
#
# Value
#
# Output is the pruned population file
#
# Further Details
#
# Pruned back are:
# 
#  (1) The YCSs for recruitment
#  (2) Catch data
#  (3) Median catch day for shifting the selectivity for the western spawning fishery
#
#--------------------------------------------------------------  

pop$current$value = current.year
pop$final$value = pop$current$value + 5 

pop$"recruitment[E]"$last_free =  current.year - 5
pop$"recruitment[E]"$YCS_years =  1975:(current.year - 2) 
pop$"recruitment[E]"$YCS =  rep(1, (current.year - 1975 -1)) 
pop$"recruitment[E]"$year_range[2] = current.year - 2

pop$"recruitment[W]"$last_free = current.year - 5
pop$"recruitment[W]"$YCS_years =  1975:(current.year - 2)
pop$"recruitment[W]"$YCS = rep(1,(current.year - 1975 -1))
pop$"recruitment[W]"$year_range[2] = current.year - 2 

pop$first_random_year$value =  current.year - 1 
#
# keep all catch data up to the current year
#

fisheries = c("Ensp1", "Ensp2", "Wnsp1", "Wnsp2", "Esp", "Wsp" )
for(fish.index in fisheries){
  
  this.fish = eval(parse(text= paste("pop$", "'fishery[", fish.index, "]'", sep="") ))
  current.index = match(current.year, as.numeric(this.fish$years))
  if(is.na(current.index))
    stop(paste("For fishery", fish.index, "current year does not exist:", current.year))
  this.fish$years =   this.fish$years[1:current.index]
  this.fish$catches = this.fish$catches[1: current.index]
  this.fish$future_years = seq(from=current.year + 1, to=current.year + 5 )
  # use the catch in the current year for future catches
  this.fish$future_catches = rep(rev(this.fish$catches)[1], 5)
  
  eval(parse(text= paste("pop$", "'fishery[", fish.index, "]' = this.fish", sep="")))
          
}

#
#  Median catch day by year for Wsp, as used in estimating annual changes in 
#  the selectivity Wspsl. The mean value was used in all years for which there
#  was catch but no Wspage data (i.e. before 1988 and in the current fishing year)
#
# The read in version had for pop$"selectivity[Wspsl]" some oddities for
# shift_years  and   shift_E (see below) which indicates the R casal
# read in function doesn't like tabs. To fix this up tabs were changed
# to spaces in Editplus for the population.csl file, and it was read 
# in again, and the problem had gone away. 
# 
# 
# $command
# [1] "selectivity"
# 
# $value
# [1] "Wspsl"
# 
# $all
# [1] "size_based" "logistic"   "50"         "10"        
# 
# $shift_years
# [1] "1972"                                                                                                                                                                                    
# [2] "1973\t1974\t1975\t1976\t1977\t1978\t1979\t1980\t1981\t1982\t1983\t1984\t1985\t1986\t1987\t1988\t1989\t1990\t1991\t1992\t1993\t1994\t1995\t1996\t1997\t1998\t1999\t2000\t2001\t2002\t2003\t2004\t2005\t2006\t2007\t2008\t2009"
# [3] "2010"                                                                                                                                                                                    
# [4] "2011"                                                                                                                                                                                    
# 
# $shift_E
# [1] "306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t306\t299\t302\t298\t301\t306\t304\t308\t307\t312\t310\t311\t309\t309\t309\t308\t309\t307\t309\t310\t307\t301\t295\t298"
# [2] "306"                                                                                                                                                        
# 
# $shift_a
# [1] "0"
# 

this.sel = pop$'selectivity[Wspsl]'

# mean value over 1988  to  current.year - 1
start.index = match("1988", this.sel$shift_years)
current.index = match(current.year, as.numeric(this.sel$shift_years))

this.sel$shift_years = this.sel$shift_years[1:current.index]

# use mean value outside the range 1998 to current.year - 1  (inclusive)
this.sel$shift_E = this.sel$shift_E[1:current.index]
mean.value = round(mean( as.numeric(this.sel$shift_E[start.index:current.index]) ))
this.sel$shift_E[1: (start.index - 1)] =  mean.value
this.sel$shift_E[current.index] =  mean.value 

pop$'selectivity[Wspsl]' = this.sel

return(pop)

}



prune.est = function(est, current.year){

#-------------------------------------------------------------
#
# This functions takes a CASAL estimation file and prunes it 
# back to a different fishing year from what it is set up for
#  
# Arguments
#
# est:    CASAL estimation file read in by extract.csl.file
#
# current.year:  year to prune the estimation file back to 
#
# Value
#
# Output is the pruned estimation file
#
# Further Details
#
# Pruned back are:
# 
#  (1) The YCSs in the estimation block
#  (2) Relative abundance observations (i.e. trawl surveys)
#  (3) Proportions and catch at age observations  
#
#--------------------------------------------------------------
 
# 1000 for the MCMCs, but set to 500 for concatenation of chains  
est$MCMC$burn_in = 1000  

# find out where the estimate blocks are
est.indices = which(substring(names(est),1,8)=="estimate")

E.YCS.index.found = FALSE
for(index in est.indices){
  this.name = est[[index]]$parameter
  if(this.name=="recruitment[E].YCS"){
    E.YCS.index = index
    E.YCS.index.found = TRUE
    break
  }
}
if(!E.YCS.index.found) stop("Estimate block for E YCSs not found")


W.YCS.index.found = FALSE
for(index in est.indices){
  this.name = est[[index]]$parameter
  if(this.name=="recruitment[W].YCS"){
    W.YCS.index = index
    W.YCS.index.found = TRUE
    break
  }
}
if(!W.YCS.index.found) stop("Estimate block for W YCSs not found")

# E YCS estimated from 1975 up to (current.year - 2)
this.est = est[[E.YCS.index]]
no.YCSs = (current.year - 2) - 1975 + 1
this.est$mu = rep(1, no.YCSs)
this.est$cv = rep(0.95, no.YCSs)
this.est$lower_bound = rep(0.06, no.YCSs)
this.est$upper_bound = rep(8.60, no.YCSs)

est[[E.YCS.index]] = this.est

# W = E, except for parameter name
est[[W.YCS.index]] = this.est
est[[W.YCS.index]]$parameter = "recruitment[W].YCS"

#
# relative abundance observations
#

all.names = names(est)
biomass.obs.names = all.names[substring(all.names,1,18)=="relative_abundance"]

for(obs.name in biomass.obs.names){
  
  this.obs = eval(parse( text=paste("est$'",obs.name,"'", sep="") ))  
  the.years = as.numeric(this.obs$years)
  invalid.years = this.obs$years[the.years>current.year]
  if(length(invalid.years>0)){
      for(invalid.year in invalid.years){
        delete.obs = paste("this.obs$'",    invalid.year, "' = NULL", sep="")
        delete.cv  = paste("this.obs$'cv_", invalid.year, "' = NULL", sep="")
        eval(parse(text=paste(text=delete.obs)))
        eval(parse(text=paste(text=delete.cv)))      
      }
  }
  this.obs$years = this.obs$years[the.years<=current.year]
  
  obs.set = paste("est$'", obs.name, "'= this.obs", sep="")
  eval(parse(text=obs.set))

}

#
# catch at observations (commerical fishery)
#
# and 
#
# proportions at observations  (trawl survey)
#

all.names = names(est)
catch.at.obs.names = all.names[substring(all.names,1,8)=="catch_at"]
prop.at.obs.names = all.names[substring(all.names,1,14)=="proportions_at"]
both.at.obs.names = c(catch.at.obs.names, prop.at.obs.names)

for(obs.name in both.at.obs.names){
  
  this.obs = eval(parse( text=paste("est$'",obs.name,"'", sep="") ))  
  the.years = as.numeric(this.obs$years)
  invalid.years = this.obs$years[the.years>current.year]
  if(length(invalid.years>0)){
      for(invalid.year in invalid.years){
        delete.obs = paste("this.obs$'",    invalid.year, "' = NULL", sep="")
        delete.cv  = paste("this.obs$'cvs_", invalid.year, "' = NULL", sep="")
        eval(parse(text=paste(text=delete.obs)))
        eval(parse(text=paste(text=delete.cv)))      
      }
  }
  this.obs$years = this.obs$year[the.years<=current.year]
  
  obs.set = paste("est$'", obs.name, "'= this.obs", sep="")
  eval(parse(text=obs.set))

}


return(est)

}    # end of function




prune.out = function(out, current.year){
  
#-------------------------------------------------------------
#
# This functions takes a CASAL output file and prunes it 
# back to a different fishing year from what it is set up for
#  
# Arguments
#
# out:    CASAL population file read in by extract.csl.file
#
# current.year:  year to prune the estimation file back to 
#
# Value
#
# Output is the pruned output file
#
# Further Details
#
# Pruned back are:
# 
#  (1) The output for the eastern spawning stock biomass
#  (2) The output for the fit to the sub-Antarctic trawl survey
#
#--------------------------------------------------------------  

this.out = out

if("abundance[ESSB]" %in% names(this.out)){
  the.years = this.out$"abundance[ESSB]"$years  
  last.index = which(the.years==as.character(current.year))
  this.out$"abundance[ESSB]"$years = the.years[1:last.index]
}

if("abundance[SAbio]" %in% names(this.out)){
  the.years = this.out$"abundance[SAbio]"$years  
  last.index = which(the.years==as.character(current.year))
  this.out$"abundance[SAbio]"$years = the.years[1:last.index]
}
 
out = this.out 
return(out)

}

  
prune.pars = function(par.file, current.year){  
  
#-----------------------------------------------------------------
#
# This function takes a CASAL input parameter file (i.e -i par.file)
# and prunes it back to a different year from which it was set up for
#
# The pruning is just of the YCS value (there are less estimated if there
# are less years in the model). 
#  
# Arguments
#
# par.file:    CASAL par file for the 2011 fishing year model
#
# current.year:  year to prune the pars file back to
#  
#
# Value
#
# Output is the pruned output file
#
# Further Details
#
# Checking in done on par.file that is contains the right number of
# YCSs to start
#
#-------------------------------------------------------------------------  
  
#
# YCSs estimated from 1975 up to (current.year - 2)
# So number of YCSs estimated = current.year - 1975 - 1  (35 YCSs in 2011)
#    - same number estimated for E & W fisheries


# E stock

t1 = extract.free.parameters.from.table(par.file)
num.free.YCSs = current.year - 1975 - 1

col.names = names(t1)
first.E.index = match("recruitment[E].YCS.1",col.names)
end.E.YCS.name = paste("recruitment[E].YCS.",num.free.YCSs,sep="")
end.E.index =match(end.E.YCS.name,col.names)
last.E.index = match("recruitment[E].YCS.35", col.names)

# only do something if there are YCSs to drop
if(end.E.index<last.E.index){
  drop.E.index = (end.E.index + 1):last.E.index
  t1 = t1[ , -drop.E.index]
}
  
col.names = names(t1)
first.W.index = match("recruitment[W].YCS.1",col.names)
end.W.YCS.name = paste("recruitment[W].YCS.",num.free.YCSs,sep="")
end.W.index =match(end.W.YCS.name, col.names)
last.W.index = match("recruitment[W].YCS.35",col.names)

# only do something if there are YCSs to drop
if(end.W.index<last.W.index){
  drop.W.index = (end.W.index + 1):last.W.index
  t1 = t1[ , -drop.W.index]
}


# make a parameter file
ld = readLines(par.file, n=2, ok=FALSE)

new.head = gsub("35", num.free.YCSs, ld[1], fixed=T)
new.vars = paste(t1, collapse=" ")
newld = ld
newld[1] = new.head
newld[2] = new.vars

return(newld)
  
}



remove.SA.from.est = function(est, remove.years){

#-------------------------------------------------------------
#
# This functions takes a CASAL estimation file and removes
# sub-Antarctic trawl survey observations (SAbio) from it
#  
# Arguments
#
# est:    CASAL estimation file read in by extract.csl.file
#
# remove.years:  fishing years for which to remove the SAbio survey
  
#
# Value
#
# Output is an estimation file with the asked for SAbio observations removed
#
# Further Details
#
# A SAbio survey in Dec 2005  is in the 2006 fishing year
#
# Code is adapted from that for prune.est  
#
#--------------------------------------------------------------
  

#
# relative abundance SA observation
#
  
obs.name = "relative_abundance[SAsumbio]"  

this.obs = eval(parse( text=paste("est$'",obs.name,"'", sep="") ))  
invalid.years = as.character(remove.years)
if(length(invalid.years>0)){
    for(invalid.year in invalid.years){
      delete.obs = paste("this.obs$'",    invalid.year, "' = NULL", sep="")
      delete.cv  = paste("this.obs$'cv_", invalid.year, "' = NULL", sep="")
      eval(parse(text=paste(text=delete.obs)))
      eval(parse(text=paste(text=delete.cv)))      
    }
}

this.obs$years = this.obs$years[!(this.obs$years %in% invalid.years)]
obs.set = paste("est$'", obs.name, "'= this.obs", sep="")
eval(parse(text=obs.set))


#
# SA proportions at observation (trawl survey)
#
 
obs.name = "proportions_at[SAsumage]"

this.obs = eval(parse( text=paste("est$'",obs.name,"'", sep="") ))  
invalid.years = as.character(remove.years)
if(length(invalid.years>0)){
    for(invalid.year in invalid.years){
      delete.obs = paste("this.obs$'",    invalid.year, "' = NULL", sep="")
      delete.cv  = paste("this.obs$'cvs_", invalid.year, "' = NULL", sep="")
      eval(parse(text=paste(text=delete.obs)))
      eval(parse(text=paste(text=delete.cv)))      
    }
}

this.obs$years = this.obs$years[!(this.obs$years %in% invalid.years)]
obs.set = paste("est$'", obs.name, "'= this.obs", sep="")
eval(parse(text=obs.set))

return(est)

}    # end of function



make.moab.script = function(new.dir, path.dir){
 
  filename = "moab.MCMC"
  cat('#!/bin/bash\n\n', file=filename)
  
  for(chain in 1:3){
  
    cat('WD="$(pwd)"\n', file=filename, append=T)
    cat(paste('RUN="~/bin/casalMprior3new -i pars -m -q -g', round(1E8*runif(1)), '"\n'), file=filename, append=T)
    cat('echo $RUN > RUN\n', file=filename, append=T) 
    cat('CMD="msub -l nodes=1 -l walltime=150:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o MCMC.log.%j -e MCMC.err.%j -S /bin/bash RUN"\n', file=filename, append=T)
    cat(paste('echo "Running MCMC',chain, 'for CASAL on MOAB in directory:" $WD\n'), file=filename, append=T)
    cat('echo -n "Job started at: " ; date\n', file=filename, append=T)
    cat('echo $RUN\n', file=filename, append=T)
    cat('COMMAND="cd $WD ; $CMD"\n', file=filename, append=T)
    cat('ssh turbine $COMMAND\n', file=filename, append=T)
    cat('sleep 0.5\n\n', file=filename, append=T)
    
  }
  
  mcmc.script.name = paste("moab.", new.dir, sep="")
  shell(paste("cp ", filename, " '", path.dir, "\\", mcmc.script.name, "'", sep=""))
  
  return(NULL)

}



## Utility funciton for calculating upper and lower bounds for a multinomial value
## "Normal Approximation Method" of the Binomial Confidence Interval:
multi_ci = function(P, N, alpha = 0.05) {
  SE = sqrt((P * (1 - P)) / N);
  half_z = abs(qnorm(alpha / 2)); ## symmetric so abs is fine
  return(list = c(U_CI = P + half_z * SE, L_CI = P - half_z * SE))
}


#' Utility extract function
#'
#' @author Dan Fu
#'
"convert.to.lines" <-
function(filename) {
  scan(filename, what = "", sep = "\n")
}


#' Utility for extract function
#'
#' @author Craig Marsh
#'
strip = function(x) {
      tmp <- unlist(strsplit(x, "\t"))
      tmp <- unlist(strsplit(tmp, " "))
      return(as.vector(paste(tmp, collapse = " ")))
} 

#' Utility extract function
#'
#' @author Dan Fu
#'
"string.to.vector.of.words" <-
function(string)
{
  temp <- unpaste(string, sep = " ")
  return(temp[temp != ""])
}

#' Utility extract function
#'
#' @author Dan Fu
#'
"unpaste" <-
function(string, sep)
{
    return(unlist(strsplit(string, split = sep)))

}


Paste  = function( ..., sep="" ) paste( ..., sep = sep )

evalit = function(x) eval(parse(text=x))

