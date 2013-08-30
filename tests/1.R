library(assertive)
library(RUnit)
library(gtools)
library(simsem)

#
# Tests that all the relevant vignettes in Version03 directory (matrix syntax) of SimSem
# can be executed with matrixPLS
#

test.simsem.Version03 <- function()
{
	DEACTIVATED('Not implemented')
	execute.simsem.tests("simsem/SupportingDocs/Examples/Version03")
}

#
# Tests that all the relevant vignettes in Version05 directory (Lavaan syntax) of SimSem
# can be executed with matrixPLS
#

test.simsem.Version05 <- function()
{
	DEACTIVATED('Not currenly in use')
	execute.simsem.tests("simsem/SupportingDocs/Examples/Version05")
}

#
# Tests that all the relevant vignettes in Version05mx directory (OpenMX syntax) of SimSem
# can be executed with matrixPLS
#

test.simsem.Version05mx <- function()
{
	DEACTIVATED('Not implemented')
	execute.simsem.tests("simsem/SupportingDocs/Examples/Version05mx")
}
#
# A helper function to execute the SimSem tests
#

execute.simsem.tests <- function(directory){
	
	vignettes <- list.files(path=directory, recursive = TRUE)

	# Locate the vignette source files based on pattern
	vignettes <- grep("(ex[0-9]+)/\\1.R", vignettes, perl = TRUE, value = TRUE)
	
	assert_is_non_empty(vignettes)
	
	vignettes <- mixedsort(vignettes)
	
	for(vignette in vignettes){
		fileName <- paste(directory,vignette,sep="/")
		print(paste("Preparing vignette from file",fileName))

		vignetteCode <- readChar(fileName, file.info(fileName)$size)
		
		# Replace calls to sim to matrixpls.sim
		
		modifiedVignetteCode <- gsub("sim\\([0-9]+(.*))\n+?","matrixpls.sim(10\\1, stopOnError = TRUE)\n",vignetteCode)

		cat(modifiedVignetteCode)
		
		print(paste("Running vignette from file",fileName))
		
		eval(parse(text = modifiedVignetteCode))
		
		# Run the code
		tryCatch({
			checkTrue(TRUE, msg = fileName)
		}, error = function(e){
			cat(modifiedVignetteCode)
			checkTrue(FALSE, msg = fileName)
		})
		
	}
}

#
# Checks that matrixPLS produces results that are identical to plspm
#

test.plspm <- function()
{
	
	library(plspm)

	# Run the example from plspm package
	
	# load dataset satisfaction
	data(satisfaction)
	# inner model matrix
	IMAG = c(0,0,0,0,0,0)
	EXPE = c(1,0,0,0,0,0)
	QUAL = c(0,1,0,0,0,0)
	VAL = c(0,1,1,0,0,0)
	SAT = c(1,1,1,1,0,0)
	LOY = c(1,0,0,0,1,0)
	sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
	# outer model list
	sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
	# vector of modes (reflective indicators)
	sat_mod = rep("A", 6)

	# apply plspm
	satpls.plspm = plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=FALSE)
	
	# apply MatrixPLS
	satpls.matrixpls = matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=FALSE)

	checkEquals(satpls.matrixpls, satpls.plspm)

	# Repeat with bootstrapping

	# apply plspm
	satpls.plspm = plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
	
	# apply MatrixPLS
	satpls.matrixpls = matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
	
	checkEquals(satpls.matrixpls, satpls.plspm)
	
}

#
# Checks that matrixPLS produces results that are identical to semPLS
#

test.semPLS <- function()
{
	DEACTIVATED('Not implemented')
	library(semPLS)

	#checkEquals()
}

#
# Checks that the PLS algorithm converges to known results
# when using population data
#

test.populationValues <- function()
{
	DEACTIVATED('Not implemented')
	#checkEquals()
}

test.plspm()
