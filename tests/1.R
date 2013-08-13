library(assertive)
library(RUnit)

#
# Tests that all the relevant vignettes in Version03 directory (matrix syntax) of SimSem
# can be executed with matrixPLS
#

test.simsem.Version03 <- function()
{
	execute.simsem.tests("simsem/SupportingDocs/Examples/Version03")
}

#
# Tests that all the relevant vignettes in Version05 directory (Lavaan syntax) of SimSem
# can be executed with matrixPLS
#

test.simsem.Version05 <- function()
{
	execute.simsem.tests("simsem/SupportingDocs/Examples/Version05")
}

#
# Tests that all the relevant vignettes in Version05mx directory (OpenMX syntax) of SimSem
# can be executed with matrixPLS
#

test.simsem.Version05mx <- function()
{
	execute.simsem.tests("simsem/SupportingDocs/Examples/Version05mx")
}
#
# A helper function to execute the SimSem tests
#

execute.simsem.tests <- function(directory){

	library(gtools)
	library(simsem)
	
	vignettes <- list.files(path=directory, recursive = TRUE)

	# Locate the vignette source files based on pattern
	vignettes <- grep("(ex[0-9]+)/\\1.R", vignettes, perl = TRUE, value = TRUE)
	
	assert_is_non_empty(vignettes)
	
	vignettes <- mixedsort(vignettes)
	
	for(vignette in vignettes[1]){
		fileName <- paste(directory,vignette,sep="/")
		print(paste("Preparing vignette from file",fileName))

		vignetteCode <- readChar(fileName, file.info(fileName)$size)
		
		# Replace calls to sim to matrixpls.sim
		
		modifiedVignetteCode <- gsub("sim\\(","matrixpls.sim(",vignetteCode)
		
		# Run the code
		checkException(eval(parse(text = modifiedVignetteCode)), msg = fileName);
		
	}
}

#
# Checks that matrixPLS produces results that are identical to plspm
#

test.plspm <- function()
{
	DEACTIVATED('Not implemented')
	library(plspm)

#	checkEquals()
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
