# lipidIntegratedAnalysis
Our lipidomics &amp; RNA-seq integrated omics analysis tool
## File system annotation
- **README.md**
- **.gitignore**
  Indicate the files that git will ignore.
- *original*
  - **procedure_LipidSearch.R** 
  for LipidSearch output lipidomics data
  - **procedure_MS_DIAL.R**
  for MS_DIAL output lipidomics data
  - **modifyFAchainHandle_procedure.R**
  for FA chain info statistics and plot(LipidSearch data procedure)
- *lipidpreproc*
  For preprocess the lipidomics data, basic statistics analysis. 
  - **LipidSearchData.R**
  Some old codes that may be used as some reference.
  - **LipidSearchData_preprocess.R**
  Pre-process the data to be used as the pre-processed input data of other scripts. It may used with the R code form "source("LipidSearchData_preprocess.R")" seen in other scripts.
  - **LipidSearchData_useMAR.R**
  uses MetaboAnalystR package to do some basic statistics and visualization.
- *HeadgroupPlot*
  For visualize the total concentrations of lipid classes.
  - **facetPlot.R**
  uses facet to do visualization.
- *FAplot*
  For visualize the concentrations of lipid "sub-classes".
  - **FAchainVisual.R**
  [originate from modifyFAchainHandle_procedure.R]
  Use FA chain info(chain length & double-bind info) to do some visulization.
  Plot the integrated info(contain both chain length & unsaturated info) of all FA info(MS1+MS2 info) to "integrated plot".
  Use facet to visualize FA chain and double-bind, 
  also use bubble plot to do more beautiful visualization.
- *testData*
  Test data. Every test data will be in the directory which have a understandble name to describe the data. The directiry will also have two parts: "input" and "output". The input data will in the "input" directory. The results files which the program generates will be in the "output" directory.
  
  
- *utilityFunc*
  Some functions that many scripts may source from. The functions can be used directly for the ultimate version in this directory.
- *dev*
  Some codes only useful for the developers. May be integrated to the ultimate version later.
  - **setwd.R**
  For debugging the code in diffrent platforms.
- *sundries*
  Some sundries. May useful for the developers.
  
  
## Notes
First use the lipsearch data to fulfill all the procedure

## Coding style
### Comment style
- **Client notes** refer to what the client should know first about the input data.
They will be written in the top of the whole scripts following a "NOTE" comment. 
Corresponding to that note, the code afterwards will comment as a "NOTE-ref".
- **Warning info** refer to some unsolved problems or some premature methods.
The comment info will start with a "!!!!!WARNING: ".
- **Client options** will be comment in a line with the following words:
"!!!Client options". 
After the comment, there are some codes that will offer optional arguments to clients or require clients to enter some input info. Some may also have description about them. 
- **Control flow WARNING** will be comment in a line with following words:
"!!!!!Control flow WARNING".
Details are shown in the "Script style" corresponding section.
- **Developer notes** refer to some info given to developers. 
They will info the developers to do some specific actions to the following codes.
The comment info will start with a "!!!!!DEV: ".
In the ultimate version, these codes may be revised in someways.
- **Source info** indicate what contents the source will offer to the current environment. 
It will be comment as "Source will offer the following contents:" and indicate the contents in the next line(s).


### Script style
- **PROGRAM EXIT!** In some situations, scripts may use "stop()" function to make program exit.
These scripts may contain the words "PROGRAM EXIT!".

eg. 
``` R
if(min(nsamples) < 3){
  stop("At least one group have no more than 2 replicates, PROGRAM EXIT!")
}
``` 
- **Control flow WARNING** The code somewhere contains a group by group analysis, which will use control flow. The code only shows the control flow style but contains wrong coding logical structure.
The control flow WARNING will indicate this code below the comment "!!!!!Control flow WARNING".
In the ultimate version, these codes should be revised to the correct structure.

- **Client options(May NOT applied to later development)** A function having client options will have the following formats:
  - The first line of the function will contain the necessary arguments.
  - The second line of the function  will contain the client-optional arguments. 
  The default parameters will also show there.
  - The third line(if have) of the function will contain the unnecessary arguments and cannot be changed.

eg.
``` R
mSet<-Normalization(mSet,    ## necessary arguments
                    rowNorm = "ProbNormT", transNorm = "NULL", scaleNorm = "AutoNorm", ref = NULL,   ## client-optional arguments
                    ratio=FALSE, ratioNum=20)   
``` 
