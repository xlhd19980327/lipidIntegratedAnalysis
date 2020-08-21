# lipidIntegratedAnalysis
Our lipidomics &amp; RNA-seq integrated omics analysis tool
## File system annotation
- **README.md**
- *original*
  - **procedure_LipidSearch.R** 
  for LipidSearch output lipidomics data
  - **procedure_MS_DIAL.R**
  for MS_DIAL output lipidomics data
  - **modifyFAchainHandle_procedure.R**
  for FA chain info statistics and plot(LipidSearch data procedure)
- *lipidpreproc*
  - ** **
  for preprocess the data 
- *tileplot*
  - **FAchainVisual.R**
  [from modifyFAchainHandle_procedure.R]
  Use FA chain info(chain length & unsaturated info) to do some visulization.
  After tidy the data, plot the integrated info(contain both chain length & unsaturated info) of all FA info(MS1+MS2 info)to "integPlot"
  
## Notes
First use the lipsearch data to fulfill all the procedure

## Coding style
### Comment style
- **Client notes** will be writen in the top of the whole scripts following a "NOTE" comment. 
Corresponding to that note, the code afterwards will comment as a "NOTE-ref".
- **Warning info** refer to some unsolved problems.
The comment info will start with a "!!!!!WARNING: ".
- **Client options** will be comment in a line with the following words:
"!!!Client options". 
After the comment, there is a function that will offer optional arguments to clients.
- **Control flow WARNING** will be comment in a line with following words:
"!!!!!Control flow WARNING".
Details are shown in the "Script style" corresponding section
### Script style
- **PROGRAM EXIT!** In some situation, scripts may use "stop()" function to make program exit.
These scripts may contain the words "PROGRAM EXIT!".

eg. 
``` R
if(min(nsamples) < 3){
  stop("At least one group have no more than 2 replicates, PROGRAM EXIT!")
}
``` 
- **Control flow WARNING** The code somewhere contains a group by group analysis, 
which will use control flow. The code only shows the control flow style but contains wrong coding logical structure.
The control flow WARNING will indicate this code below the comment "!!!!!Control flow WARNING".
In the ultimate version, these codes should be revised to the correct structure.

- **Client options** A function having client options will following formats:
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
