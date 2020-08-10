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
- **Client notes** will be writen in the top of the code following a "NOTE" comment. 
Corresponding to that note, the code afterwards will comment as a "NOTE-ref".
- **Warning info** refer to some unsolved problems.
The comment info will start with a "!!!!!WARNING: ".