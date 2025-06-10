This code contains all necessary packages and functions to replicate statistics and graphs pertaining to the manuscript entitled: "Skeletal muscle properties in long COVID and ME/CFS differ from those induced by bed rest" (Charlton et al).
Code has been tested on R Versions 4.4.1 and 4.5.0.
Typical run time estimated: 5 minutes


Prior to Running Code: 
Install latest version of R and R Studio
Manually input working directory at line 46
Manually input date at line 45
Ensure that the file "Manuscript_data_clean_230525.xlsx" is located in the working directory
If needed, adjust the output folder

To run code: 
Download excel file to working directory, input working directory and date, and run all.

Commented areas are not necessary for replication of statistics or graphs. These can be uncommented for assessing exact means and standard deviations, assessing normality of data or comparing correlations, however they have been commmented to ensure quicker running of R code. 

Output includes all graphs found in both main and supplementary figures, .csv files containing significance tests found in figures. Statistical tests related to tables can be found in lines 332-400, 2131-2160, and 5058-5088.

