# TGMI

## Instructions for running Triple Gene Mutual Interaction Algorithm (TGMI)

1. TGMI need an extensive computing resource, we recommend running TGMI under a multi-core server environment. In order to run BWERF, please install following R packages. 

- install.packages("stats")
- install.packages("MASS")
- install.packages("iterators")
- install.packages("parallel")
- install.packages("foreach")
- install.packages("doMC")
- install.packages("infotheo")

2. Open file "run_TGMI.R" and modify the parameters such as input data files (recommended to use the sample data format) and number of cores allocated.

3. Use the following command to run the R script in the background

**nohup R CMD BATCH run_TGMI.R &**
