# PTL.PRS

This R package helps users to construct multi-ethnic polygenic risk score (PRS) using transfer learning when the individual level data is not available. It can help predict PRS of small group using summary statistics from larger group resources.

This package contains three main functions: `PTL_PRS_bwes` for train, `TL_PRS_test` for test, and `pseudo_split` for pseudo splitting.

## Installation
`PTL.PRS` requires the software 'plink' as well as the following R packages:  `data.table`, `lassosum`, `Rcpp` and `parallel`. Install them by: 

```r
install.packages(c("data.table", "lassosum", "parallel", "Rcpp"), dependencies=TRUE)
```

If you have `devtools`, you can type: 
```r
install_github("bokeumcho/PTL.PRS")
```
for the latest development version. Or you can clone the latest development version here and install yourself using `devtools`. 


## Inputs of PTL_PRS
1. `Y_name`: 
The name of Y in the model, such as "LDL". Note that the Y name must corrspond to a column in the ped file.

2. `Ytype`: 
The type of Y should be either "C"(continuous) or "B"(binary).

3. `ref_file`:
The prefix of plink file of the reference panel data in the target population. Note that this file will be used for estimating LD matrix in training and validation steps.

4. **`ref_file_ps`:**
The prefix of plink file of the reference panel data in the target population. Note that this file will be used for generating pseudo-summary statistics. It can be the same file as `ref_file` but the model is more reliable with distinct file.

5. `sum_stats_file`:
The location path of effect size file. We usually can obtain this file by existing PRS methods, such as lassosum/PRS-CS. Specifically it contains the following three columns:"SNP","CHR","A1","Beta". "SNP" is the SNPID (the same format as SNPID in plink files); "CHR" is the chromosome number of each SNP; "A1" is the alternative (effect) allele; "Beta" is the effect size. 

6. `target_sumstats_file`:
The location path of summary stats file for the target population. The file requires the following columns: "SNP", "A1", "beta", and "N", where "beta" is the effect size from summary statistics, "N" is the sample size used for calculating summary statistics. **Note that it is used to generate pseudo-summary statistics.**
* option to choose whether to generate pseudo-summary stats or to provide input for separate training/validation summary stats?

7. `LDblocks`:
This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome. Currently we support three types:"EUR.hg19","AFR.hg19","ASN.hg19", corresponding to European, African and Asian populations, respectively.

8. `outfile`:
The prefix of the file location which can be used to store the final output files. Note that the function needs to save files in this directory.

9. `cluster`:
A cluster object from the parallel package for parallel computing

10. `lr_list` (optional):
A grid of learning rates used for hyper parameter tuning. If it is not specified, a default value is min( 1,10,100,1000,5000,10000 / n(SNPs) , 1).

11. `iter` (optional):
Maximum number of iterations to decide where to stop training early. If it is not specified, a default value is 200.

## Outputs of PTL_PRS
1. `best.learning.rate`: 
the learning rate we can use in order to achieve the best risk prediction.

2. `best.iteration`: 
the number of iterations we should stop in order to achieve the best risk prediction.

3. `best.beta`: 
the data frame containing three columns: "SNP","A1","beta". Note that this is the best effect size we can use to construct PRS, selected using best.learning.rate and best.iteration.  

4. `param_table`: 
the data frame containing a grid of candidates of learning rates and the number of iterations that we consider. 

## Example R script
You can use the files in example_data folder to run the example R script.

```r
library(data.table)
library(lassosum)
library(parallel)
library(Rcpp)

library(PTL.PRS)

setwd("bokeum/PTL.PRS/example_data")  ###setup your work path.

ref_file = '1kg_SAS_ref' 
ref_file_ps = '1kg_SAS_ref'

num_cores = 2
pseudo_test=FALSE
random_seed=42

sum_stats_file = 't2d_prscs_EUR.pgs' # randomly sampled 200 SNPs from t2d prscs weights
target_sumstats_file = 'tr_val.summaries' 

LDblocks="ASN.hg19"
outfile=paste0("Output_PTLPRS") 


lr_list = c(1,10,100,1000) 
iter = 100

ps=TRUE
subprop=0.9

patience = 3 
trace=TRUE

target_sumstats_train_file = NULL 
target_sumstats_val_file = NULL

out.beta = PTL_PRS_bwes(
  ref_file,
  sum_stats_file,
  target_sumstats_file,
  subprop,
  ref_file_ps,
  LDblocks,
  outfile,
  num_cores,
  target_sumstats_train_file,
  target_sumstats_val_file,
  ps,
  pseudo_test,
  random_seed,
  lr_list,
  iter,
  patience,
  trace
)

### TEST ###
plink_file = '1kg_SAS_ALL'
by_chrom = FALSE

ped_test_file='test_ped.txt'

Y_name="Y"
Ytype='B'
Covar_name="sex"

best.beta = fread(paste0(outfile,"_best.beta.txt")) #out.beta$best.beta 
best.param = fread(paste0(outfile,"_best.param.txt")) #out.beta$best.param 

TL_PRS_test(plink_file, by_chrom, ped_test_file, best.beta, best.param, Ytype, Covar_name,Y_name)
```

## Support
If there are any further questions or problems with running or installing `PTL.PRS`, please do email me at <bokeum1810@snu.ac.kr>. 
