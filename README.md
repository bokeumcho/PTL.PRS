# TLPRS

This R package helps users to construct multi-ethnic polygenic risk score (PRS) using transfer learning. It can help predict PRS of minor ancestry using summary statistics from exsiting resources, such as UK Biobank.

This package contains two main functions: TL_PRS and PTL_PRS.

## Installation
`TL-PRS` requires the software 'plink' as well as the following R packages:  `data.table`, `lassosum`, `Rcpp` and `parallel`. Install them by: 

```r
install.packages(c("data.table", "lassosum", "parallel", "Rcpp"), dependencies=TRUE)
```

If you have `devtools`, you can type: 
```r
install_github("bokeumcho/TLPRS")
```
for the latest development version. Or you can clone the latest development version here and install yourself using `devtools`. 


## Inputs of TL_PRS
1. `ped_val_file`:
The location path of ped file, where contains the information of FID, IID, outcome (Y), index of samples in the .fam file of `plink_file` (**ukb_idx**) and covariates if needed. Note that the file only requires samples for validation.

2. `Covar_name`:
A vector of names of covariates we need to adjust in the model, such as c("Sex","Age"). Note that all names must corrspond to the columns in the ped file.

3. `Y_name`: 
The name of Y in the model, such as "LDL". Note that the Y name must corrspond to a column in the ped file.

4. `Ytype`: 
The type of Y should be either "C"(continuous) or "B"(binary).

5. <span style="background-color:#fff5b1"> **`ref_file`:** </span>
The prefix of plink file of the reference panel data in the target population. Note that this file will be used for estimating LD matrix in the training step.

6. **`plink_file`:**
The prefix of plink file that includes the validation data. Note that there is no need to generate exclusive plink file for validation samples as the package selectively reads the data of samples in ped file.

7. `sum_stats_file`:
The location path of effect size file. We usually can obtain this file by existing PRS methods, such as lassosum/PRS-CS. Specifically it contains the following three columns:"SNP",**"CHR"**,"A1","Beta". "SNP" is the SNPID (the same format as SNPID in plink files); "CHR" is the chromosome number of each SNP; "A1" is the alternative (effect) allele; "Beta" is the effect size. 

8. `target_sumstats_file`:
The location path of summary stats file for the target population. The file requires the following columns: "SNP", "A1", "beta", "N", and "p", wheren "beta" is the effect size from summary statistics, "N" is the sample size used for calculating summary statistics, and "p" is p-value of the SNP. 

9. `LDblocks`:
This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome. Currently we support three types:"EUR.hg19","AFR.hg19","ASN.hg19", corresponding to European, African and Asian populations, respectively.

10. `outfile`:
The prefix of the file location which can be used to store the final output files. Note that the function needs to save files in this directory.

11. `cluster`:
A cluster object from the parallel package for parallel computing

12. **`lr_list` (optional):**
A grid of learning rates used for hyper parameter tuning. If it is not specified, a default value is min( 1,10,100,1000,5000,10000 / n(SNPs) , 1).

13. **`iter` (optional):**
Maximum number of iterations to decide where to stop training early. If it is not specified, a default value is 200.

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

## Outputs of TL_PRS
1. `best.learning.rate`: 
the learning rate we can use in order to achieve the best risk prediction.

2. `best.iteration`: 
the number of iterations we should stop in order to achieve the best risk prediction.

3. `best.beta`: 
the data frame containing three columns: "SNP","A1","beta". Note that this is the best effect size we can use to construct PRS, selected using best.learning.rate and best.iteration.  

4. `param_table`: 
the data frame containing a grid of candidates of learning rates and the number of iterations that we consider. 
