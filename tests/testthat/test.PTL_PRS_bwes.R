test_that('run_PTLPRS', {
random_seed=42

ref_file = "ExampleData_TL-PRS/ExampleData_1000G_African_train"
ref_file_ps = "ExampleData_TL-PRS/ExampleData_1000G_African_train"

sum_stats_file = 'ExampleData_TL-PRS/Beta_lassosum_usingEuropeanSummaryStatistics_modified.txt'
target_sumstats_file = 'ExampleData_TL-PRS/African_SummaryStatistics.txt'

LDblocks="EUR.hg19"
outfile=paste0("Output_test")

num_cores=1
lr_list = c(1,10) #default c(1,10,100,1000,5000,10000)
# lr_list = c(1,10,100,1000,5000)
iter = 5

ps=TRUE
subprop=0.8235 #0.7:0.15

patience = 3 #default = 3
trace=TRUE

target_sumstats_train_file = NULL #paste0('data/T2D_',comp,'_WB.summaries')
target_sumstats_val_file = NULL #paste0('data/mapp_modified/t2d_',comp,'_val_',random_seed,'.summaries')

out.beta=PTL_PRS_bwes(ref_file,sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks,outfile,num_cores, target_sumstats_train_file, target_sumstats_val_file, ps,random_seed, lr_list, iter, early_stopping, patience, trace)
})

