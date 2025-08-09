test_that('Train_PTLPRS', {
random_seed=80

base_dir = 'testdata/'
ref_file = paste0(base_dir, '1kg_SAS_ref_200')
ref_file_ps = paste0(base_dir, '1kg_SAS_ref_200')

sum_stats_file = paste0(base_dir, 't2d_prscs_EUR_200snps.pgs')
target_sumstats_file = paste0(base_dir, 'tr_val_200snps.summaries')
LDblocks="ASN.hg19"
outfile=paste0("Output_test")

num_cores=1
lr_list = c(1,10,100,1000)
iter = 5

ps=TRUE
subprop=0.8235 #0.7:0.15
pseudo_test = TRUE

patience = 3 #default = 3
trace=FALSE

target_sumstats_train_file = NULL #'Output_test_target.train.summaries' 
target_sumstats_val_file = NULL #'Output_test_target.val.summaries'

out.beta=PTL_PRS_train(ref_file,sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks,outfile,num_cores, target_sumstats_train_file, target_sumstats_val_file, ps,pseudo_test, random_seed, lr_list, iter, patience, trace)
})

test_that('Test_PTLPRS', {
    base_dir = 'testdata/'
    
    outfile=paste0("Output_test")
    
    sum_stats_file = paste0(base_dir, 't2d_prscs_EUR_200snps.pgs')
    sum_stats = data.frame(fread(sum_stats_file))
    
    target_sumstats_test_file = paste0(outfile, '_target.test.summaries') 
    sum_stats_target_test=fread(target_sumstats_test_file)

    ref_file = paste0(base_dir, '1kg_SAS_ref_200')

    best.beta = fread(paste0(outfile,"_best.beta.txt"))
    best.param = fread(paste0(outfile,"_best.param.txt"))
    
    Y_name="Y"
    Ytype='B'
    Covar_name="sex"

    outfile=paste0("Output_test")

    system.time({
        PTL_PRS_test_pseudo(best.beta, best.param, sum_stats_target_test, sum_stats, ref_file, chunk_size = dim(best.beta)[1])
    })
})