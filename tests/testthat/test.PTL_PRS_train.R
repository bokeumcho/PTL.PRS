test_that('Train_PTLPRS_multiSeeds', {
    base_dir = 'testdata/'
    cfg <- PTL_PRS_config(
            ref_file = paste0(base_dir, '1kg_SAS_ref_200'),
            ref_file_ps = paste0(base_dir, '1kg_SAS_ref_200'),
            sum_stats_file = paste0(base_dir, 't2d_prscs_EUR_200snps.pgs'),
            target_sumstats_file = paste0(base_dir, 'tr_val_200snps.summaries'),
            LDblocks="ASN.hg19",
            outfile=paste0("Output_test"),
            num_cores=1,
            lr_list = c(1,10,100),
            iter = 5,
            ps=TRUE,
            subprop=0.8235, #0.7:0.15
            pseudo_test = TRUE,
            patience = 3, #default = 3
            trace=FALSE,
            target_sumstats_train_file = NULL, #'Output_test_target.train.summaries' 
            target_sumstats_val_file = NULL #'Output_test_target.val.summaries'    
            )

    res <- PTL_PRS_run(cfg, pseudo_test=TRUE, num_repeat_ps=3, random_seeds_repeat=seq(from=10, to=30, by=10))
    print(res$summary)
})

test_that('Test_PTLPRS_multiSeeds', {
    base_dir = 'testdata/'
    
    outfile=paste0("Output_test")
    
    Y_name="Y"
    Ytype='B'
    Covar_name="sex"

    sum_stats_file = paste0(base_dir, 't2d_prscs_EUR_200snps.pgs')
    sum_stats = data.frame(fread(sum_stats_file))
    
    ref_file = paste0(base_dir, '1kg_SAS_ref_200')

    random_seeds_repeat=seq(from=10, to=30, by=10)

    for (seed in random_seeds_repeat) {
        outfile_seed = paste0(outfile, '_seed', seed)
        target_sumstats_test_file = paste0(outfile_seed,'_target.test.summaries') 
        sum_stats_target_test=fread(target_sumstats_test_file)

        best.beta = fread(paste0(outfile_seed, "_best.beta.txt"))
        best.param = fread(paste0(outfile_seed, "_best.param.txt"))

        cat('[TEST] seed:', seed, ', outfile prefix:', outfile_seed,'\n')
        PTL_PRS_test_pseudo(best.beta, best.param, sum_stats_target_test, sum_stats, ref_file, chunk_size = dim(best.beta)[1])
    }
})
