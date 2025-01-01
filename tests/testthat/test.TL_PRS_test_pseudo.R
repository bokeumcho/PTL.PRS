
# test_that('Test_PTLPRS', {
#     # test pseudo test 
#     ref_file = '/media/leelabsg-storage1/bokeum/TLPRS/1kg_phase3/SAS/1kg_SAS_ALL' 
#     outfile=paste0("/media/leelabsg-storage1/bokeum/TLPRS/output/covid/PTLPRS_SAS_prscs_pseudoTEST_genomicc_30")

#     target_sumstats_test_file = paste0(outfile, '_target.test.summaries') 
#     sum_stats_target_test=fread(target_sumstats_test_file)
    
#     sum_stats_file = paste0("/media/leelabsg-storage1/bokeum/TLPRS/data/covid/covid_prscs_genomicc.pgs")    
#     sum_stats = data.frame(fread(sum_stats_file))
    
#     best.beta = fread(paste0(outfile,"_best.beta.txt"))
#     best.param = fread(paste0(outfile,"_best.param.txt")) #out.beta$best.param 

#     PTL_PRS_test_pseudo(best.beta, best.param, sum_stats_target_test, sum_stats, ref_file)
# })

# test_that('Test_PTLPRS_list', {
#     # test pseudo test 
#     ref_file = '/media/leelabsg-storage1/bokeum/TLPRS/1kg_phase3/SAS/1kg_SAS_ALL' 
#     outfile=paste0("/media/leelabsg-storage1/bokeum/TLPRS/output/covid/PTLPRS_SAS_prscs_pseudoTEST_genomicc_30")

#     target_sumstats_test_file = paste0(outfile, '_target.test.summaries') 
#     sum_stats_target_test=fread(target_sumstats_test_file)
    
#     sum_stats_file = paste0("/media/leelabsg-storage1/bokeum/TLPRS/data/covid/covid_prscs_genomicc.pgs")    
#     sum_stats = data.frame(fread(sum_stats_file))
    
#     beta.list = fread(paste0(outfile,"_beta.list.txt"))

#     PTL_PRS_test_pseudo_list(beta.list, sum_stats_target_test, sum_stats, ref_file)
# })