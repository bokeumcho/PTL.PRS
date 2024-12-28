
test_that('Test_PTLPRS', {
plink_file = 'ExampleData_TL-PRS/ExampleData_1000G_African_test'

by_chr = FALSE

ped_test_file='ExampleData_TL-PRS/ped_test.txt'

Y_name='Y'
Ytype="C"
Covar_name="" 

outfile=paste0("Output_test")

PTL_PRS_test(plink_file, by_chr, ped_test_file, outfile, Ytype, Covar_name,Y_name)
})