
test_that('Test_PTLPRS', {
plink_file = 'ExampleData_TL-PRS/ExampleData_1000G_African_test'
ped_test_file='ExampleData_TL-PRS/ped_test.txt'

ped_test=fread(ped_test_file,header=T, fill=TRUE)

Y_name='Y'
Ytype="B"
Covar_name="" 

TL_PRS_test(plink_file, ped_test_file, best.beta, best.param, Ytype, Covar_name,Y_name)
})