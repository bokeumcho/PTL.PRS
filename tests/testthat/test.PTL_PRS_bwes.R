# test_that('load Rcpp Module', {
#   ref_file = "ExampleData_TL-PRS/ExampleData_1000G_African_train"

#   # bedFileReader <- Rcpp::Module("BedFileReader_module", PACKAGE = "PTL.PRS")
#   # BedFileReader = loadModule(module = "BedFileReader_module", TRUE)

#   # setClass("BedFileReader", representation( pointer = "externalptr" ) )
#   # BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
#   # setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
#   # setMethod("initialize","BedFileReader", function(.Object, ...) {
#   #   .Object@pointer <-.Call(BedFileReader_method("new"), ... )
#   #   .Object})

#   bedFileReader <- new(BedFileReader, paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))
#   result = try(bedFileReader$snp_index_func(), silent=TRUE)
#   Gtemp = readSomeSnp(c('rs10045497','rs10401969'), BedFileReader = bedFileReader) #ref
#   print(head(Gtemp))
# })

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
pseudo_test = FALSE

patience = 3 #default = 3
trace=TRUE

target_sumstats_train_file = NULL 
target_sumstats_val_file = NULL 

out.beta=PTL_PRS_train(ref_file,sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks,outfile,num_cores, target_sumstats_train_file, target_sumstats_val_file, ps,pseudo_test, random_seed, lr_list, iter, patience, trace)
})