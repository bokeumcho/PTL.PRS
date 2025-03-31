library(data.table)
library(lassosum)
library(Matrix)
library(ROCR)
library(parallel)

PTL_PRS_test <- function(plink_file, by_chr, ped_test_file, outfile, Ytype, Covar_name,Y_name, plink_etc=''){
    #' Test performance of PTL-PRS
    #' @export
    #' 
  ped_test = fread(ped_test_file,header=T, fill=TRUE)

    best.beta = fread(paste0(outfile,"_best.beta.txt")) #out.beta$best.beta 
    best.param = fread(paste0(outfile,"_best.param.txt")) #out.beta$best.param 

    if (by_chr) {
    BedFileReaders <- BedFileReader_prep(plink_file, by_chr, unique(best.beta$CHR), plink_etc)

    ## by CHRs -- function!!
    PRS.all <- matrix(0, nrow=dim(ped_test)[1], ncol=2) 

    # Loop through each unique value in LDblocks2[[1]]
    for (i in unique(best.beta$CHR)) {
        # Calculate the partial result
        partial_PRS <- Calculate_PRS_direct(
            ped_test, 
            best.beta[best.beta$CHR==i,1:3], 
            best.beta[best.beta$CHR==i,4:5],
            BedFileReaders[[i]],
            plink_file, by_chr, plink_etc
            )

        # print(dim(partial_PRS))
        # Add the partial result to the cumulative total
        PRS.all <- PRS.all + partial_PRS  # Adjust this operation based on how results are combined
        rm(partial_PRS)
        gc()
    }
    } else {
        BedFileReaders <- BedFileReader_prep(plink_file, by_chr)

        PRS.all <- Calculate_PRS_direct(
            ped_test, 
            best.beta[,1:3], 
            best.beta[,4:5],
            BedFileReaders,
            plink_file, by_chr
            )
    }

#   PRS.all = cbind(ped_test[,c('FID','IID')],PRS.all)
#   colnames(PRS.all)[3:4]=c("SCORESUM","SCORESUM")

  #write.table(PRS.all, file=paste0(tempfile,"_PRS_test.txt"),row.names=F,quote=F,col.names=T)

  if (Ytype=="C"){
      base_R2 = as.numeric(linear_result_generator(PRS.all[,1],ped_test,Covar_name,Y_name))
      model_R2 = as.numeric(linear_result_generator(PRS.all[,2],ped_test,Covar_name,Y_name))
  } else {
      base_R2 = as.numeric(logistic_result_generator(PRS.all[,1],ped_test,Covar_name,Y_name))
      model_R2 = as.numeric(logistic_result_generator(PRS.all[,2],ped_test,Covar_name,Y_name))
    }

    cat("===Test Results===============================================================================\n")
    cat("Baseline R2: ", base_R2, "\n Model R2: ", model_R2, '\n Best LR: ', best.param[[1]]/dim(best.beta)[1]) #, '\n best iter: ', best.params[2],'\n'

  gc()
}

PTL_PRS_test_pseudo <- function(best.beta, best.param, sum_stats_target_test, sum_stats, ref_file){
    #' Test performance of PTL-PRS using pseudo test summary
    #' @export
    #' 

    # prep inputs
    best.beta1 = best.beta[,4:5]

    colnames(sum_stats) = c('SNP', 'V2','V3','V4')

    sum_stats_target_test=merge(sum_stats,sum_stats_target_test,by="SNP", sort=F)

    if (('p' %in% colnames(sum_stats_target_test)) & (!'cor' %in% colnames(sum_stats_target_test))) {
        sum_stats_target_test$beta = as.numeric(sum_stats_target_test$beta)
        sum_stats_target_test$p = as.numeric(sum_stats_target_test$p)

        sum_stats_target_test$cor=lassosum::p2cor(p = sum_stats_target_test$p, n = median(sum_stats_target_test$N,na.rm=T), sign=sum_stats_target_test$beta)
    }

    # allele flipping (only cor, not Beta, needs flipping)
    flag=which(sum_stats_target_test$V3 !=sum_stats_target_test$A1)
    if (length(flag)>0){sum_stats_target_test$cor[flag]=-sum_stats_target_test$cor[flag]}
    # CLUMPING
    # sum(sum_stats_target$cor>=0.2)
        
    sum_stats_target_test=sum_stats_target_test[,c("SNP","V3","cor","N")]

    colnames(sum_stats_target_test)[2]="A1"; 
    gc()

    # calculate
    sum_stats_target_test1 = sum_stats_target_test[sum_stats_target_test$SNP %in% best.beta$SNP,]
    best.beta1 = best.beta1[best.beta$SNP %in% sum_stats_target_test$SNP,]

    BedFileReader_ref <- new( BedFileReader, paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))
    result = try(BedFileReader_ref$snp_index_func(), silent=TRUE)

    geno_ref = readSomeSnp(sum_stats_target_test1$SNP, BedFileReader = BedFileReader_ref)
    geno_ref = scale(geno_ref)

    betaRho.byL = t(best.beta1)%*%sum_stats_target_test1$cor
    betaG.byL = geno_ref%*%as.matrix(best.beta1)

    R2.byL = betaRho.byL^2 * median(sum_stats_target_test1$N) / colSums(betaG.byL^2) 
    cat("baseline pseudo-R2: ", R2.byL[1,], "\n model pseudo-R2: ", R2.byL[2,], '\n best lr: ', best.param[[1]]/dim(best.beta)[1]) #, '\n best iter: ', best.params[2],'\n'
}

PTL_PRS_test_pseudo_list <- function(beta.list, sum_stats_target_test, sum_stats, ref_file){
    #' Test performance of PTL-PRS using pseudo test summary
    #' @export
    #' 

    # prep inputs
    m = ncol(beta.list)
    beta.list1 = beta.list[,9:m] 

    colnames(sum_stats) = c('SNP', 'V2','V3','V4')

    sum_stats_target_test=merge(sum_stats,sum_stats_target_test,by="SNP", sort=F)

    if (('p' %in% colnames(sum_stats_target_test)) & (!'cor' %in% colnames(sum_stats_target_test))) {
        sum_stats_target_test$beta = as.numeric(sum_stats_target_test$beta)
        sum_stats_target_test$p = as.numeric(sum_stats_target_test$p)

        sum_stats_target_test$cor=lassosum::p2cor(p = sum_stats_target_test$p, n = median(sum_stats_target_test$N,na.rm=T), sign=sum_stats_target_test$beta)
    }

    # allele flipping (only cor, not Beta, needs flipping)
    flag=which(sum_stats_target_test$V3 !=sum_stats_target_test$A1)
    if (length(flag)>0){sum_stats_target_test$cor[flag]=-sum_stats_target_test$cor[flag]}
    # CLUMPING
    # sum(sum_stats_target$cor>=0.2)
        
    sum_stats_target_test=sum_stats_target_test[,c("SNP","V3","cor","N")]

    colnames(sum_stats_target_test)[2]="A1"; 
    gc()

    # calculate
    sum_stats_target_test1 = sum_stats_target_test[sum_stats_target_test$SNP %in% beta.list$SNP,]
    beta.list1 = beta.list1[beta.list$SNP %in% sum_stats_target_test$SNP,]
    # cat(dim(sum_stats_target_test1),dim(beta.list1))

    BedFileReader_ref <- new( BedFileReader, paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))
    result = try(BedFileReader_ref$snp_index_func(), silent=TRUE)

    geno_ref = readSomeSnp(sum_stats_target_test1$SNP, BedFileReader = BedFileReader_ref)
    geno_ref = scale(geno_ref)

    betaRho.byL = t(beta.list1)%*%sum_stats_target_test1$cor
    betaG.byL = geno_ref%*%as.matrix(beta.list1)

    R2.byL = betaRho.byL^2 * median(sum_stats_target_test1$N) / colSums(betaG.byL^2) 
    print(R2.byL)
}

PRStr_main_check_pv<-function(ref_file, sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks, target_sumstats_train_file, target_sumstats_val_file, ps){ # nolint
	out1=0
	if (file.exists(paste0(ref_file,".bim")) & file.exists(paste0(ref_file,".bed")) & file.exists(paste0(ref_file,".fam"))){} else {out1="The ref file doesn't exist!"}

	if (!file.exists(sum_stats_file)){out1="The summary statistic file does not exist!"} else {
		temp=fread(sum_stats_file,nrow=1)
		if (ncol(temp)==4){
			if (sum(colnames(temp) %in% c("V1","V2","V3","V4"))==4){} else{
				if (sum(colnames(temp) %in% c("SNP","CHR","A1","Beta"))==4){} else {
					out1="The structure of sum_stats_file is wrong!"
				}
			} 
		} else {
			if (ncol(temp)>4){
				if (sum(colnames(temp) %in% c("SNP","CHR","A1","Beta"))<4){ 
					out1="The structure of sum_stats_file is wrong!"
				}
			} else {
				out1="The structure of sum_stats_file is wrong!"
			}
		} 
	}

    # need pseudo summ generation
    if (!is.null(target_sumstats_file)){
        if (!file.exists(target_sumstats_file)){out1="The target summary statistic file does not exist!"} else {
            temp=fread(target_sumstats_file,nrow=1)
            if (ncol(temp)==4 | ncol(temp)==5){
                if ("cor" %in% colnames(temp)){
                    if (sum(colnames(temp) %in% c("SNP", "A1", "cor", "N"))==4){} else{
                        out1="The structure of target_sum_stats_file is wrong!"
                    }                    
                }
                else {
                    if (sum(colnames(temp) %in% c("SNP", "A1", "beta", "p", "N"))==5){} else{
                        out1="The structure of target_sum_stats_file is wrong!"
                    }            
                }
            } else {
                if (ncol(temp)>5){
                    if (sum(colnames(temp) %in% c("SNP", "A1", "cor", "N"))<4 | sum(colnames(temp) %in% c("SNP", "A1", "beta", "p", "N"))<5){ 
                        out1="The structure of target_sum_stats_file is wrong!"
                    }
                } else {
                    out1="The structure of target_sum_stats_file is wrong!"
                }
            }
        }
        if (!is.logical(ps)) {out1='ps should be boolean'}
        if (ps) {
            if (!is.numeric(subprop) | subprop<0 | subprop>1) {out1='subprop should be a number between 0 and 1'}
            if (file.exists(paste0(ref_file_ps,".bim")) & file.exists(paste0(ref_file_ps,".bed")) & file.exists(paste0(ref_file_ps,".fam"))){} else {out1="The ref file for pseudo summary (ps) generation doesn't exist!"}
        }
    }

    # no pseudo summ generation
    else {
        if (!file.exists(target_sumstats_train_file)){out1="The target summary statistic train file does not exist!"} else {
        	temp=fread(target_sumstats_train_file,nrow=1)
            if (ncol(temp)==4 | ncol(temp)==5){
                if ("cor" %in% colnames(temp)){
                    if (sum(colnames(temp) %in% c("SNP", "A1", "cor", "N"))==4){} else{
                        out1="The structure of target_sumstats_train_file is wrong!"
                    }                    
                }
                else {
                    if (sum(colnames(temp) %in% c("SNP", "A1", "beta", "p", "N"))==5){} else{
                        out1="The structure of target_sumstats_train_file is wrong!"
                    }            
                }
            } else {
                if (ncol(temp)>5){
                    if (sum(colnames(temp) %in% c("SNP", "A1", "cor", "N"))<4 | sum(colnames(temp) %in% c("SNP", "A1", "beta", "p", "N"))<5){ 
                        out1="The structure of target_sumstats_train_file is wrong!"
                    }
                } else {
                    out1="The structure of target_sumstats_train_file is wrong!"
                }
            }
        }

        if (!file.exists(target_sumstats_val_file)){out1="The target summary statistic validation file does not exist!"} else {
        	temp=fread(target_sumstats_val_file,nrow=1)
            if (ncol(temp)==4 | ncol(temp)==5){
                if ("cor" %in% colnames(temp)){
                    if (sum(colnames(temp) %in% c("SNP", "A1", "cor", "N"))==4){} else{
                        out1="The structure of target_sumstats_val_file is wrong!"
                    }                    
                }
                else {
                    if (sum(colnames(temp) %in% c("SNP", "A1", "beta", "p", "N"))==5){} else{
                        out1="The structure of target_sumstats_val_file is wrong!"
                    }            
                }
            } else {
                if (ncol(temp)>5){
                    if (sum(colnames(temp) %in% c("SNP", "A1", "cor", "N"))<4 | sum(colnames(temp) %in% c("SNP", "A1", "beta", "p", "N"))<5){ 
                        out1="The structure of target_sumstats_val_file is wrong!"
                    }
                } else {
                    out1="The structure of target_sumstats_val_file is wrong!"
                }
            }
        }
    }

	if (!LDblocks %in% c("EUR.hg19", "AFR.hg19", "ASN.hg19")) {out1="The LDblocks name is wrong!"}
    return(out1)
}

PRS_tuning_pv_byLR <- function(beta.byL, betaRho.byL, betaG.byL, lr_list, N_val){
    # R2.byL <- c()
    # for (idx in 1:length(betaRho.byL)){
    #     R2 = betaRho.byL[idx]^2 / sum(betaG.byL[,idx]^2)
    #     R2.byL <- c(R2.byL, R2)
    # }

    R2.byL = (N_val * betaRho.byL^2) / as.vector(betaG.byL %*% t(betaG.byL)) # colSums(betaG.byL^2)

  flag=which(R2.byL==max(R2.byL))[1]
#   print(R2.byL)
#   print(flag)
  out.final=list()
  out.final$best.param=lr_list[flag] / dim(beta.byL)[1]
  out.final$best.beta = as.data.frame(beta.byL)[,c(1:3,9,8+flag)] #SNP, CHR, A1, Beta2, best.beta
  colnames(out.final$best.beta)[4:5] = c('base.beta', 'best.beta')
  out.final$R2.list = R2.byL

  return(out.final)
}

BedFileReader_prep <- function(plink_file, by_chr=FALSE, chr_list=NULL, plink_etc=''){
if (by_chr){
    # Create a list to store BedFileReader objects
    BedFileReaders <- list()

    # Loop through chromosomes
    for (chr in chr_list) {
    cat(sprintf("Processing chromosome: %d\n", chr))
    
    tic <- Sys.time() # Start timing
    
    # Assuming the file naming convention is consistent and correct
    train_file <- paste0(plink_file, chr, plink_etc)

    # Attempt to create a new BedFileReader object
    tryCatch({
        BedFileReaders[[chr]] <- new(BedFileReader, 
                                    paste0(train_file, ".fam"), 
                                    paste0(train_file, ".bim"), 
                                    paste0(train_file, ".bed"))
        # Assuming snp_index_func is a method of BedFileReader
        result <- BedFileReaders[[chr]]$snp_index_func()
    }, error = function(e) {
        cat("An error occurred: ", e$message, "\n")
    })
    
    toc <- Sys.time() # End timing
    cat(sprintf("Time taken for chromosome %d: %f seconds\n", chr, as.numeric(toc - tic)))
    }
} else{    
    # Attempt to create a new BedFileReader object
    tryCatch({
        BedFileReaders <- new(BedFileReader, 
                                    paste0(plink_file, ".fam"), 
                                    paste0(plink_file, ".bim"), 
                                    paste0(plink_file, ".bed"))
        # Assuming snp_index_func is a method of BedFileReader
        result <- BedFileReaders$snp_index_func()
    }, error = function(e) {
        cat("An error occurred: ", e$message, "\n")
    })
}

return(BedFileReaders)
}

#### for TL-PRS ####
##calculate nagelkerke R2 for binary phenotypes
nagelkerke<-function (fit, null = NULL, restrictNobs = FALSE){
    TOGGLE = (class(fit)[1] == "lm" | class(fit)[1] == "gls" |
        class(fit)[1] == "lme" | class(fit)[1] == "glm" | class(fit)[1] ==
        "negbin" | class(fit)[1] == "zeroinfl" | class(fit)[1] ==
        "clm" | class(fit)[1] == "vglm" | class(fit)[1] == "betareg" |
        class(fit)[1] == "rq")
    BOGGLE = (class(fit)[1] == "nls" | class(fit)[1] == "lmerMod" |
        class(fit)[1] == "glmerMod" | class(fit)[1] == "merModLmerTest" |
        class(fit)[1] == "lmerModLmerTest" | class(fit)[1] ==
        "clmm")
    SMOGGLE = (class(fit)[1] == "lmerMod" | class(fit)[1] ==
        "glmerMod" | class(fit)[1] == "merModLmerTest" | class(fit)[1] ==
        "lmerModLmerTest" | class(fit)[1] == "vglm")
    ZOGGLE = (class(fit)[1] == "zeroinfl")
    ZOGGLE2 = (class(fit)[1] == "rq")
    NOGGLE = is.null(null)
    ERROR = "Note: For models fit with REML, these statistics are based on refitting with ML"
    ERROR2 = "None"
    if (!restrictNobs & NOGGLE & TOGGLE) {
        null = update(fit, ~1)
    }
    if (restrictNobs & NOGGLE & TOGGLE) {
        null = update(fit, ~1, data = fit$model)
    }
    if (restrictNobs & !NOGGLE) {
        null = update(null, data = fit$model)
    }
    if (NOGGLE & BOGGLE) {
        ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"
    }
    if ((!TOGGLE) & (!BOGGLE)) {
        ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"
    }
    SMOGGLE2 = (class(null)[1] == "lmerMod" | class(null)[1] ==
        "glmerMod" | class(null)[1] == "merModLmerTest" | class(null)[1] ==
        "lmerModLmerTest" | class(null)[1] == "vglm")
    Y = matrix(rep(NA, 2), ncol = 1)
    colnames(Y) = ""
    rownames(Y) = c("Model:", "Null:")
    Z = matrix(rep(NA, 3), ncol = 1)
    colnames(Z) = c("Pseudo.R.squared")
    rownames(Z) = c("McFadden", "Cox and Snell (ML)", "Nagelkerke (Cragg and Uhler)")
    X = matrix(rep(NA, 4), ncol = 4)
    colnames(X) = c("Df.diff", "LogLik.diff", "Chisq", "p.value")
    rownames(X) = ""
    U = matrix(rep(NA, 2), ncol = 1)
    colnames(U) = ""
    rownames(U) = c("Model:", "Null:")
    if (TOGGLE | BOGGLE) {
        if (!SMOGGLE) {
            Y[1] = toString(fit$call)
        }
        if (SMOGGLE) {
            Y[1] = toString(fit@call)
        }
    }
    if (TOGGLE | (BOGGLE & !NOGGLE)) {
        if (!SMOGGLE2) {
            Y[2] = toString(null$call)
        }
        if (SMOGGLE2) {
            Y[2] = toString(null@call)
        }
        if (!ZOGGLE & !ZOGGLE2) {
            N = nobs(fit)
            U[1, 1] = nobs(fit)
            U[2, 1] = nobs(null)
        }
        if (!ZOGGLE & ZOGGLE2) {
            N = length(fit$y)
            U[1, 1] = length(fit$y)
            U[2, 1] = length(null$y)
        }
        if (ZOGGLE) {
            N = fit$n
            U[1, 1] = fit$n
            U[2, 1] = null$n
        }
        if (U[1, 1] != U[2, 1]) {
            ERROR2 = "WARNING: Fitted and null models have different numbers of observations"
        }
        m = suppressWarnings(logLik(fit, REML = FALSE))[1]
        n = suppressWarnings(logLik(null, REML = FALSE))[1]
        mf = 1 - m/n
        Z[1, ] = signif(mf, digits = 6)
        cs = 1 - exp(-2/N * (m - n))
        Z[2, ] = signif(cs, digits = 6)
        nk = cs/(1 - exp(2/N * n))
        Z[3, ] = signif(nk, digits = 6)
        o = n - m
        dfm = attr(logLik(fit), "df")
        dfn = attr(logLik(null), "df")
        if (class(fit)[1] == "vglm") {
            dfm = df.residual(fit)
        }
        if (class(fit)[1] == "vglm") {
            dfn = df.residual(null)
        }
        dff = dfn - dfm
        CHI = 2 * (m - n)
        P = pchisq(CHI, abs(dff), lower.tail = FALSE)
        X[1, 1] = dff
        X[1, 2] = signif(o, digits = 5)
        X[1, 3] = signif(CHI, digits = 5)
        X[1, 4] = signif(P, digits = 5)
    }
    W = ERROR
    WW = ERROR2
    V = list(Y, Z, X, U, W, WW)
    names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null",
        "Likelihood.ratio.test", "Number.of.observations", "Messages",
        "Warnings")
    return(V)
}

##Generate logitic regression results for binary phenotype
logistic_result_generator<-function(PRS,ped,Covar_name,Y_name){
        # datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
	datatemp2 = cbind(ped,PRS); colnames(datatemp2)[ncol(datatemp2)] = 'SCORESUM'
  	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
	if (args!=""){
        	modeltemp0=glm(paste0(Y_name,"~",args), data=datatemp2, family = "binomial")
        	modeltemp=glm(paste0(Y_name,"~ PRS +",args), data=datatemp2, family = "binomial")
	} else {
        	modeltemp0=glm(paste0(Y_name,"~1"), data=datatemp2, family = "binomial")
        	modeltemp=glm(paste0(Y_name,"~PRS"), data=datatemp2, family = "binomial")
	}
	r2=nagelkerke(modeltemp,null=modeltemp0)[[2]][3]
        return(r2)
}

##Generate linear regression results for continuous phenotype
linear_result_generator<-function(PRS,ped,Covar_name,Y_name){
	# datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
	datatemp2 = cbind(ped,PRS); colnames(datatemp2)[ncol(datatemp2)] = 'SCORESUM'
    datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))

	if (args!=""){
		modeltemp0=lm(paste0(Y_name,"~",args), data=datatemp2)
		modeltemp=lm(paste0(Y_name,"~ PRS +",args), data=datatemp2)
	} else{
		modeltemp0=lm(paste0(Y_name,"~1"), data=datatemp2)
		modeltemp=lm(paste0(Y_name,"~PRS"), data=datatemp2)
	}
	#prs.coef <- summary(modeltemp)$coeff[2,]
	#prs.beta <- as.numeric(prs.coef[1])
	#prs.se <- as.numeric(prs.coef[2])
	#prs.p <- as.numeric(prs.coef[4])
	model0.r2=summary(modeltemp0)$r.squared
	model1.r2=summary(modeltemp)$r.squared
	prs.r2 = summary(modeltemp)$r.squared-model0.r2
	#list.out=cbind(class=name1, prs.beta,prs.se,prs.p, model0.r2, model1.r2, prs.r2)
	print(confint(modeltemp))
    return(prs.r2)
}

# ped_val=ped_test; beta.info=best.beta[best.beta$CHR==i,1:3]; beta.all=best.beta[best.beta$CHR==i,4:5]; BedFileReader_new = BedFileReaders[[i]]
# i=1; ped_val=ped; beta.all = beta.all[beta.info$CHR==i,]; beta.info = beta.info[beta.info$CHR==i,]
Calculate_PRS_direct <- function(ped_val, beta.info, beta.all, BedFileReader_new, plink_file, by_chr, plink_etc=''){
   if (by_chr) {
    cat('calculate PRS for CHR',beta.info$CHR[1],'\n')

    # allele-flipping
    plink_bim = fread(paste0(plink_file, beta.info$CHR[1], plink_etc, '.bim')) # plink prefix + chrom number + .bim
   } else {
    plink_bim = fread(paste0(plink_file, '.bim'))
   }

    beta.info$order = 1:nrow(beta.info)
    bim_sum_stats=merge(plink_bim, beta.info,by.x="V2",by.y="SNP",sort=FALSE)
    bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

    flag=bim_sum_stats[bim_sum_stats$V6==bim_sum_stats$A1,'order'] # flipped. when ref=ref_ps, length(flag2)=0
    if (length(flag)>0){  beta.all[flag[[1]],]=-beta.all[flag[[1]],]}


    # G_val = readSomeSnp(beta.info$SNP, ped_val$fam_idx, BedFileReader=BedFileReader_new) 
    G_val = BedFileReader_new$readSomeSnp(beta.info$SNP, ped_val$fam_idx) 
    G_val <- do.call(cbind, G_val)

    PRS.all = as.matrix(G_val) %*% as.matrix(beta.all)
    return(PRS.all)
}

Calculate_PRS_direct <- function(ped_val, beta.info, beta.all, BedFileReader_new, plink_file, by_chr, plink_etc=''){
   if (by_chr) {
    cat('calculate PRS for CHR',beta.info$CHR[1],'\n')

    # allele-flipping
    plink_bim = fread(paste0(plink_file, beta.info$CHR[1], plink_etc, '.bim')) # plink prefix + chrom number + .bim
   } else {
    plink_bim = fread(paste0(plink_file, '.bim'))
   }

    beta.info$order = 1:nrow(beta.info)
    bim_sum_stats=merge(plink_bim, beta.info,by.x="V2",by.y="SNP",sort=FALSE)
    bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

    flag=bim_sum_stats[bim_sum_stats$V6==bim_sum_stats$A1,'order'] # flipped. when ref=ref_ps, length(flag2)=0
    if (length(flag)>0){  beta.all[flag[[1]],] = -beta.all[flag[[1]],]}

    # G_val = readSomeSnp(beta.info$SNP, ped_val$fam_idx, BedFileReader=BedFileReader_new) 
    # G_val = BedFileReader_new$readSomeSnp(beta.info$SNP, ped_val$fam_idx) 
    # G_val <- do.call(cbind, G_val)

    PRS.list = BedFileReader_new$calculatePRS_mat(beta.info$SNP, beta.all, ped_val$fam_idx) 

    PRS.all <- cbind(base.beta = PRS.list[[1]], best.beta = PRS.list[[2]])

    return(PRS.all)
}

splitgenome2<-function (CHR, POS, ref.CHR, ref.breaks, details = T, right = TRUE){
    CHR <- as.character(CHR)
    ref.CHR <- as.character(ref.CHR)
    POS <- as.integer(POS)
    ref.breaks <- as.integer(ref.breaks)
    stopifnot(all(!is.na(POS)) && all(!is.na(ref.breaks)))
    stopifnot(all(POS >= 1) && all(ref.breaks >= 1))
    stopifnot(length(CHR) == length(POS))
    stopifnot(length(ref.CHR) == length(ref.breaks))
    chr <- (unique(CHR))
    chr.ref <- (unique(ref.CHR))
    included.chr <- chr %in% chr.ref
    if (!all(included.chr))
        stop("Some chromosomes were not defined in the reference. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
    levels <- character(0)
    results <- character(length(POS))
    Details <- data.frame()
    for (C in chr.ref) {
        breaks <- sort(unique(ref.breaks[ref.CHR == C]))
        if (breaks[1] > 1)
            breaks <- c(1, breaks)
        if (breaks[length(breaks)] < Inf)
            breaks <- c(breaks, Inf)
        cut <- cut(POS[CHR == C], include.lowest = T, breaks = breaks,
            right = right)
        levels <- c(levels, paste0(C, "_", levels(cut)))
        cut <- paste0(C, "_", as.character(cut))
        results[CHR == C] <- cut
        if (details) {
            df <- data.frame(chr = C, start = breaks[-length(breaks)],
                end = breaks[-1])
            Details <- rbind(Details, df)
        }
    }
    results <- factor(results, levels = levels)
    if (details) {
        Details$counts <- as.integer(table(results))
        attr(results, "details") <- Details
    }
   out.item=list(results,Details)
   return(out.item)
}