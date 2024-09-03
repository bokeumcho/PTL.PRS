library(data.table)
library(lassosum)
library(Matrix)
library(ROCR)
library(parallel)

#### COMMON for TL-PRS & PTL-PRS ####
readSomeSnp <- function(snpList, sampleList=NULL, BedFileReader) { 
    # Preallocate a list to store each SNP's data
        geno_list <- vector("list", length(snpList))

        for (i in seq_along(snpList)) {
            snpname <- snpList[i]
            oneSnp <- BedFileReader$readOneSnp(BedFileReader$findSnpIndex(snpname))
            
            if (length(sampleList)!=0) {
                valid_samples <- sampleList[sampleList != -1 & sampleList <= length(oneSnp)]
                geno_list[[i]] <- oneSnp[valid_samples+1]
            }  else {
                geno_list[[i]] <- oneSnp
            }
            # if(i%%1000==0) {cat(i,'\n')}
        }
    
    # Convert the list to a dataframe
    geno_df <- do.call(rbind, geno_list)
    
    # print('start mean imputation \n')
    
    # handling NA (mean imputation)
    geno_df1 <- apply(geno_df, 1, function(x){
      #x[is.na(x)] <- mean(x, na.rm=TRUE)
      x[which(x==9)] <- mean(x[which(x!=9)])
      return(x)
    })
    return(geno_df1) # n*p
}

TL_PRS_test <- function(plink_file, by_chr, ped_test_file, best.beta, best.param, Ytype, Covar_name,Y_name){
    #' Test performance of PTL-PRS
    #' @export
    #' 
#     sourceCpp("/media/leelabsg-storage1/bokeum/TLPRS/source/BedFileReader.cpp")
#   setClass("BedFileReader", representation( pointer = "externalptr" ) )
#   BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
#   setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
#   setMethod("initialize","BedFileReader", function(.Object, ...) {
#     .Object@pointer <-.Call(BedFileReader_method("new"), ... )
#     .Object})

  ped_test = fread(ped_test_file,header=T, fill=TRUE)

    BedFileReaders <- BedFileReader_prep(plink_file)

    if (by_chr) {
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
            plink_file, by_chr
            )

        print(dim(partial_PRS))
        # Add the partial result to the cumulative total
        PRS.all <- PRS.all + partial_PRS  # Adjust this operation based on how results are combined
        rm(partial_PRS)
        gc()
    }
    } else {
        PRS.all <- Calculate_PRS_direct(
            ped_test, 
            best.beta[,1:3], 
            best.beta[,4:5],
            BedFileReaders
            )
    }

#   PRS.all = cbind(ped_test[,c('FID','IID')],PRS.all)
#   colnames(PRS.all)[3:4]=c("SCORESUM","SCORESUM")

  #write.table(PRS.all, file=paste0(tempfile,"_PRS_test.txt"),row.names=F,quote=F,col.names=T)

  if (Ytype=="C"){
      base_R2 = as.numeric(linear_result_generator(PRS.all[,1],ped_test,Covar_name,Y_name))
      model_R2 = as.numeric(linear_result_generator(PRS.all[,2],ped_test,Covar_name,Y_name))
    cat("baseline R2: ", base_R2, "\n model R2: ", model_R2, '\n best lr: ', best.param[[1]]/dim(best.beta)[1]) #, '\n best iter: ', best.params[2],'\n'
  } else {
      base_R2 = as.numeric(logistic_result_generator(PRS.all[,1],ped_test,Covar_name,Y_name))
      model_R2 = as.numeric(logistic_result_generator(PRS.all[,2],ped_test,Covar_name,Y_name))
      
      pred = prediction(PRS.all[,1], ped_test[,2])
      perf = performance(pred, measure='tpr', x.measure='fpr')
      auc = performance(pred, measure = "auc")
      base_AUC = auc@y.values[[1]]
      
      pred = prediction(PRS.all[,2], ped_test[,2])
      perf = performance(pred, measure='tpr', x.measure='fpr')
      auc = performance(pred, measure = "auc")
      model_AUC = auc@y.values[[1]]

    cat("baseline R2: ", base_R2, "\n model R2: ", model_R2, '\n base AUC: ', base_AUC, '\n model AUC: ', model_AUC, '\n best lr: ', best.param[[1]]/dim(best.beta)[1],'\n') #, '\n best iter: ', best.params[2],'\n'
    }

  gc()
}

BedFileReader_prep <- function(plink_file, by_chr=FALSE){
if (by_chr){
    # Create a list to store BedFileReader objects
    BedFileReaders <- list()

    # Loop through chromosomes
    for (chr in 1:22) {
    cat(sprintf("Processing chromosome: %d\n", chr))
    
    tic <- Sys.time() # Start timing
    
    # Assuming the file naming convention is consistent and correct
    train_file <- paste0(plink_file, chr)
    
    # Attempt to create a new BedFileReader object
    tryCatch({
        BedFileReaders[[chr]] <- new("BedFileReader", 
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
        BedFileReaders <- new("BedFileReader", 
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
Calculate_PRS_direct <- function(ped_val, beta.info, beta.all, BedFileReader_new, plink_file, by_chr){
   if (by_chr) {
    cat('calculate PRS for CHR',beta.info$CHR[1],'\n')

    # allele-flipping
    plink_bim = fread(paste0(plink_file, beta.info$CHR[1], '_v3.bim'))
   } else {
    plink_bim = fread(plink_file, '.bim')
   }

    beta.info$order = 1:nrow(beta.info)
    bim_sum_stats=merge(plink_bim, beta.info,by.x="V2",by.y="SNP",sort=FALSE)
    bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

    flag=bim_sum_stats[bim_sum_stats$V6==bim_sum_stats$A1,'order'][[1]] # flipped. when ref=ref_ps, length(flag2)=0
    if (length(flag)>0){  beta.all[flag,]=-beta.all[flag,]}

    tic('read ukb G')
        G_val = readSomeSnp(beta.info$SNP, ped_val$fam_idx, BedFileReader=BedFileReader_new) 
    toc()

    cat(dim(G_val), dim(beta.all),'\n')

    tic('calculate beta * G')
    PRS.all = as.matrix(G_val) %*% as.matrix(beta.all)
    toc()
    return(PRS.all)
}

PRStr_main_check<-function(ped_file,Covar_name,Y_name, Ytype,ref_file,sum_stats_file,target_sumstats_file, LDblocks, lr_list, iter){
	out1=0
	if (!file.exists(ped_file)){out1="The ped file does not exist!"} else {
		temp=fread(ped_file,header=T,nrow=1)
		if (!Y_name %in% colnames(temp)){out1="Y does not exist in the ped file!"}
		if (sum(!c("FID","IID") %in% colnames(temp))>0){out1="FID and IID do not exist in the ped file!"}
		
        if (sum(Covar_name=="")<1){
			if (sum(!Covar_name %in% colnames(temp))>0){out1="The covariates do not exist in the ped file!"}	
		}
	}
	if (!Ytype %in% c("C", "B")) {out1="The Y type is wrong! Now we only support continuous(C) and binary(B) outcomes."}			
	if (file.exists(paste0(ref_file,".bim")) & file.exists(paste0(ref_file,".bed")) & file.exists(paste0(ref_file,".fam"))){} else {out1="The ref file doesn't exist!"}
	if (!file.exists(sum_stats_file)){out1="The summary statistic file does not exist!"} else {
		temp=fread(sum_stats_file,nrow=1)
		if (ncol(temp)==3){
			if (sum(colnames(temp) %in% c("V1","V2","V3"))==3){} else{
				if (sum(colnames(temp) %in% c("SNP","A1","Beta"))==3){} else {
					out1="The structure of sum_stats_file is wrong!"
				}
			} 
		} else {
			if (ncol(temp)>3){
				if (sum(colnames(temp) %in% c("SNP","A1","Beta"))<3){ 
					out1="The structure of sum_stats_file is wrong!"
				}
			} else {
				out1="The structure of sum_stats_file is wrong!"
			}
		} 
	}
	if (!LDblocks %in% c("EUR.hg19", "AFR.hg19", "ASN.hg19")) {out1="The LDblocks name is wrong!"}
	if (lr_list)
    return(out1)
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

PRStr_tuning<-function(Beta.all, ped_val_file, tempfile, Covar_name,Y_name, Ytype, lr_list, iter, plink_file, BedFileReaders){
    tic('beta list prep')
    beta.info=Beta.all[,1:3]
    for (i in 9:ncol(Beta.all)){
      sdtemp=sd(Beta.all[,i],na.rm=T)
      if (sdtemp>1){
        Beta.all[,i:ncol(Beta.all)]=1;
      }
    }
    
    beta.all=Beta.all[,9:ncol(Beta.all)]/Beta.all$sd

    # beta.all <- copy(Beta.all)  # Make a shallow copy if you need to preserve Beta.all
    # setDT(beta.all)  # Convert to data.table, if it's not already

    # beta.all[, (10:ncol(beta.all)) := lapply(.SD, function(x) x / beta.all$sd), .SDcols = 10:ncol(beta.all)]
    # beta.all = beta.all[,10:ncol(beta.all)]
    rm(Beta.all)
    toc()

    ped = fread(ped_val_file,header=T, fill=TRUE)

    tic('calculating PRS')
        ## by CHRs
        PRS.all <- matrix(0, nrow=dim(ped)[1], ncol=dim(beta.all)[2])  # Adjust this initialization based on the actual structure of your result

        # Loop through each unique value in LDblocks2[[1]]
        for (i in unique(beta.info$CHR)) {
            # Calculate the partial result
            partial_PRS <- Calculate_PRS_direct(
                ped, 
                beta.info[beta.info$CHR==i,], 
                beta.all[beta.info$CHR==i,],
                BedFileReaders[[as.numeric(i)]],
                plink_file
            )

            print(dim(partial_PRS))
            # Add the partial result to the cumulative total
            PRS.all <- PRS.all + partial_PRS  # Adjust this operation based on how results are combined
        }

    toc()

        print('start calculating R-sqaured...')

    tic('calculating R squared')
        # system.time({
        out.item=data.frame("order"=NA,"R2"=NA, "auc"=NA)
        for (flag in 3:ncol(PRS.all)){
        out.item[flag-2,1]=flag-2
        resulttemp=PRS.all[,..flag] #c(1,2,..flag)
        colnames(resulttemp)[1]="SCORESUM"
        if (Ytype=="C"){
            out.item[flag-2,2] = as.numeric(linear_result_generator(resulttemp,ped,Covar_name,Y_name))
        } else {
            out.item[flag-2,2]=as.numeric(logistic_result_generator(resulttemp,ped,Covar_name,Y_name))
            pred = prediction(PRS.all[,..flag], ped[,2])
            # perf = performance(pred, measure='tpr', x.measure='fpr')
            auc = performance(pred, measure = "auc")
            out.item[flag-2,3] = auc@y.values[[1]]
        }
        }
        # })

        gc()

    toc()

    write.table(out.item[,1:3],file=paste0(tempfile,"_R_val.txt"),row.names=F,quote=F,col.names=T)
    # out.item <- fread(paste0(tempfile,"_R_val.txt"))

    flag=which(out.item$R2==max(out.item$R2))[1]
    param_table=data.frame("lr"=c(0,rep(1/nrow(beta.all)*lr_list,each=iter)),"iter"=c(0,rep(1:iter,length(lr_list))))   
    out.final=list()
    out.final$best.params=c(param_table$lr[flag], param_table$iter[flag])
    out.final$best.beta=cbind(beta.info,"beta"=beta.all[,c(1,..flag)]) 
    out.final$param_table=param_table 

    return(out.final)
}

PRStr_tuning_es <-function(Beta.all, ped_val_file, tempfile, Covar_name,Y_name, Ytype, lr_list, iter, plink_file, patience){
    tic('beta list prep')
    beta.info=Beta.all[,1:3]
    for (i in 10:ncol(Beta.all)){
      sdtemp=sd(Beta.all[,i],na.rm=T)
      if (sdtemp>1){
        Beta.all[,i:ncol(Beta.all)]=1;
      }
    }
    
    beta.all <- copy(Beta.all)  # Make a shallow copy if you need to preserve Beta.all
    setDT(beta.all)  # Convert to data.table, if it's not already

    beta.all[, (10:ncol(beta.all)) := lapply(.SD, function(x) x / Beta.all$sd), .SDcols = 10:ncol(beta.all)]
    beta.all = beta.all[,10:ncol(beta.all)]
    rm(Beta.all)
    toc()

    ped = fread(ped_val_file,header=T, fill=TRUE)
    
    prev_R2 = 0
    cnt=0
    R2_val = c()

    for (idx in 1:ncol(beta.all)){
        beta.vec = beta.all[,..idx]
        PRS.all.vec <- rep(0, dim(ped)[1]) 

        tic('calculating PRS')
        # Loop through each unique value by chrom
        for (i in unique(beta.info$CHR)) {
            # Calculate the partial result
            partial_PRS <- Calculate_PRS_direct(
                plink_file, ped, 
                beta.info[beta.info$CHR==i,], 
                beta.vec[beta.info$CHR==i]
            )

            print(dim(partial_PRS))
            # Add the partial result to the cumulative total
            PRS.all.vec <- PRS.all.vec + partial_PRS  # Adjust this operation based on how results are combined
        }
        toc()
        
        # delete later
        if (colnames(ped)[1] == 'FID') {ped$IID = ped$FID
        } else ped$FID = ped$IID

        PRS.all = cbind(ped[,c('FID','IID')],PRS.all.vec)

        #write.table(PRS.all, file=paste0(tempfile,"_PRS_",group,".txt"),row.names=F,quote=F,col.names=T)
        print('start calculating R-sqaured...')

        tic('calculating R squared')
        colnames(PRS.all)[3]="SCORESUM"
        if (Ytype=="C"){
            curr_R2 = as.numeric(linear_result_generator(PRS.all,ped,Covar_name,Y_name))
        } else {
            curr_R2 = as.numeric(logistic_result_generator(PRS.all,ped,Covar_name,Y_name))
        }

        R2_val = c(R2_val, curr_R2)
        if (curr_R2 < prev_R2) {
            cnt = cnt + 1
        }

        if (cnt > patience) {
            write.table(R2_val,file=paste0(tempfile,"_R_val.txt"),row.names=F,quote=F,col.names=T)
            
            param_table=data.frame("lr"=c(0,rep(1/length(beta.vec)*lr_list,each=iter)),"iter"=c(0,rep(1:iter,length(lr_list))))   
            out.final=list()
            out.final$best.params=c(param_table$lr[idx], param_table$iter[idx])
            out.final$best.beta=cbind(beta.info,"beta"=cbind(beta.all[,1],beta.vec))
            out.final$param_table=param_table 
            return (out.final)
        }
        prev_R2 = curr_R2
        gc()
        toc()
    }
}

PRS_tuning_byLR <- function(beta.byL, PRS.byL,ped_val_file,Covar_name,Y_name, Ytype, lr_list){
    ped_val = fread(ped_val_file, header=T, fill=TRUE)

    R2.byL <- c()
  for (idx in 1:ncol(PRS.byL)){
      if (Ytype=="C"){
          R2 = as.numeric(linear_result_generator(PRS.byL[,idx],ped_val,Covar_name,Y_name))
      } else {
          R2 = as.numeric(logistic_result_generator(PRS.byL[,idx],ped_val,Covar_name,Y_name))
      }
      R2.byL <- c(R2.byL, R2)
  }

  flag=which(R2.byL==max(R2.byL))[1]
  out.final=list()
  out.final$best.param=lr_list[flag]
  out.final$best.beta = as.data.frame(beta.byL)[,c(1:3,9,9+flag)] #SNP, CHR, A1, Beta2, best.beta
  out.final$R2.list = R2.byL

  return(out.final)
}

readG_byCHR <- function(ped_val, geno_info2, BedFileReaders){
    G_val = matrix(,nrow=nrow(ped_val)) #initialize

    for (chr in unique(geno_info2$V1)){
        G_temp = readSomeSnp(geno_info2$V2, ped_val$fam_idx, BedFileReaders[[as.numeric(chr)]]) 
        G_val = cbind(G_val, G_temp) # n*p
    }
    G_val = G_val[,-1]

    return(G_val)
}

#### for PTL-PRS ####
PRStr_main_check_pv<-function(ref_file, sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks, target_sumstats_train_file, target_sumstats_val_file, ps){
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

PRStr_tuning_pv<-function(Beta.all, sum_stats_target_train, sum_stats_target_val, tempfile, lr_list, iter, BedFileReader){
    tic('beta list prep')
    beta.info=Beta.all[,1:3]
    for (i in 9:ncol(Beta.all)){
      sdtemp=sd(Beta.all[,i],na.rm=T)
      if (sdtemp>1){
        Beta.all[,i:ncol(Beta.all)]=1;##print("fdsf")
      }
    }
    
    beta.all=Beta.all[,9:ncol(Beta.all)]/Beta.all$sd

    # beta.all <- copy(Beta.all)  # Make a shallow copy if you need to preserve Beta.all
    # setDT(beta.all)  # Convert to data.table, if it's not already

    # beta.all[, (10:ncol(beta.all)) := lapply(.SD, function(x) x / Beta.all$sd), .SDcols = 10:ncol(beta.all)]
    # beta.all = beta.all[,10:ncol(beta.all)]
    rm(Beta.all)

    toc()

    tic('reading ref G')
    geno_ref = readSomeSnp(beta.info$SNP, BedFileReader = BedFileReader)
    toc()

    # GG = cor(as.matrix(geno_ref))
    flag=which(apply(geno_ref, 2, sd, na.rm = TRUE)!=0) # all nonzero sd

    geno_ref = scale(geno_ref[,flag])
    
    sum_stats_target_train1 = sum_stats_target_train[sum_stats_target_train$SNP %in% beta.info$SNP,]
    sum_stats_target_val1 = sum_stats_target_val[sum_stats_target_val$SNP %in% beta.info$SNP,]

    #out.item=data.frame("order"=NA)
    out.item = as.data.frame(matrix(NA, nrow=ncol(beta.all), ncol=5))
    colnames(out.item) = c('order','R_num_p_tr','R_denom_p_tr',"R_num_p_val","R_denom_p_val")

    tic('calculating pseudo R')
    # system.time({
    for (flag in 1:ncol(beta.all)){ #594:ncol(beta.all)
        # order
        out.item[flag,1]=flag 
        
        # calculate R for each beta
        out.item[flag,2:3]=Calculate_R(beta.all[,flag], sum_stats_target_train1, geno_ref)
        out.item[flag,4:5]=Calculate_R(beta.all[,flag], sum_stats_target_val1, geno_ref)
        }
    # })
    
    out.item$R2_p_tr = out.item$R_num_p_tr/out.item$R_denom_p_tr
    out.item$R2_p_val = out.item$R_num_p_val/out.item$R_denom_p_val
    
    toc()

    write.table(out.item,file=paste0(tempfile,"_R.txt"),row.names=F,quote=F,col.names=T)
    flag=which(out.item$R2_p_val==max(out.item$R2_p_val))[1]

    param_table=data.frame("lr"=c(0,rep(1/nrow(beta.all)*lr_list,each=iter)),"iter"=c(0,rep(1:iter,length(lr_list))))   
    out.final=list()
    out.final$best.params = c(param_table$lr[flag],param_table$iter[flag])
    out.final$best.beta=cbind(beta.info,"beta"=beta.all[,c(1,flag)]) 
    out.final$param_table=param_table 

    return(out.final)
}

PRStr_tuning_pv_es<-function(Beta.all, sum_stats_target_val, tempfile, lr_list, iter, patience){
    tic('beta list prep')
    beta.info=Beta.all[,1:3]
    for (i in 10:ncol(Beta.all)){
      sdtemp=sd(Beta.all[,i],na.rm=T)
      if (sdtemp>1){
        Beta.all[,i:ncol(Beta.all)]=1;##print("fdsf")
      }
    }
    
    beta.all <- copy(Beta.all)  # Make a shallow copy if you need to preserve Beta.all
    setDT(beta.all)  # Convert to data.table, if it's not already

    beta.all[, (10:ncol(beta.all)) := lapply(.SD, function(x) x / Beta.all$sd), .SDcols = 10:ncol(beta.all)]
    beta.all = beta.all[,10:ncol(beta.all)]
    rm(Beta.all)

    toc()

    tic('reading ref G')
    geno_ref = readSomeSnp(beta.info$SNP)
    toc()

    prev_R2 = 0
    cnt=0
    R2_val = c()
    
    # loop through models
    for (idx in 1:ncol(beta.all)){ 
        beta.vec = beta.all[,..idx]
        tic('calculating pseudo R')
        R_result =Calculate_R(beta.vec, sum_stats_target_val, geno_ref)
        curr_R2 = R_result[[1]]^2/R_result[[2]]
        toc()

        R2_val = c(R2_val, curr_R2)
        
        if (curr_R2 < prev_R2) {
            cnt = cnt + 1
        }

        if (cnt > patience) {
            write.table(R2_val,file=paste0(tempfile,"_pseudoR_val.txt"),row.names=F,quote=F,col.names=T)
            
            param_table=data.frame("lr"=c(0,rep(1/length(beta.vec)*lr_list,each=iter)),"iter"=c(0,rep(1:iter,length(lr_list))))   
            out.final=list()
            out.final$best.params=c(param_table$lr[idx], param_table$iter[idx])
            out.final$best.beta=cbind(beta.info,"beta"=cbind(beta.all[,1],beta.vec))
            out.final$param_table=param_table 

            return (out.final)
        }
        prev_R2 = curr_R2
        gc()
    }
}

PRS_tuning_pv_byLR <- function(beta.byL, betaRho.byL, betaG.byL, lr_list){
    # R2.byL <- c()
    # for (idx in 1:length(betaRho.byL)){
    #     R2 = betaRho.byL[idx]^2 / sum(betaG.byL[,idx]^2)
    #     R2.byL <- c(R2.byL, R2)
    # }

    R2.byL = betaRho.byL^2 / colSums(betaG.byL^2)

  flag=which(R2.byL==max(R2.byL))[1]
  out.final=list()
  out.final$best.param=lr_list[flag] / dim(beta.byL)[1]
  out.final$best.beta = as.data.frame(beta.byL)[,c(1:3,9,9+flag)] #SNP, CHR, A1, Beta2, best.beta
  out.final$R2.list = R2.byL

  return(out.final)
}