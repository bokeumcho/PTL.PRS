#' @useDynLib PTL.PRS.c, .registration = TRUE
#' @import methods Rcpp
"_PACKAGE"

# updateCoeff <- function(beta.all, GG2, geno_ref, learningrate, maxiter, betaSD, sum_stats_target_val_block, patience, trace){
#   #' Update betas in each LD block for each learning rate, using block-wise early stopping
#   #' 
#   #' @param beta.all A list of two elements: baseline betas to update for each LD block (betatemp), and their gradients (gy-GG2%*%betatemp)
#   #' @param GG2 A covariance matrix of reference genotypes (p'*p')
#   #' @param geno_ref A genotype matrix of reference panel (n*p')
#   #' @param learningrate Learning rate
#   #' @param maxiter A maximum number of iterations
#   #' @param betaSD Standard deviation of each SNP's beta
#   #' @param sum_stats_target_val_block A matrix of target summary statistics for each LD block 
#   #' @param trace Logical; if TRUE, print stopped iteration
#   #' @returns Updated betas (betatemp), numerator (R_num) and a list of PRS to calculate denominator (PRS) for validation pseudo-R
  
#     k=1

#     if (is.null(dim(beta.all))){
#       betatemp=beta.all[2]
#       u0=beta.all[1]
#     } else{
#       betatemp=beta.all[,2]
#       u0=beta.all[,1]
#     }

#     prev_R2 = 0
#     cnt=0

#     while (k<=maxiter){
#       # if (!update_new){
#       #   for (j in 1:length(betatemp)){
#       #       beta_old=betatemp[j]
#       #       betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
#       #       u0=u0-GG2[,j]*(betatemp[j]-beta_old)
#       #   }
#       # } else {
#     # diff update
#       for (j in 1:length(betatemp)){
#           beta_old=betatemp
#           betatemp[j]=(learningrate*u0[j]+beta_old[j])/ 1
#           u0[j] = u0[j] - GG2[,j] %*% (betatemp - beta_old)
#       }
#     # }

#     betatemp1 = betatemp/betaSD
#     # n_test = sum_stats_target_val_block$N[1]
#     R_num = sum(betatemp1*sum_stats_target_val_block$cor)
#     PRS = geno_ref%*%betatemp1
    
#     curr_R2 = R_num^2/sum(PRS^2)

#     if (curr_R2 < prev_R2) {
#         cnt = cnt + 1
#     }

#     if (cnt > patience) {
#         break
#     }
#     prev_R2 = curr_R2

#     k=k+1
#     } 
    
#     if(trace) cat("stopped iteration:", k, "\n") 

#     return(list(betatemp1, R_num, PRS))
# }

# #num=1; i=unique(LDblocks2[[1]])[num]; cor=bim_sum_stats[which(LDblocks2[[1]]==i),]; nsnp=nrow(bim_sum_stats);
# block_calculation_pv_es<-function(cor,num,nsnp,temp.file, lr_list, maxiter, sum_stats_target_val, patience, trace, BedFileReader){
#   #' Update betas in each LD block for a list of learning rates using block-wise early stopping
#   #' 
#   #' @param cor Source summary statistics to update for each LD block
#   #' @param num LD block number
#   #' @param nsnp A number of total SNPs used in analysis
#   #' @param temp.file Output path to save intermediate results
#   #' @param lr_list A list of learning rates
#   #' @param maxiter A maximum number of iterations
#   #' @param sum_stats_target_val A matrix of target summary statistics for validation 
#   #' @param patience A number of times to tolerate decrease of validation pseudo-R
#   #' @param trace Logical; if TRUE, print stopped iteration
#   #' @param BedFileReader A BedFileReader object to read reference Bedfile
#   #' @returns Information of updated betas (geno_info2), a list of numerators (beta_rho) and a matrix of PRS to calculate denominator (beta_g) for validation pseudo-R with all learning rates
  
#   # Gtemp = try(readSomeSnp(cor$V2, BedFileReader = BedFileReader),silent=T) #ref
  # Gtemp = try(BedFileReader$readSomeSnp(cor$V2, sampleList=integer(0)),silent=T) #ref
  # Gtemp <- do.call(cbind, Gtemp)

  # cat('read Gtemp block ',num,'\n')

  # if ("try-error" %in% class(Gtemp)){
  #   print('error while reading Bedfile')
  #   return(NULL)
  # }else{
  #   GG = cor(as.matrix(Gtemp)) 
  #   geno_info = as.data.frame(cor[,c('V2','V5')]); colnames(geno_info)=c("SNP","A1")
  #   geno_info$mean=colMeans(as.matrix(Gtemp),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))
  # }

  # list1=which(geno_info$sd==0)
  # if (length(list1)>0){
  #   geno_info=geno_info[-list1,]
  #   GG=GG[-list1,-list1]
  #   Gtemp = Gtemp[,-list1]
  # }
  # if (nrow(geno_info)==0){
  #   return(NULL)
  # } else {
  #   geno_info$order=1:nrow(geno_info)

  #   # if no sd=0 cols, same as cbind
  #   geno_info2=merge(cor[,c("V1","V2","V5","Beta2","cor")],geno_info, by.x="V2",by.y="SNP",sort=F)

  #   GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
  #   gy=geno_info2$cor
  #   betatemp=geno_info2$Beta2*geno_info2$sd
  #   u0=gy-GG2%*%betatemp
  #   beta.init=cbind(u0, betatemp)

  #   ###
  #   Gtemp3 = scale(Gtemp) 

  #   sum_stats_target_val_block = sum_stats_target_val[sum_stats_target_val$SNP %in% geno_info2$V2,]

    # beta.all = c(geno_info2$Beta2)
    # beta_rho = c(sum(geno_info2$Beta2 * sum_stats_target_val_block$cor))
    # beta_g = c(Gtemp3 %*% geno_info2$Beta2)
    # # R2.all = c() # beta_rho^2/sum(beta_g^2)

    # for (factor1 in lr_list){
    #   learningrate=1/nsnp*factor1
    #   if (learningrate>1){learningrate=1}
    #   results = updateCoeff(beta.init, GG2, Gtemp3, learningrate, maxiter, geno_info2$sd, sum_stats_target_val_block, patience, trace)
    #   beta.all=cbind(beta.all,results[[1]])
    #   beta_rho = c(beta_rho, results[[2]])
    #   beta_g = cbind(beta_g, results[[3]])
    # }
    
    # cols_to_drop <- match(c("A1", "order"), colnames(geno_info2))
    # geno_info2 <- cbind(geno_info2[, -cols_to_drop, drop = FALSE], as.data.frame(beta.all))

    # # geno_info2=cbind(geno_info2[,-c('A1','order')],as.data.frame(beta.all)) # delete A1, order
    
    # rm(beta.init, beta.all, Gtemp, Gtemp3, betatemp)
    # gc()

    # return(list(geno_info2, beta_rho, beta_g))
#   } 
# }##function end

# PRStr_calculation_pv_es<-function(sum_stats_target_train, ref_file, LDblocks, num_cores=1,temp.file, lr_list, maxiter, sum_stats_target_val,patience, trace){
#   #' Update betas for all LD blocks for a list of learning rates using block-wise early stopping
#   #' 
#   #' @param sum_stats_target_train Source summary statistics to update
#   #' @param ref_file Prefix of PLINK file of the target reference panel
#   #' @param LDblocks Character string specifying the LD blocks, defaults to "EUR.hg19".
#   #' @param num_cores A number of cores to be used in parallel computation for training, defaults to 1.
#   #' @param temp.file Output path to save intermediate results
#   #' @param lr_list A list of learning rates
#   #' @param maxiter A maximum number of iterations
#   #' @param sum_stats_target_val A matrix of target summary statistics for validation
#   #' @param patience A number of times to tolerate decrease of validation pseudo-R
#   #' @param trace Logical; if TRUE, print stopped iteration
#   #' @returns Information of updated betas (beta.byL), a list of numerators (betaRho.byL) and a matrix of PRS to calculate denominator (betaG.byL) for validation pseudo-R with all learning rates for all LD blocks

#   possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
#                          "EUR.hg38", "AFR.hg38", "ASN.hg38") 
#   if(!is.null(LDblocks)) {
#     if(is.character(LDblocks) && length(LDblocks) == 1) {
#       if(LDblocks %in% possible.LDblocks) {
#         LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
#       } else {
#         stop(paste("I cannot recognize this LDblock. Specify one of", 
#                    paste(possible.LDblocks, collapse=", ")))
#       }
#     }
#     if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
#     if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
#       if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
#         LDblocks <- as.data.frame(LDblocks)
#         stopifnot(ncol(LDblocks) == 3)
#         stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
#         LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
#       }
#   } else {
#     stop(paste0("LDblocks must be specified. Specify one of ", 
#                 paste(possible.LDblocks, collapse=", "), 
#                 ". Alternatively, give an integer vector defining the blocks, ", 
#                 "or a .bed file with three columns read as a data.frame."))
#   }
  
#   ref.bim <- fread(paste0(ref_file, ".bim"))
#   ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
#   ref.bim <- ref.bim[ref.bim$V1 %in% 1:22,]

#   ref.bim$order=1:nrow(ref.bim)
#   bim_sum_stats=merge(ref.bim, sum_stats_target_train,by.x="V2",by.y="SNP",sort=FALSE)
#   bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
  
# # remove duplicates
#   bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

#   bim_sum_stats$Beta2=NA
#   flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
#   if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
#   flag2=which(bim_sum_stats$V6==bim_sum_stats$A1) # flipped. when ref=ref_ps, length(flag2)=0
#   if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2]; bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
  
#   bim_sum_stats=bim_sum_stats[which(!is.na(bim_sum_stats$Beta2)),c("V2","V1","V5","V6","order","Beta2","cor")] 
  
#   ## match ref with sum_stats_target_val & check allele flip 
#   bim_sum_stats_val=merge(ref.bim, sum_stats_target_val, by.x="V2",by.y="SNP",sort=FALSE)
#   bim_sum_stats_val=bim_sum_stats_val[order(bim_sum_stats_val$order),]
  
#   bim_sum_stats_val = bim_sum_stats_val[!duplicated(bim_sum_stats_val$V2),]

#   flag3=which(bim_sum_stats_val$V6==bim_sum_stats_val$A1) # flipped. when ref=ref_ps, length(flag3)=0
#   if (length(flag3)>0){  bim_sum_stats_val$cor[flag3]=-bim_sum_stats_val$cor[flag3];}
  
#   # cat('PRStr_calculation_pv_es flag1,2,3: ',length(flag1),length(flag2),length(flag3),'\n')
  
#   bim_sum_stats_val=bim_sum_stats_val[,c("V2","V5","cor","N")] 
#   colnames(bim_sum_stats_val)[1:2] = c('SNP','A1')

#   ref.extract <- rep(FALSE, nrow(ref.bim))
#   ref.extract[bim_sum_stats$order] <- TRUE
  
#   if(!is.null(LDblocks)) {
#       LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
#                               POS = ref.bim$V4[ ref.extract],
#                               ref.CHR = LDblocks[,1], 
#                               ref.breaks = LDblocks[,3])
#       # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
#   } 
  
#   # assign BedFileReader object
#   BedFileReader_ref <- new( BedFileReader, paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))
#   result = try(BedFileReader_ref$snp_index_func(), silent=TRUE)

#   unique_blocks <- unique(LDblocks2[[1]])
#   blocks <- vector("list", length(unique_blocks))
#   nsnp <- nrow(bim_sum_stats)  

#   # Create a function that processes one block.
#   process_block <- function(cur_block) {
#     # Filter rows for current block.
#     cor_df <- bim_sum_stats[ LDblocks2[[1]] == cur_block, , drop = FALSE ]
    
#     # Read SNPs for current block. (If BedFileReader_ref is not thread-safe,
#     # consider creating a new instance here.)
#     Gtemp <- try(BedFileReader_ref$readSomeSnp(cor_df$V2, sampleList = integer(0)), silent = TRUE)
#     if (inherits(Gtemp, "try-error")) {
#       cat("Error while reading Bedfile in block", cur_block, "\n")
#       return(NULL)
#     }
#     Gtemp <- do.call(cbind, Gtemp)
#     cat("read Gtemp block", cur_block, "\n")
    
#     # Compute correlation matrix on the genotype matrix
#     GG <- cor(as.matrix(Gtemp))
    
#     # Prepare geno_info from cor_df
#     geno_info <- as.data.frame(cor_df[, c("V2", "V5"), drop = FALSE])
#     colnames(geno_info) <- c("SNP", "A1")
#     geno_info$mean <- colMeans(as.matrix(Gtemp), na.rm = TRUE)
#     geno_info$maf <- geno_info$mean / 2
#     geno_info$sd <- sqrt(2 * geno_info$maf * (1 - geno_info$maf))
    
#     # Remove columns with zero standard deviation.
#     nonzero_sd <- geno_info$sd != 0
#     if (!all(nonzero_sd)) {
#       geno_info <- geno_info[nonzero_sd, , drop = FALSE]
#       GG <- GG[nonzero_sd, nonzero_sd, drop = FALSE]
#       Gtemp <- Gtemp[, nonzero_sd, drop = FALSE]
#     }
#     if (nrow(geno_info) == 0) return(NULL)
    
#     geno_info$order <- seq_len(nrow(geno_info))
    
#     # Merge with additional info from cor_df.
#     geno_info2 <- merge(cor_df[, c("V1", "V2", "V5", "Beta2", "cor")],
#                         geno_info,
#                         by.x = "V2", by.y = "SNP", sort = FALSE)
    
#     # Order matrices according to geno_info2$order.
#     GG2 <- as.matrix(GG[geno_info2$order, geno_info2$order])
#     gy <- geno_info2$cor
#     betatemp <- geno_info2$Beta2 * geno_info2$sd
#     u0 <- gy - GG2 %*% betatemp
#     beta.init <- cbind(u0, betatemp)
    
#     # Standardize genotype matrix.
#     Gtemp3 <- scale(Gtemp)
    
#     # Subset sum_stats_target_val for SNPs in this block.
#     sum_stats_target_val_block <- sum_stats_target_val[sum_stats_target_val$SNP %in% geno_info2$V2, , drop = FALSE]
    
#     # Return a list with all prepared parameters.
#     list(
#       beta_all = beta.init,             # matrix with two columns (u0 and betatemp)
#       GG2 = GG2,                        # covariance matrix for the block
#       geno_ref = Gtemp3,                # standardized genotype matrix for the block
#       lr_list = lr_list / nsnp,          # adjust learning rate by nsnp
#       maxiter = maxiter,
#       sum_stats_target_val_cor = sum_stats_target_val_block$cor,
#       patience = patience,
#       trace = trace,
#       geno_info2 = geno_info2
#     )
#   }

#   # Use mclapply to process blocks in parallel.
#   # Note: adjust mc.cores according to your available cores.
#   blocks <- mclapply(unique_blocks, process_block, mc.cores = num_cores)

#   # Remove any NULL results (e.g. if an error occurred in some block).
#   # blocks <- Filter(Negate(is.null), blocks)

#   # Proceed with the RcppParallel part on the prepared blocks.
#   results.list <- block_calculation_parallel(blocks)

#   beta.byL <- do.call(rbind, lapply(results.list, function(x) x[[1]]))
#   cols_to_drop <- match(c("A1", "order"), colnames(beta.byL))
#   beta.byL <- beta.byL[, -cols_to_drop, drop = FALSE]
#   colnames(beta.byL)[1:3] <- c("SNP", "CHR", "A1")

#   betaRho.byL <- Reduce("+", lapply(results.list, function(x) x[[2]]))
#   betaG.byL <- Reduce("+", lapply(results.list, function(x) as.numeric(unlist(x[[3]]))))
#   if (is.null(dim(betaG.byL))) {
#     betaG.byL <- matrix(betaG.byL, nrow = 1)
#   }

#   return(list(beta.byL, betaRho.byL, betaG.byL))
# }

PRStr_calculation_pv_es <- function(sum_stats_target_train, ref_file, LDblocks, num_cores = 1,
                                      temp.file, lr_list, maxiter, sum_stats_target_val,
                                      patience, trace) {
  library(pryr)
  
  # cat("Memory at start:", mem_used(), "\n")
  
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),
                                                   package = "lassosum"), header = TRUE)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of", 
                   paste(possible.LDblocks, collapse = ", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[, 3] >= LDblocks[, 2]))
        LDblocks[, 1] <- as.character(sub("^chr", "", LDblocks[, 1], ignore.case = TRUE))
      }
  } else {
    stop(paste0("LDblocks must be specified. Specify one of ", 
                paste(possible.LDblocks, collapse = ", "), 
                ". Alternatively, give an integer vector defining the blocks, ", 
                "or a .bed file with three columns read as a data.frame."))
  }
  
  ref.bim <- data.table::fread(paste0(ref_file, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = TRUE))
  ref.bim <- ref.bim[ref.bim$V1 %in% 1:22, ]
  ref.bim$order <- 1:nrow(ref.bim)
  
  bim_sum_stats <- merge(ref.bim, sum_stats_target_train, by.x = "V2", by.y = "SNP", sort = FALSE)
  bim_sum_stats <- bim_sum_stats[order(bim_sum_stats$order), ]
  bim_sum_stats <- bim_sum_stats[!duplicated(bim_sum_stats$V2), ]
  
  bim_sum_stats$Beta2 <- NA
  flag1 <- which(bim_sum_stats$V5 == bim_sum_stats$A1)
  if (length(flag1) > 0) {
    bim_sum_stats$Beta2[flag1] <- bim_sum_stats$Beta[flag1]
  }
  flag2 <- which(bim_sum_stats$V6 == bim_sum_stats$A1)
  if (length(flag2) > 0) {
    bim_sum_stats$Beta2[flag2] <- -bim_sum_stats$Beta[flag2]
    bim_sum_stats$cor[flag2] <- -bim_sum_stats$cor[flag2]
  }
  
  bim_sum_stats <- bim_sum_stats[which(!is.na(bim_sum_stats$Beta2)), c("V2", "V1", "V5", "V6", "order", "Beta2", "cor")]
  
  bim_sum_stats_val <- merge(ref.bim, sum_stats_target_val, by.x = "V2", by.y = "SNP", sort = FALSE)
  bim_sum_stats_val <- bim_sum_stats_val[order(bim_sum_stats_val$order), ]
  bim_sum_stats_val <- bim_sum_stats_val[!duplicated(bim_sum_stats_val$V2), ]
  flag3 <- which(bim_sum_stats_val$V6 == bim_sum_stats_val$A1)
  if (length(flag3) > 0) {
    bim_sum_stats_val$cor[flag3] <- -bim_sum_stats_val$cor[flag3]
  }
  
  bim_sum_stats_val <- bim_sum_stats_val[, c("V2", "V5", "cor", "N")]
  colnames(bim_sum_stats_val)[1:2] <- c("SNP", "A1")
  
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  if(!is.null(LDblocks)) {
    LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ref.extract], 
                              POS = ref.bim$V4[ref.extract],
                              ref.CHR = LDblocks[, 1], 
                              ref.breaks = LDblocks[, 3])
  }
  
  BedFileReader_ref <- new(BedFileReader, paste0(ref_file, ".fam"),
                             paste0(ref_file, ".bim"), paste0(ref_file, ".bed"))
  result <- try(BedFileReader_ref$snp_index_func(), silent = TRUE)
  
  unique_blocks <- unique(LDblocks2[[1]])
  blocks <- vector("list", length(unique_blocks))
  nsnp <- nrow(bim_sum_stats)
  
  process_block <- function(cur_block) {
    cor_df <- bim_sum_stats[LDblocks2[[1]] == cur_block, , drop = FALSE]
    Gtemp <- try(BedFileReader_ref$readSomeSnp(cor_df$V2, sampleList = integer(0)), silent = TRUE)
    if (inherits(Gtemp, "try-error")) {
      cat("Error while reading Bedfile in block", cur_block, "\n")
      return(NULL)
    }
    Gtemp <- do.call(cbind, Gtemp)
    # cat("read Gtemp block", cur_block, "\n")
    GG <- cor(as.matrix(Gtemp))
    geno_info <- as.data.frame(cor_df[, c("V2", "V5"), drop = FALSE])
    colnames(geno_info) <- c("SNP", "A1")
    geno_info$mean <- colMeans(as.matrix(Gtemp), na.rm = TRUE)
    geno_info$maf <- geno_info$mean / 2
    geno_info$sd <- sqrt(2 * geno_info$maf * (1 - geno_info$maf))
    
    nonzero_sd <- geno_info$sd != 0
    if (!all(nonzero_sd)) {
      geno_info <- geno_info[nonzero_sd, , drop = FALSE]
      GG <- GG[nonzero_sd, nonzero_sd, drop = FALSE]
      Gtemp <- Gtemp[, nonzero_sd, drop = FALSE]
    }
    if (nrow(geno_info) == 0) return(NULL)
    
    geno_info$order <- seq_len(nrow(geno_info))
    geno_info2 <- merge(cor_df[, c("V1", "V2", "V5", "Beta2", "cor")],
                        geno_info,
                        by.x = "V2", by.y = "SNP", sort = FALSE)
    GG2 <- as.matrix(GG[geno_info2$order, geno_info2$order])
    gy <- geno_info2$cor
    betatemp <- geno_info2$Beta2 * geno_info2$sd
    u0 <- gy - GG2 %*% betatemp
    beta.init <- cbind(u0, betatemp)
    
    Gtemp3 <- scale(Gtemp)
    sum_stats_target_val_block <- sum_stats_target_val[sum_stats_target_val$SNP %in% geno_info2$V2, , drop = FALSE]
    
    list(
      beta_all = beta.init,
      GG2 = GG2,
      geno_ref = Gtemp3,
      lr_list = lr_list / nsnp,
      maxiter = maxiter,
      sum_stats_target_val_cor = sum_stats_target_val_block$cor,
      patience = patience,
      trace = trace,
      geno_info2 = geno_info2
    )
  }
  
  # # generate queue
  # queue <- create_input_queue()

  # for (cur_block in unique_blocks) {
  #   input <- process_block(cur_block)
  #   if (!is.null(input)) {
  #     push_input(queue, input)
  #   }
  # }

  # finish_queue(queue)

  # results.list <- block_calculation_parallel_streamed(queue, n_threads = num_cores)
# Split blocks into subgroups
block_chunks <- split(unique_blocks, cut(seq_along(unique_blocks), num_cores, labels = FALSE))

# Define a function to process one chunk
process_block_chunk <- function(block_group) {
  queue <- create_input_queue()
  
  for (cur_block in block_group) {
    input <- process_block(cur_block)
    if (!is.null(input)) {
      push_input(queue, input)
    }
  }

  finish_queue(queue)
  block_calculation_parallel_streamed(queue, n_threads = num_cores)
}

# Parallel R processes (1 per chunk)
results.chunk <- mclapply(
  block_chunks,
  FUN = process_block_chunk,
  mc.cores = length(block_chunks)
)

# Combine
results.list <- do.call(c, results.chunk)

  # cat("Memory after processing blocks before mclapply:", mem_used(), "\n")
  # blocks <- mclapply(unique_blocks, process_block, mc.cores = num_cores)
  # cat("Memory after mclapply:", mem_used(), "\n")
  
  # results.list <- block_calculation_parallel(blocks)
  # results.list <- block_calculation_parallel(bim_sum_stats, LDblocks2, sum_stats_target_val, lr_list, ref_file, unique_blocks)
  # cat("Memory after RcppParallel part:", mem_used(), "\n")
  
  beta.byL <- do.call(rbind, lapply(results.list, function(x) x[[1]]))
  cols_to_drop <- match(c("A1", "order"), colnames(beta.byL))
  beta.byL <- beta.byL[, -cols_to_drop, drop = FALSE]
  colnames(beta.byL)[1:3] <- c("SNP", "CHR", "A1")
  
  betaRho.byL <- Reduce("+", lapply(results.list, function(x) x[[2]]))
  betaG.byL <- Reduce("+", lapply(results.list, function(x) as.numeric(unlist(x[[3]]))))
  if (is.null(dim(betaG.byL))) {
    betaG.byL <- matrix(betaG.byL, nrow = 1)
  }
  
  gc()
  # cat("Memory after gc():", mem_used(), "\n")
  
  return(list(beta.byL, betaRho.byL, betaG.byL))
}

pseudo_split <- function (target_sumstats, subprop, ref_file_ps, tempfile, random_seed=42) {
  #' Pseudo-split target summary statistics into train and validation set
  #' 
  #' @param target_sumstats Target summary statistics to split
  #' @param subprop Proportion of training samples from full samples
  #' @param ref_file_ps Prefix of PLINK file of the reference panel data in the target population
  #' @param tempfile Output path to save splitted summary statistics
  #' @param random_seed A random seed number to randomly (pseudo) split target summary statistics
  #' @export
  #' @returns Splitted summary statistics for train and validation sammples
  
    ref_bim <- fread(paste0(ref_file_ps,'.bim'))

    target_sumstats_ref <- merge(target_sumstats, ref_bim, by.x='SNP', by.y='V2', sort=F)

    target_sumstats_ref = target_sumstats_ref[!duplicated(target_sumstats_ref$SNP),]
    target_sumstats_ref = target_sumstats_ref[which(target_sumstats_ref$V1 %in% 1:22),]
    ###
    if (('p' %in% colnames(target_sumstats_ref)) & (!'cor' %in% colnames(target_sumstats_ref))) {
      target_sumstats_ref$beta = as.numeric(target_sumstats_ref$beta)
      target_sumstats_ref$p = as.numeric(target_sumstats_ref$p)
      target_sumstats_ref$cor=lassosum::p2cor(p = target_sumstats_ref$p, n = median(target_sumstats_ref$N,na.rm=T), sign=target_sumstats_ref$beta)
    }

    target_sumstats_ref$cor2=NA
    flag1=which(target_sumstats_ref$V5==target_sumstats_ref$A1)
    if (length(flag1)>0){  target_sumstats_ref$cor2[flag1]=target_sumstats_ref$cor[flag1]}
    flag2=which(target_sumstats_ref$V6==target_sumstats_ref$A1)
    if (length(flag2)>0){  target_sumstats_ref$cor2[flag2]=-target_sumstats_ref$cor[flag2]}

    target_sumstats_ref=target_sumstats_ref[which(! is.na(target_sumstats_ref$cor2)),c("SNP","V5","V6","cor2","N")] #Direction -> beta
    colnames(target_sumstats_ref)[2:3] = c('A1','A2') # A1: effect, A2: ref

    BedFileReader_ps <- new( BedFileReader, paste0(ref_file_ps,".fam"), paste0(ref_file_ps,".bim"), paste0(ref_file_ps,".bed"))
    result = try(BedFileReader_ps$snp_index_func(), silent=TRUE)
    
    # X <- readSomeSnp(target_sumstats_ref$SNP, BedFileReader = BedFileReader_ps)
    X <- BedFileReader_ps$readSomeSnp(target_sumstats_ref$SNP, sampleList=integer(0))
    X <- do.call(cbind, X)

    flag=which(apply(X, 2, sd, na.rm = TRUE)!=0)

    X2 = X[,flag] #n*p

    n <- median(target_sumstats_ref$N)
    rhos <- target_sumstats_ref$cor2[flag]

    X2 <- scale(X2)

    set.seed(random_seed)
    g <- rnorm(dim(X2)[1], mean = 0, sd = 1)
    datarands <- t(as.matrix(X2)) %*% g * sqrt(1/dim(X2)[1] * (1-subprop)/subprop) 

    subrhos <- rhos + datarands*sqrt(1/n)
    restrhos <- (rhos - subrhos * subprop) / (1-subprop)

    snp_info <- target_sumstats_ref[flag,1:2]

    sub_mat <- cbind(snp_info, subrhos, rep(subprop * n, dim(snp_info)[1])); colnames(sub_mat)=c('SNP','A1','cor','N')
    rest_mat <- cbind(snp_info, restrhos, rep((1-subprop) * n, dim(snp_info)[1])); colnames(rest_mat)=c('SNP','A1','cor','N')
  
    cat(paste0('Generated pseudo summary statistics with ',dim(sub_mat)[1], ' SNPs. \n'))
   return (list(sub_mat, rest_mat))
}

# PTL_PRS_train<-function(ref_file,sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks="EUR.hg19",outfile,num_cores=1,target_sumstats_train_file=NULL, target_sumstats_val_file=NULL, ps, pseudo_test, random_seed, lr_list, iter, patience, trace=FALSE, random_seed_test=70){
# 	#' Run PTL-PRS
#   #' 
#   #' @param ref_file Prefix of PLINK file of the reference panel data in the target population.
#   #' @param sum_stats_file Location path of source PRS
#   #' @param target_sumstats_file Location path of target summary statistics 
#   #' @param subprop Proportion of training samples from full samples
#   #' @param ref_file_ps Prefix of PLINK file of the reference panel data in the target population
#   #' @param LDblocks Character string specifying the LD blocks, defaults to "EUR.hg19".
#   #' @param outfile Output path to save training results
#   #' @param num_cores A number of cores to be used in parallel computation for training, defaults to 1.
#   #' @param target_sumstats_train_file Optional; Location path of target summary statistics for train samples. Defaults to NULL.
#   #' @param target_sumstats_val_file Optional; Location path of target summary statistics for validation samples. Defaults to NULL.
#   #' @param ps Logical; if TRUE, pseudo-split target summary statistics into train and validation.
#   #' @param pseudo_test Logical; if True, pseudo-split target summary statistics into train, validation and test.
#   #' @param random_seed A random seed number to randomly (pseudo) split target summary statistics
#   #' @param lr_list A list of learning rates
#   #' @param iter A maximum number of iterations
#   #' @param patience A number of times to tolerate decrease of validation pseudo-R
#   #' @param trace Logical; if TRUE, print stopped iteration
#   #' @export
#   #' @return Returns a list containing the following elements:
# #'   \itemize{
# #'     \item \code{best.param}: The best learning rate determined by the validation.
# #'     \item \code{best.beta}: A list of the baseline and best beta coefficients found during validation. (SNP, CHR, A1, Beta2, best.beta)
# #'     \item \code{R2.list}: A list of R-squared values computed during the validation step.
# #'   }
  
#   tempfile=outfile
# 	out1=PRStr_main_check_pv(ref_file, sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks, target_sumstats_train_file, target_sumstats_val_file, ps)
# 	if (out1!=0){stop(out1)}

#   ## Rcpp source code
#   # sourceCpp("/media/leelabsg-storage1/bokeum/TLPRS/source/BedFileReader.cpp") 
#   # setClass("BedFileReader", representation( pointer = "externalptr" ) )
#   # BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
#   # setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
#   # setMethod("initialize","BedFileReader", function(.Object, ...) {
#   #   .Object@pointer <-.Call(BedFileReader_method("new"), ... )
#   #   .Object})

# 	sum_stats=data.frame(fread(sum_stats_file))
# 	if (ncol(sum_stats)==4){ 
# 		if (sum(colnames(sum_stats) %in% c("V1","V2","V3","V4"))==4){
# 			colnames(sum_stats)=c("SNP","CHR","A1","Beta")
# 		}
# 	} 
# 	sum_stats=sum_stats[,c("SNP","CHR","A1","Beta")] 
# 	# sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
# 	# write.table(sum_stats, file=sum_stats_file,col.names=T,row.names=F,quote=F)

# if (ps){
#   cat("=====================================================================================\n")
#   cat("Step 0: Generating pseudo summary statistics... \n")
#   ## pseudo summ generation
#   sum_stats_target <- fread(target_sumstats_file)
#   sum_stats_target <- merge(sum_stats,sum_stats_target,by="SNP", sort=F)  
#   sum_stats_target <- sum_stats_target[,c('SNP','A1.y','beta','p','N')]; colnames(sum_stats_target)[2]='A1'
  
#   target_sumstats <- pseudo_split(sum_stats_target, subprop, ref_file_ps, tempfile, random_seed=random_seed_test)

#   if (!pseudo_test) {
#     target_sumstats_train <- target_sumstats[[1]]
#     target_sumstats_val <- target_sumstats[[2]]
#   } else {
#     target_sumstats_test <- target_sumstats[[2]]

#     target_sumstats <- pseudo_split(target_sumstats[[1]], subprop, ref_file_ps, tempfile, random_seed)
#     target_sumstats_train <- target_sumstats[[1]]
#     target_sumstats_val <- target_sumstats[[2]]

#     write.table(target_sumstats_test, file=paste0(tempfile, '_target.test.summaries'), row.names=F,quote=F,col.names=T)
#   }

#   write.table(target_sumstats_train, file=paste0(tempfile, '_target.train.summaries'), row.names=F,quote=F,col.names=T)
#   write.table(target_sumstats_val, file=paste0(tempfile, '_target.val.summaries'), row.names=F,quote=F,col.names=T)
  
#   cat("Psuedo summary statistics written to files. \n")

#   rm(target_sumstats)
# } else {
#   target_sumstats_train <- fread(target_sumstats_train_file)
#   target_sumstats_val <- fread(target_sumstats_val_file)
# }

# for(group in c('train','val')){
#     sum_stats_target=get(paste0('target_sumstats_',group))
    
#     sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP", sort=F)
	
#     if (('p' %in% colnames(sum_stats_target)) & (!'cor' %in% colnames(sum_stats_target))) {
#       sum_stats_target$beta = as.numeric(sum_stats_target$beta)
#       sum_stats_target$p = as.numeric(sum_stats_target$p)

#       sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, n = median(sum_stats_target$N,na.rm=T), sign=sum_stats_target$beta)
#     }

#       # allele flipping (only cor, not Beta, needs flipping -> only Beta, use ref A1)
#       flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
#       if (length(flag)>0){
#         # sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]
#         sum_stats_target$Beta[flag]=-sum_stats_target$Beta[flag]
#         }
#       # CLUMPING
#       # sum(sum_stats_target$cor>=0.2)
      
#     if(group=='train'){
#       #Beta from source summ stat !!
#       sum_stats_target=sum_stats_target[,c("SNP","A1.y","Beta","cor","N")]
#     }
#     else{
#       sum_stats_target=sum_stats_target[,c("SNP","A1.y","cor","N")]
#     }
#     colnames(sum_stats_target)[2]="A1"; 
#     assign(paste0('sum_stats_target_',group), sum_stats_target)
#     gc()
#   }

#   # block-wise gradient descent
#   cat("=====================================================================================\n")
#   cat("Step 1: Running block-wise gradient descent... \n")
# 	results=PRStr_calculation_pv_es(sum_stats_target_train, ref_file, LDblocks, num_cores, temp.file=paste0(tempfile,"_step1"), lr_list, iter, sum_stats_target_val, patience,trace)
#   cat("Gradient descent completed. \n")

# 	write.table(as.data.frame(results[[1]]),file=paste0(tempfile,"_beta.list.txt"),row.names=F,quote=F,col.names=T)
#   # write.table(as.data.frame(results[[2]]),file=paste0(tempfile,"_betaRho.txt"),row.names=F,quote=F,col.names=F)
# 	# write.table(as.data.frame(results[[3]]),file=paste0(tempfile,"_betaG.txt"),row.names=F,quote=F,col.names=T)
	
#   cat("=====================================================================================\n")
#   cat("Step 2: Running hyperparameter tuning... \n")
#   out1 = PRS_tuning_pv_byLR(results[[1]], results[[2]], results[[3]], lr_list)
#   cat("Hyperparameter tuning completed. \n")

# 	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
# 	write.table(out1$best.param,file=paste0(tempfile,"_best.param.txt"),row.names=F,quote=F,col.names=T)
#   write.table(out1$R2.list,file=paste0(tempfile,"_R2_list.txt"),row.names=F,quote=F,col.names=T)
#   cat("Results have been successfully written with the prefix:", tempfile, "\n")

#   return(out1)
# }

PTL_PRS_train <- function(ref_file, sum_stats_file, target_sumstats_file,
                          subprop, ref_file_ps, LDblocks="EUR.hg19", outfile,
                          num_cores=1, target_sumstats_train_file=NULL,
                          target_sumstats_val_file=NULL, ps, pseudo_test,
                          random_seed, lr_list, iter, patience, trace=FALSE,
                          random_seed_test=70) {
  RcppParallel::setThreadOptions(numThreads = num_cores)

  tempfile=outfile
	out1=PRStr_main_check_pv(ref_file, sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks, target_sumstats_train_file, target_sumstats_val_file, ps)
	if (out1!=0){stop(out1)}

	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==4){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3","V4"))==4){
			colnames(sum_stats)=c("SNP","CHR","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","CHR","A1","Beta")] 
	# sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	# write.table(sum_stats, file=sum_stats_file,col.names=T,row.names=F,quote=F)

if (ps){
  cat("=====================================================================================\n")
  cat("Step 0: Generating pseudo summary statistics... \n")
  ## pseudo summ generation
  sum_stats_target <- fread(target_sumstats_file)
  sum_stats_target <- merge(sum_stats,sum_stats_target,by="SNP", sort=F)  
  sum_stats_target <- sum_stats_target[,c('SNP','A1.y','beta','p','N')]; colnames(sum_stats_target)[2]='A1'
  
  target_sumstats <- pseudo_split(sum_stats_target, subprop, ref_file_ps, tempfile, random_seed=random_seed_test)

  if (!pseudo_test) {
    target_sumstats_train <- target_sumstats[[1]]
    target_sumstats_val <- target_sumstats[[2]]
  } else {
    target_sumstats_test <- target_sumstats[[2]]

    target_sumstats <- pseudo_split(target_sumstats[[1]], subprop, ref_file_ps, tempfile, random_seed)
    target_sumstats_train <- target_sumstats[[1]]
    target_sumstats_val <- target_sumstats[[2]]

    write.table(target_sumstats_test, file=paste0(tempfile, '_target.test.summaries'), row.names=F,quote=F,col.names=T)
  }

  write.table(target_sumstats_train, file=paste0(tempfile, '_target.train.summaries'), row.names=F,quote=F,col.names=T)
  write.table(target_sumstats_val, file=paste0(tempfile, '_target.val.summaries'), row.names=F,quote=F,col.names=T)
  
  cat("Psuedo summary statistics written to files. \n")

  rm(target_sumstats)
} else {
  target_sumstats_train <- fread(target_sumstats_train_file)
  target_sumstats_val <- fread(target_sumstats_val_file)
}

for(group in c('train','val')){
    sum_stats_target=get(paste0('target_sumstats_',group))
    
    sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP", sort=F)
	
    if (('p' %in% colnames(sum_stats_target)) & (!'cor' %in% colnames(sum_stats_target))) {
      sum_stats_target$beta = as.numeric(sum_stats_target$beta)
      sum_stats_target$p = as.numeric(sum_stats_target$p)

      sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, n = median(sum_stats_target$N,na.rm=T), sign=sum_stats_target$beta)
    }

      # allele flipping (only cor, not Beta, needs flipping -> only Beta, use ref A1)
      flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
      if (length(flag)>0){
        # sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]
        sum_stats_target$Beta[flag]=-sum_stats_target$Beta[flag]
        }
      # CLUMPING
      # sum(sum_stats_target$cor>=0.2)
      
    if(group=='train'){
      #Beta from source summ stat !!
      sum_stats_target=sum_stats_target[,c("SNP","A1.y","Beta","cor","N")]
    }
    else{
      sum_stats_target=sum_stats_target[,c("SNP","A1.y","cor","N")]
    }
    colnames(sum_stats_target)[2]="A1"; 
    assign(paste0('sum_stats_target_',group), sum_stats_target)
    # gc()
  }


  cat("=====================================================================================\n")
  cat("Step 1: Running block-wise gradient descent... \n")
  	results=PRStr_calculation_pv_es(sum_stats_target_train, ref_file, LDblocks, num_cores, temp.file=paste0(tempfile,"_step1"), lr_list, iter, sum_stats_target_val, patience,trace)
  cat("Gradient descent completed. \n")

	write.table(as.data.frame(results[[1]]),file=paste0(tempfile,"_beta.list.txt"),row.names=F,quote=F,col.names=T)
  # write.table(as.data.frame(results[[2]]),file=paste0(tempfile,"_betaRho.txt"),row.names=F,quote=F,col.names=F)
	# write.table(as.data.frame(results[[3]]),file=paste0(tempfile,"_betaG.txt"),row.names=F,quote=F,col.names=T)
	
  cat("=====================================================================================\n")
  cat("Step 2: Running hyperparameter tuning... \n")
  out1 = PRS_tuning_pv_byLR(results[[1]], results[[2]], results[[3]], lr_list, median(sum_stats_target_val$N))
  cat("Hyperparameter tuning completed. \n")

	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
	write.table(out1$best.param,file=paste0(tempfile,"_best.param.txt"),row.names=F,quote=F,col.names=T)
  write.table(out1$R2.list,file=paste0(tempfile,"_R2_list.txt"),row.names=F,quote=F,col.names=T)
  cat("Results have been successfully written with the prefix:", tempfile, "\n")

  return(out1)
}