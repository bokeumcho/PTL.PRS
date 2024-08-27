# Same as TLPRS
block_calculation2<-function(cor,num,nsnp,temp.file, lr_list, iter){
  # temp_file=paste0(temp.file,"_block_",num)
  # write.table(cor$V2,file=temp_file,col.names=F,row.names=F,quote=F)
  # cmd = paste0("plink --bfile ",train_file," --extract ",temp_file,   " --recodeA  --out ", temp_file,"_Geno.txt")
  # system(cmd)
  
  tic('reading ref G')
  Gtemp = try(readSomeSnp(cor$V2, BedFileReader = BedFileReader),silent=T) #ref

  cat('read Gtemp block ',num,'\n')
  toc()

  if (class(Gtemp)=="try-error"){
    print('error while reading Bedfile')
    return(NULL)
  }else{
    # Gtemp2 = Gtemp[,apply(Gtemp, 2, sd, na.rm = TRUE)!=0] #4727 -> 4710 snps, n*p

    GG = cor(as.matrix(Gtemp))
    geno_info = as.data.frame(cor[,c('V2','V5')]);colnames(geno_info)=c("SNP","A1")
    # geno_info = geno_info[apply(Gtemp, 2, sd, na.rm = TRUE)!=0,]    
    geno_info$mean=colMeans(as.matrix(Gtemp),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))

    # GG=cor(as.matrix(Gtemp[,7:ncol(Gtemp)]))
    # geno_info=as.data.frame(t(sapply(colnames(Gtemp)[7:ncol(Gtemp)],   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    # geno_info$mean=colMeans(as.matrix(Gtemp[,7:ncol(Gtemp)]),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))
  }

  list1=which(geno_info$sd==0)
  if (length(list1)>0){
    geno_info=geno_info[-list1,]
    GG=GG[-list1,-list1]
  }
  if (nrow(geno_info)==0){
    return(NULL)
  } else {
    geno_info$order=1:nrow(geno_info)
    geno_info2=merge(cor[,c("V1","V2","V5","Beta2","cor")],geno_info, by.x="V2",by.y="SNP",sort=F)

    GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
    gy=geno_info2$cor
    betatemp=geno_info2$Beta2*geno_info2$sd
    u0=gy-GG2%*%betatemp
    beta.all=cbind(u0, betatemp)

    tic("calculating gradients")
    for (factor1 in lr_list){
      k=1

    if (is.null(dim(beta.all))){
      betatemp=beta.all[2]
      u0=beta.all[1]
    } else{
      betatemp=beta.all[,2]
      u0=beta.all[,1]
    }

      while (k<=iter){
        ##betanew=c()
        learningrate=1/nsnp*factor1
        if (learningrate>1){learningrate=1}
        ##print(learningrate)
        for (j in 1:length(betatemp)){
          # beta_old=betatemp[j]
          # betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
          # u0=u0-GG2[,j]*(betatemp[j]-beta_old)
          beta_old=betatemp
          betatemp[j]=(learningrate*u0[j]+beta_old[j])/ 1
          u0[j] = u0[j] - GG2[,j] %*% (betatemp - beta_old)
        }
        beta.all=cbind(beta.all,betatemp)
        k=k+1
      } 
    }
    toc()
    
    geno_info2=cbind(geno_info2,as.data.frame(beta.all)[,-1])
    return(geno_info2)
  }
}##function end

# Same as TLPRS
PRStr_calculation2<-function(sum_stats_target, ref_file, sum_stats, LDblocks, cluster=NULL,temp.file, lr_list, iter){
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of", 
                   paste(possible.LDblocks, collapse=", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
        LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
      }
  } else {
    stop(paste0("LDblocks must be specified. Specify one of ", 
                paste(possible.LDblocks, collapse=", "), 
                ". Alternatively, give an integer vector defining the blocks, ", 
                "or a .bed file with three columns read as a data.frame."))
  }
  
  ref.bim <- fread(paste0(ref_file, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
  ref.bim <- ref.bim[ref.bim$V1 %in% 1:22,]

  ref.bim$order=1:nrow(ref.bim)
  bim_sum_stats=merge(ref.bim, sum_stats_target,by.x="V2",by.y="SNP",order=FALSE)
  bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
  
# remove duplicates
  bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

  bim_sum_stats$Beta2=NA
  flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
  if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
  flag2=which(bim_sum_stats$V6==bim_sum_stats$A1)
  if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2];  bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
  
  bim_sum_stats=bim_sum_stats[which(! is.na(bim_sum_stats$Beta2)),c("V2","V1","V4","V5","V6","order","Beta2","cor")] 
  
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  if(!is.null(LDblocks)) {
      LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                              POS = ref.bim$V4[ ref.extract],
                              ref.CHR = LDblocks[,1], 
                              ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
  } 
  
  # if(is.null(cluster)) {
  # 	results.list <- lapply(unique(LDblocks2[[1]]), function(i) {
  #   		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),nsnp=nrow(bim_sum_stats),temp.file, lr_list, iter)
  # 	})
  # } else {
  # 	results.list <-  parallel::parLapplyLB(cluster,unique(LDblocks2[[1]]), function(i) {
  #   		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),nsnp=nrow(bim_sum_stats),temp.file, lr_list, iter)
  # 	})
  # }
  
    tic('block_calculation')
  results.list <- mclapply(unique(LDblocks2[[1]])[1], function(i) {
    		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),nsnp=nrow(bim_sum_stats),temp.file, lr_list, iter)
  	}, mc.cores=5)
  toc()

  results.list<-do.call("rbind", results.list)

  return(results.list)
}

# geno_ref=Gtemp3; maxiter=iter; betaSD=geno_info2$sd; sum_stats_target_val_block
updateCoeff <- function(beta.all, GG2, geno_ref, learningrate, maxiter, betaSD, sum_stats_target_val_block, trace, update_new){
    k=1

    if (is.null(dim(beta.all))){
      betatemp=beta.all[2]
      u0=beta.all[1]
    } else{
      betatemp=beta.all[,2]
      u0=beta.all[,1]
    }

    R2_list = c()
    prev_R2 = 0
    cnt=0

    while (k<=maxiter){
      if (!update_new){
        for (j in 1:length(betatemp)){
            beta_old=betatemp[j]
            betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
            # u0[j]=u0[j]-sum(GG2[,j]*(betatemp[j]-beta_old))
            u0=u0-GG2[,j]*(betatemp[j]-beta_old)
        }
      } else {
    # diff update
      for (j in 1:length(betatemp)){
          beta_old=betatemp
          betatemp[j]=(learningrate*u0[j]+beta_old[j])/ 1
          u0[j] = u0[j] - GG2[,j] %*% (betatemp - beta_old)
      }
    }

    # beta_old=betatemp
    # betatemp=(learningrate*u0+beta_old)/ 1
    # u0=u0-(GG2%*%(betatemp-beta_old))

    betatemp1 = betatemp/betaSD
    # n_test = sum_stats_target_val_block$N[1]
    R_num = sum(betatemp1*sum_stats_target_val_block$cor) #* n_test : fixed value. distribution does not change anyway..
    PRS = geno_ref%*%betatemp1

    # R_result =Calculate_R(betatemp, sum_stats_target_val, geno_ref)
    
    curr_R2 = R_num^2/sum(PRS^2)
    # cat(betatemp[1], R_num, curr_R2,'\n')

    R2_list = c(R2_list, curr_R2)

    if (curr_R2 < prev_R2) {
        cnt = cnt + 1
    }

    if (cnt > patience) {
        # write.table(R2_val,file=paste0(tempfile,"_R_vals.txt"),row.names=F,quote=F,col.names=T)
        # R2_table[which(lr_list==factor1),1]=k
        # R2_table[which(lr_list==factor1),2]=curr_R2
        break
    }
    prev_R2 = curr_R2

    k=k+1
    # betatemp=betatemp*betaSD
    } 
    
    if(trace) cat("stopped iteration:", k, "\n") 

    return(list(betatemp1, R_num, PRS, R2_list))
}

#num=1; i=unique(LDblocks2[[1]])[num]; cor=bim_sum_stats[which(LDblocks2[[1]]==i),]; nsnp=nrow(bim_sum_stats);
block_calculation_pv_es<-function(cor,num,nsnp,temp.file, lr_list, iter, sum_stats_target_val, patience, trace=FALSE, BedFileReader){
  tic('reading ref G')
  Gtemp = try(readSomeSnp(cor$V2, BedFileReader = BedFileReader),silent=T) #ref

  cat('read Gtemp block ',num,' dim(Gtemp):',dim(Gtemp),'\n')
  toc()

  if (class(Gtemp)=="try-error"){
    print('error while reading Bedfile')
    return(NULL)
  }else{
    # sdZero=which(apply(Gtemp, 2, sd, na.rm = TRUE)!=0)
    # Gtemp2 = Gtemp[,sdZero] #4727 -> 4710 snps, n*p
    GG = cor(as.matrix(Gtemp)) #Gtemp2
    # GG[is.na(GG)] <- 0

    # removedSnps <<- removedSnps + (dim(Gtemp)[2] - dim(Gtemp2)[2])
    # cat('removed snps: ',(dim(Gtemp)[2] - dim(Gtemp2)[2]),'\n')

    geno_info = as.data.frame(cor[,c('V2','V5')]);colnames(geno_info)=c("SNP","A1")
    # geno_info = geno_info[sdZero,]    
    geno_info$mean=colMeans(as.matrix(Gtemp),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))
  }

  list1=which(geno_info$sd==0)
  if (length(list1)>0){
    geno_info=geno_info[-list1,]
    GG=GG[-list1,-list1]
    Gtemp = Gtemp[,-list1]
  }
  if (nrow(geno_info)==0){
    return(NULL)
  } else {
    geno_info$order=1:nrow(geno_info)

    # if not remove sd=0 cols, just cbind
    geno_info2=merge(cor[,c("V1","V2","V5","Beta2","cor")],geno_info, by.x="V2",by.y="SNP",sort=F)

    GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
    gy=geno_info2$cor
    betatemp=geno_info2$Beta2*geno_info2$sd
    u0=gy-GG2%*%betatemp
    beta.init=cbind(u0, betatemp)

  ###
  tic('scaling G')
  Gtemp3 = scale(Gtemp) 
  # Gtemp3 <- scale(Gtemp,scale=FALSE)
  # flag=which(apply(Gtemp, 2, sd, na.rm = TRUE)!=0)
  # Gtemp3[,flag] <- scale(Gtemp3[,flag], center=FALSE)
  toc()

  sum_stats_target_val_block = sum_stats_target_val[sum_stats_target_val$SNP %in% geno_info2$V2,]

  beta.all = c(geno_info2$Beta2)
  beta_rho = c(sum(geno_info2$Beta2 * sum_stats_target_val_block$cor))
  beta_g = c(Gtemp3 %*% geno_info2$Beta2)
  # R2.all = c() # beta_rho^2/sum(beta_g^2)

    tic("calculating gradients")
    for (factor1 in lr_list){
      learningrate=1/nsnp*factor1
      if (learningrate>1){learningrate=1}
      results = updateCoeff(beta.init, GG2, Gtemp3, learningrate, iter, geno_info2$sd, sum_stats_target_val_block, trace, update_new = FALSE)
      beta.all=cbind(beta.all,results[[1]])
      beta_rho = c(beta_rho, results[[2]])
      beta_g = cbind(beta_g, results[[3]])
      # R2.all = cbind(R2.all, results[[4]])
    }

    for (factor1 in lr_list){
      learningrate=1/nsnp*factor1
      if (learningrate>1){learningrate=1}
      results = updateCoeff(beta.init, GG2, Gtemp3, learningrate, iter, geno_info2$sd, sum_stats_target_val_block, trace, update_new = TRUE)
      beta.all=cbind(beta.all,results[[1]])
      beta_rho = c(beta_rho, results[[2]])
      beta_g = cbind(beta_g, results[[3]])
      # R2.all = cbind(R2.all, results[[4]])
  }

    toc()
    
    # write.table(R2.all, paste0(tempfile, '_R_val_block',num,'.txt'), col.names=T, row.names=F, quote=F)
    geno_info2=cbind(geno_info2[,-c('A1','order')],as.data.frame(beta.all)) # delete A1, order
    
    rm(beta.init, beta.all, Gtemp, Gtemp3, betatemp)
    gc()

    return(list(geno_info2, beta_rho, beta_g))
  } 
}##function end

# temp.file=paste0(tempfile,"_step1")
PRStr_calculation_pv_es<-function(sum_stats_target_train, ref_file, sum_stats, LDblocks, num_cores=1,temp.file, lr_list, iter, sum_stats_target_val,patience, trace){
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of", 
                   paste(possible.LDblocks, collapse=", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
        LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
      }
  } else {
    stop(paste0("LDblocks must be specified. Specify one of ", 
                paste(possible.LDblocks, collapse=", "), 
                ". Alternatively, give an integer vector defining the blocks, ", 
                "or a .bed file with three columns read as a data.frame."))
  }
  
  ref.bim <- fread(paste0(ref_file, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
  ref.bim <- ref.bim[ref.bim$V1 %in% 1:22,]

  ref.bim$order=1:nrow(ref.bim)
  bim_sum_stats=merge(ref.bim, sum_stats_target_train,by.x="V2",by.y="SNP",order=FALSE)
  bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
  
# remove duplicates
  bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

  bim_sum_stats$Beta2=NA
  flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
  if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
  flag2=which(bim_sum_stats$V6==bim_sum_stats$A1) # flipped. when ref=ref_ps, length(flag2)=0
  if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2]; bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
  
  bim_sum_stats=bim_sum_stats[which(!is.na(bim_sum_stats$Beta2)),c("V2","V1","V5","V6","order","Beta2","cor")] 
  
  ## match ref with sum_stats_target_val & check allele flip 
    bim_sum_stats_val=merge(ref.bim, sum_stats_target_val, by.x="V2",by.y="SNP",order=FALSE)
  bim_sum_stats_val=bim_sum_stats_val[order(bim_sum_stats_val$order),]
  
  bim_sum_stats_val = bim_sum_stats_val[!duplicated(bim_sum_stats_val$V2),]

  flag3=which(bim_sum_stats_val$V6==bim_sum_stats_val$A1) # flipped. when ref=ref_ps, length(flag3)=0
  if (length(flag3)>0){  bim_sum_stats_val$cor[flag3]=-bim_sum_stats_val$cor[flag3];}
  
  cat('PRStr_calculation_pv_es flag1,2,3: ',length(flag1),length(flag2),length(flag3),'\n')
  
  bim_sum_stats_val=bim_sum_stats_val[,c("V2","V5","cor","N")] 
  colnames(bim_sum_stats_val)[1:2] = c('SNP','A1')

  # remove rs12952492 17_KI270857v1_alt  G  C 53433029 0.0000000000 -0.01686098 !!!!!!
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  if(!is.null(LDblocks)) {
      LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                              POS = ref.bim$V4[ ref.extract],
                              ref.CHR = LDblocks[,1], 
                              ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
  } 
  
  # assign BedFileReader object
  BedFileReader <- new( "BedFileReader", paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))
  result = try(BedFileReader$snp_index_func(), silent=TRUE)

  tic('block_calculation')
  results.list <- mclapply(unique(LDblocks2[[1]]), function(i) {
    		block_calculation_pv_es(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),nsnp=nrow(bim_sum_stats),
            temp.file, lr_list, iter,bim_sum_stats_val, patience, trace, BedFileReader)
  	}, mc.cores=num_cores)
  toc()

  # error_indices <- sapply(results.list, function(df) {all(ncol(df[[1]]) < 21)})
  # # which(error_indices)

  # beta.errors <- do.call("rbind", lapply(results.list[error_indices], function(x) cbind(x[[1]][1,1:8],t(x[[1]][,9]))))
  # colnames(beta.errors)[9] = 'betatemp'

  # beta.normals <- do.call("rbind", lapply(results.list[!error_indices], function(x) x[[1]]))
  # beta.byL = rbind(beta.errors, beta.normals)

  beta.byL <- do.call("rbind", lapply(results.list, function(x) x[[1]]))
  colnames(beta.byL)[1:3]=c("SNP","CHR","A1")

  betaRho.byL <- Reduce("+", lapply(results.list, function(x) x[[2]])) # scalar
  betaG.byL <- Reduce("+", lapply(results.list, function(x) x[[3]])) # 1*n

  return(list(beta.byL,betaRho.byL,betaG.byL))
}

# target_sumstats=sum_stats_target
pseudo_sum <- function (target_sumstats, subprop, ref_file_ps, tempfile, random_seed=42) {
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

    cat('target_sumstats_ref flag: ', length(flag2), '\n')

    target_sumstats_ref=target_sumstats_ref[which(! is.na(target_sumstats_ref$cor2)),c("SNP","V5","V6","cor2","N")] #Direction -> beta
    colnames(target_sumstats_ref)[2:3] = c('A1','A2') #A1: effect, A2: ref

    BedFileReader_ps <- new( "BedFileReader", paste0(ref_file_ps,".fam"), paste0(ref_file_ps,".bim"), paste0(ref_file_ps,".bed"))
    result = try(BedFileReader_ps$snp_index_func(), silent=TRUE)
    
    tic('reading whole reference panel')
    X <- readSomeSnp(target_sumstats_ref$SNP, BedFileReader = BedFileReader_ps)
    toc()

    cat('dimension X: ',dim(X),'\n')

    flag=which(apply(X, 2, sd, na.rm = TRUE)!=0)

    X2 = X[,flag] #n*p
    cat('dimension X2: ',dim(X2),'\n')

    n <- median(target_sumstats_ref$N)
    rhos <- target_sumstats_ref$cor2[flag]

    #sweep(X, 2, colMeans(X),'-') 
    X2 <- scale(X2)

    # X2[,flag] <- scale(X2[,flag], center=FALSE)

    set.seed(random_seed)
    g <- rnorm(dim(X2)[1], mean = 0, sd = 1)
    datarands <- t(as.matrix(X2)) %*% g * sqrt(1/dim(X2)[1] * (1-subprop)/subprop) 

    subrhos <- rhos + datarands*sqrt(1/n)
    restrhos <- (rhos - subrhos * subprop) / (1-subprop)

    # cat(dim(target_sumstats_ref))
    # print(head(subrhos))
    snp_info <- target_sumstats_ref[flag,1:2]

    sub_mat <- cbind(snp_info, subrhos, rep(subprop * n, dim(snp_info)[1])); colnames(sub_mat)=c('SNP','A1','cor','N')
    rest_mat <- cbind(snp_info, restrhos, rep((1-subprop) * n, dim(snp_info)[1])); colnames(rest_mat)=c('SNP','A1','cor','N')
  
   write.table(sub_mat, file=paste0(tempfile, '_target.train.summaries'), row.names=F,quote=F,col.names=T)
   write.table(rest_mat, file=paste0(tempfile, '_target.val.summaries'), row.names=F,quote=F,col.names=T)

   return (list(sub_mat, rest_mat))
}

PTL_PRS<-function(ref_file,sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks="EUR.hg19",outfile,cluster=NULL, target_sumstats_train_file=NULL, target_sumstats_val_file=NULL, ps,lr_list, iter, early_stopping, patience){
	tempfile=outfile
	out1=PRStr_main_check_pv(ref_file, sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks, target_sumstats_train_file, target_sumstats_val_file, ps)
	if (out1!=0){stop(out1)}

  tic('Rcpp prep - ref')
    sourceCpp("/media/leelabsg-storage1/bokeum/TLPRS/source/BedFileReader.cpp")
  setClass("BedFileReader", representation( pointer = "externalptr" ) )
  BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
  setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
  setMethod("initialize","BedFileReader", function(.Object, ...) {
    .Object@pointer <-.Call(BedFileReader_method("new"), ... )
    .Object})

  BedFileReader <<- new( "BedFileReader", paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))

  result = try(BedFileReader$snp_index_func(), silent=TRUE)
  
  toc()

	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==4){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3","V4"))==4){
			colnames(sum_stats)=c("SNP","CHR","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","CHR","A1","Beta")] 
	sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	write.table(sum_stats, file=sum_stats_file,col.names=T,row.names=F,quote=F)

if (ps){
  ## pseudo summ generation
  sum_stats_target <- fread(target_sumstats_file)
  sum_stats_target <- merge(sum_stats,sum_stats_target,by="SNP", sort=F)  
  sum_stats_target <- sum_stats_target[,c('SNP','A1.y','beta','p','N')]; colnames(sum_stats_target)[2]='A1'
  
  tic('pseudo summ generation')
  target_sumstats <- pseudo_sum(sum_stats_target, subprop, ref_file_ps, tempfile)
  toc()

  target_sumstats_train <- target_sumstats[[1]]
  target_sumstats_val <- target_sumstats[[2]]
  rm(target_sumstats)
} else{
  target_sumstats_train <- fread(target_sumstats_train_file)
  target_sumstats_val <- fread(target_sumstats_val_file)
}

tic('summ stat prep')
for(group in c('train','val')){
    sum_stats_target=get(paste0('target_sumstats_',group))
    
    sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP", sort=F)
	
    if (('p' %in% colnames(sum_stats_target)) & (!'cor' %in% colnames(sum_stats_target))) {
      sum_stats_target$beta = as.numeric(sum_stats_target$beta)
      sum_stats_target$p = as.numeric(sum_stats_target$p)

      sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, n = median(sum_stats_target$N,na.rm=T), sign=sum_stats_target$beta)
    }

      # allele flipping (only cor, not Beta, needs flipping)
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
    } else{
      sum_stats_target=sum_stats_target[,c("SNP","A1.y","cor","N")]
    }
    colnames(sum_stats_target)[2]="A1"; 
    assign(paste0('sum_stats_target_',group), sum_stats_target)

    rm(sum_stats_target); gc()
  }
toc()

  # 1699 blocks in total
  # 	tempfile=outfile
  tic('gradient descent')
	beta_list=as.data.frame(PRStr_calculation2(sum_stats_target_train, ref_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"), lr_list, iter))
	toc()

  beta_list=as.data.frame(beta_list[,-c(6,10)]) #del A1,order; idx changed due to V1
	colnames(beta_list)[1:3]=c("SNP","CHR","A1")
	
  tic('writing beta candidates')
  write.table(beta_list,file=paste0(tempfile,"_beta.candidates.txt"),row.names=F,quote=F,col.names=T)
  toc()
  #beta_list = as.data.frame(fread(paste0(tempfile,"_beta.candidates.txt")))
	
  tic('hyperparam tuning')
  if (!early_stopping){
    out1=PRStr_tuning_pv(beta_list, sum_stats_target_train, sum_stats_target_val, tempfile, lr_list, iter)
  }
  else{
    if(is.null(patience)) patience=3 #default patience
    out1=PRStr_tuning_pv_es(beta_list, sum_stats_target_val, tempfile, lr_list, iter, patience)
  }
  toc()
  
  if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))} # nolint
	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T) # nolint: line_length_linter.
	write.table(out1$best.params,file=paste0(tempfile,"_best.params.txt"),row.names=F,quote=F,col.names=T)

	return(out1)
}

PTL_PRS_bwes<-function(ref_file,sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks="EUR.hg19",outfile,num_cores=1, target_sumstats_train_file=NULL, target_sumstats_val_file=NULL, ps,random_seed, lr_list, iter, early_stopping, patience, trace=FALSE){
	tempfile=outfile
	# out1=PRStr_main_check_pv(ref_file, sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks, target_sumstats_train_file, target_sumstats_val_file, ps)
	# if (out1!=0){stop(out1)}

  tic('Rcpp preparation')
  sourceCpp("/media/leelabsg-storage1/bokeum/TLPRS/source/BedFileReader.cpp")
  setClass("BedFileReader", representation( pointer = "externalptr" ) )
  BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
  setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
  setMethod("initialize","BedFileReader", function(.Object, ...) {
    .Object@pointer <-.Call(BedFileReader_method("new"), ... )
    .Object})

  # BedFileReader <- new( "BedFileReader", paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))
  # result = try(BedFileReader$snp_index_func(), silent=TRUE)
  
  toc()

	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==4){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3","V4"))==4){
			colnames(sum_stats)=c("SNP","CHR","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","CHR","A1","Beta")] 
	sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	write.table(sum_stats, file=sum_stats_file,col.names=T,row.names=F,quote=F)

if (ps){
  ## pseudo summ generation
  sum_stats_target <- fread(target_sumstats_file)
  sum_stats_target <- merge(sum_stats,sum_stats_target,by="SNP", sort=F)  
  sum_stats_target <- sum_stats_target[,c('SNP','A1.y','beta','p','N')]; colnames(sum_stats_target)[2]='A1'
  
  # DiffRef <- ref_file!=ref_file_ps
  tic('pseudo summ generation')
  target_sumstats <- pseudo_sum(sum_stats_target, subprop, ref_file_ps, tempfile, random_seed)
  toc()

  target_sumstats_train <- target_sumstats[[1]]
  target_sumstats_val <- target_sumstats[[2]]

  #
  # target_sumstats_test <- target_sumstats[[2]]

  # target_sumstats <- pseudo_sum(target_sumstats[[1]], subprop, ref_file_ps, tempfile)
  # target_sumstats_train <- target_sumstats[[1]]
  # target_sumstats_val <- target_sumstats[[2]]
  # write.table(target_sumstats_train, file=paste0(tempfile, '_target.train.summaries'), row.names=F,quote=F,col.names=T)
  #  write.table(target_sumstats_val, file=paste0(tempfile, '_target.val.summaries'), row.names=F,quote=F,col.names=T)
  #  write.table(target_sumstats_test, file=paste0(tempfile, '_target.test.summaries'), row.names=F,quote=F,col.names=T)
  #

  rm(target_sumstats)
} else {
  target_sumstats_train <- fread(target_sumstats_train_file)
  target_sumstats_val <- fread(target_sumstats_val_file)
}

tic('summ stat prep')
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
      cat('sum_stats_target flag: ',length(flag))
      # CLUMPING
      # sum(sum_stats_target$cor>=0.2)
      
    if(group=='train'){
      #Beta from source summ stat !!
      sum_stats_target=sum_stats_target[,c("SNP","A1.y","Beta","cor","N")] #A1.x
    }
    else{
      sum_stats_target=sum_stats_target[,c("SNP","A1.y","cor","N")] #A1.x
    }
    colnames(sum_stats_target)[2]="A1"; 
    assign(paste0('sum_stats_target_',group), sum_stats_target)
    gc()
  }
toc()

  # block-wise gradient descent
  tic('gradient descent')
	results=PRStr_calculation_pv_es(sum_stats_target_train, ref_file, sum_stats, LDblocks, num_cores,temp.file=paste0(tempfile,"_step1"), lr_list, iter, sum_stats_target_val, patience,trace)
	toc()
  
  tic('writing training results')
	write.table(as.data.frame(results[[1]]),file=paste0(tempfile,"_beta.list.txt"),row.names=F,quote=F,col.names=T)
  write.table(as.data.frame(results[[2]]),file=paste0(tempfile,"_betaRho.txt"),row.names=F,quote=F,col.names=F)
	write.table(as.data.frame(results[[3]]),file=paste0(tempfile,"_betaG.txt"),row.names=F,quote=F,col.names=T)
  toc()
  #tempfile=outfile
  #beta_list = as.data.frame(fread(paste0(tempfile,"_beta.list.txt")))
	
  tic('hyperparam tuning')
  # out1=PRStr_tuning_pv(beta_list, sum_stats_target_train, sum_stats_target_val, tempfile, lr_list, iter,BedFileReader)
  out1 = PRS_tuning_pv_byLR(results[[1]], results[[2]], results[[3]], lr_list)
  toc()

  if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))}
	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
	write.table(out1$best.param,file=paste0(tempfile,"_best.param.txt"),row.names=F,quote=F,col.names=T)
  write.table(out1$R2.list,file=paste0(tempfile,"_R2_list.txt"),row.names=F,quote=F,col.names=T)
	# write.table(out1$beta.list,file=paste0(tempfile,"_beta.list.txt"),row.names=F,quote=F,col.names=T)
	return(out1)
}

