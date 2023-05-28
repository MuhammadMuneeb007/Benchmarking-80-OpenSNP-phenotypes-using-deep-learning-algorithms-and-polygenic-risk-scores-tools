args <- commandArgs(trailingOnly = TRUE)
print(args)
if (args[2]=="1"){
  #args<-c("CEU_5_1_prsData")
  
  result <-paste(".",args[1],"files",toString("test.QC.het"),sep="//")
  dat <- read.table(result, header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  result <-paste("./",args[1],"files",toString("test.valid.sample"),sep="//")
  write.table(valid[,c(1,2)], result, quote=F, row.names=F) # print FID and IID for valid samples
  result <-paste("./",args[1],"test",toString("test.bim"),sep="//")
  bim <- read.table(result)
  
  colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
  
  # Read in QCed SNPs
  result <-paste("./",args[1],"files",toString("test.QC.snplist"),sep="//")
  
  qc <- read.table(result, header = F, stringsAsFactors = F)
  # Read in the GWAS data
  result <-paste("./",args[1],"files",toString("Data.QC.gz"),sep="//")
  #result <-paste("./",args[1],"train",toString("train.assoc.fisher"),sep="//")
  
  height <-read.table(gzfile(result),
               header = T,
               stringsAsFactors = F, 
               sep="\t")
  # Change all alleles to upper case for easy comparison
  
  height$A1 <- toupper(height$A1)
  height$A2 <- toupper(height$A2)
  bim$B.A1 <- toupper(bim$B.A1)
  bim$B.A2 <- toupper(bim$B.A2)
  info <- merge(bim, height, by = c("SNP", "CHR", "BP"))
  # Filter QCed SNPs
  
  info <- info[info$SNP %in% qc$V1,]
  
  # Function for finding the complementary allele
  
  complement <- function(x) {
    switch (
      x,
      "A" = "T",
      "C" = "G",
      "T" = "A",
      "G" = "C",
      return(NA)
    )
  }
  
  # Get SNPs that have the same alleles across base and target
  info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
  # Identify SNPs that are complementary between base and target
  info$C.A1 <- sapply(info$B.A1, complement)
  info$C.A2 <- sapply(info$B.A2, complement)
  info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
  # Update the complementary alleles in the bim file
  # This allow us to match the allele in subsequent analysis
  
  complement.snps <- bim$SNP %in% info.complement$SNP
  bim[complement.snps,]$B.A1 <-
    sapply(bim[complement.snps,]$B.A1, complement)
  bim[complement.snps,]$B.A2 <-
    sapply(bim[complement.snps,]$B.A2, complement)
  
  # identify SNPs that need recoding
  info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
  # Update the recode SNPs
  recode.snps <- bim$SNP %in% info.recode$SNP
  tmp <- bim[recode.snps,]$B.A1
  bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
  bim[recode.snps,]$B.A2 <- tmp
  
  # identify SNPs that need recoding & complement
  info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
  # Update the recode + strand flip SNPs
  com.snps <- bim$SNP %in% info.crecode$SNP
  tmp <- bim[com.snps,]$B.A1
  bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
  bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
  result <-paste("./",args[1],"files",toString("test.a1"),sep="//")
  
  # Output updated bim file
  write.table(
    bim[,c("SNP", "B.A1")],
    result,
    quote = F,
    row.names = F,
    col.names = F,
    sep="\t"
  )
  mismatch <-
    bim$SNP[!(bim$SNP %in% info.match$SNP |
                bim$SNP %in% info.complement$SNP | 
                bim$SNP %in% info.recode$SNP |
                bim$SNP %in% info.crecode$SNP)]
  result <-paste("./",args[1],"files",toString("test.mismatch"),sep="//")
  
  write.table(
    mismatch,
    result,
    quote = F,
    row.names = F,
    col.names = F
  )
  
  
}
if (args[2]=="2"){
  result <-paste("./",args[1],"files",toString("Data.QC.gz"),sep="//")
  #result <-paste("./",args[1],"train",toString("train.assoc.fisher"),sep="//")
  valid <- read.table(result, header=T)
  result <-paste("./",args[1],"files",toString("test.QC.sexcheck"),sep="//")
  dat <- read.table(result, header=T)
  valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
  result <-paste("./",args[1],"files",toString("test.QC.valid"),sep="//")
  write.table(valid[,c("FID", "IID")], result, row.names=F, col.names=F, sep="\t", quote=F) 
}
if (args[2]=="3"){
  
  result <-paste("./",args[1],"files",toString("Data.QC.gz"),sep="//")
  #result <-paste("./",args[1],"train",toString("train.assoc.fisher"),sep="//")
  
  dat <- read.table(gzfile(result), header=T)
  dat$BETA <- log(dat$OR)
  result <-paste("./",args[1],"files",toString("Data.QC.Transformed"),sep="//")
  
  write.table(dat, result, quote=F, row.names=F)
  }
if (args[2]=="4"){

  result <-paste("./",args[1],"test",toString("YRI.pheno"),sep="//")
  
  p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
  phenotypes <- read.table(result, sep="\t",header=T)
  result <-paste("./",args[1],"files",toString("test.eigenvec"),sep="//")
  pcs <- read.table(result, header=F)
  xx <-ncol(pcs)-2
  colnames(pcs) <- c("FID", "IID", paste0("PC",1:xx)) 
  
  # Read in the covariates (here, it is sex)
  
  #pcs$FID <- as.character(pcs$FID)
  #pcs$IID <- as.character(pcs$FID)
  
  #pcs$FID <- paste(pcs$FID, pcs$FID,sep="_")
  #pcs$IID <- paste(pcs$IID, pcs$IID,sep="_")
  
  #result <-paste("./",args[1],"test",toString("YRI.covariate"),sep="//")
  
  #print(head(phenotypes))
  #print(head(covariate))
  #print(head(pcs))
  # Now merge the files
  
  pheno <- merge(phenotypes, pcs, by=c("FID","IID"))
  pheno$phenotype[pheno$phenotype==1]<-0
  pheno$phenotype[pheno$phenotype==2]<-1
  
  
  #pheno <- merge(merge(phenotypes, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
  
  # We can then calculate the null model (model with PRS) using a linear regression 
  # (as height is quantitative)
  
  #null.model <- lm(phenotype~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
  null.model <- glm(phenotype~., data=pheno[,!colnames(pheno)%in%c("FID","IID")],family=binomial)
  
  library(fmsb)
  # And the R2 of the null model is 
  #null.r2 <- summary(null.model)$r.squared
  null.r2 <- fmsb::NagelkerkeR2(null.model)
  #prs.coef <- summary(null.model)$coeff["SCORE",]
  #print(prs.coef)
  #exit(0)
  prs.result <- NULL
  for(i in p.threshold){
    # Go through each p-value threshold
    result <-paste("./",args[1],"files",paste0("test.",i,".profile"),sep="//")
    prs <- read.table(result, header=T)
    #prs$FID <- as.character(prs$FID)
    #prs$IID <- as.character(prs$FID)
  
    #prs$FID <- paste(prs$FID, prs$FID,sep="_")
    #prs$IID <- paste(prs$IID, prs$IID,sep="_")
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- glm(phenotype~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")],family=binomial)
    
    # model R2 is obtained as 
    model.r2 <- fmsb::NagelkerkeR2(model)
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2[[1]]-null.r2[[1]]
   
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    
    #print(summary(model$coeff["SCORE",]))
   

    prs.coef <- summary(model)$coeff["SCORE"]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    
    # We can then store the results
    #print(prs.coef)
      
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
  }
  # Best result is:
  print(prs.result[which.max(prs.result$R2),])

}
if (args[2]=="5"){
  info <- readRDS("map.rds")
 # install.packages("remotes")
  library(remotes)
  #remotes::install_github("https://github.com/privefl/bigsnpr.git",force = TRUE)
  library(bigsnpr)
  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  library(data.table)
  library(magrittr)
 
  result <-paste("./",args[1],"test","YRI.pheno",sep="//")
  phenotype <- fread(result)
  #result <-paste("./",args[1],"test","YRI.covariate",sep="//")
  #covariate <- fread(result)
  result <-paste("./",args[1],"files","test.eigenvec",sep="//")
  pcs <- fread(result)
  # rename columns
  colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
  # generate required table
  
  print(head(phenotype))
  print(head(pcs))
  pcs$FID <- as.character(pcs$FID)
  pcs$IID <- as.character(pcs$IID)
  phenotype$FID <- as.character(phenotype$FID)
  phenotype$IID <- as.character(phenotype$IID)
  #covariate$FID <- as.character(covariate$FID)
  #covariate$IID <- as.character(covariate$IID)
      
  #pcs$FID <- paste(pcs$FID, pcs$FID,sep="_")
  #pcs$IID <- paste(pcs$IID, pcs$IID,sep="_")
  pheno <- merge(phenotype, pcs ,sort = F)
  #pheno <- phenotype

  print(head(pheno))

  result <-paste("./",args[1],"files","Data.QC.gz",sep="//")
  
  sumstats <- bigreadr::fread2(result) 
  # LDpred 2 require the header to follow the exact naming
  names(sumstats) <-
    c("chr",
      "pos",
      "rsid",
      "a1",
      "a0",
      "n_eff",
      "beta_se",
      "p",
      "OR",
      "INFO",
      "MAF")

  # Transform the OR into log(OR)
  sumstats$beta <- log(sumstats$OR)
  # Filter out hapmap SNPs
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

    
    # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  
  tmp <- tempfile(tmpdir = "tmp-data")
  #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  # We want to know the ordering of samples in the bed file 
  fam.order <- NULL
  # preprocess the bed file (only need to do once for each data set)
  result <-paste("./",args[1],"files","test.QC2.bed",sep="//")
  snp_readBed(result)
  # now attach the genotype object
  result <-paste("./",args[1],"files","test.QC2.rds",sep="//")
  obj.bigSNP <- snp_attach(result)

  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  #print(map)
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # perform SNP matching
  map$a0 <- as.character(map$a0)
  map$a1 <- as.character(map$a1)

  
  print(head(map))
  print(head(sumstats))

  info_snp <- snp_match(sumstats, map)
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  #print(CHR)
  #print(POS)
  print("DONE1")
  
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  print(POS2)
  print("DONE2")
  
  # calculate LD
  info_snp$chr <- as.numeric(info_snp$chr)
  print(head(info_snp)) 

  for (chr in 1:22) {  
    
    print(chr)
    # Extract SNPs that are included in the chromosome
    #print(info_snp$chr)
    try({
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
 
   
    
    corr0 <- snp_cor(
      genotype,
      ind.col = ind.chr2,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      size = 3 / 1000
    )   
    
    if (chr == 1) {
      
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp)
    } else {
      
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
    },silent=TRUE)
  }
  # We assume the fam order is the same across different chromosomes
  fam.order <- as.data.table(obj.bigSNP$fam)
  print(head(fam.order))
  #fam.order$FID <- as.character(pheno$FID)
  #fam.order$IID <- as.character(pheno$IID)

  fam.order$affection[fam.order$affection==1]<-0
  fam.order$affection[fam.order$affection==2]<-1


  # Rename fam order
  setnames(fam.order,
          c("family.ID", "sample.ID"),
          c("FID", "IID"))
    df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  print(length(ld))
  print(length((df_beta$beta / df_beta$beta_se)^2))
  print(length((df_beta$beta / df_beta$beta_se)))
  
  print(length(df_beta$beta ))
  print(length( df_beta$beta_se))
  chi2 <- (df_beta$beta / df_beta$beta_se)^2
  print(length(chi2))
  ldsc <- snp_ldsc(   ld, 
                      length(ld), 
                      chi2 = (df_beta$beta / df_beta$beta_se)^2,
                      sample_size = df_beta$n_eff, 
                      blocks = NULL)
  h2_est <- ldsc[["h2"]]
  print("Xas") 
  print(h2_est)
  print(ldsc[["h2"]])
  print(ldsc) 
  library(fmsb)
  print((fam.order))
  pheno$FID <- as.character(pheno$FID)
  pheno$IID <- as.character(pheno$IID)
  fam.order$FID <- as.character(fam.order$FID)
  fam.order$IID <- as.character(fam.order$IID)

  #Data$genome <- NULL
  # Reformat the phenotype file such that y is of the same order as the 
  # sample ordering in the genotype file
  y <- pheno[fam.order, on = c("FID", "IID")]
  #y$i.FID <- NULL
  #y$i.IID <- NULL
  y$affection <- NULL
  y$paternal.ID <- NULL
  y$maternal.ID <- NULL
  y$sex <- NULL
  
  #y$i.IID <- NULL
  y$phenotype[y$phenotype==1]<-0
  y$phenotype[y$phenotype==2]<-1
  #y$phenotype[y$phenotype==1]<-2
  #y$phenotype[y$phenotype==0]<-1
  
  print(h2_est)

  
  #null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
 #    paste0("phenotype~+", .) %>%
#     as.formula %>%
#      glm(., data = y, family=binomial) %>%
 #     summary
 
  
 # print(summary)
 # print(null.model)
#  null.r2 <- fmsb::NagelkerkeR2(null.model)
 # print(null.r2)
  

  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
  if(is.null(obj.bigSNP)){
    result <-paste("./",args[1],"files","test.QC2.rds",sep="//")
      obj.bigSNP <- snp_attach(result)
  }
  genotype <- obj.bigSNP$genotypes
  # calculate PRS for all samples
  ind.test <- 1:nrow(genotype)
  pred_inf <- big_prodVec(    genotype,
                              beta_inf,
                              ind.row = ind.test,
                              ind.col = info_snp$`_NUM_ID_`)
  print(pred_inf)
  result <-paste("./",args[1],"result","infldpredtest.txt",sep="//")
  write.csv(pred_inf, result, sep = "\t",row.names = FALSE)
  
  # Prepare data for grid model
  p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
  h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
  grid.param <-
    expand.grid(p = p_seq,
                h2 = h2_seq,
                sparse = c(FALSE, TRUE))
  # Get adjusted beta from grid model
  beta_grid <-
    snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
  if(is.null(obj.bigSNP)){
        result <-paste("./",args[1],"files","test.QC2.rds",sep="//")
      obj.bigSNP <- snp_attach(result)
  }
  genotype <- obj.bigSNP$genotypes
  # calculate PRS for all samples
  ind.test <- 1:nrow(genotype)
  pred_grid <- big_prodMat(   genotype, 
                              beta_grid, 
                              ind.col = info_snp$`_NUM_ID_`)
  
  print(pred_grid)
  result <-paste("./",args[1],"result","gridldpredtest.txt",sep="//")
  write.csv(pred_grid, result, sep = "\t",row.names = FALSE)
  # Get adjusted beta from the auto model
  multi_auto <- snp_ldpred2_auto(
    corr,
    df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
    ncores = NCORES
  )
  beta_auto <- sapply(multi_auto, function(auto)
    auto$beta_est)
  if(is.null(obj.bigSNP)){
        result <-paste("./",args[1],"files","test.QC2.rds",sep="//")
      obj.bigSNP <- snp_attach(result)
  }
  genotype <- obj.bigSNP$genotypes
  # calculate PRS for all samples
  ind.test <- 1:nrow(genotype)
  pred_auto <-
    big_prodMat(genotype,
                beta_auto,
                ind.row = ind.test,
                ind.col = info_snp$`_NUM_ID_`)
  # scale the PRS generated from AUTO
  pred_scaled <- apply(pred_auto, 2, sd)
  final_beta_auto <-
    rowMeans(beta_auto[,
                       abs(pred_scaled -
                             median(pred_scaled)) <
                         3 * mad(pred_scaled)])
  pred_auto <-
    big_prodVec(genotype,
                final_beta_auto,
                ind.row = ind.test,
                ind.col = info_snp$`_NUM_ID_`)
  
  print(pred_auto)
  result <-paste("./",args[1],"result","autoldpredtest.txt",sep="//")
   write.csv(pred_auto, result, sep = "\t",row.names = FALSE)
}

if (args[2]=="6"){
  library(lassosum)
  library(data.table)
  library(methods)
  library(magrittr)
  library(parallel)
  cl <- makeCluster(2)
  result <- paste("./",args[1],"files","Data.QC.gz",sep="//")
  sum.stat <- result
  result <-paste("./",args[1],"files","test.QC",sep="//")
  bfile <- result
  # Read in and process the covariates
  #result <-paste("./",args[1],"test","YRI.covariate",sep="//")
  #covariate <- fread(result)
  result <-paste("./",args[1],"files","test.eigenvec",sep="//")
  result <-paste("./",args[1],"files",toString("test.eigenvec"),sep="//")
  pcs <- read.table(result, header=F)
  xx <-ncol(pcs)-2
  pcs <- fread(result) %>% setnames(., colnames(.), c("FID","IID", paste0("PC",1:xx)))
  # Need as.data.frame here as lassosum doesn't handle data.table 
  # covariates very well
  
  #pcs$FID <- as.character(pcs$FID)
  #pcs$IID <- as.character(pcs$FID)
  
  #pcs$FID <- paste(pcs$FID, pcs$FID,sep="_")
  #pcs$IID <- paste(pcs$IID, pcs$IID,sep="_")
  
  cov <-pcs
  
  # We will need the EUR.hg19 file provided by lassosum 
  # which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
  ld.file <- "EUR.hg19"
  # output prefix
  prefix <- "EUR"
  # Read in the target phenotype file
  result <-paste("./",args[1],"test","YRI.pheno",sep="//")
  #bfile <-paste("./",args[1],"train","train",sep="//")
  bfile <-paste("./",args[1],"test","test",sep="//")

  #bfile <-paste("./",args[1],"files","test.QC",sep="//")
  
  target.pheno <- fread(result)[,c("FID", "IID", "phenotype")]
  
  # Read in the summary statistics
  ss <- fread(sum.stat)
  # Remove P-value = 0, which causes problem in the transformation
  ss <- ss[!P == 0]
  # Transform the P-values into correlation
  cor <- p2cor(p = ss$P,
               n = ss$N,
               sign = log(ss$OR)
  )
  result <-paste("./",args[1],"test","test.fam",sep="//")
  
  fam <- fread(result)
  #fam$V1 <- as.character(fam$V1)
  #fam$V2 <- as.character(fam$V2)
  
  #fam$V1 <- paste(fam$V1, fam$V1,sep="_")
  #fam$V2 <- paste(fam$V2, fam$V2,sep="_")
  #fam$V6 <- paste(fam$V2, fam$V2,sep="_")
  fam$V6[fam$V6==1]<-0
  fam$V6[fam$V6==2]<-1

  fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]
  
  # Run the lassosum pipeline
  # The cluster parameter is used for multi-threading
  # You can ignore that if you do not wish to perform multi-threaded processing
  out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file, 
    cluster=cl
  )
  # Store the R2 results
  result <-paste("./",args[1],"result","YRI.result",sep="//")
  pdf(file = NULL)
dev.off()
  target.res <- validate(out, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov))
  # Get the maximum R2
  help(validate)
  result <-paste("./",args[1],"result","test.txt",sep="//")
  
  lapply(target.res[["best.pgs"]], write, result, append=TRUE, ncolumns=1000)
  r2 <- max(target.res$validation.table$value)^2
  print(r2)
}


  
