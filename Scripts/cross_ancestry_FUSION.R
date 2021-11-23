### Function to use susie for expression prediction

train_susie <- function(X,
                        Y,
                        verbose,
                        par = F,
                        n.cores = NULL,
                        seed = 1218){
    
    set.seed(seed)
    train.folds = caret::createFolds(1:length(Y),
                                     k = 10,
                                     returnTrain = T)
    
    pred = vector(mode = 'numeric',
                  length = length(Y))
    for (tr in 1:length(train.folds)){
        
        Y.tr = Y[train.folds[[tr]]]
        X.tr = X[train.folds[[tr]],]
        X.test = X[-train.folds[[tr]],]
        susie.fit = susieR::susie(Matrix::Matrix(X.tr,sparse=T),
                                  Y.tr,
                                  estimate_residual_variance = TRUE,
                                  estimate_prior_variance = FALSE,
                                  scaled_prior_variance = 0.1,
                                  verbose = verbose)
            
        B.cur = as.numeric(coef(susie.fit))[-1]
        pred[-train.folds[[tr]]] = as.numeric(X.test %*% B.cur)
        
    }
    
    r2.vec = summary(lm(Y~pred))$adj.r
    susie.fit = susieR::susie(Matrix::Matrix(X,sparse=T),
                              Y,
                              estimate_residual_variance = TRUE,
                              estimate_prior_variance = FALSE,
                              scaled_prior_variance = 0.1,
                              verbose = verbose)
    
    
    return(list(Coef = coef(susie.fit)[-1],
                R2 = r2.vec))
    
    
}

require(data.table)
require(stringr)
AFR = fread('AFR.list',header=F)
EUR = fread('EUR.list',header=F)
list_of_tissues = c("Muscle_Skeletal",
                    "Whole_Blood",
                    "Artery_Tibial",
                    "Skin_Sun_Exposed_Lower_leg",
                    "Adipose_Subcutaneous")

require(bigsnpr)
### GRAB GTEX GENOTYPES
genotypes = snp_attach(snp_readBed2('GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bed',
                                    backingfile=tempfile()))

afr.geno = snp_attach(subset(genotypes,
                             ind.row = which(genotypes$fam$family.ID %in%
                                                 AFR$V1),
                             backingfile=tempfile()))
euro.geno = snp_attach(subset(genotypes,
                              ind.row = which(genotypes$fam$family.ID %in%
                                                 EUR$V1),
                              backingfile=tempfile()))
rm(genotypes)

tissue = list_of_tissues[i]
exp = fread(paste0(tissue,'.v8.normalized_expression.bed.gz'))
cov = fread(paste0(tissue,'.v8.covariates.txt.covar'))


## Regress out covariates from expression
require(limma)
exp = as.data.frame(exp)
exp[,5:ncol(exp)] = removeBatchEffect(as.matrix(exp[,5:ncol(exp)]),
                                      covariates = as.matrix(cov[,3:ncol(cov)]))


afr.exp = exp[,colnames(exp) %in% c('#chr',
                                    'start',
                                    'end',
                                    'gene_id',
                                    str_replace_all(AFR$V1,'[.]','-'))]
eur.exp = exp[,colnames(exp) %in% c('#chr',
                                    'start',
                                    'end',
                                    'gene_id',
                                    str_replace_all(EUR$V1,'[.]','-'))]
rm(exp)

afr.geno = snp_attach(subset(afr.geno,
                             ind.row = which(str_replace_all(afr.geno$fam$family.ID,
                                                             '[.]','-') %in%
                                                 colnames(afr.exp)),
                             backingfile=tempfile()))
euro.geno = snp_attach(subset(euro.geno,
                             ind.row = which(str_replace_all(euro.geno$fam$family.ID,
                                                             '[.]','-') %in%
                                                 colnames(eur.exp)),
                             backingfile=tempfile()))
afr.exp = afr.exp[,colnames(afr.exp) %in% c('#chr',
                                    'start',
                                    'end',
                                    'gene_id',
                                    str_replace_all(afr.geno$fam$family.ID,'[.]','-'))]
eur.exp = eur.exp[,colnames(eur.exp) %in% c('#chr',
                                            'start',
                                            'end',
                                            'gene_id',
                                            str_replace_all(euro.geno$fam$family.ID,'[.]','-'))]


for (j in 1:nrow(afr.exp)){
    
    gene = afr.exp$gene_id[j]
    print(gene)
    chr = strsplit(afr.exp$`#chr`[j],'r')[[1]][2]
    start = afr.exp$start[j]
    end = afr.exp$end[j]
    afr.pheno = as.numeric(afr.exp[j,5:ncol(afr.exp)])
    euro.pheno = as.numeric(eur.exp[j,5:ncol(eur.exp)])
    
    cur.afr.geno = snp_attach(subset(afr.geno,
                                     ind.col = which(afr.geno$map$chromosome == chr &
                                                         afr.geno$map$physical.pos < end + 1e6 &
                                                         afr.geno$map$physical.pos > start - 1e6),
                                     backingfile=tempfile()))
    
    cur.euro.geno = snp_attach(subset(euro.geno,
                                      ind.col = which(euro.geno$map$chromosome == chr &
                                                          euro.geno$map$physical.pos < end + 1e6 &
                                                          euro.geno$map$physical.pos > start - 1e6),
                                      backingfile=tempfile()))
    
    keep.euro = which(apply(cur.euro.geno$genotypes[],2,sd) != 0)
    keep.afr = which(apply(cur.afr.geno$genotypes[],2,sd) != 0)
    cur.afr.geno = snp_attach(subset(afr.geno,
                                     ind.col = intersect(keep.euro,
                                                         keep.afr),
                                     backingfile=tempfile()))
    cur.euro.geno = snp_attach(subset(euro.geno,
                                      ind.col = intersect(keep.euro,
                                                          keep.afr),
                                      backingfile=tempfile()))
    
    
    
    ### ELASTIC NET
    require(glmnet)
    library(doParallel)
    registerDoParallel(6)
    enet.mod.afr = cv.glmnet(x = cur.afr.geno$genotypes[],
                             y = afr.pheno,
                             keep = T,
                             nfolds = length(afr.pheno),
                             trace.it = F,
                             grouped = F,
                             intercept = F,
                             alpha = .5,
                             parallel = T)
    afr.r2.cv.enet = summary(lm(afr.pheno~
                                    enet.mod.afr$fit.preval[,
                                                            which.min(enet.mod.afr$cvm)]))$adj.r
    
    enet.mod.euro = cv.glmnet(x = cur.euro.geno$genotypes[],
                             y = euro.pheno,
                             keep = T,
                             nfolds = length(euro.pheno),
                             trace.it = F,
                             grouped = F,
                             intercept = F,
                             alpha = .5,
                             parallel=T)
    euro.r2.cv.enet = summary(lm(euro.pheno~
                                    enet.mod.euro$fit.preval[,
                                                             which.min(enet.mod.euro$cvm)]))$adj.r
    
    
    euro_in_afr.enet = as.numeric(coef(enet.mod.euro,
                                       s='lambda.min')[-1] %*% 
                                      t(cur.afr.geno$genotypes[]))
    
    afr_in_euro.enet = as.numeric(coef(enet.mod.afr,
                                       s='lambda.min')[-1] %*% 
                                      t(cur.euro.geno$genotypes[]))
    
    euro_in_afr.enet.R2 = summary(lm(afr.pheno~
                                         euro_in_afr.enet))$adj.r
    afr_in_euro.enet.R2 = summary(lm(euro.pheno~
                                         afr_in_euro.enet))$adj.r
    
    ### LASSO
    lasso.mod.afr = cv.glmnet(x = cur.afr.geno$genotypes[],
                             y = afr.pheno,
                             keep = T,
                             nfolds = length(afr.pheno),
                             trace.it = F,
                             grouped = F,
                             intercept = F,
                             alpha = 1,
                             parallel=T)
    afr.r2.cv.lasso = summary(lm(afr.pheno~
                                    lasso.mod.afr$fit.preval[,
                                                            which.min(lasso.mod.afr$cvm)]))$adj.r
    
    lasso.mod.euro = cv.glmnet(x = cur.euro.geno$genotypes[],
                              y = euro.pheno,
                              keep = T,
                              nfolds = length(euro.pheno),
                              trace.it = F,
                              grouped = F,
                              intercept = F,
                              alpha = 1,
                              parallel=T)
    euro.r2.cv.lasso = summary(lm(euro.pheno~
                                     lasso.mod.euro$fit.preval[,
                                                              which.min(lasso.mod.euro$cvm)]))$adj.r
    
    
    euro_in_afr.lasso = as.numeric(coef(lasso.mod.euro,
                                       s='lambda.min')[-1] %*% 
                                      t(cur.afr.geno$genotypes[]))
    
    afr_in_euro.lasso = as.numeric(coef(lasso.mod.afr,
                                       s='lambda.min')[-1] %*% 
                                      t(cur.euro.geno$genotypes[]))
    
    euro_in_afr.lasso.R2 = summary(lm(afr.pheno~
                                         euro_in_afr.lasso))$adj.r
    afr_in_euro.lasso.R2 = summary(lm(euro.pheno~
                                         afr_in_euro.lasso))$adj.r
    
    ### LMM
    lmm.mod.afr = cv.glmnet(x = cur.afr.geno$genotypes[],
                            y = afr.pheno,
                            keep = T,
                            nfolds = length(afr.pheno),
                            trace.it = F,
                            grouped = F,
                            intercept = F,
                            alpha = 0,
                            parallel=T)
    afr.r2.cv.lmm = summary(lm(afr.pheno~
                                     lmm.mod.afr$fit.preval[,
                                                              which.min(lmm.mod.afr$cvm)]))$adj.r
    
    lmm.mod.euro = cv.glmnet(x = cur.euro.geno$genotypes[],
                             y = euro.pheno,
                             keep = T,
                             nfolds = length(euro.pheno),
                             trace.it = F,
                             grouped = F,
                             intercept = F,
                             alpha = 0,
                             parallel=T)
    euro.r2.cv.lmm = summary(lm(euro.pheno~
                                      lmm.mod.euro$fit.preval[,
                                                                which.min(lmm.mod.euro$cvm)]))$adj.r
    
    
    euro_in_afr.lmm = as.numeric(coef(lmm.mod.euro,
                                        s='lambda.min')[-1] %*% 
                                       t(cur.afr.geno$genotypes[]))
    
    afr_in_euro.lmm = as.numeric(coef(lmm.mod.afr,
                                        s='lambda.min')[-1] %*% 
                                       t(cur.euro.geno$genotypes[]))
    
    euro_in_afr.lmm.R2 = summary(lm(afr.pheno~
                                          euro_in_afr.lmm))$adj.r
    afr_in_euro.lmm.R2 = summary(lm(euro.pheno~
                                          afr_in_euro.lmm))$adj.r
    
    ### SUSIE
    susie.mod.afr = train_susie(X = cur.afr.geno$genotypes[],
                                Y = afr.pheno,
                                verbose = F,
                                par = F,
                                n.cores = NULL,
                                seed = 1218)
    
    susie.mod.euro = train_susie(X = cur.euro.geno$genotypes[],
                                 Y = euro.pheno,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 seed = 1218)
    
    euro_in_afr.susie = as.numeric(susie.mod.euro$Coef %*% 
                                     t(cur.afr.geno$genotypes[]))
    
    afr_in_euro.susie = as.numeric(susie.mod.afr$Coef %*% 
                                     t(cur.euro.geno$genotypes[]))
    
    euro_in_afr.susie.R2 = summary(lm(afr.pheno~
                                        euro_in_afr.susie))$adj.r
    afr_in_euro.susie.R2 = summary(lm(euro.pheno~
                                        afr_in_euro.susie))$adj.r
    
    df = data.frame(Gene = gene,
                    Tissue = tissue,
                    EUR_in_EUR = max(c(euro.r2.cv.enet,
                                       euro.r2.cv.lasso,
                                       euro.r2.cv.lmm,
                                       susie.mod.euro$R2)),
                    AFR_in_AFR = max(c(afr.r2.cv.enet,
                                       afr.r2.cv.lasso,
                                       afr.r2.cv.lmm,
                                       susie.mod.afr$R2)),
                    EUR_in_AFR = max(c(euro_in_afr.enet.R2,
                                       euro_in_afr.lasso.R2,
                                       euro_in_afr.lmm.R2,
                                       euro_in_afr.susie.R2)),
                    AFR_in_EUR = max(c(afr_in_euro.enet.R2,
                                       afr_in_euro.lasso.R2,
                                       afr_in_euro.lmm.R2,
                                       afr_in_euro.susie.R2)))
    fwrite(df,
           'GTEX_cross_ancestry.tsv',
           append=T,
           row.names=F,
           sep='\t',
           quote=F)
    
}