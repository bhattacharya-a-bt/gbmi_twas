train_susie <- function(X,
                        Y,
                        verbose,
                        par = F,
                        n.cores = NULL,
                        seed = 1218){
    
    set.seed(seed)
    train.folds = caret::createFolds(1:length(Y),
                                     k = 5,
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


### REGRESS OUT COV FROM EXP
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


require(caret)
folds = createFolds(as.factor(c(rep('AFR',nrow(afr.geno$fam)),
                                rep('EUR',nrow(euro.geno$fam)))),
                    list = T,returnTrain = F,k=5)
folds.id = ifelse(1:(nrow(afr.geno$fam) + nrow(euro.geno$fam)) %in%
                      folds[[1]],1,
                  ifelse(1:(nrow(afr.geno$fam) + nrow(euro.geno$fam)) %in%
                             folds[[2]],2,
                         ifelse(1:(nrow(afr.geno$fam) + nrow(euro.geno$fam)) %in%
                                    folds[[3]],3,
                                ifelse(1:(nrow(afr.geno$fam) + nrow(euro.geno$fam)) %in%
                                           folds[[4]],4,5))))
                                


for (j in 1:nrow(afr.exp)){
    
    res = fread('GTEX_all_ancestry.tsv')
    gene = afr.exp$gene_id[j]
    if (!gene %in% subset(res,Tissue %in% tissue)$Gene){
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
        
        
        ppp = c(afr.pheno,euro.pheno)
        gen = rbind(cur.afr.geno$genotypes[],
                    cur.euro.geno$genotypes[])
        ### ELASTIC NET
        require(glmnet)
        enet.mod = cv.glmnet(x = gen,
                             y = ppp,
                             keep = T,
                             nfolds = 5,
                             trace.it = F,
                             grouped = F,
                             foldid = folds.id,
                             intercept = F,
                             alpha = .5)
        r2.cv.enet = summary(lm(ppp~
                                    enet.mod$fit.preval[,
                                                        which.min(enet.mod$cvm)]))$adj.r
        
        ### LASSO
        lasso.mod = cv.glmnet(x = gen,
                             y = ppp,
                             keep = T,
                             nfolds = 5,
                             trace.it = F,
                             grouped = F,
                             foldid = folds.id,
                             intercept = F,
                             alpha = .5)
        r2.cv.lasso = summary(lm(ppp~
                                    lasso.mod$fit.preval[,
                                                        which.min(lasso.mod$cvm)]))$adj.r
        
        ### LMM
        lmm.mod = cv.glmnet(x = gen,
                              y = ppp,
                              keep = T,
                              nfolds = 5,
                              trace.it = F,
                              grouped = F,
                              foldid = folds.id,
                              intercept = F,
                              alpha = .5)
        r2.cv.lmm = summary(lm(ppp~
                                     lmm.mod$fit.preval[,
                                                          which.min(lmm.mod$cvm)]))$adj.r
        
        ### SUSIE
        susie.mod = train_susie(X = gen,
                                Y = ppp,
                                    verbose = F,
                                    par = F,
                                    n.cores = NULL,
                                    seed = 1218)
        
        df = data.frame(Gene = gene,
                        Tissue = tissue,
                        Full_in_Full = max(c(r2.cv.lmm,r2.cv.enet,
                                             r2.cv.lasso,susie.mod$R2)))
        fwrite(df,
               'GTEX_all_ancestry.tsv',
               append=T,
               row.names=F,
               sep='\t',
               quote=F)}
    
}