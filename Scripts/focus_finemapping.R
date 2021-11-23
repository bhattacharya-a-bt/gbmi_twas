require(data.table)
require(bigsnpr)
source('FOCUS_functions.R')
pvalue_threshold = 2.5e-6
#summaryStats = fread('GBMI_POAG_02252021_chr1.txt')
snps = 
    snp_attach('GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.rds')
### HAD TO CONVERT SHEET 1 OF THE EXCEL SHEET TO A CSV.
### ALTERNATIVELY IF THIS WORKS FOR YOU, YOU CAN USE THIS
### twas_summary = xlsx::read.xlsx('JTI_POAG_chr1.xlsx',sheetName = 'Alltissues')

tissues = unique(twas_summary$Tissues)

### THE NEXT TWO LINES ARE NECESSARY BECAUSE BREAST IS MISLABELLED
#jti_tissues = list.files(weights_folder)
#tissues = tissues[paste0('xt_',tissues,'.txt') %in% jti_tissues]

for (t in tissues[1:49]){
    jtiModel = fread(file.path(weights_folder,paste0('xt_',t,'.txt')))
    res = subset(twas_summary, pvalue < pvalue_threshold & 
                     Tissues == t)
    un_res = subset(twas_summary, pvalue >= pvalue_threshold & 
                        Tissues == t)
    if (nrow(res) == 1){
        res$PIP = 1
        res$in_cred_set = T
        un_res$PIP = 0
        un_res$in_cred_set = F
        fwrite(rbind(res,un_res),outFile,append=T,
               row.names=F,quote=F,sep='\t')
    } else {
    res = res[order(res$Chr,res$Start),]
    chr.table = table(res$Chr)
    chr.un = unique(res$Chr)
    keep.genes = c()
    for (c in chr.un){
        res.cur = subset(res,Chr == c)
        res.cur = res.cur[order(res.cur$Start),]
        for (i in 1:(nrow(res.cur)-1)){
            if (res.cur$End[i] > res.cur$Start[i+1] - 1e6){
                keep.genes = unique(c(keep.genes,
                                      c(res.cur$gene[c(i,i+1)])))
            }
        }
    }
    
    res = subset(res,gene %in% keep.genes)
    un_res = rbind(un_res,subset(res,!gene %in% keep.genes))
    un_res$PIP = ifelse(un_res$pvalue < pvalue_threshold,1,0)
    un_res$in_cred_set = ifelse(un_res$pvalue < pvalue_threshold,T,F)
    
    chr.un = unique(res$Chr)
    if (length(chr.un) >= 1){
        for (c in chr.un){
            all.snps = c()
            omega = c()
            gene = c()
            snp.chr = c()
            this.res.tot = subset(res,Chr == c)
            this.res.tot$Group = 1
            ggg = 1
            for (i in 1:(nrow(this.res.tot)-1)){
                if (this.res.tot$End[i] <= this.res.tot$Start[i+1] - 1e6){
                    ggg = ggg+1
                    this.res.tot$Group[(i+1):nrow(this.res.tot)] = ggg
                }
            }
            for (g in unique(this.res.tot$Group)){
                this.res = subset(this.res.tot,
                                  Group == g)
                this.res = this.res[order(this.res$pvalue,
                                          decreasing = F),]
                for (i in 1:nrow(this.res)){
                    Model = subset(jtiModel,
                                   gene == this.res$gene[i])
                        if (nrow(Model) == 0){
                            Model = data.frame(SNP = c(),
                                               Effect = c())
                        } else {
                            Model = data.frame(SNP = Model$rsid,
                                               Effect = Model$weight)}
                    
                    if (length(as.character(Model$SNP)) != 
                        length(as.numeric(Model$Effect))){
                        print(this.res$Gene[i])
                        print(Model)
                    }
                    all.snps = c(all.snps,
                                 as.character(Model$SNP))
                    omega = c(omega,
                              as.numeric(Model$Effect))
                    gene = c(gene,
                             rep(this.res$gene[i],nrow(Model)))
                }
                tot.df = data.frame(SNP = all.snps,
                                    Gene = gene,
                                    Effect = omega,
                                    Chromosome = c)
                model.df = as.data.frame(matrix(nrow = length(unique(all.snps)),
                                                ncol = nrow(this.res)+1))
                colnames(model.df) = c('SNP',this.res$gene)
                model.df$SNP = as.character(unique(all.snps))
                for (j in 1:nrow(this.res)){
                    print(this.res$gene[j])
                    cur.tot.df = subset(tot.df,Gene == this.res$gene[j])
                    cur.tot.df$SNP = as.character(cur.tot.df$SNP)
                    for (i in 1:nrow(model.df)){
                        w = which(cur.tot.df$SNP == model.df$SNP[i])
                        model.df[i,j+1] = ifelse(length(w) != 0,
                                                 cur.tot.df$Effect[w],
                                                 0)
                    }
                }
                model.df$Chromosome = c
                for (i in 1:nrow(model.df)){
                    rrr = subset(tot.df,SNP == model.df$SNP[i])
                    model.df$Chromosome[i] = rrr$Chromosome[1]
                }
                snp.set = subset(snps$map,
                                 marker.ID %in% model.df$SNP)
                model.df = model.df[match(snp.set$marker.ID,
                                          model.df$SNP),]
                snpMat = as.matrix(snp_attach(subset(snps,
                                                     ind.col = 
                                                         which(snps$map$marker.ID %in% 
                                                                   model.df$SNP),
                                                     backingfile=tempfile()))$genotypes[])
                rm(snp.set)
                
                require(Matrix)
                snpMat = Matrix(t(snpMat), sparse = T)
                V = tcrossprod(snpMat) / ncol(snpMat)
                if (any(model.df$Chromosome != c)){
                    V[model.df$Chromosome == c,model.df$Chromosome != c] = 0
                    V[model.df$Chromosome != c,model.df$Chromosome == c] = 0}
                Omega = Matrix(as.matrix(model.df[,-c(1,ncol(model.df))]))
                zscores = this.res$zscore
                m = length(zscores)
                omega_nonzero = which(colSums(Omega) != 0)
                genes_nonzero = this.res$gene[omega_nonzero]
                Omega = as.matrix(Omega[,omega_nonzero])
                colnames(Omega) = genes_nonzero
                if (ncol(Omega) == 1){
                    this.res$PIP = ifelse(this.res$gene %in% 
                                              colnames(Omega),
                                          1,0)
                    this.res$in_cred_set = ifelse(this.res$gene %in% 
                                                      colnames(Omega),
                                                  T,F)
                    this.res = this.res[,-which(colnames(this.res) == 'Group')]
                    fwrite(this.res,outFile,append=T,
                           row.names= F,quote=F,sep='\t')
                } else {
                wcor = estimate_cor(Omega,V,intercept=T)[[1]]
                diag(wcor) = 1
                wcor[is.na(wcor)] = 0
                swld = estimate_cor(Omega,V,intercept=T)[[2]]
                null_res = m * log(1 - 1e-3)
                marginal = m * log(1 - 1e-3)
                comb_list = list()
                for (n in 1:min(3,length(zscores))){
                    comb_list = c(comb_list,
                                  combn(1:length(zscores),n,simplify=F))
                }
                pips = rep(0,length(zscores))
                zscores = get_resid(zscores,swld,wcor)[[1]]
                for (j in 1:length(comb_list)){
                    subset = comb_list[[j]]
                    local = bayes_factor(zscores,
                                         idx_set = subset,
                                         wcor = wcor)
                    marginal = log(exp(local) + exp(marginal))
                    for (idx in subset){
                        if (pips[idx] == 0){
                            pips[idx] = local
                        } else {
                            pips[idx] = log(exp(pips[idx]) + exp(local))
                        }
                    }
                    print(pips)
                    print(marginal)
                }
                
                pips = exp(pips - marginal)
                null_res = exp(null_res - marginal)
                this.res$PIP = pips
                this.res = this.res[order(this.res$PIP,decreasing = T),]
                npost = this.res$PIP/sum(this.res$PIP)
                csum = cumsum(npost)
                this.res$in_cred_set = csum < .9
                this.res = this.res[,-which(colnames(this.res) == 'Group')]
                fwrite(this.res,outFile,append=T,
                       row.names= F,quote=F,sep='\t')}}
        }
    }
    fwrite(un_res,outFile,append=T,
           row.names=F,quote=F,sep='\t')
    }
}