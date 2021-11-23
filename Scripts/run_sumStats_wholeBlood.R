common_models = intersect(list.files('AFR'),
                          list.files('EUR'))

sumStatFile = 'BioMe.Lee.Simon.2.ASTHMA.ALL.AFR.1636.5642.SAIGE.20200730.out.munged.MAC.20.0.INFO.0.3.gz.gnomad_v3_b38_ref_afr.postGWASQC.gz'
sumStatFilePath = file.path(sumStatFolder,
                            sumStatFile)


population = ifelse(grepl('.AFR.',sumStatFile),
                    'AFR','EUR')
refFile = file.path(refFolder,
                    paste0(population,'_all_1000G.rds'))
require(rtracklayer)
require(liftOver)
ch = import.chain(file.path(refFolder,
                            'hg19ToHg38.over.chain'))

require(data.table)
require(bigsnpr)
reference = snp_attach(refFile)
sumStat = fread(sumStatFilePath)
sumStat$ChrPos = paste(sumStat$`#CHR`,
                       sumStat$POS,
                       sep=':')
outFile = file.path(outFolder,
                    paste0('TWAS_',population,'_',
                           sumStatFile,'_sumStats.tsv'))

for (m in common_models){
    print(m)
    load(file.path(population,
                   m))
    if (exists('Model_AFR')){
        
        Model = Model_AFR
        rm(Model_AFR)
        
    } 
    if (exists('Model_EUR')){
        
        Model = Model_EUR
        rm(Model_EUR)
        
    }
    
    if (nrow(Model) > 0){
        
        Model$ChrPos = paste(Model$Chromosome,
                             Model$Position,sep=':')
        sumStat.cur = subset(sumStat,
                             ChrPos %in% Model$ChrPos)
        ref.cur = snp_attach(subset(reference,
                                    ind.col = which(reference$map$marker.ID %in%
                                                        Model$SNP),
                                    backingfile = tempfile()))
        Model = Model[match(ref.cur$map$marker.ID,
                            Model$SNP),]
        ref.cur$map$ChrPos = Model$ChrPos
        
        
        theseSNPs = intersect(intersect(sumStat.cur$ChrPos,
                                        Model$ChrPos),
                              ref.cur$map$ChrPos)
        
        Model = subset(Model,ChrPos %in% theseSNPs)
        sumStat.cur = subset(sumStat.cur, ChrPos %in% theseSNPs)
        ref.cur = snp_attach(subset(ref.cur,
                                    ind.col = which(ref.cur$map$ChrPos %in%
                                                        theseSNPs),
                                    backingfile = tempfile()))
        sumStat.cur = sumStat.cur[match(ref.cur$map$ChrPos,
                                        sumStat.cur$ChrPos),]
        Model = Model[match(ref.cur$map$ChrPos,
                            Model$ChrPos),]
        
        LD = tcrossprod(t(ref.cur$genotypes[]))/(nrow(ref.cur$fam)-1)
        sumStat.cur$Z = sumStat.cur$BETA/sumStat.cur$SE
        
        Z.top = (as.vector(Model$Effect) %*% as.vector(sumStat.cur$Z))
        Z.bottom = (t(as.vector(Model$Effect)) %*% LD %*% 
                        as.vector(Model$Effect))
        Z = Z.top/sqrt(Z.bottom)
        P = 2*pnorm(-abs(Z))
        
        df = data.frame(Tissue = 'Whole Blood',
                        SummaryStats = sumStatFile,
                        Population = population,
                        Gene = strsplit(m,'[.]')[[1]][1],
                        Z = Z,
                        P = P)
        fwrite(df,
               outFile,
               sep='\t',
               append=T,
               row.names=F,
               quote=F)
        
    }
    
    
    
}