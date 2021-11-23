
list_of_tissues = unique(list.files(dataFolder))
tissue = list_of_tissues[k]
setwd(file.path(dataFolder,tissue))


if (file.exists('.RData')){
  file.remove('.RData')
}

require(bigsnpr)
require(MOSTWAS)
require(data.table)

snpObj = snp_attach(snp_readBed2('qtltools_euro.bed',
                                 backingfile = tempfile()))
mediator = data.table::fread('eqtl_euro_exp.tsv')
colnames(mediator)[1] = 'Mediator'
medLocs = data.table::fread('genelocs.tsv')
medLocs$chr = sapply(strsplit(medLocs$chr,'r'),
                     function(x) x[2])
covariates = data.table::fread('eqtl_euro_cov.tsv')
dimNumeric = nrow(covariates)


if (file.exists('trans_eqtl_tsv')){
  qtlTra = data.table::fread('trans_eqtl_tsv')
  qtMed = data.table::fread('cis_eqtl.tsv')

  qtlTra_parts = paste0('fold',c(1,2,3),'_trans_eqtl_tsv')
  qtMed_parts = paste0('fold',c(1,2,3),'_cis_eqtl.tsv')

  if (!dir.exists('DePMAModels')){
    dir.create('DePMAModels')
  }

  if (!file.exists('DePMA_done_genes.txt')){
    data.table::fwrite(data.frame(Gene = 'Hello'),
                       'DePMA_done_genes.txt',
                       quote=F,col.names=T,row.names=F,sep='\t')
  }

  runWithTemp <- function(gene){

    if (!dir.exists(paste0(gene,'_temp_bigsnp'))){
      dir.create(paste0(gene,'_temp_bigsnp'))
    }

    if (!dir.exists('DePMAModels/')){
      dir.create('DePMAModels/')
    }

    if (!paste0(gene,'.wgt.med.RData') %in% list.files('DePMAModels/')){
      if (!gene %in% data.table::fread('DePMA_done_genes.txt')$Gene){


        DePMA(geneInt = gene,
              snpObj,
              mediator,
              medLocs,
              covariates,
              cisDist = 1e6,
              qtlTra,
              qtMed,
              h2Pcutoff = .05,
              dimNumeric = dimNumeric,
              verbose = T,
              seed = 1218,
              sobel = F,
              nperms = 1000,
              k = 5,
              parallel = F,
              parType = 'no',
              prune = F,
              ldThresh = .5,
              cores = 5,
              qtlTra_parts,
              qtMed_parts,
              modelDir = 'DePMAModels/',
              tempFolder = paste0(gene,'_temp_bigsnp/'),
              R2Cutoff = 0.01)

      }
    }

    if (dir.exists(paste0(gene,'_temp_bigsnp'))){
      system(paste0('rm -r ',paste0(gene,'_temp_bigsnp/')))
    }

    if (any(grepl('core.',list.files()))){
      file.remove(list.files()[grepl('core.',list.files())])
    }

    data.table::fwrite(data.frame(Gene = gene),
                       'DePMA_done_genes.txt',
                       quote=F,append=T,row.names=F,sep='\t')

  }

  lapply(mediator$Mediator,FUN = runWithTemp)
}

