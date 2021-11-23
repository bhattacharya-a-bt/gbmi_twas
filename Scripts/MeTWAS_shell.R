dataFolder='/u/scratch/a/abtbhatt/GTEX_Gene_Gene/data'
setwd(file.path(dataFolder,tissue))

if (file.exists('.RData')){
  file.remove('.RData')
}

require(bigsnpr)
require(MOSTWAS)

snpObj = snp_attach(snp_readBed2('qtltools_euro.bed',
                                 backingfile = tempfile()))
mediator = data.table::fread('eqtl_euro_exp.tsv')
colnames(mediator)[1] = 'Mediator'
medLocs = data.table::fread('genelocs.tsv')
medLocs$chr = sapply(strsplit(medLocs$chr,'r'),
                     function(x) x[2])
covariates = data.table::fread('eqtl_euro_cov.tsv')
dimNumeric = nrow(covariates)
if (file.exists('euro_trans_med_to_gene.tsv')){
  qtlFull = data.table::fread('euro_trans_med_to_gene.tsv')
  colnames(qtlFull) = c('SNP','gene',
                        'beta','t-stat',
                        'p-value','FDR')
  if (!dir.exists('MeTWASModels')){
    dir.create('MeTWASModels')
  }

  if (!file.exists('MeTWAS_done_genes.txt')){
    data.table::fwrite(data.frame(Gene = 'Hello'),
                       'MeTWAS_done_genes.txt',
                       quote=F,col.names=T,row.names=F,sep='\t')
  }

  runWithTemp <- function(gene){

    if (!dir.exists(paste0(gene,'_temp_bigsnp'))){
      dir.create(paste0(gene,'_temp_bigsnp'))
    }

    if (!dir.exists('MeTWASModels/')){
      dir.create('MeTWASModels/')
    }

    if (!paste0(gene,'.wgt.med.RData') %in% list.files('MeTWASModels/')){
      if (!gene %in% data.table::fread('MeTWAS_done_genes.txt')$Gene){


    MOSTWAS::MeTWAS(geneInt = gene,
           snpObj = snpObj,
           mediator = mediator,
           medLocs = medLocs,
           covariates = covariates,
           dimNumeric = dimNumeric,
           qtlFull = qtlFull,
           h2Pcutoff = .05,
           numMed = 10,
           seed = 1218,
           k = 5,
           cisDist = 1e6,
           parallel = T,
           prune = F,
           ldThresh = .5,
           cores = 5,
           verbose = T,
           R2Cutoff = 0.01,
           modelDir = 'MeTWASModels/',
           tempFolder = paste0(gene,'_temp_bigsnp/'))

      }
    }

    if (dir.exists(paste0(gene,'_temp_bigsnp'))){
      system(paste0('rm -r ',paste0(gene,'_temp_bigsnp/')))
    }

    if (any(grepl('core.',list.files()))){
      file.remove(list.files()[grepl('core.',list.files())])
    }

    data.table::fwrite(data.frame(Gene = gene),
                       'MeTWAS_done_genes.txt',
                       quote=F,append=T,row.names=F,sep='\t')

  }

  lapply(mediator$Mediator,FUN = runWithTemp)
}

