#,===============================================================================
#' Description:
#' Pre-processing of poppy objects, transcript and AmpliconArchitect and 
#' AmpliconClassifier outputs to generate a RData and reproduce the figures
#'
#===============================================================================
#' Article: Cascade 
#===============================================================================
#' Author: 
#' - Thaidy Moreno-Rodriguez
#===============================================================================
#' Date:  May 2023
#===============================================================================
#****************************************************************************
#' Script to process and to generate the data for figures and downstream analysis
#****************************************************************************
################################################################################

source('functions.R') # at scripts folder
# Optional: load the environment
# This environment is the result of processing this script.
load('Reproduce_Cascade_paper_environment.RData')  # at reproduce folder

#===============================================================================
#'  Loading libraries and metadata files
#===============================================================================

library( utils )
library( dplyr )
library( data.table )
library( foreach )
library( rowr )
library( tidyr )
library( sjmisc )
library( stringr )
library( tibble )
library( biomaRt )
library( readxl )
library( GSVA )
library( ConsensusClusterPlus )
library( aardvark ) # available at David's github (https://github.com/DavidQuigley/aardvark)
library( poppy ) # available at David's github (https://github.com/DavidQuigley/poppy)
library( deconstructSigs )
library( BSgenome.Hsapiens.UCSC.hg38 )
library( StructuralVariantAnnotation )
library( TxDb.Hsapiens.UCSC.hg38.knownGene )
library( rlang )
library( CHORD )
library( mutSigExtractor )
library( GenomicRanges )
library( IRanges )
library( Gviz )
library( corrplot)
# Mutation signatures

#===============================================================================
#'  Loading cascade and wcdt poppy objects (it will take 5-10 minutes)
#===============================================================================
dir_cascade_base =("~/Marlowe/projects/human_prostate_CASCADE/")
dir_wcdt_base=("~/Marlowe/projects/WCDT_WGS_2021/variant_cnv_sv_calls/")

dir_cascade_results = paste0(dir_cascade_base,'reproduce/results/')

load(paste0(dir_cascade_results,'2022_10_31_CASCADE_poppy_results_integrated.Rdata')) # latest poppy object with copy number weighed calls and tgt data
load(paste0(dir_wcdt_base,'2022_11_04_WCDT_poppy.Rdata')) #latest poppy object with copy number weighed calls for WCDT

# Delete sample "CA093-02_WGS" from cascade poppy object due low purity and failed QC by Hartwig tools
if('CA093-02_WGS' %in% names(SO)) SO <- SO[ - which(names(SO) == "CA093-02_WGS")]

#===============================================================================
#'  Loading cascade and wcdt TPMs and counts 
#===============================================================================

# TPMs and counts produce by pseudo-alingment by Kallisto
# and processed to remove  mitochondrial and ribosomal RNAs to then 
# re-calculated the TPMs for both cohorts

load(paste0(dir_cascade_results,'tx_mod.Rdata')) #latest poppy object with copy number weighed calls and tgt data
load(paste0(dir_cascade_results,'tx_wcdt_mod.Rdata')) #latest poppy object with copy number weighed calls for WCDT

dir_cascade_metadata = paste0(dir_cascade_base,'reproduce/publication_reproduce/metadata/')

sbs.grch38 <- read.table(paste0(dir_cascade_metadata,'annotations/2021_05_12_COSMIC_signatures.txt'), header = TRUE, sep = "\t") # at metadata/annotations 
Cascade_key_annot <- read.table(paste0(dir_cascade_metadata,'sample_data/Cascade_patients_annotations.txt'), sep = "\t", header = TRUE, stringsAsFactors = FALSE) # at metadata/sample_data 
Genes_by_Cat_subset<- read.table(paste0(dir_cascade_metadata,'annotations/Genes_grouped_subset.txt'), sep = "\t", header = T, stringsAsFactors = FALSE) # at metadata/annotations 

# Obtained from https://www.cancerhotspots.org/ and manually curated
# Single residue and in-frame indel mutation hotspots identified in 24,592 tumor samples 
# by the algorithm described in [Chang et al. 2017] and [Chang et al. 2016]

Hotspot<- read.table(paste0(dir_cascade_metadata,'annotations/hotspots.txt'), sep = "\t", header = T, stringsAsFactors = FALSE) # at metadata/annotations 
AA_annot<- read.table(paste0(dir_cascade_metadata,'annotations/Aminoacid_annotation.txt'), sep = "\t", header = T, stringsAsFactors = FALSE) # at metadata/annotations 
ETS_family <- read.table(paste0(dir_cascade_metadata,'annotations/ETS.txt'), sep = "\t", header = T, stringsAsFactors = FALSE) # at metadata/annotations 

#===============================================================================
#'  Calculating mutational signatures using Cosmic v2 and WGS data
#===============================================================================

# modified the mut signature file downladed from Cosmic to match deconstructSigs 
# format, where rows are signatures, columns are trinucleotide contexts

sbs.grch38_mod<- as.data.frame(t(sbs.grch38[,4:33] ))
colnames(sbs.grch38_mod) <- sbs.grch38[,3]
  
lt_wgs_strelka_mutect<- list()
for (i in 1:length(names(SO))){
  lt_wgs_strelka_mutect[[i]] <- (SO[[i]]$strelka_mutect) 
} 
names(lt_wgs_strelka_mutect) <- gsub("_WGS", "", names(SO))
lt_wgs_strelka_mutect<-rbindlist(lt_wgs_strelka_mutect, use.names=TRUE, fill=TRUE, idcol="sample_id")
sample_id <- do.call('rbind',strsplit(as.character(names(SO)), '_' ,fixed=TRUE))[,1]
patients_id <- do.call('rbind',strsplit(as.character(names(SO)), '-' ,fixed=TRUE))[,1]

Cascade.sigs.input.strelka.mutect.wgs <- mut.to.sigs.input(lt_wgs_strelka_mutect, 
                                                           sample.id = "sample_id", 
                                                           chr = "chrom", 
                                                           pos = "pos", 
                                                           ref = "ref", 
                                                           alt = "alt", 
                                                           bsg = BSgenome.Hsapiens.UCSC.hg38)
s.m.results.wgs.G.cosmic <- vector("list", nrow(Cascade.sigs.input.strelka.mutect.wgs))
names(s.m.results.wgs.G.cosmic) <- row.names(Cascade.sigs.input.strelka.mutect.wgs)
# run the estimation of exposures for each sample and save the results in the list
for( sID in row.names(Cascade.sigs.input.strelka.mutect.wgs) ){
  s.m.results.wgs.G.cosmic[[sID]] <- whichSignatures(Cascade.sigs.input.strelka.mutect.wgs, 
                                                     sample.id=sID, 
                                                     signatures.ref= sbs.grch38_mod, #can be use also the package's signature df "signatures.cosmic" with the exact same results
                                                     tri.counts.method="genome", 
                                                     contexts.needed=TRUE) 
}

Cascade.s.m.wgs.G.cosmic <- matrix(NA, ncol = 30, nrow = 54)
cascade_weigths_List.g <- list()
for (i in 1:length(names(s.m.results.wgs.G.cosmic))) {
  weights_Cascade.g<- s.m.results.wgs.G.cosmic[[i]]$weights
  cascade_weigths_List.g[[i]]<- weights_Cascade.g
}
Cascade.s.m.wgs.G.cosmic <- do.call(rbind, cascade_weigths_List.g)
colnames(Cascade.s.m.wgs.G.cosmic) <- gsub("Signature.", "Sig_", colnames(Cascade.s.m.wgs.G.cosmic) )

#===============================================================================
#'  Calculating chromothripsis
#===============================================================================

fn_chrom_lengths = paste(dir_cascade_base, 'reproduce/metadata/HG38/HG38_chromosome_lengths.txt',sep='') # at metadata/resources
fn_centromeres = paste(dir_cascade_base, 'reproduce/metadata/HG38/HG38_centromere_loci.txt',sep='') # at metadata/resources

chrom_lengths = read.table(fn_chrom_lengths,row.names=1, stringsAsFactors=FALSE)
chrom_names = rownames(chrom_lengths)
centromeres = read.table( fn_centromeres, header=TRUE, stringsAsFactors = FALSE)

list_sv_gridss<-data.frame()
list_sv_gridss <- foreach(i=names(SO),.combine = 'rbind') %do% {SO[[i]]$list_sv_gripss}

sample_ca_ids = names(SO)

# This step takes a long time (40 -120 minutes) 
# optional load the previous results 
# ca_chromo_scores<- read.table(paste0(dir_cascade_results,'CASCADE_chromo_scores_2500000.txt'),  sep = "",header = T)

for(n in 1:length(sample_ca_ids)){
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_ca_ids[i], "\n", sep=" ")
  for(i in 1:24){
    chrom = chrom_names[i]
    cs = chromo_score( sample_id = sample_ca_ids[n], chrom=chrom )
    if( n==1 & i == 1){
      ca_chromo_scores = cs
    }else{
      ca_chromo_scores = rbind(ca_chromo_scores, cs)
    }
  }
}

ca_chromo_scores$med_CNA = round( ca_chromo_scores$med_CNA, 3 )

n_samples = length(sample_ca_ids)
chrom_maxima = matrix(0, n_samples,24)
for(i in 1:length(sample_ca_ids)){
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_ca_ids[i], "\n", sep=" ")
  sc = ca_chromo_scores[ca_chromo_scores$sample_id == sample_ca_ids[i], ]
  print( paste("Loading chromothripsis scores:", sample_ca_ids[i] ))
  for(j in 1:24){
    chrom=rownames(chrom_lengths)[j]   
    sc_chrom = sc[sc$chrom==chrom,]
    for(x in 1:50){
      if( sum( sc_chrom$n_inv>=x & sc_chrom$n_del>=10 & 
               sc_chrom$n_cna_switch>=x )>0){
        chrom_maxima[i,j]=x   
      }else{
        break
      }
    }
  }
}

# filtering steps
chrom_maxima_filtered = chrom_maxima
chrom_maxima_filtered[chrom_maxima_filtered<15]=0
list_chromo=data.frame( which(chrom_maxima_filtered>0, arr.ind = TRUE) )
list_chromo[,1] = sample_ca_ids[list_chromo[,1]]
list_chromo$col = paste("chr", list_chromo$col, sep='')
list_chromo = list_chromo[order(list_chromo[,1]),]
names(list_chromo) = c("sample_id", "chrom")
list_chromo$chrom[list_chromo$chrom=="chr23"] = "chrX"

#===============================================================================
#' Summarizing SV (Structural Variants) called by GRIDSS/GRIPSS
#===============================================================================

SV_gridss_summary = data.frame()
hm_size = 3 # detecting deletions with microhomology of at least 3 bp 
for (i in 1:length(names(SO))){
  Inversion = SO[[i]]$gn_inv
  Breakend = SO[[i]]$gn_bnd
  Insertion = SO[[i]]$gn_ins
  Duplication = SO[[i]]$gn_dup
  DEL_vcf <-  info(SO[[i]]$gridss)[info(SO[[i]]$gridss)$EVENTTYPE == "DEL",]
  DEL <- length(unique(info(SO[[i]]$gridss)[info(SO[[i]]$gridss)$EVENTTYPE == "DEL",]$EVENT))
  DEL_df = data.frame(svtype = as.character(info(SO[[i]]$gridss)[info(SO[[i]]$gridss)$EVENTTYPE == "DEL",]$EVENTTYPE), 
                      sv_id = as.character(info(SO[[i]]$gridss)[info(SO[[i]]$gridss)$EVENTTYPE == "DEL",]$EVENT),
                      hom_len = as.integer(info(SO[[i]]$gridss)[info(SO[[i]]$gridss)$EVENTTYPE == "DEL",]$HOMLEN) )
  if (dim(DEL_df)[1] == 0){
    DEL_df_mh = DEL_df
    DEL_mh = 0
    DEL_nmh = 0 
  }
  else{
    DEL_df_mh <- unique(cbind(DEL_df[ which(DEL_df$hom_len >=hm_size),], SAMPLE_ID=(SO[[i]]$sample_base)))
    DEL_mh <- sum( DEL_df_mh$hom_len >=hm_size )
    DEL_nmh = DEL - DEL_mh
  }
  
  SV_SO_ID = cbind(SO[[i]]$sample_base, DEL_nmh, DEL_mh, Inversion ,Breakend, Insertion, Duplication ) #
  SV_gridss_summary = rbind(SV_gridss_summary, SV_SO_ID)
  
}
# SV_gridss_summary[2:7] <- lapply(SV_gridss_summary[2:7], as.numeric)
rownames(SV_gridss_summary) <- SV_gridss_summary$V1
SV_gridss_summary <- SV_gridss_summary[,-1]
SV_gridss_summary <- SV_gridss_summary %>% mutate_if(is.character, as.numeric)
SV_gridss_summary <- SV_gridss_summary[Cascade_key_annot$Sample_ID,]

#===============================================================================
#' Calculating of Homologous Recombination Deficiency using (CHORD)
#' Paper: Pan-cancer landscape of homologous recombination deficiency
#' Luan Nguyen, et al. Nat Commun 11, 5584 (2020). 
#' https://www.nature.com/articles/s41467-020-19406-4
#===============================================================================
# randomForest is required by CHORD
# install.packages('randomForest')
# library( CHORD )
# library( mutSigExtractor )

# This step takes over 20 minutes
# optional load the previous results 
# SO_chord_gridss_output<- read.delim(paste0(dir_cascade_results,'SO_chord_gridss_output.txt'), sep="\t") # at results folder

## Extract contexts for all samples using gridss output
chord_output_gridss <-list()
for(i in 1:length(names(SO))){
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", names(SO)[i], "\n", sep=" ")
  if (is.data.frame(SO[[i]]$list_sv_gripss) && nrow(SO[[i]]$list_sv_gripss) !=0){
    chord_output_gridss[[i]] <- extractSigsChord(
      df.snv = SO[[i]]$strelka_mutect[,1:4],
      #df.indel = is provided on the df.snv,
      df.sv = SO[[i]]$list_sv_gripss[,c(1,10)],
      sample.name=SO[[i]]$sample_base,
      sv.caller='gridss',
      ref.genome= BSgenome.Hsapiens.UCSC.hg38  )
  }
  else{
    next
  }
}

SO_gridss_merged_contexts <- do.call(rbind, chord_output_gridss)
SO_chord_gridss_output <- chordPredict(SO_gridss_merged_contexts, do.bootstrap=T, verbose=F)

#===============================================================================
#'  Create a matrix to annotate clinical/general information to plot in Figure 1A 
#===============================================================================

ETS_fusion<-ETS_family$Member
SO_summary = data.frame()

for(i in 1:length(names(SO)) ){
    patient_id = do.call('rbind',strsplit(as.character(SO[[i]]$sample_base), '-' ,fixed=TRUE))[,1]
    sample_id = SO[[i]]$sample_base
    purity_PURPLE= as.numeric(SO[[i]]$purity*100)
    ploidy_PURPLE= as.numeric(SO[[i]]$ploidy)
    tmbPerMb = as.numeric(SO[[i]]$tmbPerMb)
    tmbPerMb_Log2 = log2(as.numeric(tmbPerMb))
    SVtmbPerMb = as.numeric(SO[[i]]$svTumorMutationalBurden)
    WG_doubling = SO[[i]]$wholeGenomeDuplication
    WG_doubling <- gsub("true", "TRUE", WG_doubling)
    WG_doubling <- gsub("false", "FALSE", WG_doubling)
    ETS = any(SO[[i]]$fusions_linx$geneStart %in% ETS_fusion | SO[[i]]$fusions_linx$geneEnd %in% ETS_fusion)
    ## match info from Cascade_key_annot file
    m1 = match.idx(Cascade_key_annot$patient_id, patient_id)
    HRD_drug = Cascade_key_annot$HRD_drug[ m1$idx.A ]
    ASI_drug = Cascade_key_annot$ASI_drug[ m1$idx.A ]
   
    # Mutation signatures
    if (sample_id %in% rownames(Cascade.s.m.wgs.G.cosmic) ){
      m2 = match.idx(rownames(Cascade.s.m.wgs.G.cosmic), sample_id)
      Sig_1 = Cascade.s.m.wgs.G.cosmic$Sig_1[ m2$idx.A ]
      Sig_3 = Cascade.s.m.wgs.G.cosmic$Sig_3[ m2$idx.A ]
      Sig_4 = Cascade.s.m.wgs.G.cosmic$Sig_4[ m2$idx.A ]
      Sig_6 = Cascade.s.m.wgs.G.cosmic$Sig_6[ m2$idx.A ]
      Sig_8 = Cascade.s.m.wgs.G.cosmic$Sig_8[ m2$idx.A ]
      Sig_15 = Cascade.s.m.wgs.G.cosmic$Sig_15[ m2$idx.A ]
      Sig_18 = Cascade.s.m.wgs.G.cosmic$Sig_18[ m2$idx.A ]
      Sig_24 = Cascade.s.m.wgs.G.cosmic$Sig_24[ m2$idx.A ]
      Sig_25 = Cascade.s.m.wgs.G.cosmic$Sig_25[ m2$idx.A ]
      Sig_26 = Cascade.s.m.wgs.G.cosmic$Sig_26[ m2$idx.A ]
      Sig_29 = Cascade.s.m.wgs.G.cosmic$Sig_29[ m2$idx.A ]
    }
    else{
      Sig_1 = NA
      Sig_3 = NA
      Sig_4 = NA
      Sig_6 = NA
      Sig_8 = NA
      Sig_15 = NA
      Sig_18 = NA
      Sig_24 = NA
      Sig_25 = NA
      Sig_26 = NA
      Sig_29 = NA
    }
    m3 = match.idx(Cascade_key_annot$Sample_ID, sample_id)
    HR_status = Cascade_key_annot$HR_alt[ m3$idx.A ]
    if (sample_id %in% do.call('rbind',strsplit(as.character(list_chromo$sample_id), '_' ,fixed=TRUE))[,1] ){
      m4 = match.idx(do.call('rbind',strsplit(as.character(list_chromo$sample_id), '_' ,fixed=TRUE))[,1] , sample_id)
      if ( list_chromo$chrom[ m4$idx.A ] %in% c(paste0("chr", c(rep(1:22), "X") ))){
        Chromothripsis = "TRUE"
      }
      else{
        Chromothripsis = "FALSE"
      }
    }
    else{
      Chromothripsis = "FALSE"
    }
    if (sample_id %in% SO_chord_gridss_output$sample) {
      HRD_p_hrd = SO_chord_gridss_output[ SO_chord_gridss_output$sample == sample_id,]$p_hrd
      hr_status = SO_chord_gridss_output[ SO_chord_gridss_output$sample == sample_id,]$hr_status
    }
    else{
      HRD_p_hrd = 0
      hr_status = "cannot_be_determined"
    }
    m5 = match.idx(Cascade_key_annot$Sample_ID, sample_id)
    Tissue = as.character(Cascade_key_annot$Location_red[ m5$idx.A ])
    summary_row<-cbind(patient_id, sample_id, purity_PURPLE, ploidy_PURPLE, 
                       tmbPerMb, tmbPerMb_Log2, WG_doubling, ETS, SVtmbPerMb, 
                       Sig_1, Sig_3, Sig_4, Sig_6, Sig_8, Sig_15, Sig_18, Sig_24,
                       Sig_25,Sig_26,Sig_29,Chromothripsis, HRD_p_hrd, hr_status,
                       HR_status, Tissue, HRD_drug, ASI_drug)
    SO_summary = rbind(SO_summary, summary_row)

}
# Added sample CA093-09 fusion manually. 
# We confirm that it was filtered out by gridss but there is evidence of the fusion sharing the same BD as the other samples
SO_summary[SO_summary$sample_id == "CA093-09",]$ETS = TRUE

SO_summary[c(3:6,9:20,22)] <- lapply(SO_summary[c(3:6,9:20,22)], as.numeric)
SO_summary<- SO_summary[match(Cascade_key_annot$Sample_ID, SO_summary$sample_id),]

#===============================================================================
#'  Creating a matrix with germline and somatic alterations for CASCADE 
#'  patients to plot in Figure 1A 
#===============================================================================

load(paste0(dir_cascade_metadata,'resources/genome.Rdata'))

# First we add chromosome to the Copy number weighed data frame from poppy 
for (i in seq_len(length(SO)) ){
  if ( all(rownames(SO[[i]]$CNA_genes) %in% genome[["gene_locs"]]$symbol) ){
    m = match.idx( rownames(SO[[i]]$CNA_genes), genome[["gene_locs"]]$symbol )
    SO[[i]]$CNA_genes$chrom = genome[["gene_locs"]]$chrom[ m$idx.B ]
  }
  else{
    print(paste0( rownames(SO[[i]]$CNA_genes), " it's not among the genes from Gencode28_refFlat ..." ))
    next
  }
}

#===============================================================================
# Set sex chromosomes and use the thresholds to call copy number gain and loss 
# and biallelic loss status usign tumor purity, tumor ploidy. 
# Copy number < 0.5 = biallelic loss. 
# Sex chromosomes, copy gain, if weighted mean copy number > (tumor ploidy x 0.8), 
# and copy loss, if weighted mean copy number < (tumor ploidy x 0.3). 
# For autosomal, copy gain,  if the weighted mean copy number > (tumor ploidy x 1.95) 
# and copy loss if the weighted mean copy number was < (tumor ploidy x 0.67).

sex.chr <- c('chrX','chrY')
cnv.matrix <- list()
for(i in 1:length(names(SO))){
  sample_id <- names(SO)[i]
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing:", sample_id, "\n", sep=" ")
  cnv <- SO[[i]]$CNA_genes[c('chrom','weighted_mean')]
  cnv$gene <- rownames(SO[[i]]$CNA_genes)
  cnv.matrix[[sample_id]] <- cnv
  cnv.matrix[[i]]$ploidy <- SO[[i]]$ploidy
  
  cnv.matrix[[i]]$cnv.def <-
    ifelse(cnv.matrix[[i]]$chrom %in% sex.chr & cnv.matrix[[i]]$weighted_mean > cnv.matrix[[i]]$ploidy*0.8,'Copy gain',
           ifelse(cnv.matrix[[i]]$chrom %out% sex.chr &cnv.matrix[[i]]$weighted_mean > cnv.matrix[[i]]$ploidy*1.95,'Copy gain',
                  ifelse(cnv.matrix[[i]]$chrom %in% sex.chr & cnv.matrix[[i]]$weighted_mean < cnv.matrix[[i]]$ploidy* 0.3,'Copy loss',
                         ifelse(cnv.matrix[[i]]$chrom %out% sex.chr & cnv.matrix[[i]]$weighted_mean < cnv.matrix[[i]]$ploidy* 0.67,'Copy loss',''))))
  cnv.matrix[[i]]$Biallelic <- ifelse(cnv.matrix[[i]]$weighted_mean < 0.5,'Biallelic','') }


# Subset by the genes we previously select 
subset.CNA <- lapply(cnv.matrix,function(x) x[x$gene%in%Genes_by_Cat_subset$Genes,])
names(subset.CNA) <- names(cnv.matrix)


mat.CNA <- matrix(ncol = 1,nrow=length(Genes_by_Cat_subset$Genes))
mat.CNA[,1] = Genes_by_Cat_subset$Genes
colnames(mat.CNA) <- 'gene'
mat.CNA <- foreach(i=names(subset.CNA),.combine = 'cbind.fill') %do% {merge(mat.CNA,subset.CNA[[i]][,c(3,5:6)],by='gene')}
rownames(mat.CNA) = mat.CNA$gene;mat.CNA <-  mat.CNA[,c(F,T,T)];colnames(mat.CNA) = rep(names(subset.CNA),each=2)

#===============================================================================
# Adding AR.enhancer to the matrix
# AR enhancer region reported in cell Quigley 2018 and Nature genetics Zhao et al.

regions<- data.frame("chrom"="chrX","start"=66895158,"end"= 66910158, "subset.CNA" = "AR_enhancer")


AR_enhancer <- data.frame()
for (i in names(SO)){
  region_weighted <- somatic_add_CNA_by_region_from_segments( SO[[i]], regions)
  region_weighted <- cbind(sample_id = SO[[i]]$sample_id, region_weighted )
  AR_enhancer <- rbind( AR_enhancer, region_weighted)
}

AR_enhancer$cnv.def <- c("")
AR_enhancer$Biallelic <- c("")
for(i in 1:dim(AR_enhancer)[1]){
  ploidy <- SO[[i]]$ploidy
  AR_enhancer[i,]$cnv.def <-
    ifelse(AR_enhancer[i,]$chrom == "chrX" & AR_enhancer[i,]$weighted_mean >= ploidy*0.8,'Copy gain',
           ifelse(AR_enhancer[i,]$chrom == "chrX" & AR_enhancer[i,]$weighted_mean <= ploidy*0.3 , 'Copy loss','') )
  AR_enhancer[i,]$Biallelic <- ifelse(AR_enhancer[i,]$weighted_mean < 0.5,'Biallelic','')           
}

AR_enhancer_mod <- foreach(i=1:nrow(AR_enhancer),.combine = 'cbind.fill') %do% {cbind( AR_enhancer[i,c(6:7)])}
colnames(AR_enhancer_mod) <- rep(AR_enhancer$sample_id,each = 2)

# Combine the Copy number list of genes of interest and the AR enhancer copy numbers
mat.CNA_mod<- rbind(mat.CNA, "AR_enhancer" = AR_enhancer_mod)

#===============================================================================
##### Mutation matrix of the long tail gene list from WGS, WES and TGT sequencing data

subset.onco <- list()
for(i in 1:length(names(SO))){
  sample_id <- names(SO)[i]
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing sample :", sample_id, "\n", sep=" ")
  
  # First using wgs data   
  wgs <- data.frame()
  
  if(nrow(SO[[i]]$strelka_mutect) != 0 & any(SO[[i]]$strelka_mutect$gene %in% Genes_by_Cat_subset$Genes) ){
    
    genes <- SO[[i]]$strelka_mutect[SO[[i]]$strelka_mutect$gene %in% Genes_by_Cat_subset$Genes ,]$gene
    pos <- SO[[i]]$strelka_mutect[SO[[i]]$strelka_mutect$gene %in% Genes_by_Cat_subset$Genes ,]$pos
    
    df_genes_wgs <- as.data.frame(cbind(genes, pos))
    df_genes_wgs<- df_genes_wgs[!duplicated(df_genes_wgs),]
    df_genes_wgs$effects <- c("")
    df_genes_wgs$aa <- c("")
    
    for (n in 1:nrow(df_genes_wgs)){
      effects_wgs <- c()
      effects_wgs <- append(effects_wgs, SO[[i]]$strelka_mutect[SO[[i]]$strelka_mutect$gene == df_genes_wgs$genes[n] & SO[[i]]$strelka_mutect$pos == df_genes_wgs$pos[n],]$effect)
      effects_wgs <- append(effects_wgs, SO[[i]]$strelka[SO[[i]]$strelka$gene == df_genes_wgs$genes[n] & SO[[i]]$strelka$pos == df_genes_wgs$pos[n],]$effect)
      effects_wgs <- append(effects_wgs, SO[[i]]$mutect[SO[[i]]$mutect$gene == df_genes_wgs$genes[n] & SO[[i]]$mutect$pos == df_genes_wgs$pos[n],]$effect)
      effects_wgs <- unique(effects_wgs)
      aa_wgs <- c()
      aa_wgs <- append(aa_wgs, SO[[i]]$strelka_mutect[SO[[i]]$strelka_mutect$gene == df_genes_wgs$genes[n] & SO[[i]]$strelka_mutect$pos == df_genes_wgs$pos[n],]$aa)
      aa_wgs <- append(aa_wgs, SO[[i]]$strelka[SO[[i]]$strelka$gene == df_genes_wgs$genes[n] & SO[[i]]$strelka$pos == df_genes_wgs$pos[n],]$aa)
      aa_wgs <- append(aa_wgs, SO[[i]]$mutect[SO[[i]]$mutect$gene == df_genes_wgs$genes[n] & SO[[i]]$mutect$pos == df_genes_wgs$pos[n],]$aa)
      aa_wgs <- unique(aa_wgs)
      
      if( any( str_detect(effects_wgs, "missense|frame|stop|splice" ) ) ) {
        df_genes_wgs$effects[n] <- paste0(effects_wgs[str_detect(effects_wgs, "missense|frame|stop|splice" )], collapse=';')
        df_genes_wgs$aa[n] <- paste0(aa_wgs, collapse=';')
        wgs <- rbind(wgs, df_genes_wgs[n,])
      }
      else{
        next
      }
    }
    
    if (nrow(wgs) != 0 ) {
      subset.onco[[sample_id]] <- rbind(subset.onco[[sample_id]],wgs)
    }
    else{
      print(paste0("None mutation meet the criteria for wgs data in patient: ", sample_id))
    }
  }
  
  # Now using exome data   
  exome_df <- data.frame()
  
  if(nrow(SO[[i]]$strelka_mutect_exome) != 0 & any(SO[[i]]$strelka_mutect_exome$gene %in% Genes_by_Cat_subset$Genes) ){
    
    genes <- SO[[i]]$strelka_mutect_exome[SO[[i]]$strelka_mutect_exome$gene %in% Genes_by_Cat_subset$Genes ,]$gene
    pos <- SO[[i]]$strelka_mutect_exome[SO[[i]]$strelka_mutect_exome$gene %in% Genes_by_Cat_subset$Genes ,]$pos
    
    df_genes_exome <- as.data.frame(cbind(genes, pos))
    df_genes_exome<- df_genes_exome[!duplicated(df_genes_exome),]
    df_genes_exome$effects <- c("")
    df_genes_exome$aa <- c("")
    
    for (n in 1:nrow(df_genes_exome)){
      effects_exome <- c()
      effects_exome <- append(effects_exome, SO[[i]]$strelka_mutect_exome[SO[[i]]$strelka_mutect_exome$gene == df_genes_exome$genes[n] & SO[[i]]$strelka_mutect_exome$pos == df_genes_exome$pos[n],]$effect)
      effects_exome <- append(effects_exome, SO[[i]]$strelka_exome[SO[[i]]$strelka_exome$gene == df_genes_exome$genes[n] & SO[[i]]$strelka_exome$pos == df_genes_exome$pos[n],]$effect)
      effects_exome <- append(effects_exome, SO[[i]]$mutect_exome[SO[[i]]$mutect_exome$gene == df_genes_exome$genes[n] & SO[[i]]$mutect_exome$pos == df_genes_exome$pos[n],]$effect)
      effects_exome <- unique(effects_exome)
      aa_exome <- c()
      aa_exome <- append(aa_exome, SO[[i]]$strelka_mutect_exome[SO[[i]]$strelka_mutect_exome$gene == df_genes_exome$genes[n] & SO[[i]]$strelka_mutect_exome$pos == df_genes_exome$pos[n],]$aa)
      aa_exome <- append(aa_exome, SO[[i]]$strelka_exome[SO[[i]]$strelka_exome$gene == df_genes_exome$genes[n] & SO[[i]]$strelka_exome$pos == df_genes_exome$pos[n],]$aa)
      aa_exome <- append(aa_exome, SO[[i]]$mutect_exome[SO[[i]]$mutect_exome$gene == df_genes_exome$genes[n] & SO[[i]]$mutect_exome$pos == df_genes_exome$pos[n],]$aa)
      aa_exome <- unique(aa_exome)
      
      if( any( str_detect(effects_exome, "missense|frame|stop|splice" ) ) ) {
        df_genes_exome$effects[n] <- paste0(effects_exome[str_detect(effects_exome, "missense|frame|stop|splice" )], collapse=';')
        df_genes_exome$aa[n] <- paste0(aa_exome, collapse=';')
        exome_df <- rbind(exome_df, df_genes_exome[n,])
      }
      else{
        next
      }
    }
    
    if (nrow(exome_df) != 0 ) {
      subset.onco[[sample_id]] <- rbind(subset.onco[[sample_id]],exome_df)
    }
    else{
      print(paste0("None mutation meet the criteria for exome data in patient: ", sample_id))
    }
  }
  
  
  # Now using tgt data   
  tgt_df <- data.frame()
  
  if(nrow(SO[[i]]$strelka_mutect_tgt) != 0 & any(SO[[i]]$strelka_mutect_tgt$gene %in% Genes_by_Cat_subset$Genes) ){
    
    genes <- SO[[i]]$strelka_mutect_tgt[SO[[i]]$strelka_mutect_tgt$gene %in% Genes_by_Cat_subset$Genes ,]$gene
    pos <- SO[[i]]$strelka_mutect_tgt[SO[[i]]$strelka_mutect_tgt$gene %in% Genes_by_Cat_subset$Genes ,]$pos
    
    df_genes_tgt <- as.data.frame(cbind(genes, pos))
    df_genes_tgt<- df_genes_tgt[!duplicated(df_genes_tgt),]
    df_genes_tgt$effects <- c("")
    df_genes_tgt$aa <- c("")
    
    for (n in 1:nrow(df_genes_tgt)){
      effects_tgt <- c()
      effects_tgt <- append(effects_tgt, SO[[i]]$strelka_mutect_tgt[SO[[i]]$strelka_mutect_tgt$gene == df_genes_tgt$genes[n] & SO[[i]]$strelka_mutect_tgt$pos == df_genes_tgt$pos[n],]$effect)
      effects_tgt <- append(effects_tgt, SO[[i]]$strelka_tgt[SO[[i]]$strelka_tgt$gene == df_genes_tgt$genes[n] & SO[[i]]$strelka_tgt$pos == df_genes_tgt$pos[n],]$effect)
      effects_tgt <- append(effects_tgt, SO[[i]]$mutect_tgt[SO[[i]]$mutect_tgt$gene == df_genes_tgt$genes[n] & SO[[i]]$mutect_tgt$pos == df_genes_tgt$pos[n],]$effect)
      effects_tgt <- unique(effects_tgt)
      aa_tgt <- c()
      aa_tgt <- append(aa_tgt, SO[[i]]$strelka_mutect_tgt[SO[[i]]$strelka_mutect_tgt$gene == df_genes_tgt$genes[n] & SO[[i]]$strelka_mutect_tgt$pos == df_genes_tgt$pos[n],]$aa)
      aa_tgt <- append(aa_tgt, SO[[i]]$strelka_tgt[SO[[i]]$strelka_tgt$gene == df_genes_tgt$genes[n] & SO[[i]]$strelka_tgt$pos == df_genes_tgt$pos[n],]$aa)
      aa_tgt <- append(aa_tgt, SO[[i]]$mutect_tgt[SO[[i]]$mutect_tgt$gene == df_genes_tgt$genes[n] & SO[[i]]$mutect_tgt$pos == df_genes_tgt$pos[n],]$aa)
      aa_tgt <- unique(aa_tgt)
      
      if( any( str_detect(effects_tgt, "missense|frame|stop|splice" ) ) ) {
        df_genes_tgt$effects[n] <- paste0(effects_tgt[str_detect(effects_tgt, "missense|frame|stop|splice" )], collapse=';')
        df_genes_tgt$aa[n] <- paste0(aa_tgt, collapse=';')
        tgt_df <- rbind(tgt_df, df_genes_tgt[n,])
      }
      else{
        next
      }
    }
    
    if (nrow(tgt_df) != 0 ) {
      subset.onco[[sample_id]] <- rbind(subset.onco[[sample_id]],tgt_df)
    }
    else{
      print(paste0("None mutation meet the criteria for tgt data in patient: ", sample_id))
    }
  }
  
  
  rownames(subset.onco[[sample_id]])<- NULL
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Sample: ", sample_id, " has been processed.", "\n", sep=" ")
}

#===============================================================================
# Adding ZFHX3 mutation to sample CA113-20, manually curated and it looks real
# but was filter out from Strelka marked as LowEVS

subset.onco[["CA113-20_WGS"]] <- rbind(subset.onco[["CA113-20_WGS"]] , cbind("genes"="ZFHX3","pos"="72800097","effects"="frameshift_variant&splice_region_variant", "aa"= "p.Val1289fs;p.Val375fs"))
subset.onco_all <- lapply(subset.onco,function(x) x[!duplicated(x),])

subset.onco_all<-lapply(subset.onco_all, function(x) 
  cbind(x, aa_mod = (strsplit(as.character(x$aa), '.' ,fixed=TRUE) 
                     %>%  gsub("p.","", .)
                     %>%  gsub('\\"',"", .)
                     %>%  gsub(";","", .)
                     %>%  gsub("c\\(\\,","", .)
                     %>%  gsub("\\)","", .)
                     %>%  gsub("[[:space:]]","", .)
                     %>%  gsub("(^\\,)","", .)
  ) ) )

for (i in 1:length(subset.onco_all) ){
  subset.onco_all[[i]]$effects_mod <- c("")
  for (n in 1:nrow(subset.onco_all[[i]])){
    if(subset.onco_all[[i]][n,]$genes %in% Hotspot$Gene){
      aa_mod_sep <- as.data.frame(cbind(str_split_1(subset.onco_all[[i]][n,]$aa_mod,","))) #,simplify = TRUE)
      aa_mod_sep_2 <- aa_mod_sep  %>% separate(V1, into = c("AA_1", "num", "AA_2"),  sep = "(?=[A-Za-z])(?<=[0-9])|(?<=[0-9])(?=[_])|(?<=[0-9])(?=[*])|(?=[0-9])(?<=[A-Za-z])" ) %>% 
        mutate(AA_1 =str_replace_all(AA_1, setNames(AA_annot$abbreviation_1, AA_annot$abbreviation_3) ) )  %>%
        mutate(AA_2 =str_replace_all(AA_2, setNames(AA_annot$abbreviation_1, AA_annot$abbreviation_3) ) )  %>%
        mutate(aa_mod2 = paste0(AA_1, num ))
      for (m in 1:nrow(aa_mod_sep_2)){
        if (paste0(subset.onco_all[[i]][n,]$genes, aa_mod_sep_2$aa_mod2[m]) %in% paste0(Hotspot$Gene,Hotspot$Residue) ){
          subset.onco_all[[i]][n,] <- subset.onco_all[[i]][n,] %>% 
            mutate(across(everything(),~ gsub("_variant|conservative_", "", .)))  %>%
            mutate(effects_mod = paste0(str_split(effects, "&", simplify = T)[, 1],"_","hotspot") )
        }
        else {
          subset.onco_all[[i]][n,] <- subset.onco_all[[i]][n,] %>% 
            mutate(across(everything(),~ gsub("_variant|conservative_", "", .)))  %>%
            mutate(effects = paste0(str_split(effects, "&", simplify = T)[, 1],"_","VUS") )
        }
      }
    }
    else {
      subset.onco_all[[i]][n,] <- subset.onco_all[[i]][n,] %>% 
        mutate(across(everything(),~ gsub("_variant|conservative_", "", .)))  %>%
        mutate(effects = paste0(str_split(effects, "&", simplify = T)[, 1],"_","VUS") )
    }
  }
}

subset.onco_all<-lapply(subset.onco_all, function(x) x %>%
                          mutate(across(everything(),~ gsub("_VUS_VUS_VUS_VUS", "_VUS", .))) %>%
                          mutate(across(everything(),~ gsub("_VUS_VUS_VUS", "_VUS", .))) %>%
                          mutate(across(everything(),~ gsub("_VUS_VUS", "_VUS", .))) %>%
                          mutate(across(everything(),~ gsub("stop_gained_VUS", "Stop gained", .))) %>%
                          mutate(across(everything(),~ gsub("stop_gained_hotspot", "Stop gained hotspot", .))) %>%
                          mutate(across(everything(),~ gsub("frameshift_VUS", "Frameshift", .))) %>%
                          mutate(across(everything(),~ gsub("inframe_deletion_VUS", "Inframe", .))) %>%
                          mutate(across(everything(),~ gsub("inframe_deletion_hotspot", "Inframe hotspot", .))) %>%
                          mutate(across(everything(),~ gsub("inframe_insertion_VUS", "Inframe", .))) %>%
                          mutate(across(everything(),~ gsub("inframe_insertion_hotspot", "Inframe hotspot", .))) %>%
                          mutate(across(everything(),~ gsub("stop_lost_VUS", "Stop lost", .))) %>% 
                          mutate(across(everything(),~ gsub("splice_region_VUS", "Splice region", .))) %>%
                          mutate(across(everything(),~ gsub("stop_", "Stop ", .))) %>%
                          mutate(across(everything(),~ gsub("missense", "Missense", .))) %>%
                          mutate(across(everything(),~ gsub("_VUS", " VUS", .))) %>%
                          mutate(across(everything(),~ gsub("_hotspot", " hotspot", .))) ) 

tmp.df <- data.frame(genes=Genes_by_Cat_subset$Genes)

hotspot_onco <- foreach(i=names(subset.onco_all)) %do% { 
  merge(tmp.df,subset.onco_all[[i]],by='genes',all.x=T ) 
}
names(hotspot_onco) <- names(subset.onco_all)
hotspot_onco <- foreach(i=names(hotspot_onco)) %do% { hotspot_onco[[i]][,c(1,3,6)] }
names(hotspot_onco) <- names(subset.onco_all)

for (i in 1:length(names(hotspot_onco))){
  hotspot_onco[[i]][is.na(hotspot_onco[[i]])] <- ""
  for (n in 1:length(Genes_by_Cat_subset$Genes)){
    if( sum(hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n]) > 1){
      if( length(unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) > 1 & all(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod == "") ){
        hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0(paste0( unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects ), collapse = ";"), ";","Multi-hit"),
                                                                                                       dim( hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1])
      }
      if( length(unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) > 1 & any(str_detect(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod, pattern = "hotspot"))  ){
        hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0( unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n]& hotspot_onco[[i]]$effects_mod != "",]$effects_mod ), ";","Multi-hit"),
                                                                                                       dim( hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1])
      }
      if( length(unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) == 1 & all(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod == "") ){
        hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0( unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects ) , ";","Multi-hit"),
                                                                                                       dim( hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1])
      }
      if (length(unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) == 1 & any(str_detect(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod, pattern = "hotspot") ) ){
        hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0( unique(hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n] & hotspot_onco[[i]]$effects_mod != "",]$effects_mod ) , ";","Multi-hit"),
                                                                                                       dim( hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1] )
      }
    }
    
    else if (sum(hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n]) == 1 & hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod == ""){
      hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <- hotspot_onco[[i]][hotspot_onco[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects
    }
    else{
      next
    }
  }
}
hotspot_onco <- foreach(i=names(hotspot_onco)) %do% { hotspot_onco[[i]][,c(1,3)] }
names(hotspot_onco) <- names(subset.onco_all)
hotspot_onco <- lapply(hotspot_onco,function(x) x[!duplicated(x),])
hotspot_onco <- as.data.frame(do.call(cbind,hotspot_onco))

rownames(hotspot_onco) = hotspot_onco[,1];hotspot_onco <-  hotspot_onco[,c(F,T)];colnames(hotspot_onco) = names(subset.onco_all)
missing_cols <- setdiff(names(SO),colnames(hotspot_onco))
hotspot_onco[ , missing_cols] <- ""
hotspot_onco <- hotspot_onco[ , names(SO)]
hotspot_onco[is.na(hotspot_onco)] <- ""

#===============================================================================
# Creating matrix for structural variants

all_sv <- list()
for(n in 1:length(names(SO))){
  sample_id <- names(SO)[n]
  sv_by_gene <- data.frame()
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_id, "\n", sep=" ")
  if ( !is_empty(SO[[n]]$list_sv_gripss) ){
    gridss_gr <- GenomicRanges::makeGRangesFromDataFrame(SO[[n]]$list_sv_gripss,
                                                         keep.extra.columns=T,
                                                         ignore.strand=T,
                                                         seqinfo=NULL,
                                                         seqnames.field=c("chrom_1"),
                                                         start.field="start_1",
                                                         end.field=c("end_1", "stop") )
    for(i in 1:length(Genes_by_Cat_subset$Genes) ){
      gene_query <- Genes_by_Cat_subset$Genes[i]
      gr_gene = GenomicRanges::GRanges(seqnames=genome$gene_locs[genome$gene_locs$symbol == gene_query,]$chrom,
                                       ranges = IRanges(start=genome$gene_locs[genome$gene_locs$symbol == gene_query,]$start,
                                                        end=genome$gene_locs[genome$gene_locs$symbol == gene_query,]$end))
      df_sv_overlapping_gene = data.frame( subsetByOverlaps( gridss_gr, gr_gene) )
      if (dim(df_sv_overlapping_gene)[1]>0 ){
        sv_by_gene <- rbind(sv_by_gene, cbind( "gene" = gene_query, "SV"= "Structural variant") )
      }
      else{
        sv_by_gene <- rbind(sv_by_gene, cbind( "gene" = gene_query, "SV"= NA ) )
      }
    }
    all_sv[[sample_id]] <- rbind(sv_by_gene)
  }
  else{
    all_sv[[sample_id]] <- cbind("gene" =Genes_by_Cat_subset$Genes, "SV"= NA)
    next
  }
}

sv_onco <- do.call(cbind,all_sv)
rownames(sv_onco) = sv_onco[,1]
sv_onco <-  sv_onco[,c(F,T)];colnames(sv_onco) = names(SO)

#===============================================================================
# Fusions also annotated as Structural variant

all_fusion <- list()
for(n in 1:length(names(SO))){
  sample_id <- names(SO)[n]
  fusion_by_gene <- data.frame()
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_id, "\n", sep=" ")
  if ( !is_empty(SO[[n]]$fusions_linx) ){
    for(i in 1:length(Genes_by_Cat_subset$Genes) ){
      for(m in 1:dim(SO[[n]]$fusions_linx)[1] ){
        if (any(SO[[n]]$fusions_linx[m,]$geneStart %in% Genes_by_Cat_subset$Genes[i] | SO[[n]]$fusions_linx[m,]$geneEnd %in% Genes_by_Cat_subset$Genes[i]) & 
            SO[[n]]$fusions_linx[m,]$geneStart != SO[[n]]$fusions_linx[m,]$geneEnd ){
          fusion_by_gene <- rbind(fusion_by_gene, cbind( "gene" = Genes_by_Cat_subset$Genes[i], "Fusion"= "Structural variant") )
        }
        else{
          next
        }
      }
      fusion_by_gene <- rbind(fusion_by_gene, cbind( "gene" = Genes_by_Cat_subset$Genes[i], "Fusion"= NA ))
    }
    all_fusion[[sample_id]] <- rbind(fusion_by_gene)
    all_fusion[[sample_id]] <- all_fusion[[sample_id]] %>% 
      distinct(gene, .keep_all = T)
  }
  else{
    all_fusion[[sample_id]] <- cbind("gene" = Genes_by_Cat_subset$Genes, "Fusion"= NA)
    next
  }
}

fusion_onco <- do.call(cbind,all_fusion)
rownames(fusion_onco) = fusion_onco[,1];fusion_onco <-  fusion_onco[,c(F,T)]
colnames(fusion_onco) = names(SO)

#===============================================================================
# Germline alteration

germline_ca <- list()
for(i in 1:length(names(SO))){
  sample_id <- names(SO)[i]
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_id, "\n", sep=" ")
  germline <- subset(SO[[i]]$germline_pathogenic,SO[[i]]$germline_pathogenic$symbol %in% Genes_by_Cat_subset$Genes)[c('symbol','effect')]
  if (dim(germline)[1]>0){
    germline$germline <- "Germline"
    germline_ca[[sample_id]]$germline <- germline
  }
  else{
    germline_ca[[sample_id]]$germline<- germline %>% add_column("germline" = NA )
  }
}

foreach(i=names(germline_ca),.combine = 'cbind.fill') %do% { germline_ca[[i]]$germline$effect[ str_detect(germline_ca[[i]]$germline$effect, "missense|frame|stop|splice" ) ] <- NA}
mat.ger.df <- matrix(ncol = 1, nrow=length(Genes_by_Cat_subset$Genes));mat.ger.df[,1] = Genes_by_Cat_subset$Genes

mat.ger.df <- foreach(i=names(germline_ca)) %do% {
  onco.ger.df <- germline_ca[[i]]$germline[c('symbol','effect','germline')] %>% group_by(effect,symbol, germline) %>% dplyr::summarise() %>% as.data.frame() 
}


mat.ger.df <- lapply(mat.ger.df,function(x) x[!duplicated(x$symbol),])
names(mat.ger.df) <- names(germline_ca);tmp.df <- data.frame(gene=Genes_by_Cat_subset$Genes)
colnames(tmp.df)<- 'symbol'
mat.ger.onco <- foreach(i=names(mat.ger.df)) %do% { 
  merge(tmp.df,mat.ger.df[[i]][,c(1,2,3)],by='symbol',all.x=T ) %>% mutate(across(everything(),~ gsub("_variant", "", .)))  %>% mutate(effect = str_split(effect, "&", simplify = T)[, 1])
}
names(mat.ger.onco) <- names(germline_ca)

mat.ger.onco <- do.call(cbind,mat.ger.onco)

rownames(mat.ger.onco) = mat.ger.onco[,1]
mat.ger.onco <-  mat.ger.onco[,c(F,T,T)]
colnames(mat.ger.onco) = rep(names(SO),each=2)

#===============================================================================
# Gathering all the alterations in the same matrix to plot with complexheatmap

gene.idx <- Genes_by_Cat_subset$Genes
mut.matrix.all <- foreach(i=names(hotspot_onco),.combine = 'cbind') %do%{
  df <- cbind(hotspot_onco[gene.idx,][i],mat.ger.onco[gene.idx,c(T,F)][i],mat.ger.onco[gene.idx,c(F,T)][i],
              mat.CNA_mod[gene.idx,c(T,F)][i],mat.CNA_mod[gene.idx,c(F,T)][i],
              sv_onco[gene.idx,][i],fusion_onco[gene.idx,][i])
  colnames(df) <- c('ONCO','EF','GL','CNA','BI','SV','FUS')
  df[df == ""]<- NA
  df <- df %>% unite('test',sep = ';',remove = T,na.rm = T)
}
colnames(mut.matrix.all) <- colnames(hotspot_onco)
colnames(mut.matrix.all)<- gsub("_WGS", "", colnames(mut.matrix.all))
mat_onco=mut.matrix.all
mat_onco <- mat_onco[, Cascade_key_annot$Sample_ID]

# Added sample CA093-09 fusion manually. 
# We confirm that it was filtered out by gridss but there is evidence of the fusion 
# sharing the same BD as the other samples

mat_onco["ETV1","CA093-09"]="Structural variant"

#===============================================================================
#'  Reversion mutations post processing and plot Figure 1B
#===============================================================================

# AARDVARK is available at https://github.com/DavidQuigley/aardvark
# For more information, visit the repository or contact: david.quigley@ucsf.edu

load( paste0(dir_cascade_results,'2023_04_20_aardvark_CA071-GERM.Rdata' )) # at results folder
fns = paste0( dir_cascade_results,'2023_04_20_aardvark_CA071-0',1:9,'.Rdata' ) # at results folder
pathogenic_mut = aardvark::Mutation("chr13", 32339657, "CTT", "C", transcript_BRCA2)

sums_gl =  summarize_candidates(reads_gl, transcript_BRCA2, pathogenic_mut)
sums_all = list()
for(i in 1:9){
  print(i)
  load( fns[i] )
  sums_all = c(sums_all, list( summarize_candidates(reads, transcript_BRCA2, pathogenic_mut) ) )
}
names(sums_all) = paste0("CA071-0",1:9)

slugs_germline = dimnames(sums_gl$summary)[[1]]

# remove reversions found in germline
for(i in 1:length(sums_all)){
  idx_del = match.idx( slugs_germline, dimnames(sums_all[[i]]$summary)[[1]] )$idx.B
  idx_keep = setdiff( 1:dim(sums_all[[i]]$summary)[1], idx_del)
  sums_all[[i]]$summary = sums_all[[i]]$summary[ idx_keep,]
}

#===============================================================================
#'  Gene expression data from CASCADE to calculate AR and NE score
#'  using GSVA library and selecting genes up tp then plot in Figure 2A 
#===============================================================================

TPM <-as.matrix(tx_mod$abundance)
Log2.TPM  <- log2( 1+TPM )
Log2.TPM_ft <- Log2.TPM[, colnames(Log2.TPM) %in% do.call('rbind',strsplit(as.character(names(SO)), '_' ,fixed=TRUE))[,1] ] #to remove samples without matched WGS-seq 
Log2.TPM.zc<- t(scale(t(Log2.TPM_ft),scale = T,center = T))
Log2_TPM_zc_gsva <- Log2.TPM.zc

#===============================================================================
# The AR signalling genes comes from: https://www.sciencedirect.com/science/article/pii/S1535610806002820
# Hieronymus, H. et al. Gene expression signature–based chemical genomic prediction 
# identifies a novel class of HSP90 pathway modulators.  Cancer Cell 10, 321–330 (2006).
# A gene expression signature of androgen stimulation was defined from gene expression 
# profiles of LNCaP cells stimulated with the synthetic androgen R1881 for 12 hr and 24 hr, 
# as compared to androgen-deprived LNCaP cells. The 27 gene signature contains both 
# androgen-induced and androgen-repressed genes.
# Here, I'm using only 21 genes that are androgen-induced, as they do in the science paper:
# Tang F. et al. Chromatin profiles classify castration-resistant prostate cancers suggesting therapeutic targets 
# https://www.science.org/doi/10.1126/science.abe1505#supplementary-materials
# Modified gene due a change in the gene name:
# BM039	-> CENPN
# TMEPAI ->	PMEPA1
# SARG ->	C1orf116


ensembl_hg38 <- useMart(biomart='ENSEMBL_MART_ENSEMBL',
                        dataset='hsapiens_gene_ensembl', verbose = T)

# Optional load the supplementary table from Hieronymus's paper
# Hieronymus_AR_Score <- readxl::read_excel(paste0(dir_cascade_metadata,'external_resources/1-t3.0-S1535610806002820-mmc.xlsx'), sheet = 1, skip = 1) #at metadata/external_resources
# genes_Hieronymus_AR_Up <- Hieronymus_AR_Score[Hieronymus_AR_Score$Direction == "up",]$Symbol
genes_Hieronymus_AR_Up[str_detect(genes_Hieronymus_AR_Up, "BM039|TMEPAI|SARG") ]

genes_Hieronymus_AR_Up <- c("FKBP5", "KLK2", "ELL2", "TMPRSS2", "ZBTB10", "CENPN", "PMEPA1", 
                            "NKX3-1", "KLK3", "EAF2", "MED28", "ABCC4", "C1orf116", "GNMT", 
                            "MPHOSPH9", "NNMT", "ADAM7", "ACSL3", "MAF", "HERC3", "PTGER4" )

AR_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=genes_Hieronymus_AR_Up, mart=ensembl_hg38, useCache = FALSE)
AR_sigID <- AR_sigID[!duplicated(AR_sigID$external_gene_name),] #remove duplicate
row.names(AR_sigID) <- AR_sigID$external_gene_name


#===============================================================================
# The NE signalling genes comes from: https://www.nature.com/articles/nm.4045
# Beltran et al: Divergent clonal evolution of castration-resistant neuroendocrine prostate cancer
# doi:10.1038/nm.4045 in Nature Med
# Modified gene due a change in the gene name:
# C7orf76 -> SEM1

genes_Beltran_NE_Up <- c("ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9",
                         "BRINP1","SEM1","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1",
                         "PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP",
                         "SYT11","AURKA","DNMT1","EZH2","MYCN")

# Optional load the supplementary table from Beltran's paper
# Beltran_NE_Score <- readxl::read_excel(paste0(dir_cascade_metadata,'external_resources/41591_2016_BFnm4045_MOESM29_ESM.xlsx'), sheet = 10, skip = 1) #at metadata/external_resources
# modify gene name C7orf76 -> SEM1
# genes_Beltran_NE_Up <- Beltran_NE_Score[Beltran_NE_Score$`RNA (CRPC-NE vs CRPC-Adeno)` == "Over-expressed",]$`HGNC ID`
# Beltran_NE_Score[Beltran_NE_Score$`HGNC ID` == "C7orf76",]$`HGNC ID` <-  "SEM1"
# genes_Beltran_NE_Up <- Beltran_NE_Score_ft[Beltran_NE_Score_ft$`RNA (CRPC-NE vs CRPC-Adeno)` == "Over-expressed",]$`HGNC ID`

NEPC_sigID <- getBM(filters= "external_gene_name",
                    attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                    values=genes_Beltran_NE_Up, mart= ensembl_hg38, useCache = FALSE)
row.names(NEPC_sigID) <- NEPC_sigID$external_gene_name

#===============================================================================
# We used GSVA to transform the gene expression measurements into enrichment scores 
# for these two gene sets

row.names(Log2_TPM_zc_gsva)[which(row.names(Log2_TPM_zc_gsva) %in% AR_sigID$external_gene_name)] <- AR_sigID[row.names(Log2_TPM_zc_gsva)[which(row.names(Log2_TPM_zc_gsva) %in% AR_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(Log2_TPM_zc_gsva)[which(row.names(Log2_TPM_zc_gsva) %in% NEPC_sigID$external_gene_name)] <- NEPC_sigID[row.names(Log2_TPM_zc_gsva)[which(row.names(Log2_TPM_zc_gsva) %in% NEPC_sigID$external_gene_name)],'ensembl_gene_id_version']
gsvaScore_all <- gsva(as.matrix(Log2_TPM_zc_gsva), list(AR_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=5, max.sz=500, method="zscore")
rownames( gsvaScore_all ) <- c("AR_score","NEPC_score")

#===============================================================================
#' Clustering & Cluster assignment of the 5 mCRPC subtypes based on JCI paper
#' from (Labrecque et al) for Figure 2A & Supp Figure 3A
#===============================================================================

P.Nelson_Sig<- read.table(paste0(dir_cascade_metadata,'annotations/Pete_Nelson_Signature_all.txt'), header = T, sep = "\t") #at metadata/annotations

#TPMs from Kallisto Gene level
sample_ca_ids_rna <- Cascade_key_annot$Sample_ID[order(Cascade_key_annot$Sample_ID)]
Log2.TPM_ft <- Log2.TPM_ft[,colnames(Log2.TPM_ft) %in% sample_ca_ids_rna]
Log2.TPM_mod.zc<- t(scale(t(Log2.TPM_ft),scale = T,center = T))
Log2.TPM_mod.zc.nelson_ft <- as.matrix(Log2.TPM_mod.zc[rownames(Log2.TPM_mod.zc) %in% P.Nelson_Sig$Gene_ID, ])
Log2.TPM_mod.zc.nelson_ft_K = sweep(Log2.TPM_mod.zc.nelson_ft,1, apply(Log2.TPM_mod.zc.nelson_ft,1,median,na.rm=T))

results_CC_nelson_ft_K_ILC = ConsensusClusterPlus(Log2.TPM_mod.zc.nelson_ft_K, # normalized expression matrix
                                                  maxK=6, # maximum k to make clusters
                                                  reps=50, # resampling
                                                  pItem=0.8, # item resampling
                                                  pFeature=1, # gene resampling
                                                  title= "geneExp_nelson_ft_K", # output title 
                                                  clusterAlg="pam", # clustering algorithm
                                                  distance="euclidean", # distance measure
                                                  seed=29,# random number 39 
                                                  innerLinkage="complete", 
                                                  plot="png") # output file format

CC_Pnelson_mod<- as.data.frame(results_CC_nelson_ft_K_ILC[[5]][["consensusClass"]][names(results_CC_nelson_ft_K_ILC[[5]][["consensusClass"]]) %in% Cascade_key_annot$Sample_ID ] )
colnames(CC_Pnelson_mod)<- "Cluster_ID"
CC_Pnelson_mod
CC_Pnelson_mod$Cluster_ID <-
  ifelse(CC_Pnelson_mod$Cluster_ID == 1, 'AR+/NE-',
         ifelse(CC_Pnelson_mod$Cluster_ID == 2, 'ARL/NE-',
                ifelse(CC_Pnelson_mod$Cluster_ID == 3, 'AR-/NE+',
                       ifelse(CC_Pnelson_mod$Cluster_ID == 4, 'AR+/NE+', 'AR-/NE-') ) ) )

#===============================================================================
#' Create a matrix with clinical/general information to plot in Figure 2A and
#' Supplementary Figure 3A
#===============================================================================

cascade_gsvaScore <- as.data.frame(t(gsvaScore_all))
RNA_Cascade = data.frame()

for(i in 1:length(colnames(Log2.TPM_ft)) ){
  
  patient_id = do.call('rbind',strsplit(as.character(colnames(Log2.TPM_ft)[i]), '-' ,fixed=TRUE))[,1]
  sample_id = colnames(Log2.TPM_ft)[i]
  ## match info from Cascade_key_annot file
  m = match.idx(Cascade_key_annot$patient_id, patient_id)
  ASI_drug = Cascade_key_annot$ASI_drug[ m$idx.A ]
  m2 = match.idx(Cascade_key_annot$Sample_ID, sample_id)
  Tissue = as.character(Cascade_key_annot$Location_red[ m2$idx.A ])

  if (sample_id %in% rownames(cascade_gsvaScore) ) {
    m3 = match.idx(rownames(cascade_gsvaScore), sample_id)
    AR_score = as.numeric(cascade_gsvaScore$AR_score[ m3$idx.A ])
    m4 = match.idx(rownames(cascade_gsvaScore), sample_id)
    NE_score = as.numeric(cascade_gsvaScore$NEPC_score[ m4$idx.A ])
  }
  else{
    AR_score = NA 
    NE_score = NA
  }
  
  if (sample_id %in% rownames(CC_Pnelson_mod) ) {
    m5 = match.idx(rownames(CC_Pnelson_mod), sample_id)
    Cluster_ID = as.character(CC_Pnelson_mod$Cluster_ID[ m5$idx.A ])
  }
  else{
    Cluster_ID = NA 
  }
  
  if (sample_id %in% gsub("_WGS","",names(SO)) ) { # Adding cnvs for Suppl Figure 3A
    m6 = match.idx( gsub("_WGS","",names(SO)), sample_id)
    AR_cnv = as.numeric(SO[[m6$idx.A]]$CNA_genes["AR",]$weighted_mean)
    AR_Log2_cnv <- as.numeric(log2(SO[[m6$idx.A]]$CNA_genes["AR",]$weighted_mean))
    MYC_cnv <- as.numeric(SO[[m6$idx.A]]$CNA_genes["MYC",]$weighted_mean)
    MYC_Log2_cnv <- as.numeric(log2(SO[[m6$idx.A]]$CNA_genes["MYC",]$weighted_mean))
    PCAT1_cnv <- as.numeric(SO[[m6$idx.A]]$CNA_genes["PCAT1",]$weighted_mean)
    PCAT1_Log2_cnv <- as.numeric(log2(SO[[m6$idx.A]]$CNA_genes["PCAT1",]$weighted_mean))
  }
  else{
    next
  }
  summary_row<-cbind( patient_id, sample_id,Cluster_ID,Tissue, ASI_drug, AR_score, 
                     NE_score, AR_cnv, AR_Log2_cnv, MYC_cnv, MYC_Log2_cnv, 
                     PCAT1_cnv, PCAT1_Log2_cnv ) 
  RNA_Cascade = rbind( RNA_Cascade, summary_row )
  
}

RNA_Cascade[c(6:13)] <- lapply(RNA_Cascade[c(6:13)], as.numeric)


idx_keep = which( P.Nelson_Sig$Gene_ID %in% rownames(Log2.TPM_mod.zc.nelson_ft) )
P.Nelson_Sig_ft = P.Nelson_Sig[idx_keep,]$Gene_ID
# Ordering the columns and rows to plot with ComplexHeatmap "Figure 2A"
Log2.TPM_mod.zc.nelson_ft <- Log2.TPM_mod.zc.nelson_ft[ P.Nelson_Sig_ft ,]
Log2.TPM_mod.zc.nelson_ft <- Log2.TPM_mod.zc.nelson_ft[, RNA_Cascade$sample_id]

#===============================================================================
#'  For figures 2 b,c,d,and e!
#'  Those plots zoom in ChrX at AR locus to plot copy numbers segments as lines
#'  to compare the different tumors from the same patients: 
#'  2b:CA071, 2c:CA104, 2d:CA113 and 2e:CA108.
#===============================================================================
# No need to pre-process any data. Use the R script Figures_2b_e.R and this
# script will reproduce all the individuals figures and save them as pdf

#===============================================================================
#'  Figure 3a: Frequency and distribution of complex event in CASCADE cohort 
#'             This data frames will be used also for Supp Fig 3A
#===============================================================================

# Use the outputs from AmpliconArchitect and AmpliconClassifier and save as RData
# object. Includes the list of files processed, the genes amplified per sample
# (all types of complex events and linear amplifications), the intervals 
# of ecDNAs and BFBs events, and the ecDNA counts per sample
# The script to gather the files is not included but is available upon request!

load(paste0(dir_cascade_results, 'cascade_complex_events.Rdata'))

# Table created for genes found amplified by ecDNA and the Pathway 
AA_genes_freq<- read.table(paste0(dir_cascade_metadata,'annotations/AA_genes_freq.txt'), sep = "\t", header = T) # at metadata/annotations

#===============================================================================
# Quantification of the events including those classified as ecDNA, BFB, Complex,
# unknown, and Linear amplifications. Linear amplifications were excluded from
# Figure 3a and Supplementary Figure 3a

cascade_events_mod <- data.frame()
samples_id = unique(beds_cascade_df$sample_id)

for (i in 1:length(samples_id) ) {
  sample_id = samples_id[i]
  features_unique <- unique( beds_cascade_df[str_detect(beds_cascade_df$sample_id, sample_id), ]$feature )
  ecDNA = sum(str_detect(features_unique, pattern = "ecDNA") )
  BFB = sum(str_detect(features_unique, pattern = "BFB") )
  Linear = sum(str_detect(features_unique, pattern = "Linear"))
  Unknown = sum(str_detect(features_unique, pattern = "unknown"))
  Complex = sum(str_detect(features_unique, pattern ="Complex"))
  n_event =  sum(ecDNA, BFB, Complex, Unknown)
  event <- cbind(sample_id, n_event, ecDNA, BFB, Complex, Unknown, Linear)
  cascade_events_mod<- rbind(cascade_events_mod,event)
}

missing_calls <- setdiff(do.call('rbind',strsplit(as.character(AA_counts_list_all_cascade$sample_id), '_' ,fixed=TRUE))[,1], cascade_events_mod$sample_id)

for (i in 1:length(missing_calls)){
  new_col = missing_calls[i]
  cascade_events_mod <- cascade_events_mod  %>% 
    mutate_at(grep("^(n_event|ecDNA|BFB|Linear|Unknown|Complex)",colnames(.)), as.numeric) %>%
    add_row( sample_id = new_col, n_event = 0, ecDNA = 0, BFB = 0, Complex = 0 , Unknown =0 , Linear =0 )
} 

cascade_events_mod<- cascade_events_mod[order(cascade_events_mod$n_event, decreasing = T),]
rownames(cascade_events_mod) <- cascade_events_mod$sample_id
genes_to_include <- c()
id_order_ca<- cascade_events_mod$sample_id

#===============================================================================
# Adding AR.enhancer amplifiedr or co-amplified by ecDNA to the matrix above.
# AR enhancer region reported in cell Quigley 2018 and Nature genetics Zhao et al.

# regions<- data.frame("chrom"="chrX","start"=66895158,"end"= 66910158, "subset.CNA" = "AR_enhancer")

# We used all the interval bed files to extract the infomation.
# we excluded also Linear amplifications
# The script to gather the files is not included but is available upon request!
# We data is included in the cascade_complex_events.Rdata

ls_events_int_ca
events_int_ca <- rbindlist(ls_events_int_ca, use.names = F, fill = F, idcol = T)
colnames(events_int_ca)<- c("sample_id","chrom", "start", "end")

events_AR_enh_ca<- annotate_enhancer_interval(events_int_ca, anno_region = regions )
events_AR_enh_ca_df <-as.data.frame(events_AR_enh_ca)

AR_enh_list_ca<- as.data.frame(cbind(gene = "AR_enhancer", str_split(events_AR_enh_ca_df$sample_id,"_", simplify = T)[, 1:3]))
colnames(AR_enh_list_ca) <- c("gene", "sample_name", "amplicon_number", "feature")

AR_enh_list_ca_mod<- AR_enh_list_ca %>% dplyr::select(sample_name, feature, gene) %>%  distinct() %>% pivot_wider(names_from = sample_name, values_from = feature)

#===============================================================================
# Creating a matrix to plot with complexheatmap's oncoPrint function
# with samples as columns and genes as rows

matrix_AA_gene_types_ca_list<- AA_gene_list_cascade %>% separate(gene, "gene") %>% mutate(feature= gsub('_[[:digit:]]+', '', feature)) %>% 
  dplyr::select(sample_name, feature, gene) %>%  distinct() %>% pivot_wider(names_from = sample_name, values_from = feature)

matrix_AA_gene_types_ca_df<- as.data.frame(matrix_AA_gene_types_ca_list)
missing_cols<- setdiff(colnames(matrix_AA_gene_types_ca_df), colnames(AR_enh_list_ca_mod))
AR_enh_list_ca_mod[missing_cols] <- ""

matrix_AA_gene_types_ca_df <- rbind(matrix_AA_gene_types_ca_df, AR_enh_list_ca_mod)
matrix_AA_gene_types_ca_df[is.na(matrix_AA_gene_types_ca_df)] = ""
rownames(matrix_AA_gene_types_ca_df) = matrix_AA_gene_types_ca_df[, 1]
matrix_AA_gene_types_ca_df = matrix_AA_gene_types_ca_df[, -1]
matrix_AA_gene_types_ca_df = t(as.matrix(matrix_AA_gene_types_ca_df))

matrix_AA_gene_types_ca_df[is.na(matrix_AA_gene_types_ca_df)] = NULL
matrix_AA_gene_types_ca_df<- matrix_AA_gene_types_ca_df[order(rownames(matrix_AA_gene_types_ca_df)),]

AA_gene_ca_freq_sub <- t(matrix_AA_gene_types_ca_df[, colnames(matrix_AA_gene_types_ca_df) %in% AA_genes_freq$Gene])
AA_gene_ca_freq_sub[1:3,1:3]
AA_gene_ca_freq_sub[is.na(AA_gene_ca_freq_sub)] = ""

for (i in 1:ncol(AA_gene_ca_freq_sub)) {
  for (n in 1:nrow(AA_gene_ca_freq_sub)){
    if(is_list(AA_gene_ca_freq_sub[n,i])){
      if( !is_empty(unlist(AA_gene_ca_freq_sub[n,i])) ){
        AA_gene_ca_freq_sub[n,i] <- unlist(AA_gene_ca_freq_sub[n,i])
      }
      else{
        AA_gene_ca_freq_sub[n,i] <- ""
        next
      }
    }
    else{
      next
    }
  }
  
  
}

missing_genes<- setdiff(AA_genes_freq$Gene, rownames(AA_gene_ca_freq_sub))
missing_genes_ca_DF<- as.data.frame(matrix(nrow= length(missing_genes), ncol = dim(AA_gene_ca_freq_sub)[2] ) )
rownames(missing_genes_ca_DF)<- missing_genes
colnames(missing_genes_ca_DF)<- colnames(AA_gene_ca_freq_sub)
AA_gene_ca_freq_sub_new <- rbind(AA_gene_ca_freq_sub, missing_genes_ca_DF)

missing_samples_ca <- setdiff(AA_counts_list_all_cascade$sample_id, colnames(AA_gene_ca_freq_sub_new))
missing_samples_ca_DF<- as.data.frame(matrix(ncol= length(missing_samples_ca), nrow = dim(AA_gene_ca_freq_sub_new)[1] ) )
colnames(missing_samples_ca_DF)<- missing_samples_ca
rownames(missing_samples_ca_DF)<- rownames(AA_gene_ca_freq_sub_new)
matrix_AA_gene_ca_freq_sub_new <- cbind(AA_gene_ca_freq_sub_new, missing_samples_ca_DF)

matrix_AA_gene_ca_freq_sub_new_mod <- apply(matrix_AA_gene_ca_freq_sub_new, 2, function(y) gsub("Linear", "", y ) )
rownames(matrix_AA_gene_ca_freq_sub_new_mod)<- rownames(matrix_AA_gene_ca_freq_sub_new)
matrix_AA_gene_ca_freq_sub_new_mod <- gsub("^NA", "", matrix_AA_gene_ca_freq_sub_new_mod)

split_ecDNA_ca_gene_list <- c()
for (i in 1:nrow(matrix_AA_gene_ca_freq_sub_new_mod)){
  
  m = match.idx(AA_genes_freq$Gene, rownames(matrix_AA_gene_ca_freq_sub_new_mod)[i])
  pathway = AA_genes_freq$Pathway[ m$idx.A ]
  split_ecDNA_ca_gene_list <- rbind(split_ecDNA_ca_gene_list, pathway)
  
}

matrix_AA_gene_ca_freq_sub_new_mod <- matrix_AA_gene_ca_freq_sub_new_mod[, id_order_ca]



#===============================================================================
# Creating a matrix of logical values for the presence/absence of any type of 
# alteration per sample. This matrix will be use for the statistics in Fig 4

all_alterations<- mat_onco %>% mutate_all(na_if,"") # Replace "" by "NA"
all_alterations <- all_alterations %>%                            
  replace(!is.na(.), TRUE) %>%  # Replace any alteration by TRUE
  replace(is.na(.), FALSE) 
all_alterations<- as.data.frame(t(all_alterations))


SO_summary_sub <- SO_summary[ match(id_order_ca , SO_summary$sample_id ),]
idx_pt = match.idx(SO_summary_sub$sample_id, rownames(all_alterations),  allow.multiple.B = F)
SO_summary_sub$TP53 = c("")
SO_summary_sub[idx_pt$idx.A,]$TP53 = all_alterations[idx_pt$idx.B, ]$TP53
SO_summary_sub$PTEN = c("")
SO_summary_sub[idx_pt$idx.A,]$PTEN = all_alterations[idx_pt$idx.B, ]$PTEN


#===============================================================================
#'  Figure 3b: Frequency and distribution of complex event in WCDT cohort
#===============================================================================
#
load(paste0(dir_cascade_results, 'wcdt_complex_events.Rdata'))

#===============================================================================
# Quantification of the events including those classified as ecDNA, BFB, Complex,
# unknown, and Linear amplifications. Linear amplifications were excluded from
# Figure 3b and Supplementary Figure 3b

wcdt_events_mod <- data.frame()
samples_id = unique(beds_wcdt_df$sample_id)

for (i in 1:length(samples_id) ) {
  sample_id = samples_id[i]
  features_unique <- unique( beds_wcdt_df[str_detect(beds_wcdt_df$sample_id, sample_id), ]$feature )
  ecDNA = sum(str_detect(features_unique, pattern = "ecDNA") )
  BFB = sum(str_detect(features_unique, pattern = "BFB") )
  Linear = sum(str_detect(features_unique, pattern = "Linear"))
  Unknown = sum(str_detect(features_unique, pattern = "unknown"))
  Complex = sum(str_detect(features_unique, pattern ="Complex"))
  n_event =  sum(ecDNA, BFB, Complex, Unknown)
  event <- cbind(sample_id, n_event, ecDNA, BFB, Complex, Unknown, Linear)
  wcdt_events_mod<- rbind(wcdt_events_mod,event)
}

missing_calls <- setdiff(do.call('rbind',strsplit(as.character(AA_counts_list_all_wcdt$sample_id), '_' ,fixed=TRUE))[,1], wcdt_events_mod$sample_id)

for (i in 1:length(missing_calls)){
  new_col = missing_calls[i]
  wcdt_events_mod <- wcdt_events_mod  %>% 
    mutate_at(grep("^(n_event|ecDNA|BFB|Linear|Unknown|Complex)",colnames(.)), as.numeric) %>%
    add_row( sample_id = new_col, n_event = 0, ecDNA = 0, BFB = 0, Complex = 0 , Unknown =0 , Linear =0 )
} 

wcdt_events_mod<- wcdt_events_mod[order(wcdt_events_mod$n_event, decreasing = T),]
rownames(wcdt_events_mod) <- wcdt_events_mod$sample_id
id_order_wcdt<- wcdt_events_mod$sample_id

#===============================================================================
# Adding AR.enhancer amplifiedr or co-amplified by ecDNA to the matrix above.
# AR enhancer region reported in cell Quigley 2018 and Nature genetics Zhao et al.

# regions<- data.frame("chrom"="chrX","start"=66895158,"end"= 66910158, "subset.CNA" = "AR_enhancer")

# We used all the interval bed files to extract the infomation.
# we excluded also Linear amplifications
# The script to gather the files is not included but is available upon request!
# We data is included in the wcdt_complex_events.Rdata

ls_events_int_wcdt
events_int_wcdt <- rbindlist(ls_events_int_wcdt, use.names = F, fill = F, idcol = T)
colnames(events_int_wcdt)<- c("sample_id","chrom", "start", "end")

events_AR_enh_wcdt<- annotate_enhancer_interval(events_int_wcdt, anno_region = regions )
events_AR_enh_wcdt_df <-as.data.frame(events_AR_enh_wcdt)

AR_enh_list_wcdt<- as.data.frame(cbind(gene = "AR_enhancer", str_split(events_AR_enh_wcdt_df$sample_id,"_", simplify = T)[, 1:3]))
colnames(AR_enh_list_wcdt) <- c("gene", "sample_name", "amplicon_number", "feature")

AR_enh_list_wcdt_mod<- AR_enh_list_wcdt %>% dplyr::select(sample_name, feature, gene) %>%  distinct() %>% pivot_wider(names_from = sample_name, values_from = feature)

#===============================================================================
# Creating a matrix to plot with complexheatmap's oncoPrint function
# with samples as columns and genes as rows

matrix_AA_gene_types_wcdt_list<- AA_gene_list_wcdt %>% separate(gene, "gene") %>% mutate(feature= gsub('_[[:digit:]]+', '', feature)) %>% 
  dplyr::select(sample_name, feature, gene) %>%  distinct() %>% pivot_wider(names_from = sample_name, values_from = feature)


matrix_AA_gene_types_wcdt_df<- as.data.frame(matrix_AA_gene_types_wcdt_list)
missing_cols<- setdiff(colnames(matrix_AA_gene_types_wcdt_df), colnames(AR_enh_list_wcdt_mod))
AR_enh_list_wcdt_mod[missing_cols] <- ""

matrix_AA_gene_types_wcdt_df <- rbind(matrix_AA_gene_types_wcdt_df, AR_enh_list_wcdt_mod)
matrix_AA_gene_types_wcdt_df[is.na(matrix_AA_gene_types_wcdt_df)] = ""
rownames(matrix_AA_gene_types_wcdt_df) = matrix_AA_gene_types_wcdt_df[, 1]
matrix_AA_gene_types_wcdt_df = matrix_AA_gene_types_wcdt_df[, -1]
matrix_AA_gene_types_wcdt_df = t(as.matrix(matrix_AA_gene_types_wcdt_df))

matrix_AA_gene_types_wcdt_df[is.na(matrix_AA_gene_types_wcdt_df)] = NULL
matrix_AA_gene_types_wcdt_df<- matrix_AA_gene_types_wcdt_df[order(rownames(matrix_AA_gene_types_wcdt_df)),]

AA_gene_wcdt_freq_sub <- t(matrix_AA_gene_types_wcdt_df[, colnames(matrix_AA_gene_types_wcdt_df) %in% AA_genes_freq$Gene])
AA_gene_wcdt_freq_sub[1:3,1:3]
AA_gene_wcdt_freq_sub[is.na(AA_gene_wcdt_freq_sub)] = ""

for (i in 1:ncol(AA_gene_wcdt_freq_sub)) {
  for (n in 1:nrow(AA_gene_wcdt_freq_sub)){
    if(is_list(AA_gene_wcdt_freq_sub[n,i])){
      if( !is_empty(unlist(AA_gene_wcdt_freq_sub[n,i])) ){
        AA_gene_wcdt_freq_sub[n,i] <- unlist(AA_gene_wcdt_freq_sub[n,i])
      }
      else{
        AA_gene_wcdt_freq_sub[n,i] <- ""
        next
      }
    }
    else{
      next
    }
  }
  
  
}

missing_genes<- setdiff(AA_genes_freq$Gene, rownames(AA_gene_wcdt_freq_sub))
missing_genes_wcdt_DF<- as.data.frame(matrix(nrow= length(missing_genes), ncol = dim(AA_gene_wcdt_freq_sub)[2] ) )
rownames(missing_genes_wcdt_DF)<- missing_genes
colnames(missing_genes_wcdt_DF)<- colnames(AA_gene_wcdt_freq_sub)
AA_gene_wcdt_freq_sub_new <- rbind(AA_gene_wcdt_freq_sub, missing_genes_wcdt_DF)

missing_samples_wcdt <- setdiff(AA_counts_list_all_wcdt$sample_id, colnames(AA_gene_wcdt_freq_sub_new))
missing_samples_wcdt_DF<- as.data.frame(matrix(ncol= length(missing_samples_wcdt), nrow = dim(AA_gene_wcdt_freq_sub_new)[1] ) )
colnames(missing_samples_wcdt_DF)<- missing_samples_wcdt
rownames(missing_samples_wcdt_DF)<- rownames(AA_gene_wcdt_freq_sub_new)
matrix_AA_gene_wcdt_freq_sub_new <- cbind(AA_gene_wcdt_freq_sub_new, missing_samples_wcdt_DF)

matrix_AA_gene_wcdt_freq_sub_new_mod <- apply(matrix_AA_gene_wcdt_freq_sub_new, 2, function(y) gsub("Linear", "", y ) )
rownames(matrix_AA_gene_wcdt_freq_sub_new_mod)<- rownames(matrix_AA_gene_wcdt_freq_sub_new)
matrix_AA_gene_wcdt_freq_sub_new_mod <- gsub("^NA", "", matrix_AA_gene_wcdt_freq_sub_new_mod)

split_ecDNA_wcdt_gene_list <- c()
for (i in 1:nrow(matrix_AA_gene_wcdt_freq_sub_new_mod)){
  
  m = match.idx(AA_genes_freq$Gene, rownames(matrix_AA_gene_wcdt_freq_sub_new_mod)[i])
  pathway = AA_genes_freq$Pathway[ m$idx.A ]
  split_ecDNA_wcdt_gene_list <- rbind(split_ecDNA_wcdt_gene_list, pathway)
  
}

matrix_AA_gene_wcdt_freq_sub_new_mod <- matrix_AA_gene_wcdt_freq_sub_new_mod[, id_order_wcdt]

#===============================================================================
# Creating a matrix of genomic alteration to then transform it to a matrix of 
# logical values for the presence/absence of any type of alteration per sample. 
# This matrix will be use for the statistics in Fig 4

# First we add chromosome to the Copy number weighed data frame from poppy 
for (i in seq_len(length(SO_WCDT_WGS)) ){
  if ( all(rownames(SO_WCDT_WGS[[i]]$CNA_genes) %in% genome[["gene_locs"]]$symbol) ){
    m = match.idx( rownames(SO_WCDT_WGS[[i]]$CNA_genes), genome[["gene_locs"]]$symbol )
    SO_WCDT_WGS[[i]]$CNA_genes$chrom = genome[["gene_locs"]]$chrom[ m$idx.B ]
  }
  else{
    print(paste0( rownames(SO_WCDT_WGS[[i]]$CNA_genes), " it's not among the genes from Gencode28_refFlat ..." ))
    next
  }
}
#===============================================================================
# Using the same threshold and setting sex chromosomes as with cascade cohort
# to call copy number gain, loss, and biallelic loss
# and biallelic loss status using tumor purity, tumor ploidy. 

sex.chr <- c('chrX','chrY')
cnv_wcdt.matrix <- list()
for(i in 1:length(names(SO_WCDT_WGS))){
  sample_id <- names(SO_WCDT_WGS)[i]
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing:", sample_id, "\n", sep=" ")
  cnv <- SO_WCDT_WGS[[i]]$CNA_genes[c('chrom','weighted_mean')]
  cnv$gene <- rownames(SO_WCDT_WGS[[i]]$CNA_genes)
  cnv_wcdt.matrix[[sample_id]] <- cnv
  cnv_wcdt.matrix[[i]]$ploidy <- SO_WCDT_WGS[[i]]$ploidy
  
  cnv_wcdt.matrix[[i]]$cnv.def <-
    ifelse(cnv_wcdt.matrix[[i]]$chrom %in% sex.chr & cnv_wcdt.matrix[[i]]$weighted_mean > cnv_wcdt.matrix[[i]]$ploidy*0.8,'Copy gain',
           ifelse(cnv_wcdt.matrix[[i]]$chrom %out% sex.chr &cnv_wcdt.matrix[[i]]$weighted_mean > cnv_wcdt.matrix[[i]]$ploidy*1.95,'Copy gain',
                  ifelse(cnv_wcdt.matrix[[i]]$chrom %in% sex.chr & cnv_wcdt.matrix[[i]]$weighted_mean < cnv_wcdt.matrix[[i]]$ploidy* 0.3,'Copy loss',
                         ifelse(cnv_wcdt.matrix[[i]]$chrom %out% sex.chr & cnv_wcdt.matrix[[i]]$weighted_mean < cnv_wcdt.matrix[[i]]$ploidy* 0.65,'Copy loss',''))))
  cnv_wcdt.matrix[[i]]$Biallelic <- ifelse(cnv_wcdt.matrix[[i]]$weighted_mean < 0.5,'Biallelic','') }


# Subset by the genes we previously select 
subset_wcdt.CNA <- lapply(cnv_wcdt.matrix,function(x) x[x$gene%in%Genes_by_Cat_subset$Genes,])
names(subset_wcdt.CNA) <- names(cnv_wcdt.matrix)


mat_wcdt_CNA <- matrix(ncol = 1,nrow=length(Genes_by_Cat_subset$Genes))
mat_wcdt_CNA[,1] = Genes_by_Cat_subset$Genes
colnames(mat_wcdt_CNA) <- 'gene'
mat_wcdt_CNA <- foreach(i=names(subset_wcdt.CNA),.combine = 'cbind.fill') %do% {merge(mat_wcdt_CNA,subset_wcdt.CNA[[i]][,c(3,5:6)],by='gene')}
rownames(mat_wcdt_CNA) = mat_wcdt_CNA$gene;mat_wcdt_CNA <-  mat_wcdt_CNA[,c(F,T,T)];colnames(mat_wcdt_CNA) = rep(names(subset_wcdt.CNA),each=2)

#===============================================================================
# Adding AR.enhancer to the matrix, as defined above
# AR enhancer region reported in cell Quigley 2018 and Nature genetics Zhao et al.

# regions<- data.frame("chrom"="chrX","start"=66895158,"end"= 66910158, "subset.CNA" = "AR_enhancer")

AR_enhancer_wcd <- data.frame()
for (i in names(SO_WCDT_WGS)){
  region_weighted <- somatic_add_CNA_by_region_from_segments( SO_WCDT_WGS[[i]], regions)
  region_weighted <- cbind(sample_id = SO_WCDT_WGS[[i]]$sample_id, region_weighted )
  AR_enhancer_wcd <- rbind( AR_enhancer_wcd, region_weighted)
}

AR_enhancer_wcd$cnv.def <- c("")
AR_enhancer_wcd$Biallelic <- c("")
for(i in 1:dim(AR_enhancer_wcd)[1]){
  ploidy <- SO_WCDT_WGS[[i]]$ploidy
  AR_enhancer_wcd[i,]$cnv.def <-
    ifelse(AR_enhancer_wcd[i,]$chrom == "chrX" & AR_enhancer_wcd[i,]$weighted_mean >= ploidy*0.8,'Copy gain',
           ifelse(AR_enhancer_wcd[i,]$chrom == "chrX" & AR_enhancer_wcd[i,]$weighted_mean <= ploidy*0.3 , 'Copy loss','') )
  AR_enhancer_wcd[i,]$Biallelic <- ifelse(AR_enhancer_wcd[i,]$weighted_mean < 0.5,'Biallelic','')           
}

AR_enhancer_wcdt_mod <- foreach(i=1:nrow(AR_enhancer_wcd),.combine = 'cbind.fill') %do% {cbind( AR_enhancer_wcd[i,c(6:7)])}
colnames(AR_enhancer_wcdt_mod) <- rep(AR_enhancer_wcd$sample_id,each = 2)

# Combine the Copy number list of genes of interest and the AR enhancer copy numbers
mat_wcdt_CNA_mod<- rbind(mat_wcdt_CNA, "AR_enhancer" = AR_enhancer_wcdt_mod)

#===============================================================================
##### Mutation matrix of the long tail gene list from WCDT WGS sequencing data

subset.wcdt.onco <- list()
for(i in 1:length(names(SO_WCDT_WGS))){
  sample_id <- names(SO_WCDT_WGS)[i]
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing sample :", sample_id, "\n", sep=" ")
  # Wgs data   
  wgs <- data.frame()
  if(nrow(SO_WCDT_WGS[[i]]$strelka_mutect) != 0 & any(SO_WCDT_WGS[[i]]$strelka_mutect$gene %in% Genes_by_Cat_subset$Genes) ){
    
    genes <- SO_WCDT_WGS[[i]]$strelka_mutect[SO_WCDT_WGS[[i]]$strelka_mutect$gene %in% Genes_by_Cat_subset$Genes ,]$gene
    pos <- SO_WCDT_WGS[[i]]$strelka_mutect[SO_WCDT_WGS[[i]]$strelka_mutect$gene %in% Genes_by_Cat_subset$Genes ,]$pos
    df_genes_wgs <- as.data.frame(cbind(genes, pos))
    df_genes_wgs<- df_genes_wgs[!duplicated(df_genes_wgs),]
    df_genes_wgs$effects <- c("")
    df_genes_wgs$aa <- c("")
    
    for (n in 1:nrow(df_genes_wgs)){
      effects_wgs <- c()
      effects_wgs <- append(effects_wgs, SO_WCDT_WGS[[i]]$strelka_mutect[SO_WCDT_WGS[[i]]$strelka_mutect$gene == df_genes_wgs$genes[n] & SO_WCDT_WGS[[i]]$strelka_mutect$pos == df_genes_wgs$pos[n],]$effect)
      effects_wgs <- append(effects_wgs, SO_WCDT_WGS[[i]]$strelka[SO_WCDT_WGS[[i]]$strelka$gene == df_genes_wgs$genes[n] & SO_WCDT_WGS[[i]]$strelka$pos == df_genes_wgs$pos[n],]$effect)
      effects_wgs <- append(effects_wgs, SO_WCDT_WGS[[i]]$mutect[SO_WCDT_WGS[[i]]$mutect$gene == df_genes_wgs$genes[n] & SO_WCDT_WGS[[i]]$mutect$pos == df_genes_wgs$pos[n],]$effect)
      effects_wgs <- unique(effects_wgs)
      aa_wgs <- c()
      aa_wgs <- append(aa_wgs, SO_WCDT_WGS[[i]]$strelka_mutect[SO_WCDT_WGS[[i]]$strelka_mutect$gene == df_genes_wgs$genes[n] & SO_WCDT_WGS[[i]]$strelka_mutect$pos == df_genes_wgs$pos[n],]$aa)
      aa_wgs <- append(aa_wgs, SO_WCDT_WGS[[i]]$strelka[SO_WCDT_WGS[[i]]$strelka$gene == df_genes_wgs$genes[n] & SO_WCDT_WGS[[i]]$strelka$pos == df_genes_wgs$pos[n],]$aa)
      aa_wgs <- append(aa_wgs, SO_WCDT_WGS[[i]]$mutect[SO_WCDT_WGS[[i]]$mutect$gene == df_genes_wgs$genes[n] & SO_WCDT_WGS[[i]]$mutect$pos == df_genes_wgs$pos[n],]$aa)
      aa_wgs <- unique(aa_wgs)
      
      if( any( str_detect(effects_wgs, "missense|frame|stop|splice" ) ) ) {
        df_genes_wgs$effects[n] <- paste0(effects_wgs[str_detect(effects_wgs, "missense|frame|stop|splice" )], collapse=';')
        df_genes_wgs$aa[n] <- paste0(aa_wgs, collapse=';')
        wgs <- rbind(wgs, df_genes_wgs[n,])
      }
      else{
        next
      }
    }
    
    if (nrow(wgs) != 0 ) {
      subset.wcdt.onco[[sample_id]] <- rbind(subset.wcdt.onco[[sample_id]],wgs)
    }
    else{
      print(paste0("None mutation meet the criteria in patient: ", sample_id))
    }
  }
  
  rownames(subset.wcdt.onco[[sample_id]])<- NULL
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Sample :", sample_id, "processed.", "\n", sep=" ")
}
subset.wcdt.onco_all <- lapply(subset.wcdt.onco,function(x) x[!duplicated(x),])

subset.wcdt.onco_all<-lapply(subset.wcdt.onco_all, function(x) 
  cbind(x, aa_mod = (strsplit(as.character(x$aa), '.' ,fixed=TRUE) 
                     %>%  gsub("p.","", .)
                     %>%  gsub('\\"',"", .)
                     %>%  gsub(";","", .)
                     %>%  gsub("c\\(\\,","", .)
                     %>%  gsub("\\)","", .)
                     %>%  gsub("[[:space:]]","", .)
                     %>%  gsub("(^\\,)","", .)

  ) ) )

for (i in 1:length(subset.wcdt.onco_all) ){
  subset.wcdt.onco_all[[i]]$effects_mod <- c("")
  for (n in 1:nrow(subset.wcdt.onco_all[[i]])){
    if(subset.wcdt.onco_all[[i]][n,]$genes %in% Hotspot$Gene){
      aa_mod_sep <- as.data.frame(cbind(str_split_1(subset.wcdt.onco_all[[i]][n,]$aa_mod,","))) #,simplify = TRUE)
      aa_mod_sep_2 <- aa_mod_sep  %>% separate(V1, into = c("AA_1", "num", "AA_2"),  sep = "(?=[A-Za-z])(?<=[0-9])|(?<=[0-9])(?=[_])|(?<=[0-9])(?=[*])|(?=[0-9])(?<=[A-Za-z])" ) %>% 
        mutate(AA_1 =str_replace_all(AA_1, setNames(AA_annot$abbreviation_1, AA_annot$abbreviation_3) ) )  %>%
        mutate(AA_2 =str_replace_all(AA_2, setNames(AA_annot$abbreviation_1, AA_annot$abbreviation_3) ) )  %>%
        mutate(aa_mod2 = paste0(AA_1, num ))
      for (m in 1:nrow(aa_mod_sep_2)){
        if (paste0(subset.wcdt.onco_all[[i]][n,]$genes, aa_mod_sep_2$aa_mod2[m]) %in% paste0(Hotspot$Gene,Hotspot$Residue) ){
          subset.wcdt.onco_all[[i]][n,] <- subset.wcdt.onco_all[[i]][n,] %>% 
            mutate(across(everything(),~ gsub("_variant|conservative_", "", .)))  %>%
            mutate(effects_mod = paste0(str_split(effects, "&", simplify = T)[, 1],"_","hotspot") )
        }
        else {
          subset.wcdt.onco_all[[i]][n,] <- subset.wcdt.onco_all[[i]][n,] %>% 
            mutate(across(everything(),~ gsub("_variant|conservative_", "", .)))  %>%
            mutate(effects = paste0(str_split(effects, "&", simplify = T)[, 1],"_","VUS") )
        }
      }
    }
    else {
      subset.wcdt.onco_all[[i]][n,] <- subset.wcdt.onco_all[[i]][n,] %>% 
        mutate(across(everything(),~ gsub("_variant|conservative_", "", .)))  %>%
        mutate(effects = paste0(str_split(effects, "&", simplify = T)[, 1],"_","VUS") )
      
    }
  }
}

subset.wcdt.onco_all<-lapply(subset.wcdt.onco_all, function(x) x %>%
                               mutate(across(everything(),~ gsub("_VUS_VUS_VUS_VUS", "_VUS", .))) %>%
                               mutate(across(everything(),~ gsub("_VUS_VUS_VUS", "_VUS", .))) %>%
                               mutate(across(everything(),~ gsub("_VUS_VUS", "_VUS", .))) %>%
                               mutate(across(everything(),~ gsub("stop_gained_VUS", "Stop gained", .))) %>%
                               mutate(across(everything(),~ gsub("stop_gained_hotspot", "Stop gained hotspot", .))) %>%
                               mutate(across(everything(),~ gsub("frameshift_VUS", "Frameshift", .))) %>%
                               mutate(across(everything(),~ gsub("frameshift", "Frameshift", .))) %>%
                               mutate(across(everything(),~ gsub("inframe_deletion_VUS", "Inframe", .))) %>%
                               mutate(across(everything(),~ gsub("inframe_deletion_hotspot", "Inframe hotspot", .))) %>%
                               mutate(across(everything(),~ gsub("inframe_insertion_VUS", "Inframe", .))) %>%
                               mutate(across(everything(),~ gsub("inframe_insertion_hotspot", "Inframe hotspot", .))) %>%
                               mutate(across(everything(),~ gsub("stop_lost_VUS", "Stop lost", .))) %>% 
                               mutate(across(everything(),~ gsub("splice_region_VUS", "Splice region", .))) %>%
                               mutate(across(everything(),~ gsub("stop_", "Stop ", .))) %>%
                               mutate(across(everything(),~ gsub("missense", "Missense", .))) %>%
                               mutate(across(everything(),~ gsub("_VUS", " VUS", .))) %>%
                               mutate(across(everything(),~ gsub("_hotspot", " hotspot", .))) ) 

tmp.df <- data.frame(genes=Genes_by_Cat_subset$Genes)

hotspot_onco_wcdt <- foreach(i=names(subset.wcdt.onco_all)) %do% { 
  merge(tmp.df,subset.wcdt.onco_all[[i]],by='genes',all.x=T ) 
}
names(hotspot_onco_wcdt) <- names(subset.wcdt.onco_all)

hotspot_onco_wcdt <- foreach(i=names(hotspot_onco_wcdt)) %do% { hotspot_onco_wcdt[[i]][,c(1,3,6)] }
names(hotspot_onco_wcdt) <- names(subset.wcdt.onco_all)

for (i in 1:length(names(hotspot_onco_wcdt))){
  hotspot_onco_wcdt[[i]][is.na(hotspot_onco_wcdt[[i]])] <- ""
  for (n in 1:length(Genes_by_Cat_subset$Genes)){
    if( sum(hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n]) > 1){
      if( length(unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) > 1 & all(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod == "") ){
        hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0(paste0( unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects ), collapse = ";"), ";","Multi-hit"),
                                                                                                                 dim( hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1])
      }
      if( length(unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) > 1 & any(str_detect(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod, pattern = "hotspot"))  ){
        hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0( unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n]& hotspot_onco_wcdt[[i]]$effects_mod != "",]$effects_mod ), ";","Multi-hit"),
                                                                                                                 dim( hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1])
      }
      if( length(unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) == 1 & all(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod == "") ){
        hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0( unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects ) , ";","Multi-hit"),
                                                                                                                 dim( hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1])
      }
      if (length(unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects)) == 1 & any(str_detect(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod, pattern = "hotspot") ) ){
        hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <-  rep( paste0( unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n] & hotspot_onco_wcdt[[i]]$effects_mod != "",]$effects_mod ) , ";","Multi-hit"),
                                                                                                                 # unique(hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n] & hotspot_onco_wcdt[[i]]$effects_mod != "",]$effects_mod ), ";","Multi-hit"),
                                                                                                                 dim( hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],])[1] )
      }
    }
    
    else if (sum(hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n]) == 1 & hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod == ""){
      hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects_mod <- hotspot_onco_wcdt[[i]][hotspot_onco_wcdt[[i]]$genes==Genes_by_Cat_subset$Genes[n],]$effects
    }
    else{
      next
    }
  }
}
hotspot_onco_wcdt <- foreach(i=names(hotspot_onco_wcdt)) %do% { hotspot_onco_wcdt[[i]][,c(1,3)] }
names(hotspot_onco_wcdt) <- names(subset.wcdt.onco_all)
hotspot_onco_wcdt <- lapply(hotspot_onco_wcdt,function(x) x[!duplicated(x),])
hotspot_onco_wcdt <- as.data.frame(do.call(cbind,hotspot_onco_wcdt))
rownames(hotspot_onco_wcdt) = hotspot_onco_wcdt[,1]
hotspot_onco_wcdt <-  hotspot_onco_wcdt[,c(F,T)];colnames(hotspot_onco_wcdt) = names(subset.wcdt.onco_all)

missing_cols <- setdiff(names(SO_WCDT_WGS),colnames(hotspot_onco_wcdt))
hotspot_onco_wcdt[ , missing_cols] <- ""
hotspot_onco_wcdt <- hotspot_onco_wcdt[ , names(SO_WCDT_WGS)]
hotspot_onco_wcdt[is.na(hotspot_onco_wcdt)] <- ""
#===============================================================================
# Creating matrix for structural variants in WCDT

all_sv_wcdt <- list()
for(n in 1:length(names(SO_WCDT_WGS))){
  sample_id <- names(SO_WCDT_WGS)[n]
  sv_by_gene <- data.frame()
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_id, "\n", sep=" ")
  if ( !is_empty(SO_WCDT_WGS[[n]]$list_sv_gripss) ){
    gridss_gr <- GenomicRanges::makeGRangesFromDataFrame(SO_WCDT_WGS[[n]]$list_sv_gripss,
                                                         keep.extra.columns=T,
                                                         ignore.strand=T,
                                                         seqinfo=NULL,
                                                         seqnames.field=c("chrom_1"),
                                                         start.field="start_1",
                                                         end.field=c("end_1", "stop") )
    
    for(i in 1:length(Genes_by_Cat_subset$Genes) ){
      gene_query <- Genes_by_Cat_subset$Genes[i]
      gr_gene = GenomicRanges::GRanges(seqnames=genome$gene_locs[genome$gene_locs$symbol == gene_query,]$chrom,
                                       ranges = IRanges(start=genome$gene_locs[genome$gene_locs$symbol == gene_query,]$start,
                                                        end=genome$gene_locs[genome$gene_locs$symbol == gene_query,]$end))
      df_sv_overlapping_gene = data.frame( subsetByOverlaps( gridss_gr, gr_gene) )
      if (dim(df_sv_overlapping_gene)[1]>0 ){
        sv_by_gene <- rbind(sv_by_gene, cbind( "gene" = gene_query, "SV"= "Structural variant") )
      }
      else{
        sv_by_gene <- rbind(sv_by_gene, cbind( "gene" = gene_query, "SV"= NA ) )
      }
    }
    all_sv_wcdt[[sample_id]] <- rbind(sv_by_gene)
  }
  else{
    all_sv_wcdt[[sample_id]] <- cbind("gene" =Genes_by_Cat_subset$Genes, "SV"= NA)
    next
  }
}

sv_onco_wcdt <- do.call(cbind,all_sv_wcdt)
rownames(sv_onco_wcdt) = sv_onco_wcdt[,1];sv_onco_wcdt <-  sv_onco_wcdt[,c(F,T)];colnames(sv_onco_wcdt) = names(SO_WCDT_WGS)

#===============================================================================
# Fusions in WCDT also annotated as Structural variant

all_fusion_wcdt <- list()

for(n in 1:length(names(SO_WCDT_WGS))){
  sample_id <- names(SO_WCDT_WGS)[n]
  fusion_by_gene <- data.frame()
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_id, "\n", sep=" ")
  if ( !is_empty(SO_WCDT_WGS[[n]]$fusions_linx) ){
    for(i in 1:length(Genes_by_Cat_subset$Genes) ){
      for(m in 1:dim(SO_WCDT_WGS[[n]]$fusions_linx)[1] ){
        if (any(SO_WCDT_WGS[[n]]$fusions_linx[m,]$geneStart %in% Genes_by_Cat_subset$Genes[i] | SO_WCDT_WGS[[n]]$fusions_linx[m,]$geneEnd %in% Genes_by_Cat_subset$Genes[i]) & 
            SO_WCDT_WGS[[n]]$fusions_linx[m,]$geneStart != SO_WCDT_WGS[[n]]$fusions_linx[m,]$geneEnd ){
          fusion_by_gene <- rbind(fusion_by_gene, cbind( "gene" = Genes_by_Cat_subset$Genes[i], "Fusion"= "Structural variant") )
        }
        else{
          next
        }
      }
      fusion_by_gene <- rbind(fusion_by_gene, cbind( "gene" = Genes_by_Cat_subset$Genes[i], "Fusion"= NA ))
    }
    all_fusion_wcdt[[sample_id]] <- rbind(fusion_by_gene)
    all_fusion_wcdt[[sample_id]] <- all_fusion_wcdt[[sample_id]] %>% 
      distinct(gene, .keep_all = T)
  }
  else{
    all_fusion_wcdt[[sample_id]] <- cbind("gene" = Genes_by_Cat_subset$Genes, "Fusion"= NA)
    next
  }
}

fusion_onco_wcdt <- do.call(cbind,all_fusion_wcdt)
rownames(fusion_onco_wcdt) = fusion_onco_wcdt[,1];fusion_onco_wcdt <-  fusion_onco_wcdt[,c(F,T)];colnames(fusion_onco_wcdt) = names(SO_WCDT_WGS)

#===============================================================================
# Gathering all the alterations in the same matrix

gene.idx <- intersect(rownames(hotspot_onco_wcdt),Genes_by_Cat_subset$Genes )
mut_wcdt_matrix_all <- foreach(i=names(hotspot_onco_wcdt),.combine = 'cbind') %do%{
  df <- cbind(hotspot_onco_wcdt[gene.idx,][i],
              mat_wcdt_CNA_mod[gene.idx,c(T,F)][i],mat_wcdt_CNA_mod[gene.idx,c(F,T)][i],
              sv_onco_wcdt[gene.idx,][i],fusion_onco_wcdt[gene.idx,][i])
  colnames(df) <- c('ONCO','CNA','BI','SV','FUS')
  df[df == ""]<- NA
  df <- df %>% unite('test',sep = ';',remove = T,na.rm = T)
}
colnames(mut_wcdt_matrix_all) <- colnames(hotspot_onco_wcdt)

#===============================================================================
#'  Calculating chromothripsis in WCDT

# Using the same files as in the Cascade section
# fn_chrom_lengths = paste(dir_cascade_base, 'reproduce/metadata/HG38/HG38_chromosome_lengths.txt',sep='')
# fn_centromeres = paste(dir_cascade_base, 'reproduce/metadata/HG38/HG38_centromere_loci.txt',sep='')
# 
# chrom_lengths = read.table(fn_chrom_lengths,row.names=1, stringsAsFactors=FALSE) # at metadata/resources
# chrom_names = rownames(chrom_lengths)
# centromeres = read.table( fn_centromeres, header=TRUE, stringsAsFactors = FALSE) # at at metadata/resources

list_sv_wcdt_gridss<-data.frame()
list_sv_wcdt_gridss <- foreach(i=names(SO_WCDT_WGS),.combine = 'rbind') %do% {SO_WCDT_WGS[[i]]$list_sv_gripss}

sample_wcdt_ids = names(SO_WCDT_WGS)

# This step takes a long time (2 - 4 hours) 
# optional load the previous results, which also take more than 20 minutes (file is 2.17 GB)
# or the Rdata with all the dataframes generated in this section
load(paste0(dir_wcdt_base,'2022_10_05_WCDT_chromothripsis.Rdata'))
# chromo_scores_wcdt<- read.table(paste0(dir_wcdt_base,'WCDT_chromo_scores_2500000.txt'),  sep = "",header = T) # at results folder


for(n in 1:length(sample_wcdt_ids)){
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_wcdt_ids[i], "\n", sep=" ")
  for(i in 1:24){
    chrom = chrom_names[i]
    cs = chromo_score( sample_id = sample_wcdt_ids[n], chrom=chrom )
    if( n==1 & i == 1){
      chromo_scores_wcdt = cs
    }else{
      chromo_scores_wcdt = rbind(chromo_scores_wcdt, cs)
    }
  }
}

chromo_scores_wcdt$med_CNA = round( chromo_scores_wcdt$med_CNA, 3 )

n_samples = length(sample_wcdt_ids)
chrom_maxima_wcdt = matrix(0, n_samples,24)
for(i in 1:length(sample_wcdt_ids)){
  cat(format(Sys.time(), "%a %b %d %X %Y"), "Processing: ", sample_wcdt_ids[i], "\n", sep=" ")
  sc = chromo_scores_wcdt[chromo_scores_wcdt$sample_id == sample_wcdt_ids[i], ]
  print( paste("Loading chromothripsis scores:", sample_wcdt_ids[i] ))
  for(j in 1:24){
    chrom=rownames(chrom_lengths)[j]   
    sc_chrom = sc[sc$chrom==chrom,]
    for(x in 1:50){
      if( sum( sc_chrom$n_inv>=x & sc_chrom$n_del>=10 & 
               sc_chrom$n_cna_switch>=x )>0){
        chrom_maxima_wcdt[i,j]=x   
      }else{
        break
      }
    }
  }
}

# filtering steps
chrom_maxima_wcdt_filtered = chrom_maxima_wcdt
chrom_maxima_wcdt_filtered[chrom_maxima_wcdt_filtered<15]=0
list_chromo_wcdt=as.data.frame( which(chrom_maxima_wcdt_filtered>0, arr.ind = TRUE) )
list_chromo_wcdt[,1] = sample_wcdt_ids[list_chromo_wcdt[,1]]
list_chromo_wcdt$col = paste("chr", list_chromo_wcdt$col, sep='')
list_chromo_wcdt = list_chromo_wcdt[order(list_chromo_wcdt[,1]),]
names(list_chromo_wcdt) = c("sample_id", "chrom")
rownames(list_chromo_wcdt)<- NULL
list_chromo_wcdt$chrom[list_chromo_wcdt$chrom=="chr23"] = "chrX"

#===============================================================================
# create DF using poppy object from WCDT cohort

SO_WCDT_summary = data.frame()

for(i in 1:length(names(SO_WCDT_WGS)) ){
  
  sample_id = as.character(SO_WCDT_WGS[[i]]$sample_base)
  purity_PURPLE= as.numeric(SO_WCDT_WGS[[i]]$purity*100)
  ploidy_PURPLE= as.numeric(SO_WCDT_WGS[[i]]$ploidy)
  MSS = SO_WCDT_WGS[[i]]$msStatus
  tmbPerMb = as.numeric(SO_WCDT_WGS[[i]]$tmbPerMb)
  tmbPerMb_Log2 = log2(as.numeric(tmbPerMb))
  SVtmbPerMb = as.numeric(SO_WCDT_WGS[[i]]$svTumorMutationalBurden)
  WG_doubling = SO_WCDT_WGS[[i]]$wholeGenomeDuplication
  WG_doubling <- gsub("true", "TRUE", WG_doubling)
  WG_doubling <- gsub("false", "FALSE", WG_doubling)

  ETS = any(SO_WCDT_WGS[[i]]$fusions_linx$geneStart %in% ETS_fusion | SO_WCDT_WGS[[i]]$fusions_linx$geneEnd %in%ETS_fusion)
  ## match info from list_chromo_wcdt df
  if(sample_id %in% list_chromo_wcdt$sample_id ) {
    m = match.idx(list_chromo_wcdt$sample_id, sample_id)
    Chromothripsis = list_chromo_wcdt$chrom[ m$idx.A ]
    Chromo = "TRUE"
  }
  else{
    Chromothripsis = NA
    Chromo = "FALSE"
  }
  ## match info from wcdt_events_mod df
  if(sample_id %in% wcdt_events_mod$sample_id ) {
    m = match.idx(wcdt_events_mod$sample_id, sample_id)
    ecDNA_counts = wcdt_events_mod$ecDNA[ m$idx.A ]
    BFB_counts = wcdt_events_mod$BFB[ m$idx.A ]
  }
  else{
    ecDNA_counts = NA
    BFB_counts = NA
  }
  summary_row<-cbind(sample_id, purity_PURPLE, ploidy_PURPLE,ETS,
                     MSS, tmbPerMb, tmbPerMb_Log2, SVtmbPerMb, 
                     WG_doubling, ecDNA_counts, BFB_counts,
                     Chromothripsis, Chromo )
  SO_WCDT_summary = rbind(SO_WCDT_summary, summary_row)
}       
SO_WCDT_summary[c(2:3,6:8,10:11)] <- lapply(SO_WCDT_summary[c(2:3,6:8,10:11)], as.numeric)

#===============================================================================
# Creating a matrix of logical values for WCDT cohort any type of alteration 
# per sample. This matrix will be use for the statistics in Fig 4

all_wcdt_altertations<- mut_wcdt_matrix_all %>% mutate_all(na_if,"") # Replace "" by "NA"
all_wcdt_altertations <- all_wcdt_altertations %>%                            
  replace(!is.na(.), TRUE) %>%  # Replace any alteration by TRUE
  replace(is.na(.), FALSE) 
all_wcdt_altertations<- as.data.frame(t(all_wcdt_altertations))

# Subset and order WCDT based on AmpliconArchitect output
SO_WCDT_summary_sub <- SO_WCDT_summary[ match(id_order_wcdt , SO_WCDT_summary$sample_id ),]
idx_pt = match.idx(SO_WCDT_summary_sub$sample_id, rownames(all_wcdt_altertations),  allow.multiple.B = F)
SO_WCDT_summary_sub$TP53 = c("")
SO_WCDT_summary_sub[idx_pt$idx.A,]$TP53 = all_wcdt_altertations[idx_pt$idx.B, ]$TP53
SO_WCDT_summary_sub$PTEN = c("")
SO_WCDT_summary_sub[idx_pt$idx.A,]$PTEN = all_wcdt_altertations[idx_pt$idx.B, ]$PTEN

#===============================================================================
#'  Figure 4A
#===============================================================================
# Code adapted from David's script to plot AR locus amplifications by ecDNA

bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
locs = read.table(paste0(dir_cascade_base,"reproduce/metadata/grch38.gencode.v32.digest.txt"), header=TRUE) # at at metadata/resources

wcdt_ecDNA_BFB <- AA_gene_list_wcdt_mod[str_detect(AA_gene_list_wcdt_mod$feature, "ecDNA|BFB"), ]
wcdt_ecDNA_BFB$is_ecdna = rep(FALSE, dim(wcdt_ecDNA_BFB)[1])
wcdt_ecDNA_BFB$is_ecdna[ grep("ecDNA", wcdt_ecDNA_BFB$feature) ] = TRUE
wcdt_ecDNA_BFB$slug = paste(wcdt_ecDNA_BFB$sample_name, wcdt_ecDNA_BFB$amplicon_number, sep="_")
tt = sort(table( wcdt_ecDNA_BFB$gene[wcdt_ecDNA_BFB$is_ecdna]))

# fEC: ecDNA gene frequency
fEC = data.frame( count=as.numeric(tt), row.names=names(tt) )
m = match.idx( dimnames(fEC)[[1]], locs$symbol)
fEC$chrom = rep("unknown", dim(fEC)[1])
fEC$pos = rep(0, dim(fEC)[1])
fEC$chrom[m$idx.A] = locs$chrom[m$idx.B]
fEC$pos[m$idx.A] = locs$start[m$idx.B]
# add a few annotations
idx = which( dimnames(fEC)[[1]]=="PRNCR1")
fEC$chrom[idx]="chr8"
fEC$pos[idx] = 127079874
idx = which( dimnames(fEC)[[1]]=="PRUNE1")
fEC$chrom[idx]="chr1"
fEC$pos[idx] = 151008449
# fboth: ecDNA and BFB gene frequency
tt = sort(table( wcdt_ecDNA_BFB$gene))
fboth = data.frame( count=as.numeric(tt), row.names=names(tt) )
m = match.idx( dimnames(fboth)[[1]], locs$symbol)
fboth$chrom = rep("unknown", dim(fboth)[1])
fboth$pos = rep(0, dim(fboth)[1])
fboth$chrom[m$idx.A] = locs$chrom[m$idx.B]
fboth$pos[m$idx.A] = locs$start[m$idx.B]
fboth = fboth[order(fboth$count, decreasing = TRUE),]

# create a genomicRanges object for ecDNA intervals beds
gr_EC = GenomicRanges::makeGRangesFromDataFrame(ecDNA_BFB_int_all_wcdt, keep.extra.columns = TRUE, seqnames.field = "chrom", start.field = "start", end.field = "end")
gr_EC_X = GenomicRanges::makeGRangesFromDataFrame(ecDNA_BFB_int_all_wcdt[ecDNA_BFB_int_all_wcdt$chrom=="chrX",], keep.extra.columns = TRUE, seqnames.field = "chrom", start.field = "start", end.field = "end")
gr_EC_8 = GenomicRanges::makeGRangesFromDataFrame(ecDNA_BFB_int_all_wcdt[ecDNA_BFB_int_all_wcdt$chrom=="chr8",], keep.extra.columns = TRUE, seqnames.field = "chrom", start.field = "start", end.field = "end")
EC_coverage_X = coverage(gr_EC_X)
EC_coverage_8 = coverage(gr_EC_8)
# unique bins for gr_EC_X and their counts as genomicRanges
ends = cumsum( runLength(EC_coverage_X$chrX ) )
starts = c(1, (ends[1:(length(ends)-1)])+1)
gr_EC_counts_X = makeGRangesFromDataFrame(
  data.frame( seqnames=rep("chrX", length(ends)),
              start=starts,
              end=ends,
              count=runValue( EC_coverage_X$chrX )),
  keep.extra.columns = TRUE
)

#===============================================================================
#'  Figure 4B
#===============================================================================
# Figure 4B take it from AmpliconArchitect output and modified the colors and the
# figure legend in Illustrator

#===============================================================================
#'  Data frame to reproduce statistics and Figure 4C
#===============================================================================
# Modified the previues df SO_summary_sub to include presence/absence of ecDNA

ecDNA_stats_cascade <- SO_summary_sub

idx_pt = match.idx(ecDNA_stats_cascade$sample_id,  cascade_events_mod$sample_id, 
                   allow.multiple.B = F)
ecDNA_stats_cascade$ecDNA = "NA"
ecDNA_stats_cascade$BFB = "NA"
ecDNA_stats_cascade$ecDNA_cat = "NA"
ecDNA_stats_cascade$BFB_cat = "NA"
ecDNA_stats_cascade[idx_pt$idx.A,]$ecDNA = cascade_events_mod[idx_pt$idx.B, ]$ecDNA
ecDNA_stats_cascade[idx_pt$idx.A,]$BFB <- cascade_events_mod[idx_pt$idx.B, ]$BFB
ecDNA_stats_cascade$ecDNA_cat<- ifelse(ecDNA_stats_cascade$ecDNA == 0,"ecDNA(-)",ifelse(ecDNA_stats_cascade$ecDNA >=1,"ecDNA(+)","NA"))
ecDNA_stats_cascade$BFB_cat<- ifelse(ecDNA_stats_cascade$BFB == 0,"BFB(-)",ifelse(ecDNA_stats_cascade$BFB >=1,"BFB(+)","NA"))

#===============================================================================
#'  Data frame to reproduce statistics and Figure 4D
#===============================================================================
# Modified the prior df SO_WCDT_summary_sub to include presence/absence of ecDNA


ecDNA_stats_wcdt <- SO_WCDT_summary_sub

idx_pt = match.idx(ecDNA_stats_wcdt$sample_id,  wcdt_events_mod$sample_id, 
                   allow.multiple.B = F)
ecDNA_stats_wcdt$ecDNA = "NA"
ecDNA_stats_wcdt$BFB = "NA"
ecDNA_stats_wcdt$ecDNA_cat = "NA"
ecDNA_stats_wcdt$BFB_cat = "NA"
ecDNA_stats_wcdt[idx_pt$idx.A,]$ecDNA = wcdt_events_mod[idx_pt$idx.B, ]$ecDNA
ecDNA_stats_wcdt[idx_pt$idx.A,]$BFB <- wcdt_events_mod[idx_pt$idx.B, ]$BFB
ecDNA_stats_wcdt$ecDNA_cat<- ifelse(ecDNA_stats_wcdt$ecDNA == 0,"ecDNA(-)",ifelse(ecDNA_stats_wcdt$ecDNA >=1,"ecDNA(+)","NA"))
ecDNA_stats_wcdt$BFB_cat<- ifelse(ecDNA_stats_wcdt$BFB == 0,"BFB(-)",ifelse(ecDNA_stats_wcdt$BFB >=1,"BFB(+)","NA"))

#===============================================================================
#'  Data frame to phylogenetic reconstructions in Suppl Figures 1 a-f
#===============================================================================
# Uses the function generate_coding_matrix_mod to create binary matrix with alterations 
# to them calculate the distances and create a tree using phangorn::NJ
# Neighbor-Joining (NJ) algorithm because allows for unequal rates of evolution.

#===============================================================================
#' Supplementary Figures 2 was generated by AARDVARK
#===============================================================================

# AARDVARK is available at https://github.com/DavidQuigley/aardvark
# For more information, visit the repository or contact: david.quigley@ucsf.edu

load( paste0(dir_cascade_results,'2023_04_20_aardvark_CA071-GERM.Rdata' ))
fns = paste0( dir_cascade_results,'2023_04_20_aardvark_CA071-0',1:9,'.Rdata' )
pathogenic_mut = aardvark::Mutation("chr13", 32339657, "CTT", "C", transcript_BRCA2)

sums_gl =  summarize_candidates(reads_gl, transcript_BRCA2, pathogenic_mut)
sums_all = list()
for(i in 1:9){
  print(i)
  load( fns[i] )
  sums_all = c(sums_all, list( summarize_candidates(reads, transcript_BRCA2, pathogenic_mut) ) )
}
names(sums_all) = paste0("CA071-0",1:9)

#===============================================================================
#' Supplementary Figures 3A: plots copy number, expression levels and distribution of ecDNA 
#'  and other complex events for AR, MYC, and PCAT1 in Cascade cohort. 
#===============================================================================
# Change order based on AR expression
id_order_ar_cascade<- colnames(Log2.TPM_ft[,order(Log2.TPM_ft["AR",], decreasing = T)])
Log2.TPM_mod.zc_order <- as.matrix(Log2.TPM_mod.zc[, match(id_order_ar_cascade, colnames(Log2.TPM_mod.zc)) ])

#===============================================================================
#' Supplementary Figures 3B: plots copy number, expression levels and distribution of ecDNA 
#'  and other complex events for AR, MYC, and PCAT1 in WCDT cohort. 
#===============================================================================
# Extract the tpms

TPM_wcdt <-as.matrix(tx_wcdt_mod$abundance)
# rename
colnames(TPM_wcdt)[41] <- "DTB-055-PRO"
colnames(TPM_wcdt)[212] <- "PR-063-BL"

id_order_ar_wcdt<- colnames(TPM_wcdt[,order(TPM_wcdt["AR",], decreasing = T)])
# Removing samples "DTB-193-BL" and "DTB-201-BL" beacuse we don't have expression data
id_order_ar_wcdt_mod <- id_order_ar_wcdt[id_order_ar_wcdt!= c("DTB-201-BL","DTB-193-BL")]
# Subset for the samples we have AA outputs
id_order_ar_wcdt_mod <- id_order_ar_wcdt_mod[id_order_ar_wcdt %in% colnames(matrix_AA_gene_wcdt_freq_sub_new_mod)]

Log2.TPM_wcdt  <- log2( 1+ TPM_wcdt )
# Remove samples without matched WGS-seq /AA output
Log2.TPM_wcdt_ft <- Log2.TPM_wcdt[, colnames(Log2.TPM_wcdt) %in% id_order_ar_wcdt_mod ] 
Log2.TPM_wcdt.zc<- t(scale(t(Log2.TPM_wcdt_ft),scale = T,center = T))

#===============================================================================
# The AR and NE signalling scores as I calculated above
# I used GSVA to transform the gene expression measurements into enrichment scores 
# for these two gene sets

Log2_TPM_wcdt_zc_gsva <- Log2.TPM_wcdt.zc
row.names(Log2_TPM_wcdt_zc_gsva)[which(row.names(Log2_TPM_wcdt_zc_gsva) %in% AR_sigID$external_gene_name)] <- AR_sigID[row.names(Log2_TPM_wcdt_zc_gsva)[which(row.names(Log2_TPM_wcdt_zc_gsva) %in% AR_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(Log2_TPM_wcdt_zc_gsva)[which(row.names(Log2_TPM_wcdt_zc_gsva) %in% NEPC_sigID$external_gene_name)] <- NEPC_sigID[row.names(Log2_TPM_wcdt_zc_gsva)[which(row.names(Log2_TPM_wcdt_zc_gsva) %in% NEPC_sigID$external_gene_name)],'ensembl_gene_id_version']
gsvaScore_wcdt_all <- gsva(as.matrix(Log2_TPM_wcdt_zc_gsva), list(AR_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=5, max.sz=500, method="zscore")
rownames( gsvaScore_wcdt_all ) <- c("AR_score","NEPC_score")

#===============================================================================
# Generate a df for WCDT copy number and AR/NE scores to plot in Supp Fig 3B

wcdt_gsvaScore <- as.data.frame(t(gsvaScore_wcdt_all))

RNA_wcdt = data.frame()
for(i in 1:length(colnames(Log2.TPM_wcdt.zc)) ){
  
  sample_id = colnames(Log2.TPM_wcdt.zc)[i]
  if (sample_id %in% rownames(wcdt_gsvaScore) ) {
    
    m1 = match.idx(rownames(wcdt_gsvaScore), sample_id)
    AR_score = as.numeric(wcdt_gsvaScore$AR_score[ m1$idx.A ])
    m2 = match.idx(rownames(wcdt_gsvaScore), sample_id)
    NE_score = as.numeric(wcdt_gsvaScore$NEPC_score[ m2$idx.A ])
  }
  else{
    AR_score = NA 
    NE_score = NA
  }
  summary_row<-cbind(sample_id, AR_score, NE_score)
  RNA_wcdt = rbind(RNA_wcdt, summary_row)
}
colnames(RNA_wcdt)

Log2.TPM_wcdt.zc_mod <- Log2.TPM_wcdt.zc[,id_order_ar_wcdt_mod]

RNA_wcdt$AR_Log2_cnv <- c("")
RNA_wcdt$MYC_Log2_cnv <- c("")
RNA_wcdt$PCAT1_Log2_cnv <- c("")
RNA_wcdt$AR_cnv <- c("")
RNA_wcdt$MYC_cnv <- c("")
RNA_wcdt$PCAT1_cnv <- c("")

for(i in 1:length(names(SO_WCDT_WGS))){
  sample_id <- names(SO_WCDT_WGS)[i]
  if (sample_id %in% RNA_wcdt$sample_id){
    RNA_wcdt[RNA_wcdt$sample_id == sample_id,]$AR_cnv <- as.numeric(SO_WCDT_WGS[[i]]$CNA_genes["AR",]$weighted_mean)
    RNA_wcdt[RNA_wcdt$sample_id == sample_id,]$AR_Log2_cnv <- as.numeric(log2(SO_WCDT_WGS[[i]]$CNA_genes["AR",]$weighted_mean))
    RNA_wcdt[RNA_wcdt$sample_id == sample_id,]$MYC_cnv <- as.numeric(SO_WCDT_WGS[[i]]$CNA_genes["MYC",]$weighted_mean)
    RNA_wcdt[RNA_wcdt$sample_id == sample_id,]$MYC_Log2_cnv <- as.numeric(log2(SO_WCDT_WGS[[i]]$CNA_genes["MYC",]$weighted_mean))
    RNA_wcdt[RNA_wcdt$sample_id == sample_id,]$PCAT1_cnv <- as.numeric(SO_WCDT_WGS[[i]]$CNA_genes["PCAT1",]$weighted_mean)
    RNA_wcdt[RNA_wcdt$sample_id == sample_id,]$PCAT1_Log2_cnv <- as.numeric(log2(SO_WCDT_WGS[[i]]$CNA_genes["PCAT1",]$weighted_mean))
  }
  else{
    next
  }
}

RNA_wcdt[c(2:9)] <- lapply(RNA_wcdt[c(2:9)], as.numeric)
RNA_wcdt<- RNA_wcdt[match(id_order_ar_wcdt_mod, RNA_wcdt$sample_id), ]


#===============================================================================
#'  Supplementary Figures 4_A-B: Using output from the script amplicon_similarity.py 
#'  from AmpliconClassifier tool
#===============================================================================
# All the AR events from all tumors from patient CA071

AR_CA071 <- AA_gene_list_cascade[str_detect(AA_gene_list_cascade$sample_name,"CA071") & str_detect(AA_gene_list_cascade$gene, "^AR$"),]

CA071_sim_score<- read.table(paste0(dir_cascade_results,'CA071_similarity_amplicon_similarity_scores.txt'), header = T, sep = "\t") # at results folder
CA071_sim_score$Amp1<- gsub(paste(paste0(sprintf("%s_WGS/AmpliconArchitect/output_data/",AR_CA071$sample_name),AR_CA071$sample_name,"_AA_OUT/|"), collapse =""),  "",CA071_sim_score$Amp1)
CA071_sim_score$Amp2<- gsub(paste(paste0(sprintf("%s_WGS/AmpliconArchitect/output_data/",AR_CA071$sample_name),AR_CA071$sample_name,"_AA_OUT/|"), collapse =""),  "",CA071_sim_score$Amp2)

CA071_sim_score <- CA071_sim_score %>%
  dplyr::mutate_at(vars(c("Amp1", "Amp2")), ~ str_replace(., "_Pre_AI", ""))

CA071_sim_score_m<- matrix(data= 1, nrow = 9, ncol = 9)
rownames(CA071_sim_score_m) <- unique(c(CA071_sim_score$Amp1,CA071_sim_score$Amp2) )
colnames(CA071_sim_score_m) <- unique(c(CA071_sim_score$Amp1,CA071_sim_score$Amp2) )

for (i in 1:nrow(CA071_sim_score)) {
  for (n in 1:nrow(CA071_sim_score_m)) {
    for (p in 1:ncol(CA071_sim_score_m)) {
      
      if ( rownames(CA071_sim_score_m)[n] != colnames(CA071_sim_score_m)[p] &
           CA071_sim_score$Amp1[i] ==  rownames(CA071_sim_score_m)[n] & 
           CA071_sim_score$Amp2[i] ==  colnames(CA071_sim_score_m)[p]) {
        CA071_sim_score_m[n, p] <- as.numeric(CA071_sim_score[CA071_sim_score$Amp1 ==  rownames(CA071_sim_score_m)[n] & CA071_sim_score$Amp2 ==  colnames(CA071_sim_score_m)[p], ]$SimilarityScore)
        CA071_sim_score_m[p, n] <- as.numeric(CA071_sim_score[CA071_sim_score$Amp1 ==  rownames(CA071_sim_score_m)[n] & CA071_sim_score$Amp2 ==  colnames(CA071_sim_score_m)[p], ]$SimilarityScore)
      }
      else{
        next
      }
    }
  }
}

CA071_jac_score_m<- matrix(data= NA, nrow = 9, ncol = 9)
rownames(CA071_jac_score_m) <- unique(c(CA071_sim_score$Amp1,CA071_sim_score$Amp2) )
colnames(CA071_jac_score_m) <- unique(c(CA071_sim_score$Amp1,CA071_sim_score$Amp2) )

for (i in 1:nrow(CA071_sim_score)) {
  for (n in 1:nrow(CA071_jac_score_m)) {
    for (p in 1:ncol(CA071_jac_score_m)) {
      if ( rownames(CA071_jac_score_m)[n] != colnames(CA071_jac_score_m)[p] &
           CA071_sim_score$Amp1[i] ==  rownames(CA071_jac_score_m)[n] & 
           CA071_sim_score$Amp2[i] ==  colnames(CA071_jac_score_m)[p]) {
        CA071_jac_score_m[n, p] <- as.numeric(CA071_sim_score[CA071_sim_score$Amp1 ==  rownames(CA071_jac_score_m)[n] & CA071_sim_score$Amp2 ==  colnames(CA071_jac_score_m)[p], ]$JaccardBreakpoint)
        CA071_jac_score_m[p, n] <- as.numeric(CA071_sim_score[CA071_sim_score$Amp1 ==  rownames(CA071_jac_score_m)[n] & CA071_sim_score$Amp2 ==  colnames(CA071_jac_score_m)[p], ]$JaccardBreakpoint)
      }
      else{
        next
      }
    }
  }
}

CA071_upper_tri <- get_upper_tri(CA071_sim_score_m)
CA071_lower_tri <- get_lower_tri(CA071_jac_score_m)
CA071_correlation_similarity <- CA071_upper_tri
for (i in 1:nrow(CA071_correlation_similarity)){
  for (n in 1:ncol(CA071_correlation_similarity)){
    if(is.na(CA071_correlation_similarity[i,n])){
      CA071_correlation_similarity[i,n]<- CA071_lower_tri[i,n]
    }
    else{
      next
    }
  }
}
CA071_correlation_similarity

rownames(CA071_correlation_similarity) <-  do.call('rbind',strsplit(as.character(rownames(CA071_correlation_similarity)), '_' ,fixed=TRUE))[,1]
colnames(CA071_correlation_similarity) <- do.call('rbind',strsplit(as.character(colnames(CA071_correlation_similarity)), '_' ,fixed=TRUE))[,1]


# The AR events from metastasis 82, 86, and 94 from patient CA108 

AR_CA108 <- AA_gene_list_cascade[str_detect(AA_gene_list_cascade$sample_name,"CA108") & str_detect(AA_gene_list_cascade$gene, "^AR$"),]
CA108_sim_score<- read.table(paste0(dir_cascade_results,'CA108_similarity_amplicon_similarity_scores.tsv'), header = T, sep = "\t") # at results folder
CA108_sim_score <- CA108_sim_score %>% 
  dplyr::mutate_all(sub, pattern = '/data1/projects/human_prostate_CASCADE/2021_03_pilot/', replacement = '')

CA108_sim_score$Amp1<- gsub(paste(paste0(sprintf("%s_WGS/AmpliconArchitect/output_data/",AR_CA108$sample_name),AR_CA108$sample_name,"_AA_OUT/|"), collapse =""),  "",CA108_sim_score$Amp1)
CA108_sim_score$Amp2<- gsub(paste(paste0(sprintf("%s_WGS/AmpliconArchitect/output_data/",AR_CA108$sample_name),AR_CA108$sample_name,"_AA_OUT/|"), collapse =""),  "",CA108_sim_score$Amp2)

CA108_sim_score <- CA108_sim_score %>%
  dplyr::mutate_at(vars(c("Amp1", "Amp2")), ~ str_replace(., "_Pre_AI", ""))

CA108_sim_score_m<- matrix(data= 1, nrow = 3, ncol = 3)
rownames(CA108_sim_score_m) <- unique(c(CA108_sim_score$Amp1,CA108_sim_score$Amp2) )
colnames(CA108_sim_score_m) <- unique(c(CA108_sim_score$Amp1,CA108_sim_score$Amp2) )

for (i in 1:nrow(CA108_sim_score)) {
  for (n in 1:nrow(CA108_sim_score_m)) {
    for (p in 1:ncol(CA108_sim_score_m)) {
      if ( rownames(CA108_sim_score_m)[n] != colnames(CA108_sim_score_m)[p] &
           CA108_sim_score$Amp1[i] ==  rownames(CA108_sim_score_m)[n] & 
           CA108_sim_score$Amp2[i] ==  colnames(CA108_sim_score_m)[p]) {
        CA108_sim_score_m[n, p] <- as.numeric(CA108_sim_score[CA108_sim_score$Amp1 ==  rownames(CA108_sim_score_m)[n] & CA108_sim_score$Amp2 ==  colnames(CA108_sim_score_m)[p], ]$SimilarityScore)
        CA108_sim_score_m[p, n] <- as.numeric(CA108_sim_score[CA108_sim_score$Amp1 ==  rownames(CA108_sim_score_m)[n] & CA108_sim_score$Amp2 ==  colnames(CA108_sim_score_m)[p], ]$SimilarityScore)
      }
      else{
        next
      }
    }
  }
}

CA108_jac_score_m<- matrix(data= NA, nrow = 3, ncol = 3)
rownames(CA108_jac_score_m) <- unique(c(CA108_sim_score$Amp1,CA108_sim_score$Amp2) )
colnames(CA108_jac_score_m) <- unique(c(CA108_sim_score$Amp1,CA108_sim_score$Amp2) )

for (i in 1:nrow(CA108_sim_score)) {
  for (n in 1:nrow(CA108_jac_score_m)) {
    for (p in 1:ncol(CA108_jac_score_m)) {
      if ( rownames(CA108_jac_score_m)[n] != colnames(CA108_jac_score_m)[p] &
           CA108_sim_score$Amp1[i] ==  rownames(CA108_jac_score_m)[n] & 
           CA108_sim_score$Amp2[i] ==  colnames(CA108_jac_score_m)[p]) {
        CA108_jac_score_m[n, p] <- as.numeric(CA108_sim_score[CA108_sim_score$Amp1 ==  rownames(CA108_jac_score_m)[n] & CA108_sim_score$Amp2 ==  colnames(CA108_jac_score_m)[p], ]$JaccardBreakpoint)
        CA108_jac_score_m[p, n] <- as.numeric(CA108_sim_score[CA108_sim_score$Amp1 ==  rownames(CA108_jac_score_m)[n] & CA108_sim_score$Amp2 ==  colnames(CA108_jac_score_m)[p], ]$JaccardBreakpoint)
      }
      else{
        next
      }
    }
  }
}

CA108_upper_tri <- get_upper_tri(CA108_sim_score_m)
CA108_lower_tri <- get_lower_tri(CA108_jac_score_m)
CA108_correlation_similarity <- CA108_upper_tri
for (i in 1:nrow(CA108_correlation_similarity)){
  for (n in 1:ncol(CA108_correlation_similarity)){
    if(is.na(CA108_correlation_similarity[i,n])){
      CA108_correlation_similarity[i,n]<- CA108_lower_tri[i,n]
      
    }
    else{
      next
    }
  }
  
}
CA108_correlation_similarity

rownames(CA108_correlation_similarity) <-  do.call('rbind',strsplit(as.character(rownames(CA108_correlation_similarity)), '_' ,fixed=TRUE))[,1]
colnames(CA108_correlation_similarity) <- do.call('rbind',strsplit(as.character(colnames(CA108_correlation_similarity)), '_' ,fixed=TRUE))[,1]

#===============================================================================
#'  Saving the environment to reproduce the plots
#===============================================================================

save.image(file='Reproduce_Cascade_paper_environment.RData')
