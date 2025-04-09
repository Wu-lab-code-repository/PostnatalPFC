library(BSgenome)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(magrittr)
library(Rcpp)
library(ggplot2)
library(viridis)
library(batchelor)
sourceCpp(code='
          #include <Rcpp.h>
          using namespace Rcpp;
          using namespace std;
          // [[Rcpp::export]]
          IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
          if(x1.size() != y1.size()){
          stop("width must equal size!");
          }
          IntegerVector x = clone(x1);
          IntegerVector y = clone(y1);
          int n = x.size();
          IntegerVector rx = seq(xmin,xmax);
          IntegerVector ry = seq(ymin,ymax);
          IntegerMatrix mat( ry.size() , rx.size() );
          int xi,yi;
          for(int i = 0; i < n; i++){
          xi = (x[i] - xmin);
          yi = (y[i] - ymin);
          if(yi >= 0 && yi < ry.size()){
          if(xi >= 0 && xi < rx.size()){
          mat( yi , xi ) = mat( yi , xi ) + 1; 
          }
          }
          }
          return mat;
          }'
)


txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
feature <- txdb %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique

insertionProfileSingles <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
  
  insertionProfileSingles_helper <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
    #Convert To Insertion Sites
    if(getInsertions){
      insertions <- c(
        GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
        GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
      )
      by <- "RG"
    }else{
      insertions <- fragments
    }
    remove(fragments)
    gc()
    
    #center the feature
    center <- unique(resize(feature, width = 1, fix = fix, ignore.strand = FALSE))
    
    #get overlaps between the feature and insertions only up to flank bp
    overlap <- DataFrame(findOverlaps(query = center, subject = insertions, maxgap = flank, ignore.strand = TRUE))
    overlap$strand <- strand(center)[overlap[,1]]
    overlap$name <- mcols(insertions)[overlap[,2],by]
    overlap <- transform(overlap, id=match(name, unique(name)))
    ids <- length(unique(overlap$name))
    
    #distance
    overlap$dist <- NA
    minus <- which(overlap$strand == "-")
    other <- which(overlap$strand != "-")
    overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(insertions[overlap[minus,2]])
    overlap$dist[other] <- start(insertions[overlap[other,2]]) - start(center[overlap[other,1]])
    
    #Insertion Mat
    profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
    colnames(profile_mat) <- unique(overlap$name)
    profile <- rowSums(profile_mat)
    
    #normalize
    profile_mat_norm <- apply(profile_mat, 2, function(x) x/max(mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]), 0.5)) #Handles low depth cells
    profile_norm <- profile/mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])
    
    #smooth
    profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
    profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)
    
    #enrichment
    max_finite <- function(x){
      suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
    }
    e_mat <- apply(profile_mat_norm_smooth, 2, function(x) max_finite(x[(flank-range):(flank+range)]))
    names(e_mat) <- colnames(profile_mat_norm_smooth)
    e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])
    
    #Summary
    df_mat <- data.frame(
      enrichment = e_mat,
      insertions = as.vector(table(mcols(insertions)[,by])[names(e_mat)]),
      insertionsWindow = as.vector(table(overlap$name)[names(e_mat)])
    )
    df_sum <- data.frame(bp = (-flank):flank, profile = profile, norm_profile = profile_norm, smooth_norm_profile = profile_norm_smooth, enrichment = e)
    rownames(df_sum) <-  NULL
    
    return(list(df = df_sum, dfall = df_mat, profileMat = profile_mat_norm, profileMatSmooth = profile_mat_norm_smooth))
  }
  
  uniqueTags <- as.character(unique(mcols(fragments)[,by]))
  splitTags <- split(uniqueTags, ceiling(seq_along(uniqueTags)/batchSize))
  
  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  batchTSS <- lapply(seq_along(splitTags), function(x){
    setTxtProgressBar(pb, round(x * 100/length(splitTags), 0))
    profilex <- insertionProfileSingles_helper(
      feature=feature, 
      fragments=fragments[which(mcols(fragments)[,by] %in% splitTags[[x]])], 
      by = by, 
      getInsertions = getInsertions,
      fix = fix, 
      flank = flank, 
      norm = norm, 
      smooth = smooth, 
      range = range
    )
    
    return(profilex)
  })
  df <- lapply(batchTSS, function(x) x$df) %>% Reduce("rbind",.)
  dfall <- lapply(batchTSS, function(x) x$dfall) %>% Reduce("rbind",.)
  profileMat <- lapply(batchTSS, function(x) x$profileMat) %>% Reduce("cbind",.)
  profileMatSmooth <- lapply(batchTSS, function(x) x$profileMatSmooth) %>% Reduce("cbind",.)
  return(list(df = df, dfall = dfall, profileMat = profileMat, profileMatSmooth = profileMatSmooth))
}

#--------------------------------------------
# Input
#--------------------------------------------

minFrags <- 100
filterFrags <- 500
filterTSS <- 1
file_fragments <- "/~/ATACdata/Y16ATAC_fragments.tsv.gz"
out_fragments <- "/~/ATACdata/fragments/Y16ATAC_fragment.rds"
name <- "Y16"
#-----------------
# Reading Fragment Files
#-----------------
message("Reading in fragment files...")
fragments <- data.frame(readr::read_tsv(file_fragments, col_names=FALSE))
posna <- grep('TRUE',is.na(fragments$X1))
#fragments=fragments[-posna,] 
fragments <- GRanges(
  seqnames = fragments[,1], 
  IRanges(fragments[,2]+1, fragments[,3]), 
  RG = fragments[,4], 
  N = fragments[,5]
)

message("Filtering Lowly Represented Cells...")
tabRG <- table(fragments$RG)
keep <- names(tabRG)[which(tabRG >= minFrags)]
fragments <- fragments[fragments$RG %in% keep,]
fragments <- sort(sortSeqlevels(fragments))

#-----------------
# TSS Profile
#-----------------
feature <- txdb %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique
tssProfile <- insertionProfileSingles(feature = feature, fragments = fragments, 
                                      getInsertions = TRUE, batchSize = 1000)
tssSingles <- tssProfile$dfall
tssSingles$uniqueFrags <- 0
tssSingles[names(tabRG),"uniqueFrags"] <- tabRG
tssSingles$cellCall <- 0
#tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags & tssSingles$enrichment >= filterTSS] <- 1

#tssSingles_raw <- tssSingles

#pos <- grep('TRUE',is.na(tssSingles$enrichment))

#tssSingles <- tssSingles[-pos,]
dim(tssSingles)
#CAUTION
pos <- grep('TRUE',tssSingles$enrichment>2)
length(pos)
pos1 <- grep('TRUE',tssSingles$uniqueFrags>filterFrags)
length(pos1)

tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags & tssSingles$enrichment >= filterTSS] <- 1


dim(tssSingles)

tssSingles <- tssSingles[complete.cases(tssSingles),]




#-----------------
# Plot Stats
#-----------------

nPass  <- sum(tssSingles$cellCall==1)
nTotal <- sum(tssSingles$uniqueFrags >= filterFrags)

#Filter
fragments <- fragments[mcols(fragments)$RG %in% rownames(tssSingles)[tssSingles$cellCall==1]]
fragments$RG <- paste0(name,"#",fragments$RG)
fragments
#Save
saveRDS(fragments, out_fragments)

#######################
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg19)
set.seed(1)
#FUNCTION

countInsertions <- function(query, fragments, by = "RG"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

seuratLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  #TF IDF LSI adapted from flyATAC
  cs <- Matrix::colSums(mat)
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    mat <- mat[head(order(Matrix::rowSums(mat),decreasing = TRUE),nFeatures),] 
  }
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + ncol(mat) / Matrix::rowSums(mat)), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  #Make Seurat Object
  message("Making Seurat Object...")
  mat <- mat[1:100,] + 1
  obj <- CreateSeuratObject(mat, project='Human_M1', min.cells=0, min.features =0)
  obj[['pca']] <- CreateDimReducObject(embeddings = matSVD,key = 'pca',assay = DefaultAssay(obj))
  
  return(obj)
}

addClusters <- function(obj, minGroupSize = 50, dims.use = seq_len(50), initialResolution = 0.8){
  #First Iteration of Find Clusters
  currentResolution <- initialResolution
  obj <- FindNeighbors(obj,dims.use=dims.use,reduction = 'pca')
  obj <- FindClusters(object = obj,  resolution = currentResolution, print.output = FALSE)
  minSize <- min(table(obj@meta.data[[paste0("res.",currentResolution)]]))
  nClust <- length(unique(paste0(obj@meta.data[[paste0("res.",currentResolution)]])))
  message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
  #If clusters are smaller than minimum group size
  while(minSize <= minGroupSize){
    obj@meta.data <- obj@meta.data[,-which(colnames(obj@meta.data)==paste0("res.",currentResolution))]
    currentResolution <- currentResolution*initialResolution
    obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, resolution = currentResolution, print.output = FALSE, force.recalc = TRUE)
    minSize <- min(table(obj@meta.data[[paste0("res.",currentResolution)]]))
    nClust <- length(unique(paste0(obj@meta.data[[paste0("res.",currentResolution)]])))
    message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
  }
  return(obj)
}

extendedPeakSet <- function(df, BSgenome = NULL, extend = 150, blacklist = NULL, nSummits = 100000){
  #Helper Functions
  readSummits <- function(file){
    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score"))))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }
  nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
    }
    if(verbose){
      message("Converging", appendLF = FALSE)
    }
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)
  }
  #Check-------
  #stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #------------
  #input
  #Time to do stuff
  chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  groups <- unique(df$groups)
  groupGRList <- GenomicRanges::GenomicRangesList(lapply(seq_along(groups), function(i){
    df_group = df[which(df$groups==groups[i]),]
    grList <- GenomicRanges::GenomicRangesList(lapply(paste0(df_group$summits), function(x){
      extended_summits <- readSummits(x) %>%
        resize(., width = 2 * extend + 1, fix = "center") %>%     
        subsetByOverlaps(.,chromSizes,type="within") %>%
        subsetByOverlaps(.,blacklist,invert=TRUE) %>%
        nonOverlappingGRanges(., by="score", decreasing=TRUE)
      extended_summits <- extended_summits[order(extended_summits$score,decreasing=TRUE)]
      if(!is.null(nSummits)){
        extended_summits <- head(extended_summits, nSummits)
      }
      mcols(extended_summits)$scoreQuantile <- trunc(rank(mcols(extended_summits)$score))/length(mcols(extended_summits)$score)
      extended_summits
    }))
    #Non Overlapping
    grNonOverlapping <- nonOverlappingGRanges(unlist(grList), by = "scoreQuantile", decreasing = TRUE)
    #Free Up Memory
    remove(grList)
    gc()
    grNonOverlapping
  }))
  grFinal <- nonOverlappingGRanges(unlist(groupGRList), by = "scoreQuantile", decreasing = TRUE)
  grFinal <- sort(sortSeqlevels(grFinal))
  return(grFinal)
}

groupSums <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}


#-------------------------------------------------------------------------------------------------
# Start
#-------------------------------------------------------------------------------------------------
fragmentFiles <- list.files("/~/ATACdata/fragments/", pattern = ".rds", full.names = TRUE)
#-------------------------------------------------------------------------------------------------
# Get Counts
#-------------------------------------------------------------------------------------------------

genomeH <- BSgenome.Hsapiens.UCSC.hg19
chromSizes <- GRanges(names(seqlengths(genomeH)), IRanges(1, seqlengths(genomeH)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
windows <- unlist(tile(chromSizes, width = 2500))
countsList <- lapply(seq_along(fragmentFiles), function(i){
  message(sprintf("%s of %s", i, length(fragmentFiles)))
  counts <- countInsertions(windows, readRDS(fragmentFiles[i]), by = "RG")[[1]]
  counts
})
mat <- lapply(countsList, function(x) x) %>% Reduce("cbind",.)

row_n=paste0(seqnames(windows),'-',ranges(windows))
head(row_n)
rownames(mat)<-row_n
#-------------------------------------------------------------------------------------------------
# Run LSI Clustering with Seurat
#-------------------------------------------------------------------------------------------------
message("Making Seurat LSI Object...")
obj <- seuratLSI(mat, nComponents = 25, nFeatures = 20000)
message("Adding Graph Clusters...")
obj <- addClusters(obj, dims.use = 2:25, minGroupSize = 200, initialResolution = 1)
clusterResults <- split(rownames(obj@meta.data), paste0("Cluster",obj@meta.data$RNA_snn_res.1))
#-------------------------------------------------------------------------------------------------
# Get Cluster Beds
#-------------------------------------------------------------------------------------------------
dirClusters <- "/~/ATACdata/Y16ATAC/LSI-Cluster-Beds/"
dir.create(dirClusters)
for(i in seq_along(fragmentFiles)){
  fragments <-readRDS(fragmentFiles[i])
  for(j in seq_along(clusterResults)){
    message(sprintf("%s of %s", j, length(clusterResults)))
    fragmentsj <- fragments[fragments$RG %in% clusterResults[[j]]]
    if(length(fragmentsj) > 0){
      out <- data.frame(
        chr = c(seqnames(fragmentsj), seqnames(fragmentsj)), 
        start = c(as.integer(start(fragmentsj) - 1), as.integer(end(fragmentsj) - 1)), 
        end = c(as.integer(start(fragmentsj)), as.integer(end(fragmentsj)))
      ) %>% readr::write_tsv(
        x = ., 
        append = TRUE, 
        path = paste0(dirClusters, paste0(names(clusterResults)[j], ".bed")), 
        col_names = FALSE)
    }
  }
}


#-------------------------------------------------------------------------------------------------
# Run MACS2
#-------------------------------------------------------------------------------------------------
dirPeaks <- "/~/ATACdata/Y16ATAC/LSI-Cluster-Peaks/"
dir.create(dirPeaks)
method <- "q"
cutoff <- 0.05
shift <- -75
extsize <- 150
genome_size <- 2.86e9
for(j in seq_along(clusterResults)){
  message(sprintf("%s of %s", j, length(clusterResults)))
  clusterBedj <- paste0(dirClusters,names(clusterResults)[j],".bed")
  cmdPeaks <- sprintf(
    "/~/ENVmacs2/bin/macs2 callpeak -g %s --name %s --treatment %s --outdir %s -B --format BED --nomodel --call-summits --nolambda --keep-dup all", 
    genome_size, 
    names(clusterResults)[j], 
    clusterBedj, 
    dirPeaks
  )
  if (!is.null(shift) & !is.null(extsize)) {
    cmdPeaks <- sprintf("%s --shift %s --extsize %s", cmdPeaks, shift, extsize)
  }
  if (tolower(method) == "p") {
    cmdPeaks <- sprintf("%s -p %s", cmdPeaks, cutoff)
  }else {
    cmdPeaks <- sprintf("%s -q %s", cmdPeaks, cutoff)
  }
  message("Running Macs2...")
  message(cmdPeaks)
  system(cmdPeaks, intern = TRUE)
}

dirPeaks
#-------------------------------------------------------------------------------------------------
# Make Non-Overlapping Peak Set
#-------------------------------------------------------------------------------------------------
df <- data.frame(
  samples = gsub("\\_summits.bed","",list.files(dirPeaks, pattern = "\\_summits.bed", full.names = FALSE)),
  groups = "Y16",
  summits = list.files(dirPeaks, pattern = "\\_summits.bed", full.names = TRUE)
)

unionPeaks <- extendedPeakSet(
  df = df,
  BSgenome = genomeH, 
  extend = 250,
  nSummits = 200000
)
unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% c(levels(seqnames(unionPeaks)))] 
unionPeaks <- keepSeqlevels(unionPeaks,c(levels(seqnames(unionPeaks))))

#Create Counts list
countsPeaksList <- lapply(seq_along(fragmentFiles), function(i){
  message(sprintf("%s of %s", i, length(fragmentFiles)))
  gc()
  countInsertions(unionPeaks, readRDS(fragmentFiles[i]), by = "RG")
})


#CountsMatrix
mat <- lapply(countsPeaksList, function(x) x[[1]]) %>% Reduce("cbind",.)
frip <- lapply(countsPeaksList, function(x) x[[2]]) %>% unlist
total <- lapply(countsPeaksList, function(x) x[[3]]) %>% unlist

se <- SummarizedExperiment(
  assays = SimpleList(counts = mat), 
  rowRanges = unionPeaks
)
rownames(se) <- paste(seqnames(se),start(se),end(se),sep="-")
colData(se)$FRIP <- frip
colData(se)$uniqueFrags <- total / 2

mat <- assay(se)
mat@x[mat@x > 0] <- 1
save(mat,file='/~/ATACdata/Y16ATACmat.Robj')

#######################
library(Signac)
library(Seurat)
library(BSgenome)
library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v75)

load(file='/~/ATACdata/Y16ATACmat.Robj')
metadata <- read.csv(file = "/~/ATACdata/Y16ATAC_singlecell.csv",header = TRUE,row.names = 1)
cell <- colnames(mat)
cell <- gsub('Y16#','',cell)
pos <- which(rownames(metadata)%in%cell)
metadata <- metadata[pos,]
chrom_assay <- CreateChromatinAssay(
  counts = mat,
  sep = c("-", "-"),
  genome = 'hg19',
  fragments = '/~/ATACdata/Y16ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200)

ATACY16 =CreateSeuratObject(counts = chrom_assay,assay = "peaks",meta.data = metadata)
ATACY16@meta.data['sample_name'] <- rownames(ATACY16@meta.data)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
Annotation(ATACY16) <- annotations
ATACY16 <- NucleosomeSignal(object = ATACY16)
ATACY16 <- TSSEnrichment(object = ATACY16, fast = FALSE)
ATACY16$pct_reads_in_peaks <- ATACY16$peak_region_fragments / ATACY16$passed_filters * 100
ATACY16$blacklist_ratio <- ATACY16$blacklist_region_fragments / ATACY16$peak_region_fragments
ATACY16$high.tss <- ifelse(ATACY16$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(ATACY16, group.by = 'high.tss') + NoLegend()
ATACY16$nucleosome_group <- ifelse(ATACY16$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = ATACY16, group.by = 'nucleosome_group')
VlnPlot(object = ATACY16,features = c('pct_reads_in_peaks', 'peak_region_fragments',
                                      'TSS.enrichment', 'nucleosome_signal',' nFeature_peaks', 'unmapped'), pt.size = 0.1,ncol = 6)

ATACY16 <- subset(x = ATACY16,
                  subset = peak_region_fragments > 350 &
                    peak_region_fragments < 35000 &
                    pct_reads_in_peaks > 20 &
                    nucleosome_signal < 4 &
                    TSS.enrichment > 2 &
                    nFeature_peaks > 4000 &
                    unmapped < 5000)

ATACY16 <- RunTFIDF(ATACY16)
ATACY16 <- FindTopFeatures(ATACY16)
ATACY16 <- RunSVD(ATACY16)
DepthCor(ATACY16,n=50)
DefaultAssay(ATACY16)<-'peaks'
ATACY16 <- RunUMAP(object = ATACY16, reduction = 'lsi', dims = 2:22)
ATACY16 <- FindNeighbors(object = ATACY16, reduction = 'lsi', dims = 2:22)
ATACY16 <- FindClusters(object = ATACY16, verbose = FALSE, algorithm = 3,resolution = c(0.5,1,2,5,10))
DimPlot(object = ATACY16, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(ATACY16)
ATACY16[['RNA']] <- CreateAssayObject(counts = gene.activities)
ATACY16 <- NormalizeData(
  object =ATACY16,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATACY16$nCount_RNA))

