######################################################################
## Creates DisplayPar objects to be used in conjunction with the     #
## gdObjects.                                                        #
######################################################################
setClass("DisplayPars", representation(pars = "environment"))

setMethod("initialize", "DisplayPars", function(.Object) {
    .Object@pars <- new.env(parent = emptyenv())
    .Object
})

##-- write the arguments of 1 into 2 return 2. 
dpConcat <- function(dp1, dp2) {
    if (is.null(dp1) && is.null(dp2))
        return(DisplayPars())
    if (is.null(dp1)) 
        return(dp2)
    if (is.null(dp2))
        return(dp1)
    
    sapply(ls(dp1@pars), function(nm) {
        setPar(dp2, nm, getPar(dp1, nm))
    })
    return(dp2)
}

DisplayPars <- function(...) {
    args <- list(...)
    dp <- new("DisplayPars")
    
    i <- match(names(args), "dp")
    if (length(i) > 0 && !is.na(i)) {
        if (class(dpOrig <- args[[i]]) == "DisplayPars") {
            args <- args[-i]
            dp <- dpOrig
        }
    }

    for (i in seq(along = args)) {
        setPar(dp, names(args)[i], args[[i]])
    }
    return(dp)
}

##
## The methods for the DisplayPars class. (they need to be in here
## because they are called in prototypes)
##
setMethod("show", "DisplayPars", function(object) {
    sapply(ls(object@pars), function(x) {
        y <- getPar(object, x)
        if (length(y) > 10)
            y <- y[1:10]
        cat(x, " = ", toString(y), "\n")
    })
})

setGeneric("setPar", def = function(obj, name, val, ...) standardGeneric("setPar"))
setMethod("setPar", "DisplayPars", function(obj, name, val) {
    assign(name, val, obj@pars)
})

setGeneric("getPar", def = function(obj, name, ...) standardGeneric("getPar"))
setMethod("getPar", "DisplayPars", function(obj, name) {
    if (!exists(name, obj@pars))
        return(NULL)
    
    get(name, obj@pars)
})

setClassUnion("dfOrNULL", c("data.frame", "NULL"))

##
## The gdObject is the parent of all GenomeGraphs objects in the system
##
setClass("gdObject", representation = representation(dp = "DisplayPars"), 
         prototype = prototype(dp = DisplayPars()))

setMethod("getPar", "gdObject", function(obj, name) {
    getPar(obj@dp, name)
})

setMethod("setPar", "gdObject", function(obj, name, val) {
    setPar(obj@dp, name, val)
})

##
## The initialize method for the gdObjects is only complex because we want
## to preserve graphical parameters in a certain order. 
##
setMethod("initialize", "gdObject", function(.Object, ...) {
    ## save the old one. 
    protDP <- .Object@dp

    ## make a new one. 
    newDP <- DisplayPars(color = "black", size = 1)
    
    sapply(ls(protDP@pars), function(nm) {
        setPar(newDP, nm, getPar(protDP, nm))
    })
    
    args <- list(...)
    ii <- which("dp" == names(args))
    
    if (length(ii) > 0) {
        argDP <- args[[ii]]
        sapply(ls(argDP@pars), function(nm) {
            setPar(newDP, nm, getPar(argDP, nm))
        })
    }

    ## set the other things then reset. 
    .Object <- callNextMethod()
    .Object@dp <- newDP
    return(.Object)
})

setClass("AnnotationTrack", contains = "gdObject",
         representation(chr = "numeric", strand = "numeric", regions = "dfOrNULL"),
         prototype(columns = c("start", "end", "feature", "group", "ID"),
                   featureColumnName = "feature",
                   dp = DisplayPars(size = 1,
                   plotId = FALSE,
                   idRotation = 0,
                   idColor = "white"))
         )

setMethod("initialize", "AnnotationTrack", function(.Object, ...) {
    .Object <- callNextMethod()

    if (!is.null(.Object@regions) && !all(.Object@columns %in% colnames(.Object@regions))) {
        stop(cat("Problem initializing AnnotationTrack need the following columns:",
                 paste(.Object@columns, collpase = ", ")), "\n")
    }
    return(.Object)
})

makeAnnotationTrack <- function(regions = NULL, chr = NULL, strand = NULL, start = NULL,
                                end = NULL, feature = NULL, group = NULL, ID = NULL,
                                dp = NULL) {
    pt <- getClass("AnnotationTrack")@prototype
    if (is.null(dp)) dp <- pt@dp

    if (missing(regions)) {
        if (missing(start) || missing(end))
            stop("Must specify either regions or start and end.")
        if (length(start) != length(end))
            stop("Start and end must be vectors of the same length")

        if (length(start) > 0) {
            if (missing(feature))
                feature <- rep("unknown", length(start))
            if (missing(group))
                group <- 1:length(start)
            if (missing(ID))
                ID <- 1:length(start)

            regions <- data.frame(start = start, end = end, feature = feature, group = group, ID = ID)
        }
    }
    if (is.null(chr))
        chr <- 0
    if (is.null(strand))
        strand <- 0

    return(new("AnnotationTrack", chr = chr, strand = strand, regions = regions, dp = dp))
}

geneBiomart <- function(id, biomart, type = "ensembl_gene_id", dp = NULL) {
   # ens <- getBM(c("structure_gene_stable_id", "structure_transcript_stable_id", "ensembl_exon_id","exon_chrom_start", "exon_chrom_end", "rank", "structure_transcript_chrom_strand", "structure_biotype"),filters = type, values = id, mart = biomart)
   ens <- getBM(c("structure_gene_stable_id", "structure_transcript_stable_id", "ensembl_exon_id","exon_chrom_start", "exon_chrom_end", "rank", "structure_transcript_chrom_strand"),filters = type, values = id, mart = biomart)
   if(!is.null(ens)){
     ens <- cbind(ens, biotype=rep("protein_coding", length(ens[,1])))
   }
    dp <- dpConcat(dp, DisplayPars(size = 1, color = "orange", plotId = FALSE, idRotation = 90, idColor = "white"))
    makeAnnotationTrack(start = ens[,4], end = ens[,5], feature = ens[,8], group = ens[,1],
                        ID = ens[,3], dp = dp)
}

geneRegionBiomart <- function(chr, start, end, strand, biomart, dp = NULL, chrFunction = function(x) x,
                              strandFunction = function(x) x) {
   # ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id", "ensembl_exon_id","exon_chrom_start","exon_chrom_end", "rank","structure_transcript_chrom_strand","structure_biotype"),filters = c("chromosome_name", "start", "end", "strand"),values = list(chrFunction(chr), start, end, strandFunction(strand)), mart = biomart)
    ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id", "ensembl_exon_id","exon_chrom_start","exon_chrom_end", "rank","structure_transcript_chrom_strand"),filters = c("chromosome_name", "start", "end", "strand"),values = list(chrFunction(chr), start, end, strandFunction(strand)), mart = biomart)

    if(!is.null(ens)){
     ens <- cbind(ens, biotype=rep("protein_coding", length(ens[,1])))
   }
    dp <- dpConcat(dp, DisplayPars(size = 1, color = "orange", plotId = FALSE, idRotation = 90, idColor = "white",
                                   C_segment = "burlywood4",
                                   D_segment = "lightblue",
                                   J_segment = "dodgerblue2",
                                   miRNA     = "cornflowerblue",
                                   miRNA_pseudogene = "cornsilk",
                                   misc_RNA         = "cornsilk3",
                                   misc_RNA_pseudogene = "cornsilk4",
                                   Mt_rRNA = "yellow",
                                   Mt_tRNA = "darkgoldenrod",
                                   Mt_tRNA_pseudogene = "darkgoldenrod1",
                                   protein_coding = "gold4",
                                   pseudogene  = "brown1",         
                                   retrotransposed = "blueviolet",
                                   rRNA = "darkolivegreen1",
                                   rRNA_pseudogene = "darkolivegreen" ,   
                                   scRNA = "darkorange",
                                   scRNA_pseudogene = "darkorange2",
                                   snoRNA = "cyan",           
                                   snoRNA_pseudogene = "cyan2",
                                   snRNA = "coral",
                                   snRNA_pseudogene = "coral3",   
                                   tRNA_pseudogene = "antiquewhite3",
                                   V_segment =  "aquamarine", dp = dp))

    makeAnnotationTrack(chr = chr, strand = strand, start = ens[,4],
                        end = ens[,5], feature = ens[,8], group =
                        ens[,1], ID = ens[,3], dp = dp)
}

###############################################
setClass("Gene", contains = "gdObject",
         representation(id = "character",
                        type = "character",
                        biomart = "Mart",
                        ens = "dfOrNULL"),
         prototype(type = "ensembl_gene_id",
                   dp = DisplayPars(size = 1,
                   color = "orange",
                   plotId = FALSE,
                   idRotation = 90,
                   idColor = "white"))
         );

setMethod("initialize", "Gene", function(.Object, ...){
    .Object <- callNextMethod()
    #.Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","structure_transcript_chrom_strand", "structure_biotype"), filters = .Object@type, values=.Object@id, mart=.Object@biomart)
 .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id",
                           "exon_chrom_start","exon_chrom_end","rank","structure_transcript_chrom_strand"), filters = .Object@type, values=.Object@id, mart=.Object@biomart)
    if(!is.null(.Object@ens)){
     .Object@ens <- cbind(.Object@ens, biotype=rep("protein_coding", length(.Object@ens[,1])))
   }
    if (is.null(.Object@ens)) {
        setPar(.Object, "size", 0)
    }
    .Object
})

makeGene <- function(id, type, biomart, dp = NULL){
 if(missing(id)) stop("Need to specify a gene identifier for creating a Gene")
  pt <- getClass("Gene")@prototype
 if (is.null(dp))
   dp <- pt@dp
 if(missing(type))
   type=pt@type
 new("Gene", id = id, type = type, biomart = biomart, dp = dp)
}
###############################################
## Why can't I just use the setClassUnion Mechanism?
setClassUnion("MartOrNULL", c("Mart", "NULL"))

setClass("GeneRegion", contains = "gdObject",
         representation(start = "numeric",
                        end = "numeric",
                        chromosome = "character",
                        strand = "character",
                        biomart = "MartOrNULL",
                        ens = "dfOrNULL"),
         prototype(biomart = NULL,
                   strand = "+",
                   dp = DisplayPars(size = 1,
                   color = "orange",
                   plotId = FALSE,
                   idRotation = 90,
                   idColor = "white",
                   "C_segment" = "burlywood4",
                   "D_segment" = "lightblue",
                   "J_segment" = "dodgerblue2",
                   "miRNA"     = "cornflowerblue",
                   "miRNA_pseudogene" = "cornsilk",
                   "misc_RNA"         = "cornsilk3",
                   "misc_RNA_pseudogene" = "cornsilk4",
                   "Mt_rRNA" = "yellow",
                   "Mt_tRNA" = "darkgoldenrod",
                   "Mt_tRNA_pseudogene" = "darkgoldenrod1",
                   "protein_coding" = "gold4",
                   "pseudogene"  = "brown1",         
                   "retrotransposed" = "blueviolet",
                   "rRNA" = "darkolivegreen1",
                   "rRNA_pseudogene" = "darkolivegreen" ,   
                   "scRNA" = "darkorange",
                   "scRNA_pseudogene" = "darkorange2",
                   "snoRNA" = "cyan",           
                   "snoRNA_pseudogene" = "cyan2",
                   "snRNA" = "coral",
                   "snRNA_pseudogene" = "coral3",   
                   "tRNA_pseudogene" = "antiquewhite3",
                   "V_segment" =  "aquamarine"        
                   ))
         );

setMethod("initialize", "GeneRegion", function(.Object,...){
    .Object <- callNextMethod()
    strand  = switch(.Object@strand, "+" = 1, "-" = -1)
    
    ##-- changing start and end positions to capture genes on the edges. 
    .Object@start <- .Object@start - 2000
    .Object@end <- .Object@end + 2000

    if (!is.null(.Object@biomart)) {
       # .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank", "structure_transcript_chrom_strand","structure_biotype"),filters=c("chromosome_name", "start", "end", "strand"), values=list(.Object@chromosome,.Object@start, .Object@end, strand), mart=.Object@biomart)
        
 .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank", "structure_transcript_chrom_strand"),filters=c("chromosome_name", "start", "end", "strand"), values=list(.Object@chromosome,.Object@start, .Object@end, strand), mart=.Object@biomart)
         if(!is.null(.Object@ens)){
     .Object@ens <- cbind(.Object@ens, biotype=rep("protein_coding", length(.Object@ens[,1])))
   }
    }
    
    if (is.null(.Object@ens)) {
        setPar(.Object, "size", 0)
    }
    .Object
})

makeGeneRegion <- function(start, end, chromosome, strand, biomart, dp = NULL){
    if(missing(start)) stop("Need to specify a start for creating a GeneRegion")
    pt <- getClass("GeneRegion")@prototype
    if (is.null(dp))
        dp <- pt@dp
    if(is.numeric(chromosome))
        chromosome = as.character(chromosome)
    new("GeneRegion", start = start, end = end, chromosome = chromosome, strand = strand ,biomart = biomart, dp = dp)
}

###########################################
setClass("Transcript", contains = "gdObject", 
         representation(id = "character",
                        type = "character",
                        transcriptSize = "numeric",
                        numOfTranscripts = "numeric",
                        biomart = "Mart",
                        ens = "dfOrNULL"),
         
         prototype(type = "ensembl_gene_id",
                   transcriptSize = 1,
                   numOfTranscripts = 0,
                   dp = DisplayPars(size = 4,
                   color = "cornflowerblue",
                   plotId = FALSE))
         );

setMethod("initialize", "Transcript", function(.Object,...){
    .Object <- callNextMethod()
    .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id",
                           "exon_chrom_start","exon_chrom_end","rank",
                           "structure_transcript_chrom_strand","structure_biotype"),
                         filters = .Object@type, values=.Object@id,mart=.Object@biomart)

    if (is.null(.Object@ens)) {
        setPar(.Object@dp, "size", 0)
        .Object@numOfTranscripts <- 0
    }
    else {
        .Object@numOfTranscripts <- length(unique(.Object@ens[,2]))
    }
    .Object
})

makeTranscript <- function(id, type, biomart, dp = NULL){
 if(missing(id)) stop("Need to specify a gene identifier for creating a Transcript")
  pt <- getClass("Transcript")@prototype
 if (is.null(dp))
   dp <- pt@dp
 if(missing(type))
   type=pt@type
 new("Transcript", id = id, type = type, biomart = biomart, dp = dp)
}
########################################
setClass("TranscriptRegion", contains = "gdObject", 
         representation(start = "numeric",
                        end = "numeric",
                        chromosome = "character",
                        biomart = "Mart",
                        ens = "data.frame"),
         prototype = prototype(dp = DisplayPars(size = 1))
         );
########################################
setClass("Ideogram", contains = "gdObject", 
         representation(chromosome = "character"),
         prototype(dp = DisplayPars(size = 1, color = "firebrickred3"))
         );

makeIdeogram <- function(chromosome, dp = NULL){
 if(missing(chromosome)) stop("Need to specify chromosome for creating an Ideogram")
 if(is.numeric(chromosome)){
   chromosome = as.character(chromosome)
 }
 if (is.null(dp))
   dp <- getClass("Ideogram")@prototype@dp
 new("Ideogram", chromosome = chromosome, dp = dp)
}
########################################
setClass("Title", contains = "gdObject", 
         representation(title = "character"),
         prototype(dp = DisplayPars(size = 1, cex = 1,
                   color = "black"))
         );

makeTitle <- function(text, cex, color, size) {
  dp <- getClass("Title")@prototype@dp
  if (!missing(cex)) setPar(dp, "cex", cex)
  if (!missing(color)) setPar(dp, "color", color)
  if (!missing(size)) setPar(dp, "size", size)
  new("Title", title = text, dp = dp)
}
#######################################
setClass("Legend", contains = "gdObject",
         representation(legend = "character"),
         prototype(dp = DisplayPars(size = 1,
                   cex = 1,
                   color = "black"))
         );

makeLegend <- function(text, fill, cex) {
  dp <- getClass("Legend")@prototype@dp
  if (!missing(cex)) setPar(dp, "cex", cex)
  if (!missing(fill)) setPar(dp, "color", fill)
  new("Legend", legend = text, dp = dp)
}
#######################################
setClass("GenomeAxis", contains = "gdObject", 
         representation(add53 = "logical",
                        add35 = "logical",
                        littleTicks = "logical"),
         prototype(add53 = FALSE,
                   add35 = FALSE,
                   littleTicks = FALSE,
                   dp = DisplayPars(size = 1, color = "black",
                     cex=1, byValue = 1, distFromAxis = 1, labelPos = "alternating"))
         );

makeGenomeAxis <- function(add53 = FALSE, add35 = FALSE, littleTicks = FALSE, dp = NULL){
 if (is.null(dp))
   dp <- getClass("GenomeAxis")@prototype@dp
 new("GenomeAxis", add53 = add53, add35 = add35, dp = dp)
}
#######################################
setClass("Segmentation", contains = "gdObject",
         representation(segments = "list",
                        segmentStart = "list",
                        segmentEnd = "list"),
         prototype(dp = DisplayPars(color = "dodgerblue",
                   lwd = 1, lty = "solid"))
         );

setClassUnion("SegmentationOrNULL", c("Segmentation", "NULL"))
setClass("Segmentable", representation(segmentation = "SegmentationOrNULL"),
         prototype(segmentation = NULL))

makeSegmentation <- function(start, end, value, dp = NULL) {
  if (!is.list(value) && !is.list(start) && !is.list(end)) {
    start <- list(start)
    end <- list(end)
    value <- list(value)
  }
  if (is.null(dp))
    dp <- getClass("RectangleOverlay")@prototype@dp
  new("Segmentation", segments = value, segmentStart = start, segmentEnd = end, dp = dp)
}
#########################################
setClass("GenericArray", contains = c("gdObject", "Segmentable"), 
         representation(intensity = "matrix",
                        probeStart = "numeric",
                        probeEnd = "numeric"),
         prototype(
                   dp = DisplayPars(color = "darkred",
                   lty = "solid",
                   pch = 16,
                   pointSize = .2,
                   lwd = 1,
                   size = 5,
                   type = "point"), segmentation = NULL)
         );

makeGenericArray <- function(intensity, probeStart, probeEnd, segmentation, dp = NULL){
 pt <- getClass("GenericArray")@prototype 
 if (is.null(dp))
   dp <- pt@dp
 if(missing(probeEnd))
    probeEnd <- pt@probeEnd
 if(missing(segmentation))
   segmentation <- pt@segmentation
 if(missing(probeStart)) stop("Need probeStart argument to know where to plot the data on the genome")
 new("GenericArray", intensity = intensity, probeStart = probeStart, probeEnd = probeEnd, dp = dp, segmentation = segmentation)
}
#########################################
setClass("ExonArray", contains = "gdObject", 
         representation(intensity = "matrix",
                        probeStart = "numeric",
                        probeEnd = "numeric",
                        probeId = "character",
                        nProbes = "numeric",
                        displayProbesets = "logical"),
         prototype(dp = DisplayPars(size = 5,
                   displayProbesets = TRUE,
                   probesetSize = 1,
                   color = "firebrick1",
                   mapColor = "dodgerblue2",
                   plotMap = TRUE,
                   lwd = 1,
                   lty = "solid",
                   probeSetLwd = 1,
                   probeSetColor = "grey"
                   ), displayProbesets = FALSE)
         );

makeExonArray <- function(intensity, probeStart, probeEnd, probeId, nProbes, displayProbesets = FALSE, dp = NULL){
  pt <- getClass("ExonArray")@prototype 
  if (is.null(dp))
    dp <- pt@dp
  if(missing(probeEnd))
    probeEnd <- pt@probeEnd
  if(missing(probeId))
    probeId <- pt@probeId
  if(missing(nProbes))
    nProbes <- pt@nProbes
 if (is.null(dp))
   dp <- getClass("ExonArray")@prototype@dp
 new("ExonArray", intensity = intensity, probeStart = probeStart, probeEnd = probeEnd,probeId = probeId,
     nProbes = nProbes, displayProbesets = displayProbesets,dp = dp)
}

###########################################
setClass("GeneModel", contains = "gdObject", 
         representation(exonStart = "numeric",
                        exonEnd = "numeric",
                        chromosome = "numeric"),
         prototype(dp = DisplayPars(color = "darkgreen",
                   size = 1))
         );

makeGeneModel <- function(start, end, chromosome, dp = NULL){
 if (is.null(dp))
   dp <- getClass("GeneModel")@prototype@dp
 new("GeneModel", exonStart = start, exonEnd = end, dp = dp)
}
##########################################
setClass("BaseTrack", contains = c("gdObject", "Segmentable"), 
         representation(base = "numeric",
                        value = "numeric",
                        strand = "character"),
         prototype(strand = "+",
                   dp = DisplayPars(size = 5,
                   color = "orange",
                   lty = "solid",
                   type = "p", 
                   lwd = 1), segmentation = NULL)
         );

makeBaseTrack <- function(base, value, strand, segmentation, dp = NULL){
 pt <- getClass("BaseTrack")@prototype 
 if (is.null(dp))
   dp <- pt@dp
 if(missing(strand))
    strand <- pt@strand
 if(missing(segmentation))
   segmentation <- pt@segmentation
 if(missing(base)) stop("Need base argument to know the base positions to plot the data on the genome")
 if(missing(value)) stop("Need value argument")
 
 new("BaseTrack", base = base, value = value, strand = strand, dp = dp, segmentation = segmentation)
}
##############################################
setClass("MappedRead", contains = "gdObject", 
         representation(start = "numeric",
                        end = "numeric",
                        strand = "character",
                        chromosome = "character"),
         prototype(dp = DisplayPars(size = 5,
                   color = "orange",
                   lty = "solid",
                   lwd = 1))
         );




