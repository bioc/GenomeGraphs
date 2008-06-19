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
                   idColor = "white")))

setMethod("initialize", "AnnotationTrack", function(.Object, ...) {
    .Object <- callNextMethod()

    if (!all(.Object@columns %in% colnames(.Object@regions))) {
        stop(cat("problem initializing AnnotationTrack need the following columns:",
                 paste(.Object@columns, collpase = ", ")), "\n")
    }
    return(.Object)
})

geneBiomart <- function(id, biomart, type = "ensembl_gene_id", dp = NULL) {
    ens <- getBM(c("structure_gene_stable_id", "structure_transcript_stable_id", "structure_exon_stable_id",
                   "structure_exon_chrom_start", "structure_exon_chrom_end", "structure_exon_rank",
                   "structure_transcript_chrom_strand", "structure_biotype"),
                 filters = type, values = id, mart = biomart)
    
    dp <- dpConcat(dp, DisplayPars(size = 1, color = "orange", plotId = FALSE, idRotation = 90, idColor = "white"))
    new("AnnotationTrack", chr = 0, strand = -1,
        regions = data.frame(start = ens[,4], end = ens[,5], feature = ens[,8], group = ens[,1], ID = ens[,3]),
        dp = dp)
}

geneRegionBiomart <- function(chr, start, end, strand, biomart, dp = NULL, chrFunction = function(x) x,
                              strandFunction = function(x) x) {
    ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id", "structure_exon_stable_id",
                   "structure_exon_chrom_start","structure_exon_chrom_end", "structure_exon_rank",
                   "structure_transcript_chrom_strand","structure_biotype"),
                 filters = c("chromosome_name", "start", "end", "strand"),
                 values = list(chrFunction(chr), start, end, strandFunction(strand)), mart = biomart)
    
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
    
    new("AnnotationTrack", chr = chr, strand = strand,
        regions = data.frame(start = ens[,4], end = ens[,5], feature = ens[,8], group = ens[,1], ID = ens[,3]),
        dp = dp)
}

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
    .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","structure_exon_stable_id",
                           "structure_exon_chrom_start","structure_exon_chrom_end","structure_exon_rank",
                           "structure_transcript_chrom_strand", "structure_biotype"),
                         filters = .Object@type, values=.Object@id, mart=.Object@biomart)

    if (is.null(.Object@ens)) {
        setPar(.Object, "size", 0)
    }
    .Object
})

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
        .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id",
                               "structure_exon_stable_id","structure_exon_chrom_start","structure_exon_chrom_end",
                               "structure_exon_rank", "structure_transcript_chrom_strand","structure_biotype"),
                             filters=c("chromosome_name", "start", "end", "strand"),
                             values=list(.Object@chromosome,.Object@start, .Object@end, strand),
                             mart=.Object@biomart)
    }
    
    if (is.null(.Object@ens)) {
        setPar(.Object, "size", 0)
    }
    .Object
})

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
    .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","structure_exon_stable_id",
                           "structure_exon_chrom_start","structure_exon_chrom_end","structure_exon_rank",
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

setClass("TranscriptRegion", contains = "gdObject", 
         representation(start = "numeric",
                        end = "numeric",
                        chromosome = "character",
                        biomart = "Mart",
                        ens = "data.frame"),
         prototype = prototype(dp = DisplayPars(size = 1))
         );

setClass("Ideogram", contains = "gdObject", 
         representation(chromosome = "character"),
         prototype(dp = DisplayPars(size = 1, color = "firebrickred3"))
         );

setClass("Title", contains = "gdObject", 
         representation(title = "character"),
         prototype(dp = DisplayPars(size = 1, cex = 1,
                   color = "black"))
         );

setClass("Legend", contains = "gdObject",
         representation(legend = "character"),
         prototype(dp = DisplayPars(size = 1,
                   cex = 1,
                   color = "black"))
         );

setClass("GenomeAxis", contains = "gdObject", 
         representation(add53 = "logical",
                        add35 = "logical",
                        littleTicks = "logical"),
         prototype(add53 = FALSE,
                   add35 = FALSE,
                   littleTicks = FALSE,
                   dp = DisplayPars(size = 1, color = "black"))
         );

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

setClass("GenericArray", contains = c("gdObject", "Segmentable"), 
         representation(intensity = "matrix",
                        probeStart = "numeric",
                        probeEnd = "numeric"),
         prototype(dp = DisplayPars(color = "darkred",
                   lty = "solid",
                   pch = 16,
                   pointSize = .2,
                   lwd = 1,
                   size = 5,
                   type = "point"), segmentation = NULL)
         );

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
                   ))
         );

setClass("GeneModel", contains = "gdObject", 
         representation(exonStart = "numeric",
                        exonEnd = "numeric",
                        chromosome = "numeric"),
         prototype(dp = DisplayPars(color = "darkgreen",
                   size = 1))
         );


setClass("BaseTrack", contains = c("gdObject", "Segmentable"), 
         representation(base = "numeric",
                        value = "numeric",
                        strand = "character"),
         prototype(strand = "+",
                   dp = DisplayPars(size = 5,
                   color = "orange",
                   lty = "solid",
                   lwd = 1), segmentation = NULL)
         );

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

setClassUnion("numericOrNull", c("numeric", "NULL"))
setClass("HighlightRegion", contains = "gdObject",
         representation(start = "numeric",
                        end = "numeric",
                        region = "numericOrNull",
                        coords = "character"),
         prototype(region = NULL,
                   coords = "genomic",
                   dp = DisplayPars(color = "black",  alpha = 0, lwd = 1, lty = "solid")
                   ));

