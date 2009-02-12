setMethod("show",signature(object="Gene"),
  function(object){
    len = min(5, length(object@ens[,1]))
    res = paste("Object of class 'Gene':\n ID:",object@id,"\n Type:",object@type,"\n Exons in Ensembl: \n", sep="")
    cat(res)
    print(object@ens[1:len,])
    res = paste("\n There are ",length(object@ens[,1])-len," more rows", sep="")
    cat(res)
})

setMethod("show",signature(object="Transcript"),
  function(object){
    len = min(5, length(object@ens[,1]))
    res = paste("Object of class 'Transcript':\n ID:",object@id,"\n Type:",object@type,"\n Exons in Ensembl: \n", sep="")
    cat(res)
    print(object@ens[1:len,])
    res = paste("\n There are ",length(object@ens[,1])-len," more rows", sep="")
    cat(res)
})

setMethod("show",signature(object="GeneRegion"),
  function(object) {
    len <- min(5, length(object@ens[,1]))
    res <- paste("Object of class 'GeneRegion':\n Start:",object@start,"\n End:",
                 object@end,"\n Chromosome: ",object@chromosome,"\n Exons in Ensembl: \n", sep="")
    cat(res)
    print(object@ens[1:len,])
    res <- paste("\n There are ", length(object@ens[,1]) - len, " more rows", sep="")
    cat(res)
})

setMethod("show",signature(object="TranscriptRegion"),
  function(object){
    len <- min(5, length(object@ens[,1]))
    res <- paste("Object of class 'TranscriptRegion':\n Start:",object@start,"\n End:",
                 object@end,"\n Chromosome: ",object@chromosome,"\n Exons in Ensembl: \n", sep="")
    cat(res)
    print(object@ens[1:len,])
    res <- paste("\n There are ",length(object@ens[,1])-len," more rows", sep="")
    cat(res)
})

setMethod("show",signature(object="GenericArray"),
  function(object){
    len <- min(5, length(object@intensity[,1]))
    cat("Object of class 'GenericArray':\n ProbeStart: \n")
    print(object@probeStart[1:len])
    cat("Intensity: \n")
    print(object@intensity[1:len,])
    res = paste("\n There are ",length(object@intensity[,1])-len," more rows", sep="")
    cat(res)
    show(object@dp)
})

setMethod("show",signature(object="BaseTrack"),
  function(object){
    len <- min(5, length(object@base))
    cat("Object of class 'BaseTrack':\n base position: \n")
    print(object@base[1:len])
    cat("Values: \n")
    print(object@value[1:len])
    res = paste("\n There are ", length(object@base)-len," more rows", sep="")
    cat(res)
    show(object@dp)
})

setMethod("show",signature(object="ExonArray"),
  function(object){
    len = min(5, length(object@intensity[,1]))
    cat("Object of class 'ExonArray':\n ProbeStart: \n")
    print(object@probeStart[1:len])
    cat("ProbeEnd: \n")
    print(object@probeEnd[1:len])
    cat("Intensity: \n")
    print(object@intensity[1:len,])
    res = paste("\n There are ",length(object@intensity[,1])-len," more rows", sep="")
    cat(res)
    show(object@dp)
})

setGeneric("getColor",def=function(obj,...) standardGeneric("getColor"))
setMethod("getColor", signature("gdObject"), function(obj) {
    getPar(obj@dp, "color")
})

setGeneric("getBiotypeColor",def=function(obj, biotype,...) standardGeneric("getBiotypeColor"))
setMethod("getBiotypeColor", signature("gdObject"), function(obj,biotype) {
    getPar(obj@dp, biotype)
})

setGeneric("getSize",def=function(obj,...)standardGeneric("getSize"))
setMethod("getSize",signature("gdObject"),function(obj) {
    getPar(obj@dp, "size")
})

setGeneric("getID", def= function(obj,...) standardGeneric("getID"))
setMethod("getID", signature("Gene"), function(obj) obj@id)
setMethod("getID", signature("Transcript"), function(obj) obj@id)

##
## XXX: this is being used in two different ways
##
setGeneric("getType", def=function(obj,...) standardGeneric("getType"))
setMethod("getType", signature("Gene"),function(obj) obj@type)
setMethod("getType", signature("Transcript"),function(obj) obj@type)
setMethod("getType", signature("GenericArray"),function(obj) {
    getPar(obj@dp, "type")
})

setGeneric("getExonModel",def=function(obj,...)standardGeneric("getExonModel"))
setMethod("getExonModel",signature("Gene"),function(obj) obj@ens)
setMethod("getExonModel",signature("Transcript"),function(obj) obj@ens)
setMethod("getExonModel",signature("GeneRegion"),function(obj) obj@ens)
setMethod("getExonModel",signature("TranscriptRegion"),function(obj) obj@ens)

setGeneric("getBiomart",def=function(obj,...)standardGeneric("getBiomart"))
setMethod("getBiomart",signature("Gene"),function(obj) obj@biomart)
setMethod("getBiomart",signature("Transcript"),function(obj) obj@biomart)
setMethod("getBiomart",signature("GeneRegion"),function(obj) obj@biomart)
setMethod("getBiomart",signature("TranscriptRegion"),function(obj) obj@biomart)

setGeneric("getChromosome",def=function(obj,...)standardGeneric("getChromosome"))
setMethod("getChromosome",signature("GeneRegion"),function(obj) obj@chromosome)
setMethod("getChromosome",signature("TranscriptRegion"),function(obj) obj@chromosome)
setMethod("getChromosome",signature("Ideogram"),function(obj) obj@chromosome)


setGeneric("getTranscriptSize",def=function(obj,...)standardGeneric("getTranscriptSize"))
setMethod("getTranscriptSize",signature("Transcript"),function(obj) obj@transcriptSize)
setMethod("getTranscriptSize",signature("TranscriptRegion"),function(obj) obj@transcriptSize)

setGeneric("getNumOfTranscripts",def=function(obj,...)standardGeneric("getNumOfTranscripts"))
setMethod("getNumOfTranscripts",signature("Transcript"),function(obj) obj@numOfTranscripts)

setGeneric("getAdd53",def=function(obj,...)standardGeneric("getAdd53"))
setMethod("getAdd53",signature("GenomeAxis"),function(obj) obj@add53)

setGeneric("getAdd35",def=function(obj,...)standardGeneric("getAdd35"))
setMethod("getAdd35",signature("GenomeAxis"),function(obj) obj@add35)

setGeneric("plotLittleTicks", def=function(obj, ...) standardGeneric("plotLittleTicks"))
setMethod("plotLittleTicks", signature("GenomeAxis"), function(obj) obj@littleTicks)

setGeneric("getTitle",def=function(obj,...)standardGeneric("getTitle"))
setMethod("getTitle",signature("Title"),function(obj) obj@title)

setGeneric("getLegend",def=function(obj,...) standardGeneric("getLegend"))
setMethod("getLegend",signature("Legend"),function(obj) obj@legend)

setGeneric("getCex",def=function(obj,...) standardGeneric("getCex"))
setMethod("getCex", signature("gdObject"),function(obj) {
    getPar(obj@dp, "cex")
})

setGeneric("getPch",def = function(obj,...) standardGeneric("getPch"))
setMethod("getPch", signature("gdObject"), function(obj) {
    getPar(obj@dp, "pch")
})
setGeneric("getPointSize",def = function(obj,...) standardGeneric("getPointSize"))
setMethod("getPointSize", signature("gdObject"), function(obj) {
    getPar(obj@dp, "pointSize")
})

setGeneric("getLwd", def = function(obj,...) standardGeneric("getLwd"))
setMethod("getLwd", signature("gdObject"), function(obj) {
    getPar(obj@dp, "lwd")
})

setGeneric("getLty",def=function(obj,...)standardGeneric("getLty"))
setMethod("getLty", signature("gdObject"), function(obj) {
    getPar(obj@dp, "lty")
})

setGeneric("getProbeId",def=function(obj,...)standardGeneric("getProbeId"))
setMethod("getProbeId",signature("ExonArray"),function(obj) obj@probeId)

setGeneric("getProbeStart",def=function(obj,...)standardGeneric("getProbeStart"))
setMethod("getProbeStart",signature("ExonArray"),function(obj) obj@probeStart)
setMethod("getProbeStart",signature("GenericArray"),function(obj) obj@probeStart)

setGeneric("getExonStart",def=function(obj,...)standardGeneric("getExonStart"))
setMethod("getExonStart",signature("GeneModel"),function(obj) obj@exonStart)

setGeneric("getExonEnd",def=function(obj,...)standardGeneric("getExonEnd"))
setMethod("getExonEnd",signature("GeneModel"),function(obj) obj@exonEnd)

setGeneric("getProbeEnd",def=function(obj,...)standardGeneric("getProbeEnd"))
setMethod("getProbeEnd",signature("ExonArray"),function(obj) obj@probeEnd)
setMethod("getProbeEnd",signature("GenericArray"),function(obj) obj@probeEnd)

setGeneric("getIntensity",def = function(obj,...)standardGeneric("getIntensity"))
setMethod("getIntensity",signature("ExonArray"),function(obj) obj@intensity)
setMethod("getIntensity",signature("GenericArray"),function(obj) obj@intensity)

setGeneric("getSegmentation", def=function(obj,...) standardGeneric("getSegmentation"))
setMethod("getSegmentation", signature("Segmentable"), function(obj) {
    obj@segmentation
})

setGeneric("getSegments", def = function(obj, ...) standardGeneric("getSegments"))
setGeneric("getSegmentStart",def=function(obj,...)standardGeneric("getSegmentStart"))
setGeneric("getSegmentEnd",def=function(obj,...)standardGeneric("getSegmentEnd"))

setMethod("getSegments", signature("Segmentation"),function(obj) obj@segments)
setMethod("getSegmentStart",signature("Segmentation"),function(obj) obj@segmentStart)
setMethod("getSegmentEnd",signature("Segmentation"),function(obj) obj@segmentEnd)

setGeneric("getNprobes",def=function(obj,...)standardGeneric("getNprobes"))
setMethod("getNprobes",signature("ExonArray"),function(obj) obj@nProbes)

setGeneric("displayProbesets",def=function(obj,...)standardGeneric("displayProbesets"))
setMethod("displayProbesets",signature("ExonArray"),function(obj) obj@displayProbesets)

setGeneric("getMapColor",def=function(obj,...)standardGeneric("getMapColor"))
setMethod("getMapColor",signature("ExonArray"),function(obj) {
    getPar(obj@dp, "mapColor")
})

setGeneric("getProbeSetLwd",def=function(obj,...)standardGeneric("getProbeSetLwd"))
setMethod("getProbeSetLwd",signature("ExonArray"),function(obj) {
    getPar(obj@dp, "probeSetLwd")
})

setGeneric("getProbeSetColor",def=function(obj,...)standardGeneric("getProbeSetColor"))
setMethod("getProbeSetColor",signature("ExonArray"),function(obj) {
    getPar(obj@dp, "probeSetColor")
})

setGeneric("getBase",def=function(obj,...)standardGeneric("getBase"))
setMethod("getBase",signature("BaseTrack"),function(obj) obj@base)

setGeneric("getBaseValue", def=function(obj,...)standardGeneric("getBaseValue"))
setMethod("getBaseValue", signature("BaseTrack"), function(obj) obj@value)

setGeneric("getStart",def=function(obj,...)standardGeneric("getStart"))
setMethod("getStart",signature("MappedRead"),function(obj) obj@start)

setGeneric("getEnd",def=function(obj,...)standardGeneric("getEnd"))
setMethod("getEnd",signature("MappedRead"),function(obj) obj@end)

setGeneric("getStrand",def=function(obj,...)standardGeneric("getStrand"))
setMethod("getStrand",signature("MappedRead"),function(obj) obj@strand)

setGeneric("getPlotId",def=function(obj,...)standardGeneric("getPlotId"))

setMethod("getPlotId",signature("Transcript"),function(obj){
  getPar(obj@dp,"plotId")
})
setMethod("getPlotId",signature("Gene"),function(obj){
    getPar(obj@dp,"plotId")
})
setMethod("getPlotId",signature("GeneRegion"),function(obj){
    getPar(obj@dp,"plotId")
})
setMethod("getPlotId",signature("AnnotationTrack"),function(obj){
    getPar(obj@dp,"plotId")
})

setGeneric("getPlotMap",def=function(obj,...)standardGeneric("getPlotMap"))
setMethod("getPlotMap",signature("ExonArray"),function(obj){
    getPar(obj@dp,"plotMap")
})

##
## I could imagine that this gets refactored into gdObject, however for
## now leave it like this. 
##
setGeneric("getGenomicRange", def=function(obj,...) standardGeneric("getGenomicRange"))
setMethod("getGenomicRange",signature("gdObject"),function(obj){
    return(c(NA, NA))
})

setMethod("getGenomicRange", signature("GeneRegion"), function(obj){
    return(c(obj@start, obj@end))
})
setMethod("getGenomicRange", signature("TranscriptRegion"), function(obj){
    return(c(obj@start, obj@end))
})

.getRangeFromENS <- function(obj) {
    if (!is.null(obj@ens)) {
        c(min(obj@ens[,"exon_chrom_start"]), max(obj@ens[,"exon_chrom_end"]))
    }
    else {
        callNextMethod()
    }
}
setMethod("getGenomicRange", signature("Gene"), function(obj){
    .getRangeFromENS(obj)
})
setMethod("getGenomicRange", signature("Transcript"), function(obj){
    .getRangeFromENS(obj)
})
setMethod("getGenomicRange", signature("BaseTrack"), function(obj){
    c(min(obj@base), max(obj@base))
})
setMethod("getGenomicRange", signature("BaseTrack"), function(obj){
    c(min(obj@base), max(obj@base))
})
setMethod("getGenomicRange", signature("GenericArray"), function(obj){
    c(min(obj@probeStart), ifelse(length(obj@probeEnd) > 0, max(obj@probeEnd), max(obj@probeStart)))
})
setMethod("getGenomicRange", signature("ExonArray"), function(obj){
    c(min(obj@probeStart), ifelse(length(obj@probeEnd) > 0, max(obj@probeEnd), max(obj@probeStart)))
})

##
## The drawGD methods. 
##
setGeneric("drawGD", def=function(gdObject, minBase, maxBase, vpPosition, ...) standardGeneric("drawGD"))

setMethod("drawGD", signature("AnnotationTrack"), function(gdObject, minBase, maxBase, vpPosition) {
    regions <- gdObject@regions

    if (is.null(regions) || nrow(regions) == 0) {
        return(NULL)
    }
    pushViewport(dataViewport(xData=c(minBase, maxBase), extension = 0,
                              clip = TRUE, yscale = c(0, 40),
                              layout.pos.col = 1, layout.pos.row = vpPosition))
    
    for(i in seq(length = nrow(regions))) {
        color <- getPar(gdObject, as.character(regions[i, gdObject@featureColumnName]))
        if (is.null(color))
            color <- getPar(gdObject, "defaultFeatureColor")
        
        grid.rect(regions[i, "end"], 5, width = regions[i, "end"] - regions[i, "start"],
                  height = 30, gp = gpar(col = "black", fill = color),
                  default.units = "native", just = c("right", "bottom"))

        if (getPlotId(gdObject)) {
            rot <- getPar(gdObject, "idRotation")
            col <- getPar(gdObject, "idColor")

            if (!is.null(nm <- regions[i, "ID"])) {
                grid.text(nm, (regions[i, "start"] + regions[i, "end"])/2, 15,
                          rot = rot, gp = gpar(col=col, cex = getCex(gdObject)),
                          default.units = "native", just = c("center", "center"))
            }
        }
    }
    groups <- split(regions, regions[, "group"])
    
    lapply(Filter(function(x) nrow(x) > 1, groups), function(group) {
        ord <- order(group[,"start"])
        group <- group[ord, ]

        for (i in seq(length = (nrow(group) - 1))) {
            start <- group[i,"end"]
            end <- group[i + 1, "start"]
            mid <- start + (end - start)/4
            pp <- 15

            color <- getPar(gdObject, as.character(regions[i, gdObject@featureColumnName]))
            if (is.null(color))
                color <- getPar(gdObject, "defaultFeatureColor")
            
            grid.lines(c(start, mid), c(20, pp), default.units = "native",
                       gp = gpar(col = color, just = c("right", "bottom")))
            grid.lines(c(mid, end), c(pp, 20), default.units = "native",
                       gp = gpar(col = color, just = c("right", "bottom")))
        }
    })
    popViewport(1)
})

################################
# Plots Ensembl Gene models    #
################################
.drawGene <- function(gdObject, minBase, maxBase, vpPosition) {
    ens <- getExonModel(gdObject)
    if (is.null(ens)) {
        warning("No genes in gene region.")
        return(NULL)
    }
    pushViewport(dataViewport(xData=c(minBase, maxBase), extension = 0,
                              clip = TRUE, yscale = c(0, 40),
                              layout.pos.col=1, layout.pos.row = vpPosition))
    
    if(ens[1,6] == 1) { pp = 25 }
    else { pp = 15 }

    for(i in seq(along=ens[,1])) {
 
        color <- getColor(gdObject)
       
        if (!is.null(getBiotypeColor(gdObject, as.character(ens[i,8]))))
            color <- getBiotypeColor(gdObject, as.character(ens[i,8]))
        grid.rect(ens[i, 5], 5, width = ens[i, 5] - ens[i, 4], height = 30,
                  gp=gpar(col = "black", fill = color), default.units="native",
                  just = c("right", "bottom"))
     
        if (getPlotId(gdObject)) {
            rot <- getPar(gdObject, "idRotation")
            col <- getPar(gdObject, "idColor")
            grid.text(ens[i,1], (ens[i,4] + ens[i, 5])/2, 15, rot = rot, gp = gpar(col=col, cex = getCex(gdObject)),
                      default.units = "native", just = c("center", "center"))
        }
    }
    
    genes = unique(ens[,1])
    for(g in seq(along=genes)){
        color = getColor(gdObject)
        exons = ens[ens[,1]==genes[g],-c(1,2,3,6)]
        ord = order(exons[,1])
        exons = exons[ord,]
        if(!is.null(getBiotypeColor(gdObject,as.character(exons[1,4])))) color = getBiotypeColor(gdObject,as.character(exons[1,4]))
        for(j in seq(along=exons[,1])){
            if(j < length(exons[,1])){
                grid.lines(c(exons[j,2],exons[j,2]+((exons[j+1,1] - exons[j,2])/4)),c(20,pp),
                           default.units = "native", gp=gpar(col = color, just = c("right", "bottom")))
                grid.lines(c(exons[j,2]+((exons[j+1,1] - exons[j,2])/4), exons[j+1,1]),c(pp,20),
                           default.units = "native", gp=gpar(col = color, just = c("right", "bottom")))
            }
        }
    }
    popViewport(1)
}

setMethod("drawGD", signature("Gene"), .drawGene)
setMethod("drawGD", signature("GeneRegion"), .drawGene)

#############################
#Plots Ensembl transcripts  #
#############################

setMethod("drawGD", signature("Transcript"), function(gdObject, minBase, maxBase, vpPosition) {
    ens = getExonModel(gdObject)
    
    if (is.null(ens)) {
        warning("No transcripts in transcript region.")
        return(NULL)
    }
    
    size = getSize(gdObject) / max(getNumOfTranscripts(gdObject),getSize(gdObject)/getTranscriptSize(gdObject))
    vplayout = rep(size,getNumOfTranscripts(gdObject))
    color = getColor(gdObject)
    pushViewport(viewport(layout=grid.layout(length(vplayout), 1, height=vplayout),
                          layout.pos.col=1, layout.pos.row = vpPosition))

    transcripts = unique(ens[,2])
    for(i in seq(along = transcripts)){
        trans = ens[ens[,2] == transcripts[i],]
        pushViewport(dataViewport(xData=c(minBase,maxBase), yData=c(0,50), extension=0,
                                  layout.pos.col=1, layout.pos.row = i))
        if(getPlotId(gdObject)){
          grid.text(label=formatC(trans[1,2], format="d"), x = minBase, y = 30, just = c("left", "bottom"),
                    gp = gpar(cex = getCex(gdObject)), default.units = "native")
        } 
        for(i in seq(along=trans[,1])){
          grid.rect(trans[i,4],5,width=trans[i,4]-trans[i,5],height=25,gp=gpar(col = "black",fill = color),
                    default.units="native", just=c("right","bottom"))
        }
        if(trans[1,6] == 1){ pp = 16}
        else{ pp = 8}
        exons = unique(trans[,-2])
        ord = order(exons[,3])
        exons = exons[ord,]
        for(j in seq(along=exons[,1])){
            if(j < length(exons[,1])){
                grid.lines(c(exons[j,4],exons[j,4]+((exons[j+1,3] - exons[j,4])/4)),c(12,pp),
                           default.units = "native",gp=gpar(col = color))
                grid.lines(c(exons[j,4]+((exons[j+1,3] - exons[j,4])/4), exons[j+1,3]),c(pp,12),
                           default.units = "native", gp=gpar(col = color))
            }
        }
        popViewport()
    }
    popViewport()
})

#################################
#Plots a non-Ensembl gene model #
#################################
setMethod("drawGD", signature("GeneModel"), function(gdObject, minBase, maxBase, vpPosition) {
    model = cbind(getExonStart(gdObject),getExonEnd(gdObject))
    pushViewport(dataViewport(xData=c(minBase,maxBase), yscale=c(0,40), extension=0,
                              layout.pos.col=1, layout.pos.row=vpPosition))
    col = getColor(gdObject)
    for(i in seq(along=model[,1])){
        grid.rect(model[i,1],5,width=model[i,1]-model[i,2],height=30,
                  gp = gpar(col = "black",fill = col), default.units="native", just = c("right", "bottom"))
    }
    for(j in seq(along=model[,1])){
        if(j < length(model[,1])){
            grid.lines(c(model[j,2],model[j+1,1]),c(20,20),
                       default.units = "native",gp=gpar(col = col,just = c("right", "bottom")))
        }
    }  
    popViewport()
})
          
##############################
#Plots an exon module        #
##############################
exonModule <- function(position,  col = "burlywood3"){
    pushViewport(dataViewport(xData=c(0,20), yData=c(0,20), extension=0,
                              layout.pos.col = position, layout.pos.row=1))
    grid.rect(5,5,width=10,height=10,gp=gpar(col = col,fill = col), default.units="native",
              just=c("left","bottom"))
    grid.lines(c(0,5),c(10,10), default.units = "native",gp=gpar(col = col))
    grid.lines(c(15,20),c(10,10), default.units = "native",gp=gpar(col = col))
}

########################################
#Plots a gene/transcript module        #
########################################
transcriptModule <- function(vplayout){
    pushViewport(viewport(layout=grid.layout(1,length(vplayout), width=vplayout)))
    for(i in seq(along = vplayout)){
        exonModule(position = i)
    }
}

##############################
#GenericArray                #
##############################
setMethod("drawGD", signature("GenericArray"), function(gdObject, minBase, maxBase, vpPosition) {
  intensity = getIntensity(gdObject)

  ## abstract this operation for all functions :
  ylim <- getPar(gdObject, "ylim")
  xlim <- getPar(gdObject, "xlim")
  if (is.null(xlim)) xlim <- c(minBase, maxBase)
  if (is.null(ylim)) ylim <- range(intensity, na.rm=TRUE)
  
  pushViewport(dataViewport(xData = xlim, yData = intensity, extension = 0,
                            layout.pos.col=1, layout.pos.row = vpPosition, yscale = ylim))
  lwd = getLwd(gdObject)
  lty = getLty(gdObject)
  color = getColor(gdObject)


  whProbes <- getProbeStart(gdObject) >= minBase & getProbeStart(gdObject) < maxBase
  
  if(length(getProbeEnd(gdObject)) > 0) {
    probepos=cbind(getProbeStart(gdObject), getProbeEnd(gdObject))
    for(s in seq(along=intensity[1,])){
      for(p in seq(along=intensity[,1])){
        grid.lines(c(probepos[p,1],probepos[p,2]), c(intensity[p,s],intensity[p,s]),
                   default.units = "native", gp = gpar(col=color[1], lwd = lwd, lty = lty))
      }
    }
  }
  else{
    if(getType(gdObject) == "line"){
      probeStart = getProbeStart(gdObject) 
      ord = order(probeStart)
      probepos = probeStart[ord]
      intensity = intensity[ord,]
      for(p in seq(along=intensity[1,])){
        lwdInd = 1
        colInd = 1
        ltyInd = 1
        if(length(color) == length(intensity[1,])) colInd = p
        if(length(lwd) == length(intensity[1,])) lwdInd = p
        if(length(lty) == length(intensity[1,])) ltyInd = p
        grid.lines(probepos, intensity[,p], default.units = "native",
                   gp = gpar(col=color[colInd], lwd = lwd[lwdInd], lty = lty[ltyInd]))
      }
    }
    else{
      probeStart = getProbeStart(gdObject)
      pSize = getPointSize(gdObject)
      pch = getPch(gdObject)
      for(p in seq(along=intensity[1,])){
        pSizeInd = 1
        colInd = 1
        pchInd = 1
        if(length(color) == ncol(intensity)) colInd = p
        if(length(pSize) == ncol(intensity)) pSizeInd = p
        if(length(pch) == ncol(intensity)) pchInd = p
        
        grid.points(probeStart[whProbes], intensity[whProbes, p], default.units = "native",
                    gp = gpar(col=color[colInd]),size = unit(pSize[pSizeInd],"char"),
                    pch = pch[pchInd])
      }
    }
  }
    
  sObj <- getSegmentation(gdObject)
  if (!is.null(sObj)) {
    .drawSegments(sObj, minBase, maxBase)
  }
  
  grid.yaxis()
  popViewport(1)
})

setMethod("drawGD", signature("Segmentation"), function(gdObject, minBase, maxBase, vpPosition) {
  segments <- getSegments(gdObject)
  
  ylim <- getPar(gdObject, "ylim")
  xlim <- getPar(gdObject, "xlim")
  if (is.null(xlim)) xlim <- c(minBase, maxBase)
  if (is.null(ylim)) ylim <- range(segments, na.rm=TRUE)
  
  pushViewport(dataViewport(xData = xlim, yData = intensity, extension = 0,
                            layout.pos.col=1, layout.pos.row = vpPosition, yscale = ylim))
  .drawSegments(gdObject, minBase, maxBase)
  grid.yaxis(gp=gpar(cex=getPar(gdObject,"cex.axis")))
  popViewport(1)
})

.drawSegments <- function(sObj, minBase, maxBase) {
    segments <- getSegments(sObj)
    segmentStart <- getSegmentStart(sObj)
    segmentEnd <- getSegmentEnd(sObj)

    if(length(segments) > 0) {
        colVal <- getColor(sObj) 
        lwdVal <- getLwd(sObj)   
        ltyVal <- getLty(sObj)   

        ## this just fixes up a real annoyance of having to make things a list. 
        if (length(segments) == 1 && !is.list(colVal))
          colVal <- list(colVal)
        if (length(segments) == 1 && !is.list(ltyVal))
          ltyVal <- list(ltyVal)
        if (length(segments) == 1 && !is.list(lwdVal))
          lwdVal <- list(lwdVal)
        
        ## make the parameters lists of values equal to the number of segments
        lwdVal <- rep(lwdVal,length=length(segments)) #EP
        ltyVal <- rep(ltyVal,length=length(segments)) #EP
        colVal <- rep(colVal,length=length(segments)) #EP
        for(k in seq(along = segments)) {
            whDraw <- !((segmentStart[[k]] > maxBase) |
                        (segmentEnd[[k]] < minBase))
            currLwd<-rep(lwdVal[[k]],length=length(segments[[k]]))
            currLty<-rep(ltyVal[[k]],length=length(segments[[k]]))
            currCol<-rep(colVal[[k]],length=length(segments[[k]]))
            for(s in seq(along=segments[[k]])){
                if (whDraw[s]) {
                    ss <- segmentStart[[k]][s]
                    ee <- segmentEnd[[k]][s]
                    ss <- if (ss < minBase) minBase else ss
                    ee <- if (ee > maxBase) maxBase else ee
                    
                    grid.lines(c(ss,ee), c(segments[[k]][s], segments[[k]][s]), default.units = "native",
                               gp = gpar(col=currCol[s], lwd = currLwd[s], lty =currLty[s] )) #EP
                }
            }
        }
    }
}

setMethod("drawGD", signature("BaseTrack"), function(gdObject, minBase, maxBase, vpPosition) {
    baseValue <- getBaseValue(gdObject)

    ylim <- getPar(gdObject, "ylim")
    xlim <- getPar(gdObject, "xlim")
    if (is.null(xlim)) xlim <- c(minBase, maxBase)
    if (is.null(ylim)) ylim <- range(baseValue, na.rm = TRUE, finite = TRUE)

    drawAxis <- getPar(gdObject, "drawAxis")
    if (is.null(drawAxis)) drawAxis <- TRUE

    if (drawAxis) {
        pushViewport(dataViewport(xData = xlim, yData = ylim, extension = 0,
                                  layout.pos.col = 1, layout.pos.row = vpPosition))
    } else {
        ## XXX: THESE ARE TOTAL HACKS FOR NOW.
        if(is.na(ylim) || diff(ylim) == 0)
            ylim <- c(0,.1)
        if(any(!is.finite(ylim))) {
            ylim <- c(0,.1)
        }
        pushViewport(dataViewport(xData = xlim, yscale = ylim, extension = 0, clip = TRUE,
                                  layout.pos.col = 1, layout.pos.row = vpPosition))
    }

    ## here i probably want to vectorize these two.
    lwd <- getPar(gdObject, "lwd")
    lty <- getPar(gdObject, "lty")
    pty <- getPar(gdObject, "type")
    pos <- getBase(gdObject)
    col <- if (length(color <- getPar(gdObject, "color")) == length(pos)) {
        color
    }
    else {
        rep(color, length(pos))[1:length(pos)]
    }
    whBase <- (pos > minBase & pos < maxBase)
    col <- col[whBase]
    baseValue <- baseValue[whBase]
    pos <- pos[whBase]
    
    if (sum(whBase) > 0) {
        dP <- function() {
            grid.points(pos, baseValue, default.units = "native", gp = gpar(col=col),
                        size = unit(lwd, "char"), pch = 16)
        }
        dVL <- function() {
            mapply(function(a, b, c) {
                grid.lines(x = c(a, a),
                           y = c(min(ylim), b),
                           default.units = "native", gp = gpar(col=c, lwd = unit(lwd, "char")))
            }, pos, baseValue, col)
        }
        
        if (pty == "p") {
            dP()
        } else if (pty == "h") {
            dVL()
        }
        grid.yaxis()
    }

    sObj <- getSegmentation(gdObject)
    if (!is.null(sObj)) {
      .drawSegments(sObj, minBase, maxBase)
    }
    popViewport()
})

##############################################
#
##############################################
setMethod("drawGD", signature("MappedRead"), function(gdObject, minBase, maxBase, vpPosition) {

    pushViewport(dataViewport(xData = c(minBase, maxBase), yData=c(1:100), extension=0,
                              layout.pos.col=1, layout.pos.row = vpPosition))
    col = getColor(gdObject)
    lwd = getLwd(gdObject)
    lty = getLty(gdObject)  
    start = getStart(gdObject)
    end = getEnd(gdObject)
    strand = getStrand(gdObject)
    y = rep(1,length(start)) 
    yl = 1
    for(i in seq(along=start)){
      gp = gpar(col = "blue")
     
     if(i>1){
       if(start[i] <= end[i-1]){
         yl = yl+1
         y[i]=yl
       }
       else{
         yl = 1
       }
     }
    }
    maxy = max(y)
    for(i in seq(along=start)){
      gp = gpar(col = "blue")
      if(strand[i] == "+") gp = gpar(col = "red")
      grid.lines(c(start[i], end[i]), c((y[i]/maxy)*100,(y[i]/maxy)*100), gp = gp, default.units = "native")
    }
  popViewport()
})



##############################################
#Plots exon array data on the probe level    #
##############################################
setMethod("drawGD", signature("ExonArray"), function(gdObject, minBase, maxBase, vpPosition) {
    exonData = getIntensity(gdObject)
    xloc = seq(1:length(exonData[,1]))
    probeLocs = rowMeans(cbind(getProbeStart(gdObject), getProbeEnd(gdObject)))
    vplayout = c("exon"=5,"probemap"=1)
    
    pushViewport(viewport(layout=grid.layout(length(vplayout), 1, height=vplayout),
                          layout.pos.col=1, layout.pos.row = vpPosition))
    pushViewport(dataViewport(xData = xloc, yData=exonData, extension=0, layout.pos.col=1,
                              layout.pos.row = which(names(vplayout)=="exon")))
    color = getColor(gdObject)
    lwd = getLwd(gdObject)
    lty = getLty(gdObject)
    nprobes = getNprobes(gdObject)
  
    for(p in seq(along=exonData[1,])){
        lwdInd = 1
        colInd = 1
        ltyInd = 1
        if(length(color) == length(exonData[1,])) colInd = p
        if(length(lwd) == length(exonData[1,])) lwdInd = p
        if(length(lty) == length(exonData[1,])) ltyInd = p
        grid.lines(xloc, exonData[,p], default.units = "native", gp = gpar(col=color[colInd],
                                                                 lwd = lwd[lwdInd], lty = lty[ltyInd]))
    }
    prloc = 0
    probeSetColor = getProbeSetColor(gdObject)
    probeSetLwd = getProbeSetLwd(gdObject)
    for(i in seq(along=nprobes)){
        prloc = prloc + nprobes[i]
        colInd = 1
        lwdInd = 1 
        if(length(probeSetColor) == length(nprobes)) colInd = i
        if(length(probeSetLwd) == length(nprobes)) lwdInd = i

        if(i < length(nprobes)) {
            grid.lines(c(prloc + 0.5, prloc + 0.5), c(min(exonData),max(exonData)),
                       default.units = "native", gp = gpar(col = probeSetColor[colInd], lwd = probeSetLwd[lwdInd]))
        }
        else{
            grid.lines(c(prloc, prloc), c(min(exonData),max(exonData)),
                       default.units = "native", gp = gpar(col = probeSetColor[colInd], lwd = probeSetLwd[lwdInd]))
        }
        probeLabels = getProbeId(gdObject)
        probeStart = getProbeStart(gdObject)
        
        cex = getCex(gdObject)
        if(displayProbesets(gdObject)) {
          grid.text(label = formatC(probeLabels[i], format="d"), x=prloc - nprobes[i] + 0.5 + nprobes[i]/2 ,
                    y = min(exonData)-1,gp = gpar(cex=cex), rot = 90 ,default.units="native")
        }
    }
    grid.yaxis()
    popViewport()
    if(getPlotMap(gdObject)){
      pushViewport(dataViewport(xData = c(minBase, maxBase), yData=c(0,20), extension=0,
                                layout.pos.col=1, layout.pos.row=which(names(vplayout)=="probemap")))
      connectStep = (maxBase-minBase)/sum(nprobes)
      connectStart = connectStep/2
      mapColor = getMapColor(gdObject)
      for(i in seq(along=probeLocs)){
        colInd = 1
        if(length(mapColor) == length(probeLocs)) colInd = i
        if(i > 1) {
          grid.lines(c(probeLocs[i], (minBase + (sum(nprobes[1:(i-1)]) * connectStep) +
                                      ((nprobes[i]/2)*connectStep))), c(0,15), default.units = "native",
                     gp = gpar(col = mapColor[colInd]))
        }
        else {
          grid.lines(c(probeLocs[i], minBase + ((nprobes[i]/2)*connectStep)), c(0,15),
                     default.units = "native", gp = gpar(col = mapColor[colInd]))
        }
      }
      popViewport()
    }
    popViewport()
})

#######################################
#Add an ideogram to a plot and        #
#possibly highlight region of interest#
#######################################
setMethod("drawGD", signature("Ideogram"), function(gdObject, minBase, maxBase, vpPosition) {
    chromosome = getChromosome(gdObject)
    data("ideogram", package="GenomeGraphs")
    ideo = ideogramTab[ideogramTab[,1] == chromosome,]
    nocol=FALSE
    if(sum(is.na(ideo[,9])) == length(ideo[,9])){nocol=TRUE}
    chromL = max(ideo[,7])
    pushViewport(dataViewport(xData=c(0,max(ideo[,7])), yData=c(0,6), extension=0,
                              layout.pos.col=1, layout.pos.row=vpPosition))
    grid.rect(minBase,1.5,width=maxBase-minBase,height=3,gp=gpar(col = "darkred"),
              default.units="native", just=c("left","bottom"))
    if(!nocol) {
        pal = colorRampPalette(c("white","black"))(100)
        ideo[is.na(ideo[,9]),9] = 1
        for(i in seq(along = ideo[,1])){
            col = pal[ideo[i,9]]   
            grid.rect(x=ideo[i,6],y=2,width=ideo[i,7]-ideo[i,6],height=2,
                      gp=gpar(col = col, fill=col), default.units="native", just=c("left","bottom"))
        }
    }
    else {
        col=rep("white",length(ideo[,2]))
        col[ideo[,8] == "gpos"] = "black"
        for(i in seq(along = ideo[,1])) {   
            grid.rect(x=ideo[i,6],y=2,width=ideo[i,7]-ideo[i,6],height=2,
                      gp=gpar(col = col[i], fill=col[i]), default.units="native", just=c("left","bottom"))
        }
    }
    center = ideo[ideo[,2]=="cen",]
    grid.lines(c(0, center[7]-500000), c(4,4),default.units = "native")
    grid.lines(c(0, center[7]-500000), c(2,2),default.units = "native")
    grid.lines(c(0,0), c(2,4),default.units = "native")
    grid.lines(c(chromL,chromL), c(2,4),default.units = "native")
    grid.lines(c(center[7]+500000,chromL), c(4,4),default.units = "native")
    grid.lines(c(center[7]+500000, chromL), c(2,2),default.units = "native")
    grid.lines(c(center[7]-500000, center[7]+500000), c(4,2),default.units = "native")
    grid.lines(c(center[7]-500000, center[7]+500000), c(2,4),default.units = "native")
    popViewport()
})


## #################################
## #Adding a custom viewport       #
## #################################
## customViewport = function(gdObject){
##   pushViewport(gdObject)
##   popViewport()
## }

#################################
#Adding a title                 #
#################################
setMethod("drawGD", signature("Title"), function(gdObject, minBase, maxBase, vpPosition) {
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpPosition))
    grid.text(getTitle(gdObject), gp = gpar(col = getColor(gdObject), cex = getCex(gdObject)))
    popViewport()
})


####################################################
#Adding an axis containing the genomic coordinates #
####################################################
###################################
## Code from tilingArray package ##
###################################
.ticks <- function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/3, 10)
  fl = floor(lz)
  if( lz-fl > log(5, 10))
    fl = fl +  log(5, 10)
  tw = round(10^fl)
  i0 = ceiling(rx[1]/tw)
  i1 = floor(rx[2]/tw)
  seq(i0, i1)*tw
}



setMethod("drawGD", signature("GenomeAxis"), function(gdObject, minBase, maxBase, vpPosition) {
    pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(-1, 1), extension = 0,
                              layout.pos.col = 1, layout.pos.row = vpPosition))
    grid.lines(c(minBase, maxBase), c(0, 0), default.units = "native")

    ## here we plot the top level ticks
    tck <- .ticks(c(minBase, maxBase))

    ## reformat if byValue != 1
    tckText <- tck
    formatValue <- "d" #to be passed to formatC
    byValue <- getPar(gdObject, "byValue")

    if(byValue != 1) {
        tckText <- tckText/byValue
        formatValue <- "g"
    }

    y <- getPar(gdObject,"distFromAxis")*rep(c(.4, -.4), length(tck))[1:length(tck)]
    labelPos <- getPar(gdObject,"labelPos")
    y <- switch(labelPos,"alternating" = y,"revAlternating" = -y,"above" = abs(y),"below"= -abs(y))
    
    grid.text(label = formatC(tckText, format = formatValue), x = tck, y = y, just = c("centre", "centre"),
              gp = gpar(cex=getPar(gdObject,"cex")), default.units = "native")
    
    if (plotLittleTicks(gdObject)) {
      if (!(minBase %in% tck))
        tck <- c(minBase, tck)
      if (!(maxBase %in% tck))
        tck <- c(tck, maxBase)
      
      if (mean(diff(tck)) > diff(range(tck))/10) {
        for (j in 1:(length(tck) - 1)) {
          stcks <- .ticks(tck[j:(j+1)])
          stcks <- stcks[c(-1, -length(stcks))]
          
          grid.segments(x0 = stcks, x1 = stcks, y0 = -0.1, y1 = 0.1,  default.units = "native")
          
          y <- rep(c(.15, -.15), length(stcks))[1:length(stcks)]
          grid.text(label=formatC(stcks, format="d"), x = stcks, y = y,
                    just = c("centre", "centre"), gp = gpar(cex=.65), default.units = "native")
        }
      }
    }
    
    if(getAdd53(gdObject)){
        grid.text(label=formatC("5'", format="d"), x = minBase, y = 0.3,
                  just = c("centre", "centre"), gp = gpar(cex=.75), default.units = "native")
        grid.text(label=formatC("3'", format="d"), x = maxBase , y = 0.3,
                  just = c("centre", "centre"), gp = gpar(cex=.75), default.units = "native")
    }
    if(getAdd35(gdObject)){
        grid.text(label=formatC("3'", format="d"), x = minBase , y = -0.3,
                  just = c("centre", "centre"), gp = gpar(cex=.75), default.units = "native")
        grid.text(label=formatC("5'", format="d"), x = maxBase , y = -0.3,
                  just = c("centre", "centre"), gp = gpar(cex=.75), default.units = "native")
    }

    grid.segments(x0 = tck, x1 = tck, y0 = -0.25, y1 = 0.25,  default.units = "native")
    popViewport()
})
          
####################
# Adding a legend  #
####################
setMethod("drawGD", signature("Legend"), function(gdObject, minBase, maxBase, vpPosition) {
  pushViewport(dataViewport(xData=c(minBase,maxBase), yscale=c(0,10), extension=0,
                            layout.pos.col=1, layout.pos.row=vpPosition)) 
  color <- getPar(gdObject,"color")
  cex <- getPar(gdObject,"cex")
  legend <- gdObject@legend
  ntext <- length(legend)
  color <- rep(color,ntext)
  cex <- rep(cex,ntext)
  step <- (maxBase - minBase) /(ntext+2) 
  pos <- minBase + step
  
  for(i in 1:ntext) {
    grid.text(label = formatC(legend[i], format="d"), x= pos, y = 4,
                  just = c("center","top"), default.units="native",gp=gpar(cex=cex,lineheight=1)) # add gpar stuff here
    
    if(!is.na(color[i])) {
      grid.rect(pos,5,width=step/3,height=4,
                gp=gpar(col = "black", fill=color[i]), default.units="native",
                just=c("center","bottom"))
    }
    pos <- pos + step
  }
  popViewport()
})

setGeneric("showDisplayOptions", def = function(obj, ...) standardGeneric("showDisplayOptions"))
setMethod("showDisplayOptions", signature("character"), function(obj) {
  getClass(obj)@prototype@dp
})
setMethod("showDisplayOptions", signature("gdObject"), function(obj) {
  obj@dp
})

          





