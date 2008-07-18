#############################
## Genomic Dataset plot    ##
#############################
gdPlot <- function(gdObjects, minBase = NA, maxBase = NA, overlays = NULL,
                   labelColor = "black", labelCex = 1, labelRot = 90) {
    grid.newpage()

    if (class(gdObjects) != "list")
        gdObjects <- list(gdObjects)

    if (missing(minBase) || missing(maxBase)) {
        gr <- sapply(gdObjects, getGenomicRange)
        minBase <- ifelse(all(is.na(gr[1,])), NA, min(gr[1,], na.rm = TRUE))
        maxBase <- ifelse(all(is.na(gr[2,])), NA, max(gr[2,], na.rm = TRUE))
    }
    if (is.na(minBase) || is.na(maxBase))
        stop("Need to define a suitable minBase and maxBase; cannot determine this from the objects.")

    vplayout <- sapply(gdObjects, getSize)
    names(vplayout) <- as.character(seq(along = gdObjects))
    
    ## don't draw things that are of size 0.
    mask <- vplayout > 0
    vplayout <- vplayout[mask]
    gdObjects <- gdObjects[mask]

    ## this is so they follow recycling norms.
    labelColor <- rep(labelColor, length(gdObjects))
    labelCex <- rep(labelCex, length(gdObjects))
    labelRot <- rep(labelRot, length(gdObjects))
    
    if (!is.null(names(gdObjects))) {
        pushViewport(viewport(layout = grid.layout(1, 2, width = c(0.10, 0.9)),
                              width = .90, height = .95))
        pushViewport(viewport(layout=grid.layout(length(vplayout), 1, height=vplayout),
                              layout.pos.col = 1, layout.pos.row = 1))

        for(i in seq(along=gdObjects)) {
            nm <- names(gdObjects)[i]
            if (nm != "") {
                pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
                pushViewport(viewport(0, .5, width = .5, height = .5))
                grid.text(label = formatC(nm, format = "d"),
                          just = "centre", rot = labelRot[i], 
                          default.units = "native", gp = gpar(col = labelColor[i], cex = labelCex[i]))
                popViewport(2)
              }
          }
        popViewport()
        
        pushViewport(viewport(layout = grid.layout(length(vplayout), 1, height = vplayout),
                              layout.pos.col = 2, layout.pos.row = 1))
    }
    else{
        ## here you make the region smaller so that the axis fit. 
        pushViewport(viewport(layout = grid.layout(length(vplayout), 1, height = vplayout),
                              width = .85, height = .95))
    }
    
    for (i in seq(along = gdObjects)) {
        drawGD(gdObjects[[i]], minBase, maxBase, i)
    }
    
    if (!is.null(overlays)) {
        if (!is.list(overlays))
            overlays <- list(overlays)

        lapply(overlays, function(overlay) {
            drawOverlay(overlay, minBase, maxBase, vplayout)
        })
    }
    ##-- pop to the top of the viewport stack.
    popViewport(0)
}

