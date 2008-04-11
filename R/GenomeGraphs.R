#############################
## Genomic Dataset plot    ##
#############################
gdPlot <- function(gdObjects, minBase = NA, maxBase = NA, highlightRegions = NULL) {
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
        
        ## should there be an argument for 'text rotation?'
        grid.text(label = formatC(nm, format = "d"),
                  just = "centre", rot = 90, 
                  default.units = "native")
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
  
  if (!is.null(highlightRegions)) {
    if (!is.list(highlightRegions))
      highlightRegions <- list(highlightRegions)

    lapply(highlightRegions, function(highlightRegion) {
      .drawHighlightRegion(highlightRegion, minBase, maxBase, vplayout)
    })
  }
  ##-- pop to the top of the viewport stack.
  popViewport(0)
}

.drawHighlightRegion <- function(hr, minBase, maxBase, vplayout) {
  switch(hr@coords,
         "genomic" = {
           ss <- (hr@start - minBase)/(maxBase - minBase) 
           ee <- (hr@end - minBase)/(maxBase - minBase) 
           if (is.null(hr@region)) {
             y0 <- 0
             height <- 1
           }
           else {
             region <- hr@region
             y0 <- 1 - sum(vplayout[1:region[2]])/sum(vplayout)
             height <- sum(vplayout[region[1]:region[2]])/sum(vplayout)
           }
         },
         "absolute" = {
           ss <- hr@start
           ee <- hr@end
           y0 <- hr@region[1]
           height <- hr@region[2] - hr@region[1]
         },
         stop("Unknown coordinate specification in HighlightRegion."))
  
  color <- rgb(t(col2rgb(getColor(hr))/255), alpha = getPar(hr, "alpha"))
  grid.rect(ss, y0, width = (ee - ss), height = height, gp = gpar(fill=color, lwd = getLwd(hr), lty = getLty(hr)),
            just=c("left", "bottom"))
}

###################################
## Code from tilingArray package ##
###################################
ticks <- function(x){
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
