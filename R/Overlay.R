setClassUnion("numericOrNull", c("numeric", "NULL"))

setClass("Overlay", contains = "gdObject")

setClass("RectangleOverlay", contains = "Overlay",
         representation(start = "numeric",
                        end = "numeric",
                        region = "numericOrNull",
                        coords = "character"),
         prototype(region = NULL,
                   coords = "genomic",
                   dp = DisplayPars(color = "black",
                                    fill = "black",
                                    alpha = 0, lwd = 1, lty = "solid"))
         );

setClass("TextOverlay", contains = "Overlay",
         representation(text = "character",
                        xpos = "numeric",
                        ypos = "numeric",
                        region = "numericOrNull",
                        coords="character"),
         prototype(region = NULL,
                   coords = "genomic", 
                   dp = DisplayPars(color = "black", cex=1,just=c(0.5,0.5)))
         );


setGeneric("drawOverlay", def = function(obj, ...) standardGeneric("drawOverlay"))
setMethod("drawOverlay", signature("RectangleOverlay"), function(obj, minBase, maxBase, vplayout) {
  switch(obj@coords,
         "genomic" = {
           ss <- (obj@start - minBase)/(maxBase - minBase) 
           ee <- (obj@end - minBase)/(maxBase - minBase) 
           if (is.null(obj@region)) {
             y0 <- 0
             height <- 1
           }
           else {
             region <- obj@region
             y0 <- 1 - sum(vplayout[1:region[2]])/sum(vplayout)
             height <- sum(vplayout[region[1]:region[2]])/sum(vplayout)
           }
         },
           "absolute" = {
             ss <- obj@start
             ee <- obj@end
             y0 <- obj@region[1]
             height <- obj@region[2] - obj@region[1]
           },
         stop("Unknown coordinate specification in HighlightRegion."))
  
  fill <- rgb(t(col2rgb(getPar(obj, "fill"))/255), alpha = getPar(obj, "alpha"))
  color <- getPar(obj, "color")
  
  grid.rect(ss, y0, width = (ee - ss), height = height, gp = gpar(fill=fill, col = color, lwd = getLwd(obj),
                                                          lty = getLty(obj)),
            just=c("left", "bottom"))
})

setMethod("drawOverlay", signature("TextOverlay"), function(obj, minBase, maxBase, vplayout) {
            switch(obj@coords,
                   "genomic" = {
                     xab <- (obj@xpos - minBase)/(maxBase - minBase) 
                   },
                   "absolute"={
                     xab<-obj@xpos 
                   })
            region <- obj@region
            if(!is.null(region)){ #when user-given y is absolute position relative to within a plot
              ymin <- 1 - sum(vplayout[1:region[2]])/sum(vplayout) #absolute coordintes of bottom of specified region
              ymax <- ymin+sum(vplayout[region[1]:region[2]])/sum(vplayout)#absolute coordinates of the top of specified region
              yab<-ymin+obj@ypos*(ymax-ymin)
            }
            else
              yab<-obj@ypos

            grid.text(obj@text,x=xab, y=yab, just=getPar(obj,"just"),
                      gp = gpar(col=getColor(obj), cex=getPar(obj,"cex"), lineheight=getPar(obj,"lineheight")))
          })


##
## Public API
##
makeRectangleOverlay <- function(start, end, region = NULL, coords = c("genomic", "absolute"), 
                             dp = NULL) {
  coords <- match.arg(coords)
  if (is.null(dp))
    dp <- getClass("RectangleOverlay")@prototype@dp
  new("RectangleOverlay", start = start, end = end, region = region, coords = coords,
      dp = dp)
}

makeTextOverlay <- function(text, xpos, ypos, region = NULL, coords = c("genomic", "absolute"),
                        dp = NULL) {
  coords <- match.arg(coords)
  if (is.null(dp))
    dp <- getClass("TextOverlay")@prototype@dp
  new("TextOverlay", text = text, xpos = xpos, ypos = ypos, region = region, coords = coords,
      dp = dp)
}
