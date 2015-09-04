intpolat.visual <- function(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter") {
  
  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("Spatial.crossVald", "a1950", target, "ps", sep=".")), 
    color=TRUE, 
    paper=size
  )
    b <- xyplot( mse ~ span | factor(month)*factor(year)
      , data = rst
      , subset = year == 1950
      , xlab = list(label = "Span", cex = 1.5)
      , ylab = list(label = "MSE", cex = 1.5)
      , pch = 16
      , cex = 0.5
      , layout = c(4,3)
      , panel = function(x, y, ...) {
          panel.xyplot(x, y, ...)
        }
    )
    print(b)
  dev.off()
  
}