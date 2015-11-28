a1950Status <- function(sample=100) {

  rst <- rhread("/ln/tongx/Spatial/tmp/tmax/a1950/bymonth", max=sample)
  us.map <- map('state', plot = FALSE, fill = TRUE)

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste(par$dataset, "a1950", "status", "ps", sep=".")), 
    color = TRUE, 
    paper = "letter"
  )
  for(i in 1:sample) {
    b <- xyplot( lat ~ lon
      , data = rst[[i]][[2]]
      , groups = factor(as.numeric(!is.na(resp)))
      , xlab = list(label="Longitude", cex = 1.5)
      , ylab = list(label="Latitude", cex = 1.5)
      , sub = paste(rst[[i]][[1]][2], rst[[i]][[1]][1])
      , key = list(
          type = "p", 
          text = list(label=c("missing","valid")),  
          points = list(cex=1, pch=16, col=col[1:2]), 
          columns = 2
        )
      , pch = 16
      , cex = 0.4
      , scales = list(
          y = list(cex = 1.2),x = list(cex = 1.2)
        )
      , panel = function(x,y,...) {
          panel.polygon(us.map$x,us.map$y)   
          panel.xyplot(x,y,...)
          panel.text(x=-123, y=26, labels=paste("Missing:", sum(is.na(rst[[i]][[2]]$resp))), adj=c(0,0))
          panel.text(x=-123, y=25.5, labels=paste("Valid:    ", sum(!is.na(rst[[i]][[2]]$resp))), adj=c(0,1))
        }
    )
    print(b)
  }
  dev.off()

}