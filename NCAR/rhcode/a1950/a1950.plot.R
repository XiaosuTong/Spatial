intpolat.visual <- function(size = "letter", surf, lim=0) {
  
  rst1 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", surf, "1", "MSE"))[[1]][[2]]
  rst2 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", surf, "2", "MSE"))[[1]][[2]]
  rst <- rbind(rst1, rst2)
  rst$degree <- rep(c(1,2), each = nrow(rst1))

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("Spatial.crossVald", "a1950", par$dataset, surf, "ps", sep=".")), 
    color=TRUE, 
    paper=size
  )
    b <- xyplot( mse ~ as.numeric(span) | factor(month, levels=month.abb)*factor(year)
      , data = arrange(rst, span)
      , subset = as.numeric(span) >= lim
      , group = degree
      , xlab = list(label = "Span", cex = 1.5)
      , ylab = list(label = "MSE", cex = 1.5)
      , scale = list(x= list(cex = 1.2), y=list(cex=1.2, tick.number = 4, relation="free"))
      , key = list(
          text=list(label=c("elev degree=1","elev degree=2")), 
          lines=list(lwd=1.5, type=c("l","l"), col=col[1:2]), 
          columns=2
        )
      , layout = c(4,3)
      , panel = function(x, y, ...) {
          panel.xyplot(x, y, pch=16, type="b", cex=0.5, ...)
        }
    )
    print(b)
  dev.off()

  sub <- subset(rst, span %in% c(0.005, 0.008,0.01,0.015, 0.025))
  
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "span" ,"ps", sep=".")), 
    color=TRUE, 
    paper=size
  )
  b <- qqmath(~mse|as.factor(span)
    , data = sub
    , group = degree
    , subset = mse<=20
    , dist = qunif
    , cex = 0.5
    , auto.key = T
    , layout = c(5, 1)
    , scale = list(y=list(relation="free"))
  )
  print(b)
  dev.off()
  
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "degree" ,"ps", sep=".")), 
    color=TRUE, 
    paper=size
  )
  b <- qqmath(~mse|as.factor(degree)
    , data = sub
    , group = span
    , subset = mse<=20& span!=0.005
    , dist = qunif
    , cex = 0.5
    , auto.key = T
    , layout = c(2, 1)
  )
  print(b)
  dev.off()
  
}