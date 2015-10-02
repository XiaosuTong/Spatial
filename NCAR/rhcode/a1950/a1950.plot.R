intpolat.visual <- function(size = "letter", surf, SPsize, check=NULL) {
  
  rst1 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", surf, "1", "MSE"))[[1]][[2]]
  rst2 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", surf, "2", "MSE"))[[1]][[2]]
  rst <- rbind(rst1, rst2)
  rst$degree <- rep(c(1,2), each = nrow(rst1))

  if (is.null(check)) {

    if(SPsize == "median") {
      sub <- subset(rst, span %in% c(0.008, 0.01, 0.015, 0.02, 0.04))
    } else if (SPsize == "large") {
      sub <- subset(rst, span %in% seq(0.06, 0.1, 0.01))
    } else if (SPsize == "small"){
      sub <- subset(rst, span %in% seq(0.005, 0.009, 0.001))
    }

    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "span", SPsize, "ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse|as.factor(span)
      , data = sub
      , subset = mse <= 4
      , group = degree
      , dist = qunif
      , cex = 0.5
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
        text = list(label=c("degree=1","degree=2")),
        lines = list(pch=1, cex=1, type="p", col=col[1:2]), 
        columns = 2
        )
      , layout = c(length(unique(sub$span)), 1)
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,4,0.5), v=seq(0,1,0.25), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()
    
    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "degree", SPsize ,"ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse|as.factor(degree)
      , data = sub
      , group = span
      , subset = mse <= 4
      , dist = qunif
      , cex = 0.5
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
        text = list(label=paste("span=", sort(unique(sub$span)), sep="")),
        lines = list(pch=1, cex=1, type="p", col=col[1:length(unique(sub$span))]), 
        columns = length(unique(sub$span))
        )
      , layout = c(2, 1)
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,4,0.5), v=seq(0,1,0.2), col="lightgray")
          panel.qqmath(x,...)
        }  

    )
    print(b)
    dev.off()

  } else {

    sub <- subset(rst, (span == check[[1]][1] & degree == check[[1]][2]) | (span == check[[2]][1] & degree == check[[2]][2]))
    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "span", "check", "ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse
      , data = sub
      , group = degree
      , dist = qunif
      , cex = 0.5
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
        text = list(label=c("degree=1, span=0.008","degree=2, span=0.015")),
        lines = list(pch=1, cex=1, type="p", col=col[1:2]), 
        columns = 2
        )
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,4,0.5), v=seq(0,1,0.2), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

  } 
  
}