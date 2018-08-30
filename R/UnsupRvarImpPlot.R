#' Unsupervised random forest variable importance 
#
#' @description Extract variable importance from \code{\link{RFdist}} if importance = TRUE  and plot for both the 
#' "Synthetic and "True" data   
#
#' @param x object of class \code{\link{RFdist}} 
#' @param sort, sort variable importance by variables in the orginal data. 
#' @param top  number of top variables to plot. 
#' @param plot (logical) plot or simply return variable importance as a matrix ? 
#' @return matrix if plot = FALSE

#' @import ggplot2 
#' @export
### plot variable importance 
UnsupRFvarImpPlot <- function(x, sort=TRUE, top=min(30, nrow(x$UnsupRFvarimp)), 
                  plot = TRUE) {
#suppressMessages(requireNamespace(package="ggplot2", quietly = TRUE))
  
imp <- as.data.frame(x$UnsupRFvarimp)[, c(1,2)] 
Variable = rownames(imp)
Importance = imp[, "True.Data"] 
x1 <-  data.frame(Variable, Importance)
x1$class <- "True.Data"
Variable = rownames(imp)
Importance = imp[, "Synthetic.Data"]
x2 <-  data.frame(Variable, Importance) 
x2$class <- "Synthetic.Data"

if(sort){
ix <- order(x1$Importance, decreasing = TRUE)
x1 <- x1[ix, ][1:top, ]
x2 <- x2[ix, ][1:top, ]
} else {
x1 <- x1[1:top, ]
x2 <- x2[1:top, ]
}
tab <- rbind(x1, x2) 
if(plot){
tab$Variable <- factor(tab$Variable, levels = unique(tab$Variable))

rg <- range(tab$Importance)
  ticks <- round(seq(rg[1]-0.02, rg[2]+ 0.02, 0.03), 3)
  pp <- ggplot( ) +
    geom_point(data = tab, aes(x = Variable, y = Importance, colour = class),
     position = position_dodge(width = 0.02), size = 3) + 
    geom_hline(yintercept = 0) +
    xlab("") + ylab("")  + ggtitle("Variable Importance") + labs(colour = "Type") + scale_y_continuous(breaks = ticks, labels=function(x)x*100) + 
    coord_flip() +  scale_colour_manual(values = c("darkgreen",  "darkred")) + theme_bw() + 
    theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), axis.title.x=element_text(size=14,face="bold"), 
          axis.title.y=element_text(size=14,face="bold"), 
          axis.text.x = element_text(size = 15, face="bold",colour = "gray10"),
          legend.text = element_text(size=14,face="bold"), 
          legend.title = element_text(size=14,face="bold"),
          plot.title = element_text(size=15, face="bold.italic", hjust = 0.5),
          axis.text.y = element_text(size = 15, face="bold",colour = "gray10"))  
  
 print(pp)
 } else return(tab) 
    }





