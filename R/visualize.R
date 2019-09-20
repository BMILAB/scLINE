#' Visualization
#' @usage visualize(data, label)
#' @param data a matrix in which each row represents a cell
#' @param label true labels of the cell corresponding to  data
#'
#' @return Visualization image of the data
#' @export
#' @description Visualize the low-dimensional vectors obtained by scLINE in 2D
#' @examples
#' visualize(low_mat$cell_low,Usoskin$label)
visualize<-function(data,label){
  pca<-prcomp(data)
  temp1<-as.data.frame(predict(pca))
  score<-data.frame(PC1=temp1$PC1,PC2=temp1$PC2,lab=label)
  ggplot(score,aes(x=PC1,y=PC2,color=lab))+
    scale_color_brewer(palette = "Set3")+
    theme_bw(base_size = 12)+
    theme(
      # panel.border = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          # axis.line = element_line(colour = "black",size = 1.5),
          # axis.title.x = element_text(size = 24,color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
          # axis.title.y = element_text(size = 24,color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
          # axis.text = element_text(size = 24,face = "bold",colour = 'black'),
          plot.title = element_text(hjust = 0.5,face = "bold"),
          # legend.position = c(0.7,0.8),
          legend.text =element_text(size = 12,face="bold") ,
          legend.title=element_blank())+
    geom_point(size=4)+
    xlab("PC1")+ylab("PC2")+
    labs(title="Visualization",hjust=0.3)
}
