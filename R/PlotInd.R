PlotInd=function(objet, ellipse = TRUE, legend=TRUE){
  if(class(objet)!="PLSDA"){
    stop('First argument must be of class PLSDA')
  }

  ggplot(objet$Classe.score, aes(Comp.1,Comp.2, colour = Classe)) +
    geom_point() +
    if(ellipse==TRUE){
      stat_ellipse(aes(x=Comp.1, y=Comp.2,color=Classe),type = "norm")
    } else if(legend==FALSE){
      theme(legend.position='none')
    }
}
