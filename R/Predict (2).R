install.packages("nnet")
library(nnet)

predict=function(objet,newdata,type="class"){

  if(class(objet)!="PLSDA"){
    stop('First argument must be of class PLSDA')
  }

  SplitX <- function(x){
    l_fact <- c()

    for(i in 1:ncol(x)){
      if(is.factor(x[,i])==TRUE){
        l_fact <- cbind(l_fact, i)
      }
    }

    if(length(l_fact)>1 | is.null(l_fact)==T){
      print("Votre Data Frame ne doit poss√©der qu'une seule variable de type factorielle !")
    } else {
      instance=list()
      instance$y <- x[,l_fact]
      instance$x <- x[,-l_fact]
      return(instance)
    }
  }
  dtest = SplitX(newdata)

  #Centrer et rÈduire par la moyenne et la variance de la base d'entrainement
  dscale=t(apply(dtest$x,1,function(x){(x-objet$mean.var)/objet$sd.var}))

  score_logit <- as.matrix(dscale) %*% t(objet$fonction_logit[,-ncol(objet$fonction_logit)])

  score_logit <- as.matrix(apply(score_logit,1,function(x){
    x + objet$fonction_logit[,5]
  }))

  #SoftMax
  proba <-t(apply(t(score_logit),1,function(x){exp(x)/sum(exp(x))}))
  #Pr?diction de la classe
  Classe<-sapply(apply(proba,1,which.max),function(x){
    x=colnames(proba)[x]})

  MC <- table(dtest$y,Classe)

  inst=list()
  inst$predict # coefficients rattach?s ? chaque variable par axes factorielles
  inst$variates # coefficients dans les deux axes factorielles
  inst$B.hat # table avec les coefficients entre les variables et les classes ? pr?dire
  inst$classe=Classe # classes pr?dites en fonction de plusieurs crit?res (mahalanobis.dist,centroids.dist,max.dist)
  inst$proba=proba
  inst$MatrixConfusion=MC
  instance$x.score=objet$x.score
  instance$Classe.x.score=data.frame(Classe,x.score)
  class(inst)="PLSDA"
  return(inst)
}
"""
  if(type=="class"){
    return(instance[-5])
  } else if(type=="posterior"){
    return(instance[-4])
  }
"""

