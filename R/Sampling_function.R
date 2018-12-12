
sampling <- function(x, pourcentage=0.7){
  # Controle de la bonne entr?e des valeurs
  if (!class(x)=="data.frame" & !class(x)=="matrix"){
    stop("Classe des donn?es innapropri?e")
  }
  if (is.numeric(pourcentage)==FALSE | pourcentage>1 | pourcentage<0){
    stop("Pourcentage incorrecte")
  }
  smp = sample(1:nrow(x), nrow(x)*pourcentage)
  instance = list()
  instance$train = x[smp,]
  instance$test = x[-smp,]
  return(instance)
}
data=sampling(x)
