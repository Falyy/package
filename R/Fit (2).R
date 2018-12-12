install.packages("dummies", dependencies = T)
install.packages("plsdepot", dependencies = T)
library("dummies")
library("plsdepot")
library("nnet")
library("ggplot2")


fit=function(formula,data,ncomp=2, select=F) {
  if (!class(formula)=="formula"){
    stop('First argument must be of class formula')
  }
  dtrain <- model.frame(formula=formula, data=data)

SplitX <- function(x){
  l_fact <- c()

  for(i in 1:ncol(x)){
    if(is.factor(x[,i])==TRUE){
      l_fact <- cbind(l_fact, i)
    }
  }

  if(length(l_fact)>1 | is.null(l_fact)==T){
    print("Votre Data Frame ne doit posséder qu'une seule variable de type factorielle !")
  } else {
    instance=list()
    instance$y <- x[,l_fact]
    instance$x <- x[,-l_fact]
    return(instance)
  }
}
  dtrain = SplitX(dtrain)

  if (any(is.na(dtrain$y))){
    stop("\n NA trouvés, traitement impossible !")
  }

  if(select == T){

    #Initialisation
    passage = 1 #Décris le nombre de passage dans la boucle repeat
    repeat{

      if(passage > 1){
        x_ini = as.matrix(dtrain$x[,-test_vip])
        x= as.matrix(scale(dtrain$x[,-test_vip]))
      }else{
        x_ini = as.matrix(dtrain$x)
        x= as.matrix(scale(dtrain$x))
      }

      y <- as.matrix(dtrain$y)
      y_dum = scale(dummies::dummy(y))
      y_dum_ini = y_dum
      p=ncol(x)
      n = nrow(x)
      k = ncol(y_dum)
      rg = qr(x)$rank

      #Initialisation matrices
      w= matrix(0,p,rg)
      t= matrix(0,n,rg)
      u= matrix(0,n,rg)
      q= matrix(0,k,rg)
      P= matrix(0,p,rg)
      C= matrix(0,1,rg)
      # X= list()
      # Y= list()

      #Matrices à rendre en sortie : Initialisation
      x.score = matrix(0, n, rg)
      y.score = matrix(0, n, rg)
      x.loads = matrix(0, p, rg)
      y.loads = matrix(0, k, rg)

      for(r in 1 : rg){
        #Préparation données initiales
        u[,r] = y_dum[,1]
        w_ancien = rep(1, p)
        i = 1

        w_ancien = t(x)%*% u[,r]/sum(u[,r]^2)
        w_ancien = w_ancien/sqrt(sum(w_ancien^2)) #Normalisation

        #Partie itérative jusqu'a convergence de w
        while( w[,r] != w_ancien) {

          if(i != 1){ # hors premier passage
            w_ancien = w[,r]# on recupere l'ancienne valeur de w
          }

          t[,r] = x %*% w_ancien
          q[,r] = t(y_dum) %*% t[,r]/sum(t[,r]^2)
          u[,r] = y_dum %*% q[,r]/sum(q[,r]^2)

          w[,r] = t(x)%*% u[,r]/sum(u[,r]^2)
          w[,r] = w[,r]/sqrt(sum(w[,r]^2))

          # c.new = t(Y.old) %*% t.new/sum(t.new^2)
          # u.new = Y.old %*% c.new/sum(c.new^2)

          i = i+1
          if(i > 200){
            print("Convergence de l'algorithme")
            break
          }

        }

        P[,r] = t(x) %*% t[,r]/sum(t[,r]^2)
        x = x - t[,r] %*% t(P[,r])
        C[,r] = (t(t[,r]) %*% u[,r]) / (t(t[,r]) %*% t[,r])
        y_dum = y_dum -  t[,r] %*% t(q[,r])

        x.score[, r] = t[,r]
        x.loads[, r] = P[,r]
        y.score[, r] = u[,r]
        y.loads[, r] = q[,r]

      }
      passage = passage +1

      #correlation entre les x/y et les scores de x
      cor.xt = cor(x_ini, x.score[,1:ncomp])^2
      cor.yt = cor(y_dum_ini, x.score[,1:ncomp])^2
      # Corrélations moyennes par colonnes
      rx_mean = apply(cor.xt,2,mean)
      ry_mean = apply(cor.yt,2,mean)

      EV = cbind(rx_mean, cumsum(rx_mean), ry_mean, cumsum(ry_mean))
      #Création matrice carrée: de taille rang de X
      m_correlation = matrix(0, ncomp, ncomp)
      for (m in 1:ncomp){
        m_correlation[1:ncomp,m] = c(ry_mean[1:m], rep(0, ncomp-m))
      }
      VIP = sqrt((w[,1:ncomp]^2) %*% m_correlation %*% diag(ncol(x_ini)/sum(ry_mean), ncomp, ncomp))
      rownames(VIP)=colnames(x_ini)

      barplot(sort(VIP[,1], decreasing =T), xlab=rownames(VIP), ylab="VIP", col = c(3), main = c("Importance de chaque variable graphe : ", passage-1), ylim= c(0,1.5))
      abline(h=0.8, col="red",lwd=3, lty=2)

      test_vip = which(VIP[,1] < 0.8)

      #Si l'objet est vide on sort : signifie qu'il n'y a plus de variables à enlever
      if(length(test_vip) == 0){
        break()
      }

    }
  }else{ #Selection de variable à F
    #Initialisation

    x_ini = as.matrix(dtrain$x)
    x= as.matrix(scale(dtrain$x))

    y <- as.matrix(dtrain$y)
    y_dum = scale(dummies::dummy(y))
    y_dum_ini = y_dum
    p=ncol(x)
    n = nrow(x)
    k = ncol(y_dum)
    rg = qr(x)$rank

    #Initialisation matrices
    w= matrix(0,p,rg)
    t= matrix(0,n,rg)
    u= matrix(0,n,rg)
    q= matrix(0,k,rg)
    P= matrix(0,p,rg)
    C= matrix(0,1,rg)

    #Matrices à rendre en sortie : Initialisation
    x.score = matrix(0, n, rg)
    y.score = matrix(0, n, rg)
    x.loads = matrix(0, p, rg)
    y.loads = matrix(0, k, rg)

    for(r in 1 : rg){
      #Préparation données initiales premier passage boucle while
      u[,r] = y_dum[,1]
      w_ancien = rep(1, p)
      i = 1

      w_ancien = t(x)%*% u[,r]/sum(u[,r]^2)
      w_ancien = w_ancien/sqrt(sum(w_ancien^2)) #Normalisation

      #Partie itérative jusqu'a convergence de w
      while( w[,r] != w_ancien) {

        if(i != 1){ # hors premier passage
          w_ancien = w[,r]# on recupere l'ancienne valeur de w
        }

        t[,r] = x %*% w_ancien
        q[,r] = t(y_dum) %*% t[,r]/sum(t[,r]^2)
        u[,r] = y_dum %*% q[,r]/sum(q[,r]^2)

        w[,r] = t(x)%*% u[,r]/sum(u[,r]^2)
        w[,r] = w[,r]/sqrt(sum(w[,r]^2))

        # c.new = t(Y.old) %*% t.new/sum(t.new^2)
        # u.new = Y.old %*% c.new/sum(c.new^2)

        i = i+1
        if(i > 200){
          print("Convergence de l'algorithme")
          break}

      }

      P[,r] = t(x) %*% t[,r]/sum(t[,r]^2)
      x = x - t[,r] %*% t(P[,r])
      C[,r] = (t(t[,r]) %*% u[,r]) / (t(t[,r]) %*% t[,r])
      y_dum = y_dum -  t[,r] %*% t(q[,r])

      x.score[, r] = t[,r]
      x.loads[, r] = P[,r]
      y.score[, r] = u[,r]
      y.loads[, r] = q[,r]

    }

    #correlation entre les x/y et les scores de x
    cor.xt = cor(x_ini, x.score[,1:ncomp])^2
    cor.yt = cor(y_dum_ini, x.score[,1:ncomp])^2
    # Corrélations moyennes par colonnes
    rx_mean = apply(cor.xt,2,mean)
    ry_mean = apply(cor.yt,2,mean)

    EV = cbind(rx_mean, cumsum(rx_mean), ry_mean, cumsum(ry_mean))
    #Création matrice carrée: de taille rang de X
    m_correlation = matrix(0, ncomp, ncomp)
    for (m in 1:ncomp){
      m_correlation[1:ncomp,m] = c(ry_mean[1:m], rep(0, ncomp-m))
    }

    VIP = sqrt((w[,1:ncomp]^2) %*% m_correlation %*% diag(ncol(x_ini)/sum(ry_mean), ncomp, ncomp))
    rownames(VIP)=colnames(x_ini)

    barplot(sort(VIP[,1], decreasing =T), xlab=rownames(VIP), ylab="VIP", col = c(3), main = "Importance de chaque variable graphe" , ylim= c(0,1.5))
    abline(h=0.8, col="red",lwd=3, lty=2)


  }

  x.score=x.score[,1:ncomp]
  y.score=y.score[,1:ncomp]
  x.loads=x.loads[,1:ncomp]
  y.loads=y.loads[,1:ncomp]

  rr=rep("Comp-", ncomp)
  ss=seq(1,ncomp)
  compnames=paste(rr,ss, sep="")

  rownames(x.loads)=colnames(dtrain$x)
  colnames(x.loads)=compnames
  colnames(x.score)=compnames
  colnames(y.score)=compnames
  colnames(y.score)=compnames
  colnames(y.loads)=compnames
  rownames(y.loads)=colnames(y_dum)

  # Predict sur les donn�es d'entrainement
  reg <- multinom(formula = dtrain$y ~ x.score) # R?gression logisitique

  #score=apply(score,2,function(x){(x-mean(x))/sd(x)})

  coef <- as.matrix(t(coef(reg))) # coefficients sur les axes factorielles
  intercept = as.matrix(coef[1,], ncol=1) # extraction du coefficient fixe
  colnames(intercept) <-"intercept"  #
  coef_fonction = t(as.matrix(x.loads %*% coef[-1,])) # projection des coefficients sur les 4 axes originaux

  # Ajout de la troisi?me classe dans le tableau var/classe
  classe <- as.matrix(t(c(rep(0,ncol(coef_fonction)),1)))
  fonction_logit <- rbind(cbind(coef_fonction, intercept),classe) # objet fonction logit doit �tre transmis
  row.names(fonction_logit)[nrow(fonction_logit)]=levels(dtrain$y)[1]

  #V?rifier (train scale ?)
  score_logit <- as.matrix(scale(dtrain$x)) %*% t(fonction_logit[,-ncol(fonction_logit)])

  #score_logit <- as.matrix(scale(data$test)) %*% t(fonction_logit[,-ncol(fonction_logit)])

  score_logit <- as.matrix(apply(score_logit,1,function(x){
    x + fonction_logit[,5]
  }))

  #SoftMax
  proba <-t(apply(t(score_logit),1,function(x){exp(x)/sum(exp(x))}))
  #Pr?diction de la classe
  Classe<-sapply(apply(proba,1,which.max),function(x){
    x=colnames(proba)[x]})

  #
  MC <- table(dtrain$y,Classe)

  mean.var=apply(dtrain$x,2,mean)
  sd.var=apply(dtrain$x,2,sd)

  instance=list()
  instance$predict # coefficients rattach?s ? chaque variable par axes factorielles
  instance$variates # coefficients dans les deux axes factorielles
  instance$B.hat # table avec les coefficients entre les variables et les classes ? pr?dire
  instance$Classe=Classe # classes pr?dites en fonction de plusieurs crit?res (mahalanobis.dist,centroids.dist,max.dist)
  instance$proba=proba
  instance$x.score=x.score
  instance$y.score=y.score
  instance$x.loads=x.loads
  instance$y.loads=y.loads
  instance$MatrixConfusion=MC
  instance$mean.var=mean.var
  instance$sd.var=sd.var
  instance$fonction_logit=fonction_logit
  instance$Classe.score=data.frame(Classe,x.score)
  class(instance)="PLSDA"

  return(instance)
}

#   colnames(classe_pred)=c("Comp1","Comp2","Classe")


 """
  y_indic <- as.data.frame(dummy(d$y))

  # les rajouter en sortie
  mean.var=apply(d$x,2,mean)
  sd.var=apply(d$x,2,sd)

  d$x <- as.data.frame(scale(as.data.frame(d$x), center = T, scale = T))

  instance=list()
  class(instance)="PLSDA"
  instance$x.scores=p$x.scores
  instance$x.loads=p$x.loads
  instance$y.scores=p$y.scores
  instance$y.loads=p$y.loads
  instance$expvar=p$expvar
  instance$VIP=p$VIP
  return(instance)
}
