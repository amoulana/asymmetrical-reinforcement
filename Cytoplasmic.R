library(parallel)

calc.A <- function(pop){
  return(pop[1]+pop[5]+(pop[2]+pop[3]+pop[6]+pop[7])/2)
}

calc.a <- function(pop){
  return(pop[4]+pop[8]+(pop[2]+pop[3]+pop[6]+pop[7])/2)
}

calc.C <- function(pop){
  return(sum(pop[1:4]))
}

calc.c <- function(pop){
  return(sum(pop[5:8]))
}

random.mating <- function(A,a,C,c){
  local.ad <- c(A^2,A*a,A*a,a^2,A^2,A*a,A*a,a^2)
  return(c(C,C,C,C,c,c,c,c)*local.ad)
}

#CAA,CAa,CaA,Caa,cAA,cAa,caA,caa
#CAABB,CAaBB,CaaBB,... (32)
selection <- function(pop1,pop2,sel1,sel2,sel.comp){
  sel1 <- c(0,0,0,0,sel1,sel1,sel1,sel1)
  sel2 <- c(sel2,sel2,sel2,sel2,00,0,0)
  s1Avg <- sum((1-sel1)*(1-sel.comp)*pop1)
  s2Avg <- sum((1-sel2)*(1-sel.comp)*pop2)
  pop1s <- (1-sel1)*(1-sel.comp)*pop1/s1Avg
  pop2s <- (1-sel2)*(1-sel.comp)*pop2/s2Avg
  return(list(pop1s,pop2s))
}

migration<- function(pop1, pop2,m12,m21){
  sum1m <- sum(pop1*(1-m12) + pop2*m21)
  sum2m <- sum(pop1*m12 + pop2*(1-m21))
  pop1m <- (pop1*(1-m12) + pop2*m21)/sum1m
  pop2m <- (pop1*m12 + pop2*(1-m21))/sum2m
  return(list(pop1m,pop2m))
}

mating<-function(pop1,pop2){
  A1 <- calc.A(pop1)
  a1 <- calc.a(pop1)
  A2 <- calc.A(pop2)
  a2 <- calc.a(pop2)
  C1 <- calc.C(pop1)
  c1 <- calc.c(pop1)
  C2 <- calc.C(pop2)
  c2 <- calc.c(pop2)
  pop1t <- random.mating(A1,a1,C1,c1)
  pop2t <- random.mating(A2,a2,C2,c2)
  return(list(pop1t,pop2t))
}

generation<-function(pop1,pop2,sel1,sel2,sel.comp,m12,m21){
  selected <- selection(pop1,pop2,sel1,sel2,sel.comp)
  mig <- migration(selected[[1]],selected[[2]],m12,m21);
  mat <- mating(mig[[1]],mig[[2]]);
  delta.A1 <- calc.A(mat[[1]])-calc.A(pop1)
  delta.A2 <- calc.A(mat[[2]])-calc.A(pop2)
  A1 <- calc.A(mat[[1]])
  A2 <- calc.A(mat[[2]])
  return(list(mat, delta.A1, delta.A2, A1, A2))
}

equilibrium<-function(pop1,pop2,sel1,sel2,sel.comp,m12,m21){
  pop1x <- pop1
  pop2x <- pop2
  arraydp <- generation(pop1x,pop2x,sel1,sel2,sel.comp,m12,m21)[[3]];
  if(arraydp > 0.00001){
    vec.x <- generation(pop1x,pop2x,sel1,sel2,sel.comp,m12,m21)[[1]]
    equilibrium(vec.x[[1]],vec.x[[2]],sel1,sel2,sel.comp,m12,m21)
  }
  else{
    return(generation(pop1x,pop2x,sel1,sel2,sel.comp,m12,m21)[[4]])
  }
}

selectionAfter<-function(genos1, genos2, sel1, sel2, sel.comp){
  avg1 <- c(sum(genos1[1:8])/sum(genos1),sum(genos1[9:16])/sum(genos1),
            sum(genos1[17:24])/sum(genos1),sum(genos1[25:232])/sum(genos1))
  avg2 <- c(sum(genos2[1:8])/sum(genos2),sum(genos2[9:16])/sum(genos2),
            sum(genos2[17:24])/sum(genos2),sum(genos2[25:232])/sum(genos2))
  genos.BB <- selection(genos1[1:8],genos2[1:8],sel1,sel2,sel.comp)
  genos.Bb <- selection(genos1[9:16],genos2[9:16],sel1,sel2,sel.comp)
  genos.bB <- selection(genos1[17:24],genos2[17:24],sel1,sel2,sel.comp)
  genos.bb <- selection(genos1[25:32],genos2[25:32],sel1,sel2,sel.comp)
  genos1.s <- c(genos.BB[1]/avg1[1],genos.Bb[1]/avg1[2],genos.bB[1]/avg1[2],genos.bb[1]/avg1[3])
  genos2.s <- c(genos.BB[2]/avg2[1],genos.Bb[2]/avg2[2],genos.bB[2]/avg2[2],genos.bb[2]/avg2[3])
  return(list(genos1.s, genos2.s))
}

migrationAfter<-function(genos1, genos2, m12, m21){
  genos1mig <- genos1*(1-m12)+genos2*m21
  genos2mig <- genos2*(1-m21)+genos1*m12
  sum1 = sum(genos1mig)
  sum2 = sum(genos2mig)
  genos1m = genos1mig/sum1
  genos2m = genos2mig/sum2
  return(list(genos1m, genos2m))
}

gametes<-function(geno,r){
  arrayg<-rep(0,4)
  geno<-geno %% 16
  if(geno==1){arrayg<-c(1., 0., 0., 0.)}
  else if(geno==2){arrayg<-c(0.5, 0.5, 0., 0.)} 
  else if(geno==3){arrayg<-c(0.5, 0.5, 0., 0.)} 
  else if(geno==4){arrayg<-c(0., 1., 0., 0.)} 
  else if(geno==5){arrayg<-c(0.5, 0., 0.5, 0.)}
  else if(geno==6){arrayg<-c((1 - r)/2, r/2, r/2, (1 - r)/2)} 
  else if(geno==7){arrayg<-c(r/2, (1 - r)/2, (1 - r)/2, r/2)} 
  else if(geno==8){arrayg<-c(0., 0.5, 0., 0.5)} 
  else if(geno==9){arrayg<-c(0.5, 0., 0.5, 0.)} 
  else if(geno==10){arrayg<-c(r/2, (1 - r)/2, (1 - r)/2, r/2)} 
  else if(geno==11){arrayg<-c((1 - r)/2, r/2, r/2, (1 - r)/2)} 
  else if(geno==12){arrayg<-c(0., 0.5, 0., 0.5)} 
  else if(geno==13){arrayg<-c(0., 0., 1., 0.)} 
  else if(geno==14){arrayg<-c(0., 0., 0.5, 0.5)} 
  else if(geno==15){arrayg<-c(0., 0., 0.5, 0.5)} 
  else if(geno==0){arrayg<-c(0., 0., 0., 1.)}
  return(arrayg)
}

matingResult<-function(genos, genof, genom, adv, r){
  coefficient <- c(1., 1., 1., 1.)
  genof <- genof %% 16
  genom <- genom %% 16
  ifelse(genof==0,genof<-16,genof<-genof)
  ifelse(genom==0,genom<-16,genom<-genof)
  if(genof == 1 | genof == 5 | genof == 9 | genof == 13){ 
    coefficient <- adv;
  }
  if(genom < 5) {maleAallele <- 1} 
  else if(genom < 9){maleAallele <- 2} 
  else if(genom < 13){maleAallele <- 3}
  else{maleAallele <- 4};
  
  freqMating <- genos[genof]*genos[genom]*coefficient[maleAallele];
  gametef <- gametes(genof, r);
  gametem <- gametes(genom, r);
  
  arraym <- freqMating*c(gametef[1]*gametem[1], 
                         gametef[1]*gametem[2], gametef[2]*gametem[1], 
                         gametef[2]*gametem[2], gametef[1]*gametem[3], 
                         gametef[1]*gametem[4], gametef[2]*gametem[3], 
                         gametef[2]*gametem[4], gametef[3]*gametem[1], 
                         gametef[3]*gametem[2], gametef[4]*gametem[1], 
                         gametef[4]*gametem[2], gametef[3]*gametem[3], 
                         gametef[3]*gametem[4], gametef[4]*gametem[3], 
                         gametef[4]*gametem[4]);
  return(arraym)
}


matingAfter<-function(genos1, genos2, d1, d2, r, model, cost){
  #The following calculates the coefficient of mating. 
  #Only when the female genotype has BB is this coefficient used; 
  #otherwise coefficient is {1,1,1,1}. 
  #sumM calculated the denominator. 
  #adv1 and adv2 calculate the coefficients in each population*)
  
  sumM1<-function(genos){ 
    genos <- genos[1:16]+genos[17:32]
    return(
      sum(genos[2:4])*(1 + d1) + 
        sum(genos[6:8])*(1 + d2) +
        sum(genos[10:12]) +
        sum(genos[14:16]) +
        genos[1]*(1 + d1)*(1 - cost) +
        genos[5]*(1 + d2)*(1 - cost) +
        genos[9]*(1 - cost) +
        genos[13]*(1 - cost)
    )
  }
  
  sumM2<-function(genos){ 
    genos <- genos[1:16]+genos[17:32]
    return(
      sum(genos[2:4]) + 
        sum(genos[6:8])*(1 + d2) + 
        sum(genos[10:12]) + 
        sum(genos[14:16])*(1 + d1) + 
        genos[13]*(1 + d1)*(1 - cost) +
        genos[5]*(1 + d2)*(1 - cost) +
        genos[9]*(1 - cost) + 
        genos[1]*(1 - cost)
    )
  }
  
  if(model == 1){
    adv1 <- c((1 + d1)*(1 - cost)/sumM1(genos1), 
              (1 + d2)*(1 - cost)/sumM1(genos1), 
              1*(1 - cost)/sumM1(genos1), 
              1*(1 - cost)/sumM1(genos1))
    adv2 <- c((1 + d1)*(1 - cost)/sumM1(genos2),
              (1 + d2)*(1 - cost)/sumM1(genos2), 
              1*(1 - cost)/sumM1(genos2), 
              1*(1 - cost)/sumM1(genos2))
  }
  else{
    adv1 <- c(1*(1 - cost)/sumM2(genos1), 
              (1 + d2)*(1 - cost)/sumM2(genos1), 
              1*(1 - cost)/sumM2(genos1), 
              (1 + d1)*(1 - cost)/sumM2(genos1))
    adv2 <- c(1*(1 - cost)/sumM2(genos2), 
              (1 + d2)*(1 - cost)/sumM2(genos2), 
              1*(1 - cost)/sumM2(genos2), 
              (1 + d1)*(1 - cost)/sumM2(genos2))
  }
  #There are 32x32 possible matings and we calculate the genotype
  #frequencies resulted from each mating then sum up all results for
  #post-mating genotype frequencies. Somehow, 
  #double sum does not work on Mathematica*)
  
  X.C <- c(1:16)
  X.c <- c(17:32)
  Y <- c(1:32)
  
  newgenos1.C <- rep(0,16)
  newgenos1.c <- rep(0,16)
  newgenos2.C <- rep(0,16)
  newgenos2.c <- rep(0,16)
  
  for(i in X.C){
    for(j in Y){
      newgenos1.C <- newgenos1.C + matingResult(genos1,i,j,adv1,r)
      newgenos2.C <- newgenos2.C + matingResult(genos2,i,j,adv2,r)
    }
  }
  
  for(i in X.c){
    for(j in Y){
      newgenos1.c <- newgenos1.c + matingResult(genos1,i,j,adv1,r)
      newgenos2.c <- newgenos2.c + matingResult(genos2,i,j,adv2,r)
    }
  }
  
  newgenos1 <- c(newgenos1.C,newgenos1.c)
  newgenos2 <- c(newgenos2.C,newgenos2.c)
  return(list(newgenos1, newgenos2))
}

generationAfter<-function(genos1, genos2, m12, m21, sel1, sel2, sel.comp, cost, d1, d2, r, model){
  genos1s <- selectionAfter(genos1, genos2, sel1, sel2, sel.comp)[[1]]
  genos2s <- selectionAfter(genos1, genos2, sel1, sel2, sel.comp)[[2]]
  genos1m <- migrationAfter(genos1s, genos2s, m12, m21)[[1]];
  genos2m <- migrationAfter(genos1s, genos2s, m12, m21)[[2]]
  return(matingAfter(genos1m, genos2m, d1, d2, r, model, cost))
}

B.Allele<- function(genos1x){ 
  return(
    genos1x[1] + genos1x[5] + genos1x[9] + genos1x[13] + genos1x[2]/2 + 
      genos1x[3]/2 + genos1x[6]/2 + genos1x[7]/2 + genos1x[10]/2 + 
      genos1x[11]/2 + genos1x[14]/2 + genos1x[15]/2+
      genos1x[1+16] + genos1x[5+16] + genos1x[9+16] + genos1x[13+16] + genos1x[2+16]/2 + 
      genos1x[3+16]/2 + genos1x[6+16]/2 + genos1x[7+16]/2 + genos1x[10+16]/2 + 
      genos1x[11+16]/2 + genos1x[14+16]/2 + genos1x[15+16]/2
  )
}
