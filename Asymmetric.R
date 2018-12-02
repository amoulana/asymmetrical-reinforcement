library(parallel)
selection <- function(AA1, Aa1, aA1, aa1, AA2, Aa2, aA2, aa2, sAa, saA,sAA, saa){
  s1Avg <- (1 - sAA)*AA1 + (1 - saa)*aa1 + (1 - sAa)*Aa1 + (1 - saA)*aA1 
  s2Avg <- (1 - sAA)*AA2 + (1 - saa)*aa2 + (1 - sAa)*Aa2 + (1 - saA)*aA2
  AA1s <- ((1 - sAA)*AA1)/s1Avg
  Aa1s <- (Aa1*(1 - sAa))/s1Avg
  aA1s <- (aA1*(1 - saA))/s1Avg
  aa1s <- ((1 - saa)*aa1)/s1Avg
  AA2s <- ((1 - sAA)*AA2)/s2Avg
  Aa2s <- (Aa2*(1 - sAa))/s2Avg
  aA2s <- (aA2*(1 - saA))/s2Avg
  aa2s <- ((1 - saa)*aa2)/s1Avg
  return(c(AA1s, Aa1s, aA1s, aa1s, AA2s, Aa2s, aA2s, aa2s))
}

migration<- function(AA1, Aa1, aA1, aa1, AA2, Aa2, aA2, aa2, m12,m21){
  AA1m <- AA1*(1 - m12) + AA2*m21 
  Aa1m <- Aa1*(1 - m12) + Aa2*m21 
  aA1m <- aA1*(1 - m12) + aA2*m21 
  aa1m <- aa1*(1 - m12) + aa2*m21 
  AA2m <- AA2*(1 - m21) + AA1*m12
  Aa2m <- Aa2*(1 - m21) + Aa1*m12
  aA2m <- aA2*(1 - m21) + aA1*m12
  aa2m <- aa2*(1 - m21) + aa1*m12
  sum1m <- AA1m + aA1m + Aa1m + aa1m
  sum2m <- AA2m + aA2m + Aa2m + aa2m
  AA1m <- AA1m/sum1m
  Aa1m <- Aa1m/sum1m
  aA1m <- aA1m/sum1m
  aa1m <- aa1m/sum1m;
  AA2m <- AA2m/sum2m; Aa2m <- Aa2m/sum2m; aA2m <- aA2m/sum2m; 
  aa2m <- aa2m/sum2m;
  return(c(AA1m, Aa1m, aA1m, aa1m, AA2m, Aa2m, aA2m, aa2m))
}

mating<-function(AA1, Aa1, aA1, aa1, AA2, Aa2, aA2, aa2){
  p1 <- AA1 + Aa1/2 + aA1/2; q1 <- aa1 + Aa1/2 + aA1/2; 
  p2 <- AA2 + Aa2/2 + aA2/2; q2 <- aa2 + Aa2/2 + aA2/2; 
  AA1t <- p1^2; 
  Aa1t <- p1*q1;
  aA1t <- p1*q1; aa1t <- q1^2; AA2t <- p2^2; Aa2t <- p2*q2;
  aA2t <- p2*q2; aa2t <- q2^2;
  return(c(AA1t, Aa1t, aA1t, aa1t, AA2t, Aa2t, aA2t, aa2t))
}

generation<-function(AA1, Aa1, aA1, aa1, AA2, Aa2, aA2, aa2,sAa,saA,m12,m21,sAA,saa){
   sel <- selection(AA1, Aa1, aA1, aa1, AA2, Aa2, aA2, aa2, sAa, saA, sAA, saa)
   mig <- migration(sel[1], sel[2], sel[3], sel[4], sel[5],sel[6], sel[7], sel[8], m12, m21);
   mat <- mating(mig[1], mig[2], mig[3], mig[4], mig[5],mig[6], mig[7], mig[8]);
   deltap1 <- mat[1] + ((mat[2] + mat[3])/2) - AA1 - ((Aa1 + aA1)/2);
   deltap2 <- mat[5] + ((mat[6] + mat[7])/2) - AA2 - ((Aa2 + aA2)/2);
   pone <- mat[1] + ((mat[2] + mat[3])/2);
   ptwo <- mat[5] + ((mat[6] + mat[7])/2);
   return(list(mat, deltap1, deltap2, pone, ptwo))
}

equilibrium<-function(AA1, Aa1, aA1, aa1, AA2, Aa2, aA2, aa2,sAa,saA,m12,m21,sAA,saa){
  AA1x <- AA1;
  Aa1x <- Aa1;
  aA1x <- aA1;
  aa1x <- aa1;
  AA2x <- AA2;
  Aa2x <- Aa2;
  aA2x <- aA2;
  aa2x <- aa2;
  arraydp <- generation(AA1x, Aa1x, aA1x, aa1x, AA2x, Aa2x, aA2x, aa2x, sAa, 
               saA, m12, m21, sAA, saa)[[3]];
  if(arraydp > 0.00001){
    vec.x <- generation(AA1x, Aa1x, aA1x, aa1x, AA2x, Aa2x, aA2x, aa2x, sAa, saA, m12, m21, sAA, saa)[[1]]
    equilibrium(vec.x[1],vec.x[2],vec.x[3],vec.x[4],vec.x[5],vec.x[6],vec.x[7],vec.x[8], sAa, 
                    saA, m12, m21, sAA, saa)
  }
  else{
    return(generation(AA1x, Aa1x, aA1x, aa1x, AA2x, Aa2x, aA2x, aa2x, sAa, saA,
                      m12, m21, sAA, saa)[[4]])
  }
}

selectionAfter<-function(genos1, genos2, sAa, saA, sAA, saa, sBB){
  s1Avg <- sum(genos1) - sAA*sum(genos1[1:4]) -sAa*sum(genos1[5:8]) -
    saA*sum(genos1[9:12]) - saa*sum(genos1[13:16])
  s2Avg <- sum(genos2) - sAA*sum(genos2[1:4]) -sAa*sum(genos2[5:8]) -
    saA*sum(genos2[9:12]) - saa*sum(genos2[13:16])
  genos1s <- rep(0, 16);
  genos2s <- rep(0, 16);
  scoff <- c(sAA, sAA, sAA, sAA, sAa, sAa, sAa, sAa, saA, saA, saA, saA,
               saa, saa, saa, saa);
  for(i in 1:16){ 
  genos1s[i] = (genos1[i]*(1 - scoff[i]))/s1Avg; 
  genos2s[i] = (genos2[i]*(1 - scoff[i]))/s2Avg
  }
  return(list(genos1s, genos2s))
}

migrationAfter<-function(genos1, genos2, m12, m21){ 
  genos1mig = rep(0,16);
  genos2mig = rep(0,16);
  genos1m = rep(0,16);
  genos2m = rep(0,16);
  for(i in 1:16){
    genos1mig[i] = genos1[i]*(1 - m12) + genos2[i]*m21; 
    genos2mig[i] = genos2[i]*(1 - m21) + genos1[i]*m12;
  }
  sum1 = sum(genos1mig);
  sum2 = sum(genos2mig);
  for(i in 1:16){
    genos1m[i] = genos1mig[i]/sum1; 
    genos2m[i] = genos2mig[i]/sum2;
  }
  return(list(genos1m, genos2m))
}

gametes<-function(geno,r){
  arrayg<-rep(0,4);
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
  else if(geno==16){arrayg<-c(0., 0., 0., 1.)}
  return(arrayg)
}

matingResult<-function(genos, genof, genom, adv, r){
  coefficient <- c(1., 1., 1., 1.)
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
  #There are 16x16 possible matings and we calculate the genotype
  #frequencies resulted from each mating then sum up all results for
  #post-mating genotype frequencies. Somehow, 
  #double sum does not work on Mathematica*)
  
  newgenos1<-rep(0,16)
  newgenos2<-rep(0,16)
  for(i in 1:16){
    for(j in 1:16){
      newgenos1 <- newgenos1 + matingResult(genos1,i,j,adv1,r)
      newgenos2 <- newgenos2 + matingResult(genos2,i,j,adv2,r)
    }
  }
  return(list(newgenos1, newgenos2))
}

generationAfter<-function(genos1, genos2, m12, m21, sAa, saA, sAA, saa, 
                sBB, cost, d1, d2, r, model){
  genos1s <- selectionAfter(genos1, genos2, sAa, saA, sAA, saa, sBB)[[1]]
  genos2s <- selectionAfter(genos1, genos2, sAa, saA, sAA, saa, sBB)[[2]]
  genos1m <- migrationAfter(genos1s, genos2s, m12, m21)[[1]];
  genos2m <- migrationAfter(genos1s, genos2s, m12, m21)[[2]]
  return(matingAfter(genos1m, genos2m, d1, d2, r, model, cost))
}

bAllele<- function(genos1x){ 
  return(
    genos1x[1] + genos1x[5] + genos1x[9] + genos1x[13] + genos1x[2]/2 + 
    genos1x[3]/2 + genos1x[6]/2 + genos1x[7]/2 + genos1x[10]/2 + 
    genos1x[11]/2 + genos1x[14]/2 + genos1x[15]/2
  )
}

after<-function(m12, m21, sAa, saA, sAA, saa, sBB, cost, d1, d2, r, b, model, background){
  p1<-equilibrium(1., 0., 0., 0., 0., 0., 0., 1., sAa, saA, m12, m21, sAA, saa);
  q1 = 1 - p1;
  p2 = q1; q2 = p1;
  
  #genos1 is a vector representing the initial equilibrial
  #frequency of genotypes in population 1 after the introduction of
  #alleles B.
  if(background == 1){ 
   genos1 <- c(p1*p1*b*b, p1*p1*b*(1 - b), p1*p1*(1 - b)*b, 
     p1*p1*(1 - b)*(1 - b), p1*q1*b*b, p1*q1*b*(1 - b), 
     p1*q1*(1 - b)*b, p1*q1*(1 - b)*(1 - b), q1*p1*b*b, 
     q1*p1*b*(1 - b), q1*p1*(1 - b)*b, q1*p1*(1 - b)*(1 - b), 
     q1*q1*b*b, q1*q1*b*(1 - b), q1*q1*(1 - b)*b, 
     q1*q1*(1 - b)*(1 - b));
   genos2 <- c(0., 0., 0., p2*p2, 0., 0., 0., p2*q2, 0., 0., 0., q2*p2,
     0., 0., 0., q2*q2)
  }
  else{
    genos1 <- c(0., 0., 0., p1*p1, 0., 0., 0., p1*q1, 0., 0., 0., q1*p1,
       0., 0., 0., q1*q1);
    genos2 <- c(p2*p2*b*b, p2*p2*b*(1 - b), p2*p2*(1 - b)*b, 
     p2*p2*(1 - b)*(1 - b), p2*q2*b*b, p2*q2*b*(1 - b), 
     p2*q2*(1 - b)*b, p2*q2*(1 - b)*(1 - b), q2*p2*b*b, 
     q2*p2*b*(1 - b), q2*p2*(1 - b)*b, q2*p2*(1 - b)*(1 - b), 
     q2*q2*b*b, q2*q2*b*(1 - b), q2*q2*(1 - b)*b, 
     q2*q2*(1 - b)*(1 - b));
  }
  
  nestHelp<-function(genos1x, genos2x){
    pop1Binit<-bAllele(genos1x);
    pop2Binit<-bAllele(genos2x);
    newGenos1<- generationAfter(genos1x, genos2x, m12, m21, sAa, saA, sAA, saa, 
                                 sBB, cost, d1, d2, r, model)[[1]] 
    newGenos2 <- generationAfter(genos1x, genos2x, m12, m21, sAa, saA, sAA, saa, 
                      sBB, cost, d1, d2, r, model)[[2]]
    dB1 <- abs(pop1Binit - bAllele(newGenos1))
    dB2 <- abs(pop2Binit - bAllele(newGenos2))
    return(list(newGenos1, newGenos2, dB1, dB2))
  }
  genos1x <- genos1
  genos2x <- genos2
  lil <- 1
  db1<-1
  db2<-1
  if(p1 != 0 & p2 != 0){
  #CHANGE LIL LIMIT LATER!!
    while((lil<1e3) & (db1 > 0) & (db2 > 0)){
      result <- nestHelp(genos1x,genos2x)
      genos1x <- result[[1]]
      genos2x <- result[[2]]
      db1 <- result[[3]]
      db2 <- result[[4]]
      lil <- lil+1
    }
  }
  A1 <- 0
  A2 <- 0
  B1 <- 0
  B2 <- 0
  for(i in 1:4){
    A1 <- A1 + genos1x[i] + genos1x[4+i]/2 + genos1x[8+i]/2
    A2 <- A2 + genos2x[i] + genos2x[4+i]/2 + genos2x[8+i]/2
    B1 <- B1 + genos1x[4*(i-1)+1] + genos1x[4*(i-1)+2]/2 + genos1x[4*(i-1)+3]/2
    B2 <- B2 + genos2x[4*(i-1)+1] + genos2x[4*(i-1)+2]/2 + genos2x[4*(i-1)+3]/2
  }
  return(c(m12,m21,d1,cost,sAa,saA,p1,p2,A1,A2,B1,B2,lil))
}

coreAfter <- function(i){
  print(i)
  AA <- i[1]
  aa <- i[2]
  cc <- i[3]
  sell1 <- i[4]
  sell2 <- i[5]
  return(after(AA, aa, sell1, sell2, 0., 0., 0., cc, 1.0, 0., 0.05, 0.01, 1, 1))
}

#0.03, 0.03, AA, aa, 0., 0., 0., cc, 1, 0., 0.05, 0.01, 1, 1 (for the sel run)

#CHANGE THIS
coeff.1 <- seq(0.,0.15,0.025)
coeff.2 <- seq(0.,0.15,0.025)
selcost <- seq(0.,0.3,0.01)
select1 <- seq(0.1,0.9,0.1)
select2 <- seq(0.1,0.9,0.1)

sel <- list()
loy <- 1
for(i in coeff.1){
  for (j in coeff.2){
    for (k in selcost){
      for (l in select1){
        for (m in select2){
          sel[[loy]] <- c(i,j,k,l,m)
          loy<-loy+1
        }
      }
    }
  }
}

results <- mclapply(sel, coreAfter, mc.cores = 75)
results.df <- as.data.frame(do.call(rbind, results))
colnames(results.df) <- c('m_12','m_21','d','cost','s_Aa','s_aA','A1_init','A2_init','A_1','A_2','B_1','B_2','gen')
write.table(results.df, file = "Asymmetric_migration_melange.txt", sep = "\t",
            row.names = FALSE)
