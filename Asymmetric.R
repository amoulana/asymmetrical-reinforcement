library(parallel)

#We defne a couple of cuntions. Calc.x functions calculate the frequency of allele x in the population
#of interest.
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

#Assuming Hardy-Weinberg, then we define a function
#to calculate the genotype frequency based on allele freq
random.mating <- function(A,a,C,c){
  local.ad <- c(A^2,A*a,A*a,a^2,A^2,A*a,A*a,a^2)
  #freqs of C and c just stay constant 
  return(c(C,C,C,C,c,c,c,c)*local.ad)
}

#The sequence of genotype in a vector: without B allele-- CAA,CAa,CaA,Caa,cAA,cAa,caA,caa
#when there's B allele-- CAABB,CAaBB,CaaBB,... (32)

#First, the selection function (before the intro of B allele).
#pop1 and pop2: the genotype freq vector in pop 1 and pop 2
#sel 1 and 2: selection coefficient for allele C/c being in the other population (local adaptation).
#Enter a NUMERIC
#sel.comp: selection VECTOR due to incompatibility between cytoplasmic and nuclear allele. This must
#be in the form c(0,w,w,x,y,z,z,0)
selection <- function(pop1,pop2,sel1,sel2,sel.comp){
  sel1 <- c(0,0,0,0,sel1,sel1,sel1,sel1)
  sel2 <- c(sel2,sel2,sel2,sel2,0,0,0,0)
  s1Avg <- sum((1-sel1)*(1-sel.comp)*pop1)
  s2Avg <- sum((1-sel2)*(1-sel.comp)*pop2)
  pop1s <- (1-sel1)*(1-sel.comp)*pop1/s1Avg
  pop2s <- (1-sel2)*(1-sel.comp)*pop2/s2Avg
  
  #return a list of two VECTORS: pop1 and pop2 after selection
  return(list(pop1s,pop2s))
}

#migration function with input pop1, pop2, and NUMERICAL m12 and m21 which are mig coeffs
#migration coeffs are constant across genotypes in a population
migration<- function(pop1, pop2,m12,m21){
  #calculate the sum of freqs in the populations after migration.
  #m12: mig from pop 1 to pop2, and vice cersa (m12 fraction of pop 1 replaced by pop2)
  #after migration the freq in pop1 = freq staying * pop1 + freq coming * pop2
  pop1m <- (pop1*(1-m12) + pop2*m12)
  pop2m <- (pop1*m21 + pop2*(1-m21))
  return(list(pop1m,pop2m))
}

#mating fuction with input pop1 and pop2
mating<-function(pop1,pop2){
  #calculate allele freqs
  A1 <- calc.A(pop1)
  a1 <- calc.a(pop1)
  A2 <- calc.A(pop2)
  a2 <- calc.a(pop2)
  C1 <- calc.C(pop1)
  c1 <- calc.c(pop1)
  C2 <- calc.C(pop2)
  c2 <- calc.c(pop2)
  
  #apply random mating functions to get after-mating freqs of pop1 and pop2
  pop1t <- random.mating(A1,a1,C1,c1)
  pop2t <- random.mating(A2,a2,C2,c2)
  return(list(pop1t,pop2t))
}

#sum-up function with input pop1, pop2 being the initial geno freqs before selection
generation<-function(pop1,pop2,sel1,sel2,sel.comp,m12,m21){
  selected <- selection(pop1,pop2,sel1,sel2,sel.comp)
  mig <- migration(selected[[1]],selected[[2]],m12,m21)
  mat <- mating(mig[[1]],mig[[2]])
  
  #calculate change in allele A freqs
  delta.A1 <- calc.A(mat[[1]])-calc.A(pop1)
  delta.A2 <- calc.A(mat[[2]])-calc.A(pop2)
  A1 <- calc.A(mat[[1]])
  A2 <- calc.A(mat[[2]])
  
  #return a list of the resulting geno freq, change, and final allele freq
  return(list(mat, delta.A1, delta.A2, A1, A2))
}

#a function that runs function GENERATION up until an equilibrium
equilibrium<-function(pop1,pop2,sel1,sel2,sel.comp,m12,m21){
  pop1x <- pop1
  pop2x <- pop2
  #calculate change in allele A freq
  arraydp <- generation(pop1x,pop2x,sel1,sel2,sel.comp,m12,m21)[[3]];
  #continue if change 
  if(arraydp > 0.00001){
    vec.x <- generation(pop1x,pop2x,sel1,sel2,sel.comp,m12,m21)[[1]]
    equilibrium(vec.x[[1]],vec.x[[2]],sel1,sel2,sel.comp,m12,m21)
  }
  
  #stop and return if no change
  else{
    return(generation(pop1x,pop2x,sel1,sel2,sel.comp,m12,m21))
  }
}

#A similar function to selection, but now we have genos1 and genos2
#These VECTORS are similar to pop1 and pop2, but have length of 32 instead of 8, adding the B allele
#other inputs are similar
selectionAfter<-function(genos1, genos2, sel1, sel2, sel.comp){
  #for each B/b genotype, find the sum of the freqs stored in avg1 and avg2 
  avg1 <- c(sum(genos1[1:8])/sum(genos1),sum(genos1[9:16])/sum(genos1),
            sum(genos1[17:24])/sum(genos1),sum(genos1[25:232])/sum(genos1))
  avg2 <- c(sum(genos2[1:8])/sum(genos2),sum(genos2[9:16])/sum(genos2),
            sum(genos2[17:24])/sum(genos2),sum(genos2[25:232])/sum(genos2))
  #perform SELECTION function on each of B/b genotype. 
  #Note that since selection is independent of B/b alleles,
  #we can separate this function into all four possible B/b genotypes, and operate the function on each
  genos.BB <- selection(genos1[1:8],genos2[1:8],sel1,sel2,sel.comp)
  genos.Bb <- selection(genos1[9:16],genos2[9:16],sel1,sel2,sel.comp)
  genos.bB <- selection(genos1[17:24],genos2[17:24],sel1,sel2,sel.comp)
  genos.bb <- selection(genos1[25:32],genos2[25:32],sel1,sel2,sel.comp)

  #reconvene and normalize the freqs by the sum of freqs initially for each B/b genotype
  genos1.s <- c(genos.BB[[1]]/avg1[[1]],genos.Bb[[1]]/avg1[[2]],
                genos.bB[[1]]/avg1[[2]],genos.bb[[1]]/avg1[[3]])
  genos2.s <- c(genos.BB[[2]]/avg2[[1]],genos.Bb[[2]]/avg2[[2]],
                genos.bB[[2]]/avg2[[2]],genos.bb[[2]]/avg2[[3]])
  return(list(genos1.s, genos2.s))
}

#migration is really similar to the previous function, since m12 and m21 are independent of the genos
migrationAfter<-function(genos1, genos2, m12, m21){
  genos1mig <- genos1*(1-m12)+genos2*m12
  genos2mig <- genos2*(1-m21)+genos1*m21
  sum1 = sum(genos1mig)
  sum2 = sum(genos2mig)
  genos1m = genos1mig/sum1
  genos2m = genos2mig/sum2
  return(list(genos1m, genos2m))
}

#create a function that calculates the frequency of certain gametes after recombination.
#gametes: AB, Ab, aB, ab
#input geno denotes genotype of interest: 1=AABB, 2=AABb, etc. up until 16. r = recombination coeff
#This repeats twice up until 32 due to C/c allele addition.
gametes<-function(geno,r){
  arrayg<-rep(0,4)
  #Ignore C/c alleles for gamete production
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

#mating result given the genos (genotype freqs), genof (type of genotype for female), genom (male)
#adv is a VECTOR that gives the frequency of mating when the BB genof mates (with mating preference)
matingResult<-function(genos, genof, genom, adv, r){
  #set freq mating to normal
  coefficient <- c(1., 1., 1., 1.)
  genof <- genof %% 16
  genom <- genom %% 16
  ifelse(genof==0,genof<-16,genof<-genof)
  ifelse(genom==0,genom<-16,genom<-genof)
  
  #change freq mating (coefficient) to adv when genof=___BB
  if(genof == 1 | genof == 5 | genof == 9 | genof == 13){ 
    coefficient <- adv;
  }
  
  #calculate frequency of mating now:
  if(genom < 5) {maleAallele <- 1} 
  else if(genom < 9){maleAallele <- 2} 
  else if(genom < 13){maleAallele <- 3}
  else{maleAallele <- 4};
  freqMating <- genos[genof]*genos[genom]*coefficient[maleAallele];
  
  #calculate gamete frequency for male and female
  gametef <- gametes(genof, r);
  gametem <- gametes(genom, r);
  
  #calculate the resulting genotype of the mating
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

#the ACTUAL mating function
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
  
  #ALWAYS use model 1 where BB prefers AA not aa (as in model 2)
  #adv gies the relative mating frequency for each A/a genotype. incorporate cost as well
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
  #post-mating genotype frequencies. C mothers produce C offsprings and c mothers produce c offsprings
  #due to maternal inhertiance
  
  X.C <- c(1:16)
  X.c <- c(17:32)
  Y <- c(1:32)
  
  newgenos1.C <- rep(0,16)
  newgenos1.c <- rep(0,16)
  newgenos2.C <- rep(0,16)
  newgenos2.c <- rep(0,16)
  
  #i is for genof, j for genom. Calculate this for C genotypes. Sum all over possible genof/genom.
  #Note that C genotypes can only come from C mothers
  for(i in X.C){
    for(j in Y){
      newgenos1.C <- newgenos1.C + matingResult(genos1,i,j,adv1,r)
      newgenos2.C <- newgenos2.C + matingResult(genos2,i,j,adv2,r)
    }
  }
  
  #Now for c
  for(i in X.c){
    for(j in Y){
      newgenos1.c <- newgenos1.c + matingResult(genos1,i,j,adv1,r)
      newgenos2.c <- newgenos2.c + matingResult(genos2,i,j,adv2,r)
    }
  }
  
  #combine between the two
  newgenos1 <- c(newgenos1.C,newgenos1.c)
  newgenos2 <- c(newgenos2.C,newgenos2.c)
  return(list(newgenos1, newgenos2))
}

#combine all functions
generationAfter<-function(genos1, genos2, m12, m21, sel1, sel2, sel.comp, cost, d1, d2, r, model){
  genos1s <- selectionAfter(genos1, genos2, sel1, sel2, sel.comp)[[1]]
  genos2s <- selectionAfter(genos1, genos2, sel1, sel2, sel.comp)[[2]]
  genos1m <- migrationAfter(genos1s, genos2s, m12, m21)[[1]];
  genos2m <- migrationAfter(genos1s, genos2s, m12, m21)[[2]]
  return(matingAfter(genos1m, genos2m, d1, d2, r, model, cost))
}

#create a function to calculate B allele freq.
B.Allele<- function(genos1x){ 
  return(
    #some +16 to account for c genotypes
    genos1x[1] + genos1x[5] + genos1x[9] + genos1x[13] + genos1x[2]/2 + 
      genos1x[3]/2 + genos1x[6]/2 + genos1x[7]/2 + genos1x[10]/2 + 
      genos1x[11]/2 + genos1x[14]/2 + genos1x[15]/2+
      genos1x[1+16] + genos1x[5+16] + genos1x[9+16] + genos1x[13+16] + genos1x[2+16]/2 + 
      genos1x[3+16]/2 + genos1x[6+16]/2 + genos1x[7+16]/2 + genos1x[10+16]/2 + 
      genos1x[11+16]/2 + genos1x[14+16]/2 + genos1x[15+16]/2
  )
}

#create a function to combine the initial model (up until requilibirum), introduce B allele, and run
#simulation up until N generations
after<-function(m12, m21, sel1, sel2, sel.comp, cost, d1, d2, r, b,model,background,N){
  
  #find equilibrium freqs.
  pop1.init <- c(1, 0, 0, 0, 0, 0, 0, 0)
  pop2.init <- c(0, 0, 0, 0, 0, 0, 0, 1)
  eq.list<-equilibrium(pop1.init,pop2.init,sel1,sel2,sel.comp,m12,m21)
  pops.eq <- eq.list[[1]]
  pop1.eq.before <- pops.eq[[1]]
  pop2.eq.before <- pops.eq[[2]]
  
  #expand vector to 32-element
  pop1.eq <- rep(pop1.eq.before,each=4)
  pop2.eq <- rep(pop2.eq.before,each=4)
  
  #genos1 is a vector representing the initial equilibrial
  #frequency of genotypes in population 1 after the introduction of
  #alleles B.
  
  b.intros <- rep(c(b^2,b*(1-b),b*(1-b),(1-b)^2),8)
  
  if(background == 1){ 
    genos1.eq <- pop1.eq * b.intros
    genos2.eq <- pop2.eq
  }
  else{
    genos1.eq <- pop1.eq 
    genos2.eq <- pop2.eq * b.intros
  }
  
  #create a function that calculates the change in allele B frequency (not necessary, but just in case)
  nestHelp<-function(genos1x, genos2x){
    pop1Binit<-B.Allele(genos1x);
    pop2Binit<-B.Allele(genos2x);
    newGenos1<- generationAfter(genos1x, genos2x, m12, m21, 
                                sel1, sel2, sel.comp, cost, d1, d2, r, model)[[1]] 
    newGenos2 <- generationAfter(genos1x, genos2x, m12, m21, 
                                 sel1, sel2, sel.comp, cost, d1, d2, r, model)[[2]]
    dB1 <- abs(pop1Binit - B.Allele(newGenos1))
    dB2 <- abs(pop2Binit - B.Allele(newGenos2))
    return(list(newGenos1, newGenos2, dB1, dB2))
  }
  
  #run the simulation up until certain N generations
  N.gen <- 1
  #initialize change in allele B freq to be greater than 0 first
  db1<-1
  db2<-1
  
  #only proceed if at eq, both alleles A and a maintained
  if(calc.A(pop1.eq.before) != 0 & calc.A(pop2.eq.before) != 0){
    #run sim
    while((N.gen<N) & (db1 > 0) & (db2 > 0)){
      result <- nestHelp(genos1.eq,genos2.eq)
      genos1.eq <- result[[1]]
      genos2.eq <- result[[2]]
      db1 <- result[[3]]
      db2 <- result[[4]]
      N.gen <- N.gen+1
    }
  }
  
  #calculate frequencies of alleles A and B in pop 1 and 2
  A1 <- 0
  A2 <- 0
  B1 <- 0
  B2 <- 0
  for(i in 1:8){
    A1 <- A1 + genos1.eq[i] + genos1.eq[4+i]/2 + genos1.eq[8+i]/2
    A2 <- A2 + genos2.eq[i] + genos2.eq[4+i]/2 + genos2.eq[8+i]/2
    B1 <- B1 + genos1.eq[4*(i-1)+1] + genos1.eq[4*(i-1)+2]/2 + genos1.eq[4*(i-1)+3]/2
    B2 <- B2 + genos2.eq[4*(i-1)+1] + genos2.eq[4*(i-1)+2]/2 + genos2.eq[4*(i-1)+3]/2
  }
  
  #return all the params
  return(c(m12,m21,d1,cost,
           sel1,sel2,
           sel.comp[2],sel.comp[3],sel.comp[4],sel.comp[5],
           calc.A(pop1.eq.before),calc.A(pop2.eq.before),
           A1,A2,B1,B2,N.gen))
}

#create a function to run multiple simulations simultaneously with an input of a vector length 11
coreAfter <- function(i){
  sim.m12 <- i[1]
  sim.m21 <- i[2]
  sim.sel1 <- i[3]
  sim.sel2 <- i[4]
  sim.sel.comp <- c(0,i[5],i[5],i[6],i[7],i[8],i[8],0)
  sim.d <- i[9]
  sim.r <- i[10]
  sim.cost <- i[11]
  return(after(sim.m12,sim.m21,sim.sel1, sim.sel2, sim.sel.comp, sim.cost,sim.d,0,sim.r,0.01,1,1,100))
}

#create all possible combinations of params (in this model there are 11!)
seq.m12 <- seq(0.1,0.1,0.0)
seq.m21 <- seq(0.1,0.1,0.0)
seq.sel1 <- seq(0.1,0.1,0.0)
seq.sel2 <- seq(0.1,0.1,0.0)
seq.sel.comp1 <- seq(0.1,0.1,0.0)
seq.sel.comp2 <- seq(0.1,0.1,0.0)
seq.sel.comp3 <- seq(0.1,0.1,0.0)
seq.sel.comp4 <- seq(0.1,0.1,0.0)
seq.d <- seq(0.1,0.1,0.0)
seq.r <- seq(0.1,0.1,0.0)
seq.cost <- seq(0.1,0.1,0.0)
sim.list.0 <- expand.grid(seq.m12,seq.m21,seq.sel1,seq.sel2,
                        seq.sel.comp1,seq.sel.comp2,seq.sel.comp3,seq.sel.comp4,
                        seq.d,seq.r,seq.cost)
sim.list <- as.list(as.data.frame(t(sim.list.0)))

#run multiple simulations
results <- mclapply(sim.list, coreAfter, mc.cores = 1)
results.df <- as.data.frame(do.call(rbind, results))
colnames(results.df) <- c('m_12','m_21','d','cost',
                          's_c1','s_C2',
                          's_CAa','s_Caa','s_cAA','s_cAa',
                          'A1_eq','A2_eq',
                          'A1_fin','A2_fin','B1_fin','B2_fin','numb_gen')
write.table(results.df, file = "Asymmetric_migration_cytoplasmic.txt", sep = "\t",
            row.names = FALSE)
