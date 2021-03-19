## Calculate H-score
#-------------------------------------------------------------------------
#Copyright (c) 2019 Antoine Allier
#-------------------------------------------------------------------------

H_UC_Bridging<-function(candidate.marker.genos.i,
                                 candidate.haplotype,
                                 GeneBank,
                                 GeneBank.genotype,
                                 GeneBank.Phenotype,
                                 p.sel=20,
                                 predictions.out,
                                 parents.information=parent.selections.list,
                                 trainingpanel=list(geno=TP.genos.i,pheno=TP.phenos.i),
                                 map=map,
                                 rep.iter,
                                 genome){
  
  Geno<-rbind(candidate.marker.genos.i, GeneBank.genotype)
  
  
  # Number of iterations and burn-in should be increased for proper inference
  nIter <- 1000
  burnIn <- 100
  thin <- 2
  # Run the Bayesian Ridge Regression
  WGR <- BGLR(y = rbind( trainingpanel$pheno,GeneBank.Phenotype$mean.pheno.values),
              ETA = list(list(X = rbind(trainingpanel$geno, GeneBank.genotype),model = "BRR",saveEffects = T)),
              nIter = nIter,burnIn = burnIn,thin = thin, saveAt=paste0("Temp/rep",rep.iter,"/"),verbose=F)
  # Load .bin file including MCMC samples of marker effects
  B <- readBinMat(paste0("Temp/rep",rep.iter,"/ETA_1_b.bin"))
  # Posterior means of marker effects
  Bhat <- colMeans(B)
  # Genomic estimated breeding values
  GEBV <- Geno%*%Bhat
  
  
  #get the values from the prediction model
  #GEBV <- predictions.out$GEBV
  
  # Assuming the Elites are the 10 best lines
  PopE <- names(GEBV[order(GEBV,decreasing = TRUE),])[seq(p.sel/2)]

  PopD <- sample(setdiff(rownames(GEBV),PopE),(dim(GEBV)[1]-p.sel/2),replace = FALSE)
  
 
  
  
  
  
  
  
  
  
  
  
  
  ###########
  ###########
  Haplo<-rbind(as.matrix(do.call("rbind", candidate.haplotype)),
               as.matrix(do.call("rbind",GeneBank)))
  
  positions<-c()
  begin=0
  for (i in seq(7)){
    positions<-c(positions, genome[[i]]@pos.add.qtl$ID+begin)
    begin<-begin+genome[[i]]@num.snp.chr+genome[[i]]@num.add.qtl.chr
  }
  
  Haplo<-Haplo[,-positions]
  ObjectHEBV <- GetHEBVmat(ObjectGeno = Haplo,
                           ObjectBeta = Bhat ,
                           ObjectMap = map,
                           WindowSize = WindowSizeUI,
                           StepSize = StepSizeUI)
  
  d1<-dim(ObjectHEBV$HEBV)
  HEBV<-c()
  for(i in seq(1,d1[1],by=2)){
    ind.hebv<-rbind(ObjectHEBV$HEBV[i,],ObjectHEBV$HEBV[i+1,])
    temp.hebv<-c()
    for (j in seq(d1[2])){
      temp.hebv<-c(temp.hebv,max(ind.hebv[,j]))
    }
    HEBV<-rbind(HEBV, temp.hebv)
  }
  
  line.names <- row.names(ObjectHEBV$HEBV) %>%
    str_replace(pattern = "\\.[0-9]$", replacement = "") %>%
    unique()
  
  rownames(HEBV)<-line.names
  colnames(HEBV)<- colnames(ObjectHEBV$HEBV)
  
  ObjectHEBV<-list(HEBV=HEBV,POSITION=ObjectHEBV$POSITION)
  ########
  ########
  
  
  
  SelDonor <- c()
  
  HSelDonor <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,# elites only
                         ListPopE = PopE,
                         ListPopD = NULL)
  
  invisible(lapply(1:(p.sel/2),function(iter){
    tmp <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,
                    ListPopE = c(SelDonor,PopE),# PopE + selected donors
                    ListPopD = PopD[!PopD%in%c(SelDonor)] # yet unselected donors
    )
    #Increment with selected donor
    SelDonor <<- c(SelDonor,tmp[order(tmp$H,decreasing =TRUE),]$LINE[1])
    # Increment with H criterion value of elites + selected donors
    HSelDonor <<- rbind(HSelDonor,tmp[order(tmp$H,decreasing =TRUE),][1,])
  }))
  
  # Assuming we have selected the donors:
  SelDonors4UC <- HSelDonor$LINE[seq(2,p.sel/2+1)]
  # For each donor predict UC for all possible crosses
  ResUC_pre <- do.call(rbind,lapply(SelDonors4UC,function(Donor){
    return(do.call(rbind,lapply(PopE,function(Recipient){
      UCtmp <- UC(Bmat = B,
                  ObjectGenoP1 = Geno[Donor,],
                  ObjectGenoP2 = Geno[Recipient,],
                  ObjectMap = map,
                  Psel = 0.05,
                  h = 1)
      return(data.frame(DONOR = Donor,
                        RECIPIENT = Recipient,
                        UC = UCtmp))
    })))
  }))
  
  cros<-data.frame(ResUC_pre %>%
                   group_by(DONOR)%>%
                   top_n(n = 1)) 
  
  
  
  parent.sel<-parents.information$parent.selections.i
  cros.block<-parents.information$crossing.block.i
  
  block<-rbind(cbind(as.character(cros$DONOR), as.character(cros$RECIPIENT)))
  cros.block[seq((100-(p.sel))/2+1,50),]<-block
  
  parent.sel$lines.sel<-c(parent.sel$lines.sel,block[,1],block[,2])
  
  parent.sel$value.sel<-rbind(parent.sel$value.sel, as.matrix(GeneBank.Phenotype$mean.pheno.values[cros$DONOR,]),as.matrix(GEBV[cros$RECIPIENT,]))
  
  
  
  
  
  
  

  
  
  return(list(parent.selections.i=parent.sel,crossing.block.i=cros.block))
  
}
