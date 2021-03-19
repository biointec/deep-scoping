library(dplyr)
library(ggplot2)
library(egg)
library(gridExtra)

#map<-"own_results"
#map<-"data/Models"
map<-"data/Genenbank_grootte"
  
if (file.exists(paste0(map,"/Results"))){
  unlink(paste0(map,"/Results"), recursive=TRUE)
}

files<-dir(map, full.names = T)
names<-dir(map, full.names = F)
results<-c()
level<-c()
fac<-c()
dir.create(paste0(map,"/Results"))


iter<-0
startbc<-c(5,5,5,5,5)
endbc<-c(50,50,50,50,50)
select.exp<-c(1,2,3,4,5)
#
#select.exp<-c(76,88,68,80,274,275,456)
#exp.name<-c("Deep scoping method (80/10/10 BC5)", "Deep scoping method (80/10/10 BC20)", "Deep scoping method (50/10/40 BC5)", "Deep scoping method (50/10/40 BC20)","HUC method (80/20/00 BC5)","HUC method (80/20/00 BC20)","Truncation selection")
#startbc<-c(5,20,5,20,5,20,1)
#select.exp<-c(263,68,396)
#endbc<-c(50,50,50,50,50,50,20)
exp.name<-c("Size 0200", "Size 1000", "Size 2000" , "Size 3000","Size 5000")

for(expir in c(select.exp)){
  iter<-iter+1
  load(files[expir])
  a<-names[expir]
  n.set<-length(experiment.sub.results)
  
  acc_fixed<-c()
  for(s in seq(n.set)){
    set <- experiment.sub.results[[s]]
    repnames<-names(set)
    if(!is.null(set)  & (class(set) != "try-error")){
    nrep<-length(set)
    for(r in seq(nrep)){ 
      rep<-set[[r]]
      rname<-as.numeric(substr(repnames[r],4,nchar(repnames[r])))
      if (!any(is.na(rep))){
      genome<-rep$genome
      sim<-rep$sim.results
      
     
      out<-c()    
      geneticValue<-c()
      geneticValue.top<-c()
      
      
     # for(br.cycle in seq(metadata$n.cycles)){
      for(br.cycle in seq(startbc[iter],endbc[iter])){
      #  if(br.cycle<=startbc[iter]){
        #  gen<-rep_baseline$sim.results[[br.cycle]]
      #  }else{
          cycle.name<-paste0("cycle",br.cycle)
          gen<-sim[[cycle.name]]
      #  }
        
         
          candidate.af.i<-gen$geno.summary.stats$candidate.af
          geneticValue<-cbind(geneticValue, gen$candidate.values$mu.g/gen$QTL_info[['maximum']])
          geneticValue.top<-cbind(geneticValue.top, mean(sort(gen$candidate.values$geno.values,decreasing = TRUE)[1:10])/gen$QTL_info[['maximum']])

        
         
          #out<-cbind(out, gen$QTL_top)
          #if(is.null(out)){
          #  out<-cbind(out, gen$QTL_info)
          #}
          out<-cbind(out, gen$QTL_info)
          
         
      }

      Data<-data.frame(t(out))
      
      Data$cycle<-seq(startbc[iter],endbc[iter])
     
      # nodige data in 1 column zetten met factor for ggplot legend
      #correct the breeding cycle for the preselected results
 
      cycle<-Data$cycle#+((iter-1)*5)
      
      m1<-cbind(cycle, Data$max_reachable/Data$maximum, t(geneticValue.top), iter)
      
     
      m<-rbind(m1)
      m<-data.frame(m)
      colnames(m)<-c("cycle","Maximum_Reachable_GV","Top_GV","Legend")
      results<-rbind(results, m)
      }
     }
   }
  }
  level<-c(level, toString(iter))
  
  fac<-c(fac,exp.name[iter])

}#Close the for expir loop
results$Legend<-factor(results$Legend, levels=seq(length(fac)),labels = fac)


#bereken de gemiddelde per generatie en de std om deze dan te plotten
sum.result<-data.frame(cycle = seq(metadata$n.cycles))
  

result.cycle<-c()
 
for(b.c in seq(min(startbc),max(endbc))){
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(Maximum_Reachable_GV), SD = sd(Maximum_Reachable_GV))->temp.cycle
    temp.cycle$Type<-"Maximum Reachable GV"
    temp.cycle$cycle<-b.c
    result.cycle<-rbind(result.cycle, temp.cycle)
    
    
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(Top_GV), SD = sd(Top_GV))->temp.cycle
    temp.cycle$Type<-"Top-10 GV"
    temp.cycle$cycle<-b.c
    result.cycle<-rbind(result.cycle, temp.cycle)
    
}

result.cycle$Type<-paste0(result.cycle$Type," (",result.cycle$Legend,")")


save("result.cycle",file="GeneBank_size.RData")

title<-"Greedy parental selection with prebreeding"
f1<-ggplot(data=result.cycle, aes(x=cycle, y=GV, colour=Type,linetype=Type)) +
    geom_line(size=1.2)  +
    #geom_errorbar(aes(ymin=GV-SD, ymax=GV+SD), width=.2,position=position_dodge(.9))+ 
    scale_linetype_manual(values=c(rep("twodash",length(select.exp)),rep("solid",length(select.exp))))+
    #scale_colour_manual(values=rep(seq(length(select.exp)),2))+
    #scale_colour_manual(values=rep(c("#b2182b","#f4a582","#fddbc7","#d1e5f0","#92c5de","#2166ac","black"),2))+
    scale_colour_manual(values=rep(c("#b2182b","#f4a582","#d1e5f0","#2166ac","#2166ac"),2))+
    ylab('Genetic Value (GV)')+
    xlab('Breeding Cycle')+
    theme_bw()+
    scale_y_continuous(breaks = seq(0.0, 1, 0.2),expand=c(0,0),limits=c(0.4,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
    scale_x_continuous(breaks = seq(0, max(endbc), 5),expand=c(0,0),limits=c(5,max(endbc)),sec.axis = dup_axis(name = NULL, labels = NULL))+
    theme(legend.title =  element_blank(),text = element_text(size=12) , legend.position= c(0.70,0.20),legend.direction = "vertical")+
    guides(linetype=guide_legend(ncol=2))+
    theme(axis.text = element_text(colour="black"))+
    theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
    theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"),legend.text = element_text(size=12))+
    theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
    theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
    theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
    
  
  


titel<-"Effects of prebreeding on a greedy selection"
filename<-paste0(map,"/Results/",titel,".eps")
ggsave(filename, height=8 ,width= 12, f1,device="eps")

filename<-paste0(map,"/Results/",titel,".png")
ggsave(filename, height=8 ,width= 12, f1,device="png")

filename<-paste0(map,"/Results/",titel,".pdf")
ggsave(filename, height=8 ,width= 12, f1,device="pdf")


