# Exploratory sensitivity checks only; does not modify manuscript models.
suppressPackageStartupMessages({library(lmerTest);library(emmeans)})
out<-'output/reviews/biomass_model_audit'
o<-readRDS(file.path(out,'objects.rds'));b<-o$data
results<-list(); anovas<-list(); outliers<-list(); k<-0
for(sp in unique(b$Species))for(tp in c('T1','TF')) {
 d<-droplevels(subset(b,Species==sp & Time.Point==tp))
 for(clade in c(FALSE,TRUE))for(structure in c('parent+tank','parent','tank','none')) {
 fixed<-if(clade)'Treatment*Clade'else'Treatment'
 random<-switch(structure,'parent+tank'='+(1|Parent.ID)+(1|Tank.num)',parent='+(1|Parent.ID)',tank='+(1|Tank.num)',none='')
 form<-as.formula(paste('AFDWmg.cm2~',fixed,random))
 f<-if(structure=='none')lm(form,d)else lmer(form,d)
 q<-as.data.frame(pairs(emmeans(f,~Treatment),reverse=TRUE,adjust='none'))
 k<-k+1;results[[k]]<-data.frame(Species=sp,Time.Point=tp,clade_adjusted=clade,random=structure,n=nrow(d),singular=if(structure=='none')NA else isSingular(f),q)
 if(structure=='parent+tank') {
  if(clade) {
   a<-as.data.frame(anova(f));a$term<-rownames(a);a$Species<-sp;a$Time.Point<-tp;anovas[[paste(sp,tp)]]<-a
  } else {
   r<-residuals(f)/sigma(f);rank<-order(abs(r),decreasing=TRUE)[1:2]
   outliers[[paste(sp,tp)]]<-data.frame(d[rank,c('Frag.ID','Parent.ID','Species','Time.Point','Treatment','Tank.num','AFDWmg.cm2')],scaled_residual=r[rank])
  }
 }
 }
}
write.csv(do.call(rbind,results),file.path(out,'random_effect_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,anovas),file.path(out,'time_specific_clade_anova.csv'),row.names=FALSE)
write.csv(do.call(rbind,outliers),file.path(out,'largest_residual_observations.csv'),row.names=FALSE)
# Compare raw/log scales for P. acuta descriptively, without selecting by p-value.
pacu<-subset(b,Species=='Pocillopora acuta' & Time.Point!='T0')
checks<-list()
png(file.path(out,'pacu_raw_log_diagnostics.png'),width=1600,height=1600,res=150);par(mfrow=c(4,2),mar=c(4,4,3,1))
for(tp in c('T1','TF'))for(response in c('AFDWmg.cm2','log(AFDWmg.cm2)')) {
 d<-droplevels(subset(pacu,Time.Point==tp));f<-lmer(as.formula(paste(response,'~Treatment+(1|Parent.ID)+(1|Tank.num)')),d)
 r<-residuals(f);plot(fitted(f),r,main=paste(tp,response),xlab='Fitted',ylab='Residual');abline(h=0,lty=2);qqnorm(r,main='Residual Q-Q');qqline(r)
 q<-as.data.frame(pairs(emmeans(f,~Treatment),reverse=TRUE,adjust='none'));checks[[paste(tp,response)]]<-data.frame(Time.Point=tp,response=response,singular=isSingular(f),q)
}
dev.off()
write.csv(do.call(rbind,checks),file.path(out,'pacu_raw_log_sensitivity.csv'),row.names=FALSE)
# Record source-data exclusions after the existing metadata filters.
m<-o$environment$meta
write.csv(m[!m$Frag.ID %in% b$Frag.ID,],file.path(out,'metadata_without_biomass.csv'),row.names=FALSE)
cat('Sensitivity checks complete\n')
