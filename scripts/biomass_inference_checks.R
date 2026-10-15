suppressPackageStartupMessages({library(lmerTest);library(emmeans)})
out<-'output/reviews/biomass_model_audit';o<-readRDS(file.path(out,'objects.rds'))
results<-list();k<-0
for(id in names(o$fits)) {
 f<-o$fits[[id]]; v<-all.vars(formula(f))
 if(!('Treatment'%in%v)||any(c('Time.Point','Species')%in%v))next
 for(method in c('satterthwaite','kenward-roger')) {
  q<-as.data.frame(pairs(emmeans(f,~Treatment,lmer.df=method),reverse=TRUE,adjust='none'))
  k<-k+1;results[[k]]<-data.frame(id=id,method=method,q)
 }
}
write.csv(do.call(rbind,results),file.path(out,'degrees_of_freedom_comparison.csv'),row.names=FALSE)
f<-o$fits[['L486_Biomass_Pvar_model.T1']]
if(is.null(f)) {
 candidates<-o$fits[grepl('Biomass_Pvar_model.T1$',names(o$fits))]
 f<-candidates[[which(vapply(candidates,function(f)'Clade'%in%all.vars(formula(f)),logical(1)))]]
}
q<-pairs(emmeans(f,~Treatment|Clade,lmer.df='satterthwaite'),reverse=TRUE)
write.csv(as.data.frame(summary(q,by=NULL,adjust='holm')),file.path(out,'pvar_T1_within_clade_treatment_Holm.csv'),row.names=FALSE)
e<-o$environment
for(n in c('meta','AFDW','Slurry','SA','biomass2'))cat(n,'N181 present:', 'N181'%in%e[[n]]$Frag.ID,'\n')
