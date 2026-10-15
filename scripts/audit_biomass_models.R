# Run from the repository root. Audit the notebook's existing biomass models
# without changing model specifications or overwriting manuscript figures.
suppressPackageStartupMessages({library(lmerTest);library(emmeans)})
out <- 'output/reviews/biomass_model_audit'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readLines('ACE21_physiology.Rmd'); y <- rep('',length(x)); active <- FALSE
for(i in seq_along(x)) {
 if(startsWith(x[i],'```{r')){active <- !grepl('eval=FALSE',x[i],fixed=TRUE);next}
 if(startsWith(x[i],'```')){active <- FALSE;next}
 if(active)y[i] <- x[i]
}
z <- parse(text=y,keep.source=TRUE); e<-new.env(parent=globalenv())
fits<-list(); records<-list(); warnings<-list(); anovas<-list()
end<-grep('^PLOTS$',x)[1]
pdf(tempfile(fileext='.pdf'))
for(i in seq_along(z)) {
 line<-attr(z,'srcref')[[i]][1]; if(line>=end)break
 expr<-z[[i]]
 withCallingHandlers(invisible(eval(expr,e)),warning=function(w){warnings[[length(warnings)+1]]<<-data.frame(line=line,message=conditionMessage(w));invokeRestart('muffleWarning')})
 if(is.call(expr) && is.symbol(expr[[1]]) && as.character(expr[[1]]) %in% c('<-','=') && is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name('lmer'))) {
  name<-as.character(expr[[2]]); fit<-e[[name]]; id<-paste0('L',line,'_',name);fits[[id]]<-fit
  vc<-as.data.frame(VarCorr(fit)); drops<-attr(getME(fit,'X'),'col.dropped')
  records[[id]]<-data.frame(id=id,line=line,name=name,formula=paste(deparse(formula(fit)),collapse=' '),n=nobs(fit),singular=isSingular(fit),dropped_columns=length(drops),parent_variance=if('Parent.ID'%in%vc$grp)vc$vcov[vc$grp=='Parent.ID']else NA_real_,tank_variance=if('Tank.num'%in%vc$grp)vc$vcov[vc$grp=='Tank.num']else NA_real_,residual_variance=sigma(fit)^2,convergence=paste(fit@optinfo$conv$lme4$messages,collapse='; '))
  a<-as.data.frame(anova(fit));a$term<-rownames(a);a$id<-id;rownames(a)<-NULL;anovas[[id]]<-a
 }
}
dev.off()
write.csv(do.call(rbind,records),file.path(out,'model_inventory.csv'),row.names=FALSE)
write.csv(do.call(rbind,anovas),file.path(out,'anova_terms.csv'),row.names=FALSE)
if(length(warnings))write.csv(do.call(rbind,warnings),file.path(out,'warnings.csv'),row.names=FALSE)
b<-e$biomass2
write.csv(b,file.path(out,'biomass_input.csv'),row.names=FALSE)
counts<-aggregate(Frag.ID~Species+Time.Point+Treatment,b,length);names(counts)[4]<-'n'
write.csv(counts,file.path(out,'sample_counts.csv'),row.names=FALSE)
write.csv(as.data.frame(with(b,table(Time.Point,Treatment,Tank.num))),file.path(out,'tank_allocation.csv'),row.names=FALSE)
# Retain each treatment-only time-specific fit separately from later clade fits.
selected<-names(fits)[vapply(fits,function(f){v<-all.vars(formula(f)); 'Treatment'%in%v && !any(c('Clade','Time.Point','Species')%in%v)},logical(1))]
contrasts<-list(); diagnostics<-list()
png(file.path(out,'treatment_model_diagnostics.png'),width=1800,height=2400,res=150)
par(mfrow=c(8,2),mar=c(3,4,2,1))
for(id in selected) {
 f<-fits[[id]]; r<-residuals(f); fv<-fitted(f)
 plot(fv,r,xlab='Fitted biomass',ylab='Residual',main=id,cex.main=.65);abline(h=0,lty=2)
 qqnorm(r,main='Residual Q-Q',cex.main=.8);qqline(r)
 q<-as.data.frame(pairs(emmeans(f,~Treatment),reverse=TRUE,adjust='none'));q$id<-id;contrasts[[id]]<-q
 diagnostics[[id]]<-data.frame(id=id,residual_min=min(r),residual_max=max(r),max_abs_scaled_residual=max(abs(r/sigma(f))),residual_sd_ambient=sd(r[model.frame(f)$Treatment=='Ambient']),residual_sd_high=sd(r[model.frame(f)$Treatment=='High']))
}
dev.off()
write.csv(do.call(rbind,contrasts),file.path(out,'treatment_contrasts_unadjusted.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'residual_summary.csv'),row.names=FALSE)
saveRDS(list(fits=fits,data=b,environment=as.list(e)),file.path(out,'objects.rds'))
cat('AUDITED',length(fits),'fits;',length(selected),'time-specific treatment-only models\n')
cat('Rows',nrow(b),'unique fragments',length(unique(b$Frag.ID)),'duplicate fragment IDs',sum(duplicated(b$Frag.ID)),'missing biomass',sum(is.na(b$AFDWmg.cm2)),'\n')
print(do.call(rbind,records)[,c('id','n','singular','dropped_columns','parent_variance','tank_variance')],row.names=FALSE)
