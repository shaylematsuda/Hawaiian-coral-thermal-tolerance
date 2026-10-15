suppressPackageStartupMessages({library(lmerTest); library(emmeans)})
out <- 'output/reviews/t0_trait_audit'; dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- readLines('ACE21_physiology.Rmd'); z <- rep('', length(x)); active <- FALSE
for(i in seq_along(x)) { if(startsWith(x[i],'```{r')) { active <- !grepl('eval=FALSE',x[i],fixed=TRUE); next }; if(startsWith(x[i],'```')) { active <- FALSE; next }; if(active) z[i] <- x[i] }
end <- grep('^PLOTS$', x)[1]; z[(end+1):length(z)] <- ''
e <- new.env(parent=globalenv()); fits <- list(); rows <- list(); exprs <- parse(text=z, keep.source=TRUE)
for(i in seq_along(exprs)) {
 expr <- exprs[[i]]
 tryCatch({ eval(expr,e) }, error=function(err) NULL)
 if(is.call(expr) && is.symbol(expr[[1]]) && as.character(expr[[1]]) %in% c('<-','=') && is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name('lmer'))) {
  nm <- as.character(expr[[2]]); if(grepl('T0_model',nm,fixed=TRUE) && exists(nm,e)) { f<-e[[nm]]; a<-as.data.frame(anova(f)); a$term<-rownames(a); a$id<-nm; rows[[length(rows)+1]]<-a; if(grepl('PAM.*T0_model',nm)) saveRDS(f,file.path(out,'PAM_T0_model.rds')) }
 }
}
if(length(rows)) write.csv(do.call(rbind,rows),file.path(out,'t0_anova.csv'),row.names=FALSE)
cat('Captured',length(rows),'T0 model result tables\n')
