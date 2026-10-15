suppressPackageStartupMessages(library(lmerTest))
out <- 'output/reviews/gross_photosynthesis_t0_audit'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta <- read.csv('Data/Master_Metadata_20190402.csv')
metaPR <- meta[c('Frag.ID','Parent.ID','Treatment','Species','Tank.num')]
pr <- subset(read.csv('Data/respirometry_processed_Final.csv'), Time.Point=='T0')
d <- merge(pr,metaPR,all.y=TRUE); d <- d[!is.na(d$PR),]
drop <- list(c('N94',2),c('N211',4),c('N209',2),c('N197',4),c('N181',8),c('N104',6),c('422',4),c('3442',4),c('284',1),c('284',11),c('1556',2),c('1279',1),c('1226',8),c('1144',4),c('N92',4))
for(z in drop)d<-d[!(d$Frag.ID==z[1]&d$round==as.numeric(z[2])),]
d<-d[!(d$Frag.ID=='1736'&d$channel==8),]
d$Photo<-d$Photo*60*1000;d$Resp<-d$Resp*60*1000;d$GrossP<-d$Photo-d$Resp
d$PR_real <- d$GrossP / (-d$Resp)
fit<-lmer(PR_real~Treatment*Species+(1|Tank.num),d)
write.csv(as.data.frame(anova(fit)),file.path(out,'pr_ratio_t0_anova.csv'))
saveRDS(fit,file.path(out,'pr_ratio_t0_model.rds'))
print(anova(fit))
