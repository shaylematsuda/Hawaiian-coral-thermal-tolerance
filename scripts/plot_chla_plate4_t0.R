suppressPackageStartupMessages({library(ggplot2)})
plate <- read.csv('Data/plate4.csv', check.names=FALSE)
plate$Chla_ug_ml <- rowMeans(plate[c('Rep1','Rep2')], na.rm=TRUE)
meta <- read.csv('Data/Master_Metadata_20190402.csv', check.names=FALSE)
d <- merge(plate, meta, by='Frag.ID')
d <- subset(d, Time.Point=='T0')
d$Species <- factor(gsub('_',' ',d$Species), levels=c('Montipora capitata','Porites compressa','Pavona varians','Pocillopora acuta'))
d$Treatment <- factor(d$Treatment, levels=c('Ambient','High'))
p <- ggplot(d, aes(Species, Chla_ug_ml, colour=Treatment)) +
  geom_boxplot(aes(group=interaction(Species,Treatment)), position=position_dodge(.65), width=.55, outlier.shape=NA, alpha=.2) +
  geom_jitter(position=position_jitterdodge(dodge.width=.65,jitter.width=.08), size=2.8) +
  labs(title='Plate 4 chlorophyll measurements at T0', x=NULL, y='Chlorophyll (µg mL⁻¹)', colour='Treatment') +
  theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=20,hjust=1))
dir.create('output/reviews/chla_t0_audit',recursive=TRUE,showWarnings=FALSE)
ggsave('output/reviews/chla_t0_audit/plate4_T0_chlorophyll.png',p,width=9,height=5,dpi=220)
write.csv(d,'output/reviews/chla_t0_audit/plate4_T0_chlorophyll_values.csv',row.names=FALSE)
