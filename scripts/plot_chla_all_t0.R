suppressPackageStartupMessages(library(ggplot2))
d <- read.csv('output/reviews/chla_t0_audit/chla_t0_input.csv')
p4 <- read.csv('Data/plate4.csv')$Frag.ID
d$Plate <- ifelse(d$Frag.ID %in% p4, 'Plate 4', 'Other plates')
d$Species <- factor(gsub('_',' ',d$Species), levels=c('Montipora capitata','Porites compressa','Pavona varians','Pocillopora acuta'))
d$Treatment <- factor(d$Treatment,levels=c('Ambient','High'))
p <- ggplot(d,aes(Species,ChlaSA,colour=Treatment,shape=Plate)) +
 geom_boxplot(aes(group=interaction(Species,Treatment)),position=position_dodge(.65),width=.55,outlier.shape=NA,alpha=.15) +
 geom_jitter(position=position_jitterdodge(dodge.width=.65,jitter.width=.08),size=2.7) +
 scale_shape_manual(values=c('Other plates'=16,'Plate 4'=17)) +
 labs(title='T0 chlorophyll by eventual treatment assignment',x=NULL,y='Chlorophyll per surface area',colour='Treatment',shape='Measurement source') +
 theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave('output/reviews/chla_t0_audit/chla_all_T0_plate4_triangles.png',p,width=9,height=5,dpi=220)
