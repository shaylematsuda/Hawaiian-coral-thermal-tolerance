suppressPackageStartupMessages(library(ggplot2))
d <- model.frame(readRDS('output/reviews/protein_t0_audit/protein_t0_model.rds'))
d$Species <- factor(gsub('_',' ',d$Species), levels=c('Montipora capitata','Porites compressa','Pavona varians','Pocillopora acuta'))
d$Treatment <- factor(d$Treatment,levels=c('Ambient','High'))
p <- ggplot(d,aes(Species,mg.mm2_sqrt,colour=Treatment)) +
 geom_boxplot(aes(group=interaction(Species,Treatment)),position=position_dodge(.65),width=.55,outlier.shape=NA,alpha=.15) +
 geom_jitter(position=position_jitterdodge(dodge.width=.65,jitter.width=.08),size=2.7) +
 labs(title='T0 protein by eventual treatment assignment',x=NULL,y='Protein (square-root transformed)',colour='Treatment') + theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave('output/reviews/protein_t0_audit/protein_all_T0.png',p,width=9,height=5,dpi=220)
