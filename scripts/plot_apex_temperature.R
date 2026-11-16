# Approved mean ± SE temperature plot. Run from the project root.
# Figure 1B owner: this script refreshes derived inputs, applies the documented
# exclusions, prints the treatment means ± SE, and exports the panel for assembly.
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
source("scripts/figure1_workflow_helpers.R")
refresh_apex_records()
out<-'Data/APEX/derived'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(file.path(out,'agreeing_records.csv'))
d$Date<-as.POSIXct(d$Date,tz='Pacific/Honolulu')
# Exclusions confirmed by the investigator: these probes were removed for
# cleaning and returned to the tanks. They are logger-handling readings,
# not tank water temperatures. Exclude only these six exact tank/timestamp pairs.
cleaning <- data.frame(
 Tank=c(2,2,4,3,1,5),
 Date=as.POSIXct(c('2017-04-17 17:00:00','2017-04-17 17:15:00',
                  '2017-04-17 16:45:00','2017-04-18 11:30:00',
                  '2017-04-18 11:15:00','2017-04-19 18:45:00'),tz='Pacific/Honolulu'))
excluded_cleaning <- inner_join(d,cleaning,by=c('Tank','Date'))
stopifnot(nrow(excluded_cleaning)==6)
write.csv(excluded_cleaning,file.path(out,'sampling_window_cleaning_exclusions.csv'),row.names=FALSE)
d<-anti_join(d,cleaning,by=c('Tank','Date'))
start<-as.POSIXct('2017-03-20 00:00:00',tz='Pacific/Honolulu')
end<-as.POSIXct('2017-05-18 00:00:00',tz='Pacific/Honolulu')
# The investigator confirmed May 1-15 data for tanks 1-4 are not real tank
# temperatures (logger removal). Drop that interval from the plot only.
# Leave raw measurements untouched and save all excluded source rows for audit.
may_removed <- d %>% filter(Tank %in% 1:4,
 Date>=as.POSIXct('2017-05-01',tz='Pacific/Honolulu'),
 Date<as.POSIXct('2017-05-16',tz='Pacific/Honolulu'))
write.csv(may_removed,file.path(out,'logger_removal_May1_15_tanks1_4.csv'),row.names=FALSE)
# Bounds follow the first/last dated PAM sampling files; end is exclusive so
# May 17 is included in full. No records are extrapolated through missing periods.
d<-d %>% filter(Date>=start,Date<end,
 !(Tank %in% 1:4 & Date>=as.POSIXct('2017-05-01',tz='Pacific/Honolulu') & Date<as.POSIXct('2017-05-16',tz='Pacific/Honolulu'))) %>%
 mutate(Treatment=if_else(Tank %in% c(1,4,5,8,9),'Ambient','Elevated')) %>%
 arrange(Tank,Date) %>% group_by(Tank) %>%
 mutate(Segment=cumsum(c(TRUE,diff(as.numeric(Date))>900))) %>% ungroup()
# Independent replicates are tanks, not 15-minute measurements. Compute a
# separate treatment mean and SE across available tanks at every timestamp.
# N varies after exclusions; SE is undefined with only one tank and is not drawn.
means <- d %>% group_by(Date,Treatment) %>%
 summarise(Mean=mean(Temperature),N=n_distinct(Tank),
           SE=if(n_distinct(Tank)>1) sd(Temperature)/sqrt(n_distinct(Tank)) else NA_real_,
           .groups='drop') %>% arrange(Treatment,Date) %>% group_by(Treatment) %>%
 mutate(Segment=cumsum(c(TRUE,diff(as.numeric(Date))>900))) %>% ungroup()
# Elapsed days are a display coordinate only: day 0 is the first sampling date.
# Dates remain available in the saved data and all averaging/exclusions are unchanged.
means$Day <- as.numeric(difftime(means$Date,start,units="days"))
plot_end_day <- as.numeric(difftime(end,start,units="days"))
write.csv(means,file.path(out,'APEX_sampling_window_treatment_mean_SE.csv'),row.names=FALSE)
ribbons <- means %>% filter(!is.na(SE)) %>% group_by(Treatment) %>%
 mutate(Segment=cumsum(c(TRUE,diff(as.numeric(Date))>900))) %>% ungroup()
apex_temperature_plot<-ggplot(means,aes(Day,Mean,color=Treatment,fill=Treatment,
                         group=interaction(Treatment,Segment)))+
 geom_ribbon(data=ribbons,aes(ymin=Mean-SE,ymax=Mean+SE),alpha=.18,color=NA,na.rm=TRUE)+
 geom_line(linewidth=.65)+
 scale_color_manual(values=c(Ambient='#2166AC',Elevated='#D73027'),name=NULL)+
 scale_fill_manual(values=c(Ambient='#2166AC',Elevated='#D73027'),name=NULL)+
 scale_x_continuous(limits=c(0,plot_end_day),breaks=c(0,10,20,30,40,50,58),
                    expand=expansion(mult=c(0,0)))+
 labs(title='Mean tank temperature by treatment',subtitle='Mean ± 1 SE across available tanks',
      x='Days since first sampling',y=expression(Temperature~(degree*C)),
      caption='Means calculated across tanks at each timestamp; no averaging across time. Available tank counts vary.\nMay 1-15 tank 1-4 records and six invalid readings on April 17-19 omitted. No SE ribbon where fewer than two tanks remain.')+
 theme_classic(base_size=12)+theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=9),plot.title=element_text(face='bold'),axis.text.x=element_text(size=10))

print(apex_temperature_plot)
ggsave("Figures/Fig1B_temperature_mean_SE.pdf", apex_temperature_plot,
       width=11,height=5.5,bg="white")
ggsave("Figures/Fig1B_temperature_mean_SE.png", apex_temperature_plot,
       width=11,height=5.5,dpi=300,bg="white")
saveRDS(apex_temperature_plot,"Figures/Fig1B_temperature_mean_SE_plot.rds")

# Export a compact version for import by the Figure 1 assembly. Identical data,
# means, SE, and exclusions; only typography and non-panel captions differ.
Fig1B <- apex_temperature_plot +
 labs(title=NULL,subtitle=NULL,caption=NULL) +
 scale_x_continuous(limits=c(0,plot_end_day),breaks=c(0,20,40,58),
                    expand=expansion(mult=c(0,0))) +
 theme(text=element_text(size=10),axis.text.x=element_text(size=8,hjust=1),
       axis.title=element_text(size=10),legend.text=element_text(size=9),
       plot.margin=margin(5,5,5,5))
dir.create("Figures/assets/Fig1",recursive=TRUE,showWarnings=FALSE)
ggsave("Figures/assets/Fig1/Fig1B_temperature.pdf",Fig1B,width=5,height=3.75,bg="white")
saveRDS(Fig1B,"Figures/assets/Fig1/Fig1B_temperature_plot.rds")
