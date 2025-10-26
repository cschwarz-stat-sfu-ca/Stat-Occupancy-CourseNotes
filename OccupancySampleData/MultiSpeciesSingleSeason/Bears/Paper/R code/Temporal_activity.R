##calculate activity overlaps

library(overlap)

setwd("C:/")

temporal_events <- read.csv("Temporal_events.csv")

##create radians time column

temporal_events$time_rad<-temporal_events$time.decimal*2*pi

##split into difference species dfs

gbear<-temporal_events[temporal_events$Species=="Grizzly Bear",]
bbear<-temporal_events[temporal_events$Species=="Black Bear",]
motorised<-temporal_events[temporal_events$Species=="motorised",]
non_motorised<-temporal_events[temporal_events$Species=="non-motorised",]
all_rec<-temporal_events[temporal_events$Species=="motorised"|temporal_events$Species=="non-motorised",]

##bootstrap resamples

gbear_bstrap<-resample(gbear[,8],10000)
bbear_bstrap<-resample(bbear[,8],10000)
rec_bstrap<-resample(all_rec[,8],10000)

##calculate coefficients of overlap

gbear_bbear_est<-overlapEst(gbear[,8],bbear[,8])
overlapPlot(gbear[,8],bbear[,8],olapcol = 'light blue',rug=T,main="")
legend('topleft',c("Grizzly","Black"),lty=c(1,2),col=c(1,4),bty='n')

gbear_bbear_bstrap<-bootEst(gbear_bstrap,bbear_bstrap)
bootCI(gbear_bbear_est[1],gbear_bbear_bstrap[,1],conf=0.95)

gbear_motor_est<-overlapEst(gbear[,8],motorised[,8])
overlapPlot(gbear[,8],motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8)
legend('topleft',c("Grizzly","Motorised"),lty=c(1,2),col=c(1,4),bty='n')

bbear_motor_est<-overlapEst(bbear[,8],motorised[,8],adjust=c(NA,0.8,NA))
overlapPlot(bbear[,8],motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8)
legend('topleft',c("Black","Motorised"),lty=c(1,2),col=c(1,4),bty='n')

gbear_nmotor_est<-overlapEst(gbear[,8],non_motorised[,8],adjust=c(NA,0.8,NA))
overlapPlot(gbear[,8],non_motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8)
legend('topleft',c("Grizzly","Non-motorised"),lty=c(1,2),col=c(1,4),bty='n')

bbear_nmotor_est<-overlapEst(bbear[,8],non_motorised[,8],adjust=c(NA,0.8,NA))
overlapPlot(bbear[,8],non_motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8)
legend('topleft',c("Grizzly","Non-motorised"),lty=c(1,2),col=c(1,4),bty='n')

gbear_rec_est<-overlapEst(gbear[,8],all_rec[,8])
bbear_rec_est<-overlapEst(bbear[,8],all_rec[,8])

gbear_rec_bstrap<-bootEst(gbear_bstrap,rec_bstrap)
bootCI(gbear_rec_est[1],gbear_rec_bstrap[,1],conf=0.95)

bbear_rec_bstrap<-bootEst(bbear_bstrap,rec_bstrap)
bootCI(bbear_rec_est[1],bbear_rec_bstrap[,1],conf=0.95)

overlapPlot(bbear[,8],all_rec[,8],olapcol = 'light blue',rug=T,main="")
legend('topleft',c("Black","Recreation"),lty=c(1,2),col=c(1,4),bty='n')

overlapPlot(gbear[,8],all_rec[,8],olapcol = 'light blue',rug=T,main="")
legend('topleft',c("Grizzly","Recreation"),lty=c(1,2),col=c(1,4),bty='n')


##split datasets by locations where activity is present and absent

#motorised

gbear_motor<-gbear[gbear$Location%in%motorised$Location,]
gbear_nomotor<-gbear[!gbear$Location%in%motorised$Location,]

overlapPlot(gbear_motor[,8],gbear_nomotor[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")

bbear_motor<-bbear[bbear$Location%in%motorised$Location,]
bbear_nomotor<-bbear[!bbear$Location%in%motorised$Location,]

overlapPlot(bbear_motor[,8],bbear_nomotor[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")

#non_motorised

gbear_non_motor<-gbear[gbear$Location%in%non_motorised$Location,]
gbear_nonon_motor<-gbear[!gbear$Location%in%non_motorised$Location,]

overlapPlot(gbear_non_motor[,8],gbear_nonon_motor[,8])

bbear_non_motor<-bbear[bbear$Location%in%non_motorised$Location,]
bbear_nonon_motor<-bbear[!bbear$Location%in%non_motorised$Location,]

overlapPlot(bbear_non_motor[,8],bbear_nonon_motor[,8])

#both forms of rec

gbear_norec<-gbear[!gbear$Location%in%non_motorised$Location&!gbear$Location%in%motorised$Location,]
gbear_rec<-gbear[gbear$Location%in%non_motorised$Location&gbear$Location%in%motorised$Location,]
# gbear_rec2<-gbear[gbear$Location%in%non_motorised$Location|gbear$Location%in%motorised$Location,]

overlapPlot(gbear_norec[,8],all_rec[,8],rug=T)
gbear_norec_est<-overlapEst(gbear_norec[,8],all_rec[,8])
overlapPlot(gbear_rec[,8],all_rec[,8],rug=T)
gbear_rec_est<-overlapEst(gbear_rec[,8],all_rec[,8])

gbear_bstrap_rec<-resample(gbear_rec[,8],10000)
gbear_bstrap_norec<-resample(gbear_norec[,8],10000)

gbear_rec_bstrap<-bootEst(gbear_bstrap_rec,rec_bstrap)
gbear_norec_bstrap<-bootEst(gbear_bstrap_norec,rec_bstrap)

bootCI(gbear_rec_est[1],gbear_rec_bstrap[,1],conf=0.95)

bootCI(gbear_norec_est[1],gbear_norec_bstrap[,1],conf=0.95)

#both forms bbear

bbear_norec<-bbear[!bbear$Location%in%non_motorised$Location&!bbear$Location%in%motorised$Location,]
bbear_rec<-bbear[bbear$Location%in%non_motorised$Location&bbear$Location%in%motorised$Location,]

overlapPlot(bbear_norec[,8],all_rec[,8],rug=T)
bbear_norec_est<-overlapEst(bbear_norec[,8],all_rec[,8],adjust=c(0.8,NA,NA))
overlapPlot(bbear_rec[,8],all_rec[,8],rug=T)
bbear_rec_est<-overlapEst(bbear_rec[,8],all_rec[,8],adjust=c(0.8,NA,NA))

bbear_bstrap_rec<-resample(bbear_rec[,8],10000)
bbear_bstrap_norec<-resample(bbear_norec[,8],10000)

bbear_rec_bstrap<-bootEst(bbear_bstrap_rec,rec_bstrap)
bbear_norec_bstrap<-bootEst(bbear_bstrap_norec,rec_bstrap)

bootCI(bbear_rec_est[1],bbear_rec_bstrap[,1],conf=0.95)

bootCI(bbear_norec_est[1],bbear_norec_bstrap[,1],conf=0.95)


##bears

gbear_bbear<-gbear[gbear$Location%in%bbear$Location,]
gbear_nobbear<-gbear[!gbear$Location%in%bbear$Location,]

bbear_gbear<-bbear[bbear$Location%in%gbear$Location,]
bbear_nogbear<-bbear[!bbear$Location%in%gbear$Location,]

overlapPlot(gbear_bbear[,8],bbear[,8])
overlapPlot(gbear_nobbear[,8],bbear[,8])

overlapPlot(bbear_gbear[,8],gbear[,8])
overlapPlot(bbear_nogbear[,8],gbear[,8])

gbear_bb<-overlapEst(gbear_bbear[,8],bbear[,8])
gbear_nobb<-overlapEst(gbear_nobbear[,8],bbear[,8])

bbear_gb<-overlapEst(bbear_gbear[,8],gbear[,8])
bbear_nogb<-overlapEst(bbear_nogbear[,8],gbear[,8])

gbear_bstrap_bb<-resample(gbear_bbear[,8],10000)
gbear_bstrap_nobb<-resample(gbear_nobbear[,8],10000)

bbear_bstrap_gb<-resample(bbear_gbear[,8],10000)
bbear_bstrap_nogb<-resample(bbear_nogbear[,8],10000)

gbear_bstrap<-resample(gbear[,8],10000)
bbear_bstrap<-resample(bbear[,8],10000)

gbear_bb_bstrap<-bootEst(gbear_bstrap_bb,bbear_bstrap)
gbear_nobb_bstrap<-bootEst(gbear_bstrap_nobb,bbear_bstrap)

bbear_gb_bstrap<-bootEst(bbear_bstrap_gb,gbear_bstrap)
bbear_nogb_bstrap<-bootEst(bbear_bstrap_nogb,gbear_bstrap)

bootCI(gbear_bb[1],gbear_bb_bstrap[,1],conf=0.95)
bootCI(gbear_nobb[1],gbear_nobb_bstrap[,1],conf=0.95)

bootCI(bbear_gb[1],bbear_gb_bstrap[,1],conf=0.95)
bootCI(bbear_nogb[1],bbear_nogb_bstrap[,1],conf=0.95)

#motorised

gbear_motor_present<-overlapEst(gbear_motor[,8],motorised[,8])
overlapPlot(gbear_motor[,8],motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Grizzly","Motorised"),lty=c(1,2),col=c(1,4),bty='n')

gbear_motor_absent<-overlapEst(gbear_nomotor[,8],motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(gbear_nomotor[,8],motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Grizzly","Motorised"),lty=c(1,2),col=c(1,4),bty='n')

bbear_motor_present<-overlapEst(bbear_motor[,8],motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(bbear_motor[,8],motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Black","Motorised"),lty=c(1,2),col=c(1,4),bty='n')

bbear_motor_absent<-overlapEst(bbear_nomotor[,8],motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(bbear_nomotor[,8],motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Black","Motorised"),lty=c(1,2),col=c(1,4),bty='n')

#non-motorised

gbear_non_motor_present<-overlapEst(gbear_non_motor[,8],non_motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(gbear_non_motor[,8],non_motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Grizzly","Non_motorised"),lty=c(1,2),col=c(1,4),bty='n')

gbear_non_motor_absent<-overlapEst(gbear_nonon_motor[,8],non_motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(gbear_nonon_motor[,8],non_motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Grizzly","Non_motorised"),lty=c(1,2),col=c(1,4),bty='n')

bbear_non_motor_present<-overlapEst(bbear_non_motor[,8],non_motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(bbear_non_motor[,8],non_motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Black","Non_motorised"),lty=c(1,2),col=c(1,4),bty='n')

bbear_non_motor_absent<-overlapEst(bbear_nonon_motor[,8],non_motorised[,8],adjust=c(0.8,NA,NA))
overlapPlot(bbear_nonon_motor[,8],non_motorised[,8],olapcol = 'light blue',rug=T,adjust=0.8,main="")
legend('topleft',c("Black","Non_motorised"),lty=c(1,2),col=c(1,4),bty='n')

##bootstrap for confidence intervals

##gbear vs bbear

gbear_bstrap<-resample(gbear[,8],10000)
bbear_bstrap<-resample(bbear[,8],10000)

gbear_bbear_bstrap<-bootEst(gbear_bstrap,bbear_bstrap,adjust=c(0.8,NA,NA))
bootCI(gbear_bbear_est[1],gbear_bbear_bstrap[,1],conf=0.95)

##motorised with and without

#with

gbear_bstrap_motor<-resample(gbear_motor[,8],10000)
bbear_bstrap_motor<-resample(bbear_motor[,8],10000)
motorised_bstrap<-resample(motorised[,8],10000)

boot_gbear_motor<-bootEst(gbear_bstrap_motor,motorised_bstrap,n.grid=512)
bootCI(gbear_motor_present[1],boot_gbear_motor[,1],conf=0.95)

boot_bbear_motor<-bootEst(bbear_bstrap_motor,motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(bbear_motor_present[1],boot_bbear_motor[,1],conf=0.95)

#without

gbear_bstrap_nomotor<-resample(gbear_nomotor[,8],10000)
bbear_bstrap_nomotor<-resample(bbear_nomotor[,8],10000)
# motorised_bstrap<-resample(motorised[,8],10000)

boot_gbear_nomotor<-bootEst(gbear_bstrap_nomotor,motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(gbear_motor_absent[1],boot_gbear_nomotor[,1],conf=0.95)

boot_bbear_nomotor<-bootEst(bbear_bstrap_nomotor,motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(bbear_motor_absent[1],boot_bbear_nomotor[,1],conf=0.95)

##non_motorised with and without

#with

gbear_bstrap_non_motor<-resample(gbear_non_motor[,8],10000)
bbear_bstrap_non_motor<-resample(bbear_non_motor[,8],10000)
non_motorised_bstrap<-resample(non_motorised[,8],10000)

boot_gbear_non_motor<-bootEst(gbear_bstrap_non_motor,non_motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(0.4427099,boot_gbear_non_motor[,1],conf=0.95)

boot_bbear_non_motor<-bootEst(bbear_bstrap_non_motor,non_motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(0.6381969,boot_bbear_non_motor[,1],conf=0.95)

#without

gbear_bstrap_nonon_motor<-resample(gbear_nonon_motor[,8],10000)
bbear_bstrap_nonon_motor<-resample(bbear_nonon_motor[,8],10000)
# non_motorised_bstrap<-resample(non_motorised[,8],10000)

boot_gbear_nonon_motor<-bootEst(gbear_bstrap_nonon_motor,non_motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(0.4944207,boot_gbear_nonon_motor[,1],conf=0.95)

boot_bbear_nonon_motor<-bootEst(bbear_bstrap_nonon_motor,non_motorised_bstrap,adjust=c(0.8,NA,NA))
bootCI(0.5923677,boot_bbear_nonon_motor[,1],conf=0.95)



