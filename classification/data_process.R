#Hek293t
hek293t=read.csv('hek293t.epiotrt',sep='\t',header = FALSE)
hek293t=hek293t[,c(-3,-4,-5,-6,-8,-9,-10,-11)]
colnames(hek293t)=c('sgRNAid','sgRNA','sequence','label')
hek293t$cell='hek293t'
hek293t=hek293t[order(hek293t$sgRNAid),]
write.csv(hek293t,'hek293t_classification.csv')
#Hek293t_sgRNA
hek293t_sgRNA=hek293t[,c(1,2)]
hek293t_sgRNA=unique(hek293t_sgRNA)
write.csv(hek293t_sgRNA,'hek293t_sgRNA.csv')

#K562
k562=read.csv('k562.epiotrt',sep='\t',header = FALSE)
k562=k562[,c(-3,-4,-5,-6,-8,-9,-10,-11)]
colnames(k562)=c('sgRNAid','sgRNA','sequence','label')
k562$cell='k562'
k562=k562[order(k562$sgRNAid),]
write.csv(data,'off_classification.csv')
#K562_sgRNA
k562_sgRNA=k562[,c(1,2)]
k562_sgRNA=unique(k562_sgRNA)
write.csv(k562_sgRNA,'k562_sgRNA.csv')
#Hek293t_off_position_information
hek293t_off_position_information=read.csv('sequence_matches.csv',sep=',',header = TRUE)
hek293t_off_position_information=hek293t_off_position_information[hek293t_off_position_information$unique=='yes',]
write.csv(hek293t_off_position_information,'Unique_hek293t_off_position_information.csv')
