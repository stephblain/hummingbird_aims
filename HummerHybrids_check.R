#quick check to make sure that the ref panel's AIMs make sense

library(tidyverse); library(data.table)

setwd("/scratch/user/sblain/hummer_aims/run_hmm_autosomes/")

hybs<-fread("HummerHybrids_posterior_states",sep="\t",header=T)
srrIDs<-colnames(hybs)[startsWith(colnames(hybs),"SRR")]
hybsRef<-hybs%>%select(1:2,all_of(srrIDs))

SAMPLES<-c("rubythroated_SRR12247318", "blackchinned_SRR12247328", "blackchinned_SRR12247340", "blackchinned_SRR12247342", "blackchinned_SRR12247344", "blackchinned_SRR12247345", "blackchinned_SRR12247346", "blackchinned_SRR12247347", "blackchinned_SRR12247348", "rubythroated_SRR12247319", "rubythroated_SRR12247320", "rubythroated_SRR12247322", "rubythroated_SRR12247325", "rubythroated_SRR12247326", "rubythroated_SRR12247327", "rubythroated_SRR12247324", "blackchinned_SRR12247329", "rubythroated_SRR12247323")
spIDs<-data.frame(species=str_split_i(SAMPLES,"_",1),
           reference=str_split_i(SAMPLES,"_",2))%>%
  filter(reference%in%srrIDs)

hybsRef<-hybsRef%>%
  mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%select(-CHROM,-POS)%>%
  pivot_longer(cols=c(-ChromPos),names_to = "reference",values_to="genotype")%>%
  left_join(spIDs)

highNAs<-hybsRef%>%group_by(ChromPos)%>%
  summarise(NAs=sum(is.na(genotype))/18)%>%
  filter(NAs>0.25)%>%pull(ChromPos)

hybsRef<-hybsRef%>%filter(!ChromPos%in%highNAs)

hybsRef%>%group_by(species)%>%
  summarise(species_mean=mean(genotype,na.rm=T))

parhmm<-fread("../ref_panel/HummerRefPanel.f50.parhmm")
colnames(parhmm)<-c("CHROM","POS","refSp1","altSp1","refSp2","altSp2")
parhmm<-parhmm%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%select(-CHROM,-POS)%>%
  filter(ChromPos%in%unique(hybsRef$ChromPos))

parhmm<-parhmm%>%mutate(refSp=if_else(refSp1/(refSp1+altSp1)>refSp2/(refSp2+altSp2),"sp1","sp2"))

hybsRef<-hybsRef%>%left_join(parhmm%>%select(ChromPos,refSp))

hybsRef%>%group_by(species,refSp)%>%
  summarise(species_mean=mean(genotype,na.rm=T))
