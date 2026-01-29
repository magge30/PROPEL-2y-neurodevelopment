library(SparseMCMM)

# generalliz: phyloseqobject = your phyloseq object
pseq = subset_samples(phyloseqobject, Language != "NA")
pseq <- microbiome::transform(pseq,"compositional")

#define arguments  
data<-as(sample_data(pseq),"data.frame")
data$Supplementation<-as.numeric(data$L.reuteri)
data$Supplementation<-gsub(1,0,data$Supplementation)
data$Supplementation<-gsub(2,1,data$Supplementation)
Treatment<-as.numeric(data$Supplementation)
otu.com<-t(as(otu_table(pseq),"matrix"))
outcome<-sample_data(pseq)$Language #must be continuous

# run analysis
res1<-SparseMCMM(Treatment,otu.com,outcome,covariates=NULL,covariate.fix=NULL,num.per=100,bootstrap=10)



