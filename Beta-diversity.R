# VST normalization
dds = phyloseq_to_deseq2(phyloseqobject, ~outcome)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
phyloseqobject.VST = phyloseqobject

# NMDS
otu<-as(otu_table(phyloseqobject.VST),"matrix")
otu<-as.data.frame(otu)
dat<-as(sample_data(phyloseqobject.VST),"data.frame")

nmds<- ordinate(phyloseqobject.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds, display = "sites"))
scrs <- cbind(scrs,Language = dat$outcome, Location=dat$Location) 

ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Language, shape=Location)) + scale_shape_manual(values = c(1, 16)) +
  coord_fixed() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("chocolate1","slateblue2")) +
  theme(text = element_text(size = 5)) +
  ggtitle("") 


# ANOSIM
set.seed(123)
RunAnosim<-function(x){
  sup = get_variable(x,"outcome")
  asv = as(otu_table(x),"matrix")
  ano<-anosim(as.data.frame(asv),grouping=sup,distance="euclidean",permutations=999)
  res<-data.frame(ano$signif,ano$statistic)
  return(res)
}

RunAnosim(phyloseqobject.VST)

# Adonis
library(vegan)
# remove missing values form multilingual variable
ps_clean <- subset_samples(p3.multi, !is.na(multilingual)) # change p2.multi, p3.multi, p4.multi
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

pseq.rel <- microbiome::transform(ps_clean, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

permanova <- adonis2(t(otu) ~  multilingual + outcome,
                     data = meta, permutations=999, method = "bray",by="margin")
permanova["Pr(>F)"]

#check homogeneity
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$outcome))