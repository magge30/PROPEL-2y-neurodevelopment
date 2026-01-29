#++++++++++++++++++++++#
# Function for outcome dichotomisation
#++++++++++++++++++++++#
# ps= phyloseq object
run2v2.Language<-function(ps){
  temp = subset_samples(ps, Language_2v2 == "Normal" | Language_2v2 == "Impairment")
  temp<-prune_taxa(taxa_sums(temp) > 0, temp)
  temp <- temp %>%
    rename_sample_data(outcome = Language_2v2)
  return(temp)
}
run2v2.NDI<-function(ps){
  temp = subset_samples(ps, NDI_2v2 == "Normal" | NDI_2v2 == "Impairment")
  temp<-prune_taxa(taxa_sums(temp) > 0, temp)
  temp <- temp %>%
    rename_sample_data(outcome = NDI_2v2)
  return(temp)
}
run2v2.Motor<-function(ps){
  temp = subset_samples(ps, Motor_2v2 == "Normal" | Motor_2v2 == "Impairment")
  temp<-prune_taxa(taxa_sums(temp) > 0, temp)
  temp <- temp %>%
    rename_sample_data(outcome = Motor_2v2)
  return(temp)
}
run2v2.Cognition<-function(ps){
  temp = subset_samples(ps, Cognition_2v2 == "Normal" | Cognition_2v2 == "Impairment")
  temp<-prune_taxa(taxa_sums(temp) > 0, temp)
  temp <- temp %>%
    rename_sample_data(outcome = Cognition_2v2)
  return(temp)
}

# ex: 
p1.lan<-run2v2.Language(p1)

#++++++++++++++++++++++#
# function for tax_glom
#++++++++++++++++++++++#
run_taxa_levels<-function(phyloseqobject,taxalevel){
  phy.glom<-tax_glom(phyloseqobject,taxrank=taxalevel,NArm=TRUE)
  taxa_names(phy.glom)<-tax_table(phy.glom)[,taxalevel]
  return(phy.glom)
}

#ex:
p1.gen<-run_taxa_levels(p1,"Genus")
