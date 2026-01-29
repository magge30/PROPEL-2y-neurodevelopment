library(coda4microbiome)

# Taxnonomical level
ps.gen<-run_taxa_levels(ps,"Genus")

# prepare data timepoint to days
sample_data(ps.gen)$Timepoint_days<-gsub("1v",7,sample_data(ps.gen)$Timepoint)
sample_data(ps.gen)$Timepoint_days<-gsub("2v",14,sample_data(ps.gen)$Timepoint_days)
sample_data(ps.gen)$Timepoint_days<-gsub("3v",21,sample_data(ps.gen)$Timepoint_days)
sample_data(ps.gen)$Timepoint_days<-gsub("4v",28,sample_data(ps.gen)$Timepoint_days)
sample_data(ps.gen)$Subject<-gsub("\\-.*","",sample_data(ps.gen)$SampleID)

# Example with language
# subset by level
v2<-run2v2.Language(ps.gen)

# function summary for filter samples with at elast 3 timepoints
runFilterData<-function(phyObject,obs,n){
  # Define arguments
  x=t(as(otu_table(phyObject),"matrix")) 
  x_time = as.numeric(sample_data(phyObject)$Timepoint_days)
  subject_id=sample_data(phyObject)$Subject
  ini_time=7
  end_time=28
  metadata=sample_data(phyObject)
  temp<-as(tax_table(phyObject),"matrix")
  temp<-as.data.frame(temp)
  taxanames=temp$Genus
  
  # run function filter_longitudinal min observations: 3 of 4 timepoints
  data_filtered<- filter_longitudinal(
    x,
    taxanames,
    x_time,
    subject_id,
    metadata,
    ini_time,
    end_time,
    percent_indv = n, #percentage of TAXA with more than min_obs observations
    min_obs = obs
  )
  
  # extract data
  genPre_filt.asv<-data_filtered$`filtered abundance matrix`
  genPre_filt.tax<-data_filtered$`filtered taxa names`
  genPre_filt.dat<-data_filtered$`filtered metadata`
  return(list(asv=genPre_filt.asv,taxanames=genPre_filt.tax,metadata=genPre_filt.dat))
} 

# filter function
ps.filtered<-runFilterData(v2,3,0.1) 

# Generalize arguments
asv=ps.filtered$asv
metadata=ps.filtered$metadata
taxanames=ps.filtered$taxanames

# number of samples for analysis
length(unique(metadata$Subject))
metadata_df<-as(metadata,"data.frame")
metadata_df %>%
  group_by(outcome) %>%
  summarise(unique_subject_chars = n_distinct(Subject))

# Define Arguments
x=asv[, which(colSums(asv) != 0)] # microbiome abundance
taxanames<-as.data.frame(x) %>% names()
x_time = as.numeric(metadata$Timepoint_days)   # observation times
subject_id = as.numeric(as.factor(metadata$Subject))     # subject id
y=as.factor(metadata$outcome) # CHANGE: variable 
ini_time = 7
end_time = 28

# Run function
set.seed(1234)
res_long.min <-coda_glmnet_longitudinal(x,y, x_time, subject_id, ini_time, end_time, lambda="lambda.min",nfolds=6,showPlots=TRUE)


#plots
res_long.min$`signature plot`
res_long.min$`trajectories plot`
res_long.min$`predictions plot`

res_long.min$`apparent AUC`
res_long.min$`mean cv-AUC`
res_long.min$`sd cv-AUC`

# Extract data
res_long.min$taxa.num
res_long.min$taxa.name
res_long.min$`log-contrast coefficients`

# taxa in group1
coef<-res_long.min$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
taxanames[positives[op]]

# taxa in group2
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
taxanames[negatives[on]]


## MODIFIED PLOTS
plot_signature_modified<-function (vars, coeff, showPlots = TRUE, varnames = NULL) 
{
  Language <- ifelse(coeff > 0, 1, 0)
  Language <- factor(Language, levels = c(0, 1), labels = c("Impairment", 
                                                            "Normal"))
  df <- data.frame(vars, coefficient = round(coeff, digits = 2), 
                   Language)
  if (!is.null(varnames)) {
    df$vars <- varnames
  }
  L <- ggpubr::ggbarplot(df, x = "vars", y = "coefficient", color = "Language",  palette = c("chocolate1","slateblue2"), 
                         fill = "Language", sort.val = "asc", orientation = "horiz", 
                         position = ggplot2::position_dodge(), label = TRUE, lab.vjust = 0.2, 
                         lab.hjust = 0.5, ylab = FALSE)
  if (showPlots == TRUE) {
    print(L)
  }
  return(L)
}

plot_signature_modified(vars=res_long.min$taxa.name, coeff=res_long.min$`log-contrast coefficients`, showPlots = TRUE)

# save plots for MS
p_signature_language<-plot_signature_modified(vars=res_long.min$taxa.name, coeff=res_long.min$`log-contrast coefficients`, showPlots = TRUE)
p_trajectories_language<-res_long.min$`trajectories plot`
p_predictions_language<-res_long.min$`predictions plot`