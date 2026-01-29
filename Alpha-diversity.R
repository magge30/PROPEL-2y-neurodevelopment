# Diversity Fucntion
runAlphaDiversity<-function(phyloseqObject){
  tse<-convertFromPhyloseq(phyloseqObject)
  colData(tse)$outcome = factor(colData(tse)$outcome, levels = c("Normal", "Impairment"))
  tse <- mia::estimateRichness(tse, 
                               assay_name = "counts", 
                               index = "observed", 
                               name="Richness")
  tse <- mia::estimateDiversity(tse, 
                                assay_name = "counts",
                                index = c("shannon"),
                                name = c("Diversity"))
  tse <- mia::estimateEvenness(tse, 
                               index = c("pielou"), 
                               name = c("Evenness"))
  
}

#+++++++++++++#
# glm()
tse<-runAlphaDiversity(phyloseqObject) 
df<-as.data.frame(colData(tse))
df$outcome_num <- recode(df$outcome, "Normal" = "0", "Impairment" = "1")

# Index: replace for Richness, Diversity and Evenness
# apply glm() model according to timepoing, outcome, and confounding factors
summary(glm(outcome_num ~ Index, data = df,family = binomial(link = "logit")))
summary(glm(outcome_num ~ Index + Location, data = df,family = binomial(link = "logit")))
summary(glmer(outcome_num ~ Index+ (1|Location), data = df,family = binomial()))

#+++++++++++++#
# normality
RunShapiro<-function(tseobject){
  df <- as.data.frame(colData(tseobject)) # change data set
  res.var<-df %>% shapiro_test(Richness, Diversity,Evenness)
  res.group<-df %>%
    group_by(outcome) %>%
    shapiro_test(Richness, Diversity,Evenness)
  return(list(res.var,res.group))
  
}
RunNormalityPlots<-function(tseobject,variable){
  df <- as.data.frame(colData(tseobject))
  var=df[[variable]]
  # Density plot
  p1<-ggdensity(var, fill = "lightgray")
  # QQ plot
  p2<-ggqqplot(var)
  list(plot(p1), plot(p2))
}

#Statistics (t-test and wilcox)
runStatistics<-function(treeSE){
  df <- as.data.frame(colData(treeSE))
  pvalue_wilcox<-sapply(df[,c("Diversity","Richness","Evenness")], function(x) wilcox.test(x~df$outcome)$p.value)
  pvalue_t.test<-sapply(df[,c("Diversity","Richness","Evenness")], function(x) t.test(x~df$outcome)$p.value)
  cbind(pvalue_t.test,pvalue_wilcox)
}

# p.aj 
library(stats)
p<-c() # vector with p-values
p.adjust.methods ="BH"
p.adjust(p, method = p.adjust.methods, n = length(p))

#+++++++++++++#
# Mediation
library(mediation)

# generalize
tse<-runAlphaDiversity(phyloseqObject) 
df <- as.data.frame(colData(tse)) #
sel<- df %>% select(L.reuteri,Diversity,Richness,Evenness,Language)  

# defining variables
X <-sel$L.reuteri
Y <-sel$Language 
M <-sel$Evenness # change Diversity, Richness, Evenness

# defining models
model.m <- glm(M ~ X, family=poisson) #a path
model.y <- glm(Y ~ X + M,family=poisson) # 

#model.m <- lm(M ~ X) #a path
#model.y <- lm(Y ~ X + M) # 
# mediate
set.seed(123)
mediation_results <- mediate(model.m = model.m,
                             model.y = model.y,
                             sims=5000,
                             boot = TRUE,
                             mediator = "M",
                             treat = "X")
summary(mediation_results)

#+++++++++++++#
# Mediation wiht confounding factors
tse<-runAlphaDiversity(phyloseqObject) 
df <- as.data.frame(colData(tse)) # 
sel<- df %>% select(L.reuteri,Diversity,Richness,Evenness,Language,Gestalder,Rok,Location,BPD,Apg10,FVg)
#con<- df %>% select(Gestalder,Rok,Location,BPD,Apg10,FVg)

# defining variables
X <-sel$L.reuteri
Y <-sel$Language 
M <-sel$Richness # change Diversity, Richness, Evenness
confounder<-sel$Location #Gestalder,Rok,Location,BPD,Apg10,FVg

# defining models
b_path=glm(Y~M,family=poisson())
summary(b_path)

model.m <- lm(M ~ X) #a path
model.y <- lm(Y ~ X + M) # 

# mediate
set.seed(123)
mediation_results <- mediate(model.m = model.m,
                             model.y = model.y,
                             sims=5000,
                             boot = TRUE,
                             mediator = "M",
                             covariates = confounder,
                             treat = "X")
summary(mediation_results)