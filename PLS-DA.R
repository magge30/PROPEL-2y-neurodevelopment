library(base)
library(mixomics)

# PLS-DA - summary
extract_ASV_metadata<-function(physeqObject){
  OTU1 = as(otu_table(physeqObject), "matrix")
  # transpose if necessary
  if(taxa_are_rows(physeqObject)){OTU1 <- t(OTU1)}
  # Coerce to data.frame
  asv = as.data.frame(OTU1)
  dat = as(sample_data(physeqObject), "data.frame")
  return(list(asv=asv,dat=dat))
}

#-------------------#
# example with language week 2
out<-run2v2.Language(p1)  # language normal+ mild vs rest

out <- microbiome::transform(out,"clr")
title<-"1w"

temp<-extract_ASV_metadata(out)
X<-temp$asv
Y<-temp$dat$outcome
shape_location<-temp$dat$Location

# Define shapes manually based on group
pch_vals <- ifelse(temp$dat$Location == "0", 1, 16)

res.plsda <- plsda(X, Y, ncomp = 2)
plotIndiv(res.plsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = FALSE, title = title,
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

#AUC ROC curve
res.plsda$auc
auc.splsda = auroc(res.plsda, roc.comp = 2, print = FALSE) # AUROC for all two components

#ASV
contribution<-plotLoadings(res.plsda, comp=1, contrib = 'max', method = 'median',title = title)
TAX = as(tax_table(p1), "data.frame")

#The function selectVar() outputs the variables selected for a given component and their loading values (ranked in decreasing absolute values).
head(selectVar(res.plsda, comp = 2)$value)

df<-merge(contribution,TAX, by="row.names", all.x=TRUE)

plotIndiv(res.plsda, ind.names = FALSE, legend=TRUE,
               ellipse = TRUE, star = FALSE, title = title,
               X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2', col=c("chocolate1","slateblue2"),pch=pch_vals,cex=3)

#+++++++++++++++++++#
# sPLS-DPA

set.seed(30) 
perf.final.plsda.srbct <- perf(res.plsda, validation = 'Mfold', folds = 5, 
                               progressBar = FALSE,  # Set to TRUE to track progress
                               nrepeat = 50)         # We suggest nrepeat = 50

plot(perf.final.plsda.srbct, sd = TRUE, legend.position = 'horizontal')
perf.final.plsda.srbct$error.rate$BER[, 'max.dist']
perf.final.plsda.srbct$error.rate.class$max.dist

# background prediction
background.max <- background.predict(res.plsda, 
                                     comp.predicted = 2,
                                     dist = 'max.dist') 
plotIndiv(res.plsda, comp = 1:2, group = Y,
          ind.names = FALSE, title = 'Maximum distance',
          legend = TRUE,  background = background.max)


# Number of variables to select
# Grid of possible keepX values that will be tested for each comp
list.keepX <- c(1:10,  seq(20, 100, 10))

# I add one more
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 3, validation = 'Mfold', 
                                 folds = 5, dist = 'max.dist', 
                                 test.keepX = list.keepX, nrepeat = 10)

# Just a head of the classification error rate per keepX (in rows) and comp
head(tune.splsda.srbct$error.rate)

plot(tune.splsda.srbct, sd = TRUE)

# The optimal number of components according to our one-sided t-tests
tune.splsda.srbct$choice.ncomp$ncomp

# The optimal keepX parameter according to minimal error rate
tune.splsda.srbct$choice.keepX

# Final model performance
# Optimal number of components based on t-tests on the error rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 

# Optimal number of variables to select
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  
select.keepX

# I have to add ncom=2 for the plotting
splsda.srbct <- mixOmics::splsda(X, Y, ncomp = 2, keepX = select.keepX) 

my.colors<-c("chocolate1","slateblue2")
plotLoadings(splsda.srbct, comp=1, contrib = 'max', method = 'median', legend.color =my.colors)


# plot
plotIndiv(splsda.srbct, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = FALSE, title ="",
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2',col=c("chocolate1","slateblue2"),pch=16,cex=1,style = 'graphics')
