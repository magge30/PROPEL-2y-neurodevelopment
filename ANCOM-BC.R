library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)

# convert phylo to tse object
tse<-makeTreeSEFromPhyloseq(phyloseqobject) ## change p1.lan, p2.lan, p3.lan, p4.lan
colData(tse)$outcome = factor(colData(tse)$outcome, levels = c("Normal", "Impairment"))
levels(colData(tse)$outcome)

######
set.seed(123)

#group = NULL because it is only a two group categorical. But then I don't get the structural zeros. 
output = ancombc2(data = tse, 
                  assay_name = "counts", 
                  tax_level = "Genus",
                  fix_formula = "L.reuteri + outcome", 
                  rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, # 10%
                  lib_cut = 0, 
                  s0_perc = 0.05,
                  group = NULL, 
                  struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = FALSE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 10))

res_prim = output$res

res_prim %>%
  dplyr::select(taxon, contains("outcome")) %>%
  filter(q_outcomeImpairment < 0.05 & passed_ss_outcomeImpairment == TRUE) %>%
  arrange(q_outcomeImpairment) %>%
  head() %>%
  knitr::kable()

# Results for supplementation:
res_prim %>%
  dplyr::select(taxon, contains("L.reuteri1")) %>%
  filter(q_L.reuteri1 < 0.05 & passed_ss_L.reuteri1 == TRUE) %>%
  arrange(q_L.reuteri1) %>%
  head() %>%
  knitr::kable()