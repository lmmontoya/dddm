library(dplyr)
library(SL.ODTR)
library(SuperLearner)
set.seed(4490)
#### 0a. Load data ####

load("/Users/linamontoya/Dropbox/Interventions - Working Data (1)/Lina/process data/analysis-data-interventions.RData")
data = data[!is.na(data$Y_nonminorarrest0to2yr),]
data = data[data$W_ResearchSite == "San Francisco",]
data$W_ResearchSite = NULL

W = data[,grep("W_", colnames(data))]
A = ifelse(data$A_Condition == "Intervention", 1, 0)
Y = data$Y_nonminorarrest0to2yr


#### 0b. Helper funs and libraries ####

EM_unadj_fun = function(i, W, A, Y) {
  toreturn = list()
  mod = W[,i]
  toreturn$EM = glm(Y ~ A*mod, family = "binomial")
  main = glm(Y ~ A+mod, family = "binomial")
  marg = glm(Y ~ mod, family = "binomial")
  toreturn$marg = summary(marg)$coefficients#[grep(":", rownames(summary$coefficients)),]
  toreturn$LRT_pval = as.numeric(na.omit(anova(toreturn$EM, main, test = "LRT")$`Pr(>Chi)`))
  toreturn$levels_W = levels(mod)
  toreturn$levels_A = levels(A)
  return(toreturn)
}

for (i in colnames(W)[grep("W_", colnames(W))]) {
  eval(parse(text=paste0("SL.QAW.",i,"<-function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
    if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  fit.glm <- glm(Y ~ A*", i, ", data = X, family = family, weights = obsWeights,
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = \"response\")
  fit <- list(object = fit.glm)
  class(fit) <- \"SL.glm\"
  out <- list(pred = pred, fit = fit)
  return(out)
  }")))
}

for (i in colnames(W)[grep("W_", colnames(W))]) {
  eval(parse(text=paste0("SL.blip.",i,"<-function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
    if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  fit.glm <- glm(Y ~ ", i, ", data = X, family = family, weights = obsWeights,
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = \"response\")
  fit <- list(object = fit.glm)
  class(fit) <- \"SL.glm\"
  out <- list(pred = pred, fit = fit)
  return(out)
  }")))
}

SL.library = c("SL.mean", "SL.earth", "SL.rpart", "SL.glm", "SL.randomForest", "SL.bayesglm", "SL.stepAIC") #SL.nnet, SL.svm (removed bc crash)
QAW.SL.library = c(ls()[grep("SL.QAW.", ls())], "SL.mean", "SL.glm")
blip.SL.library = c(ls()[grep("SL.blip.", ls())], SL.library)

#### 1. Subgroup analysis ####

# EM analysis
results_subgroup = lapply(1:ncol(W), EM_unadj_fun, W = W, A = A, Y = Y)
names(results_subgroup) = colnames(W)

save(results_subgroup, file = "paper dddm/2.results/1.subgroup_results.RData")

sig_Ws = names(results_subgroup)[unlist(lapply(1:length(results_subgroup), function(i) results_subgroup[[i]]$LRT_pval) < 0.1)]
sig_W_df = data[,sig_Ws]
tobreak = (sapply(sig_W_df, class) == "numeric") & sapply(sig_W_df, function(x) length(unique(x))) > 4
sig_W_df[,tobreak] = sapply(data.frame(sig_W_df[,tobreak]), function(x) paste0("Q", cut(x, breaks = quantile(x, (0:2)/2, na.rm = T), F, T)))
colnames(sig_W_df) = names(tobreak)
sig_W_df$A = data$A_Condition
sig_W_df$Y = data$Y_nonminorarrest0to2yr

# dtr for subgroup
dtr_subgroup_fun = function(sig_W) {
  
  print(sig_W)
  
  cate_df = sig_W_df %>%
    group_by(A, sig_W_df[,sig_W]) %>%
    summarize(mean(Y))
  cate_df = data.frame(y = cate_df[cate_df$A == "Intervention",]$`mean(Y)` - cate_df[cate_df$A == "Control",]$`mean(Y)`,
                       x = levels(factor(sig_W_df[,sig_W])))
  cate_df$rule = ifelse(cate_df$y < 0, 1, 0)
  
  dtr <- data.frame(cate_df$rule[match(sig_W_df[,sig_W], cate_df$x)])
  colnames(dtr) = paste0("subgroupdtr_", sig_W)
  
  return(dtr)

}

dtrs_subgroup = do.call('cbind', lapply(sig_Ws, dtr_subgroup_fun))
dtrs_subgroup$subgroupdtr_W_all = as.numeric(rowSums(dtrs_subgroup) > 0)

#### 2. Risk score ####
Y0 = Y[A == 0]
W0 = W[A == 0,]

results_SL = SuperLearner(Y = Y0, X = W0, SL.library = blip.SL.library, family = "binomial")
results_riskscore = predict(results_SL, newdata = W)$pred

# dtr for risk score
rho.grid = seq(0.1, 0.9, by = 0.1)
dtrs_riskscore = do.call('cbind', lapply(rho.grid, function(rho) ifelse(results_riskscore > rho, 1, 0)))
colnames(dtrs_riskscore) = paste0("rho_", rho.grid)

# dtr for LSI
results_LSI = odtr(W = W, V = W, A = A, Y = 1-Y, QAW.SL.library = QAW.SL.library, g.SL.library = "SL.mean", blip.SL.library = "SL.blip.W_LSIR_TotalScore15.9", risk.type = "CV MSE", metalearner = "blip")
dtr_LSI = results_LSI$dopt

save(results_SL, results_riskscore, results_LSI, file = "paper dddm/2.results/2.riskscore_results.RData")

#### 3. Optimal rule and values- full covariate set ####
dtrs = data.frame(treatall = 1, treatnone = 0, dtrs_riskscore, dtr_LSI, dtrs_subgroup)
W$W_Antisocial_IntentTotal_subgrp = dtrs_subgroup$subgroupdtr_W_Antisocial_IntentTotal

for (i in colnames(W)[grep("W_", colnames(W))]) {
  eval(parse(text=paste0("SL.QAW.",i,"<-function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
    if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  fit.glm <- glm(Y ~ A*", i, ", data = X, family = family, weights = obsWeights,
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = \"response\")
  fit <- list(object = fit.glm)
  class(fit) <- \"SL.glm\"
  out <- list(pred = pred, fit = fit)
  return(out)
  }")))
}

for (i in colnames(W)[grep("W_", colnames(W))]) {
  eval(parse(text=paste0("SL.blip.",i,"<-function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
    if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  fit.glm <- glm(Y ~ ", i, ", data = X, family = family, weights = obsWeights,
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = \"response\")
  fit <- list(object = fit.glm)
  class(fit) <- \"SL.glm\"
  out <- list(pred = pred, fit = fit)
  return(out)
  }")))
}

SL.library = c("SL.mean", "SL.earth", "SL.rpart", "SL.glm", "SL.randomForest", "SL.bayesglm", "SL.stepAIC") #SL.nnet, SL.svm (removed bc crash)
QAW.SL.library = c(ls()[grep("SL.QAW.", ls())], "SL.mean", "SL.glm")
blip.SL.library = c(ls()[grep("SL.blip.", ls())], SL.library)

set.seed(4490)
EYdopt_results = EYdopt(W = W, V = W, A = A, Y = 1-Y, QAW.SL.library = QAW.SL.library, blip.SL.library = blip.SL.library, risk.type = "CV MSE", contrast = dtrs, VFolds = 20)
save(EYdopt_results, file = "paper dddm/2.results/3.odtr_and_value_results.RData")


