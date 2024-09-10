library(sva)
library(dplyr)
library(DASC)
library(magrittr)
library(ggplot2)
library(ggfortify)
library(proBatch)
library(parallel)
library(Harman)
library(here)

# This script runs the case study for the hidden batch correction decision tree in Figure 4.

# Set seet for reproducibility
set.seed(64318418)

# Set working directory
wd <- here()

#### Load and prepare the data ----
multipro <- read.csv(paste0(wd, "/Thinking Points for Effective Batch Correction on Biomedical Data/DIA_protein_matrix.csv")) # Load the data
exprs <- multipro[,2:37] # Expression matrix
batch_f <- as.factor(rep(c(1,2), each=18)) # Batch factor
class_f <- as.factor(rep(c(1,2), each=9, times=2)) # Class factor
pdata <- data.frame(batch=batch_f, class=class_f) # Set the metadata

# Simulate hidden batch effects
nhidden <- 4 # Number of hidden batches
hidden_f <- as.factor(sample(1:nhidden, 36, replace=T)) # Assign hidden batch labels

# Check if the hidden batch factor is confounded with class or batch
table(hidden_f, batch_f); table(hidden_f, class_f) # No confounding
pdata$hidden <- hidden_f

# PCA visualization
autoplot(prcomp(t(na.omit(exprs)), scale=F, center=T), data=as.data.frame(pdata), color="hidden")

# Simulate batch effects by adding mean and variance shifts in hidden batches
batch_means <- matrix(rnorm(nrow(exprs) * nhidden, 0, 0.5), nrow(exprs), nhidden) #  hidden batch effect means
batch_matrix <- matrix(0, nrow(exprs), ncol(exprs))
for (i in seq_len(nrow(exprs))) { # Create  
  for (j in seq_len(ncol(exprs))) {
    k <- hidden_f[j]
    batch_matrix[i, j] <- rnorm(1, batch_means[i, k], 0.1)
  }
}
exprs2 <- exprs + batch_matrix

# Batch effect correction by combat
rownames(exprs2) <- multipro$prot # Set rownames to protein names
exprs3 <- na.omit(exprs2) # Remove proteins with missing values
bc_exprs3 <- sva::ComBat(exprs3, batch=batch_f, mod=NULL) # Apply ComBat

#### Level 1 of tree: Are there still hidden/unreported sources of variance? (Yes) ----
# Visualisations
autoplot(prcomp(t(bc_exprs3), scale=F, center=T), data=as.data.frame(pdata), color="class")
autoplot(prcomp(t(bc_exprs3), scale=F, center=T), data=as.data.frame(pdata), color="batch")
# true hidden factor
autoplot(prcomp(t(bc_exprs3), scale=F, center=T), data=as.data.frame(pdata), color="hidden")

## Choosing K clusters
# DASC to detect batch effects
#' lambda = hyperparameter
#' rank = hyperparameter for the number of batches
#' nrun = number of Semi-NMF iterations

## Parameter tuning. Takes around 10-15 minutes to run
out_list <- list()
dispersionsl1 <- c()
dispersionsl2 <- c()
dispersionsl3 <- c()
starttotal = Sys.time()
for (i in 1:9){
  start = Sys.time()
  dasc_iter1 <- DASC::DASC(bc_exprs3, pdata, factor=pdata$class, lambda=0.1, rank=(i+1), nrun=50)
  dasc_iter2 <- DASC::DASC(bc_exprs3, pdata, factor=pdata$class, lambda=0.01, rank=(i+1), nrun=50)
  dasc_iter3 <- DASC::DASC(bc_exprs3, pdata, factor=pdata$class, lambda=0.001, rank=(i+1), nrun=50)
  
  out_list[i][1] <- dasc_iter1[[1]]
  out_list[i][2] <- dasc_iter2[[1]]
  out_list[i][3] <- dasc_iter3[[1]]
  
  dispersionsl1 <- c(dispersionsl1, dasc_iter1[[1]]$dispersion)
  dispersionsl2 <- c(dispersionsl2, dasc_iter2[[1]]$dispersion) 
  dispersionsl3 <- c(dispersionsl3, dasc_iter3[[1]]$dispersion)
  end = Sys.time()
  print(paste0("Iteration ",i," time taken: ",end-start))
}
endtotal = Sys.time()
print(paste0("Total optimization time taken: ",endtotal-starttotal))
dispersions <- reshape2::melt(data.frame(lambda0.1=dispersionsl1, lambda0.01=dispersionsl2, lambda0.001=dispersionsl3))
dispersions$x <- 1:9
ggplot(dispersions, aes(x=x, y=value, color=variable)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:9, labels=c(2:10)) + 
  labs(y="Dispersion", x="No. of batches", color="Lambda", title="DASC parameter tuning")

## Apply DASC with the obtained parameters
dasc_out <- DASC::DASC(bc_exprs3, NULL, factor=pdata$class, lambda=0.01, rank=4, nrun=50)
# Visualize the labels
autoplot(prcomp(t(bc_exprs3), scale=F, center=T), data=data.frame(batch=dasc_out[[1]]$class, class=class_f), fill="batch", color="black", shape="class", size=3) +
  scale_shape_manual(values=c(21,24)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  labs(title="DASC labels", fill="hidden")
# Compare it with the true hidden labels
autoplot(prcomp(t(bc_exprs3), scale=F, center=T), data=as.data.frame(pdata), fill="hidden", color="black", shape="class", size=3) +
  labs(title="True labels") +
  scale_shape_manual(values=c(21,24)) +
  guides(fill = guide_legend(override.aes = list(shape=21)))

# Manuscript Figure 4C
autoplot(prcomp(t(bc_exprs3), scale=F, center=T), data=data.frame(batch=dasc_out[[1]]$class, class=class_f), fill="batch", color="black", shape="class", size=3) +
  scale_shape_manual(values=c(21,24)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  labs(title="Known batch corrected dataset", fill="hidden \n batch")

#### Level 2 of tree: Does it account for large variance? (Yes) ----
pcs <- prcomp(t(bc_exprs3), scale=F, center=T)
# pcs <- prcomp(t(na.omit(exprs)), scale=F, center=T)

pcregression_res <- kBET::pcRegression(pcs, batch=dasc_out[[1]]$class)$maxVar
(pcregression_res <- sum(pcregression_res)/100)
# 0.399


#### Level 3 of tree: Can we explain the hidden batch factor? (Here, we assume yes) ----

#### Level 4 of tree: Are the hidden batches confounded with the factor of interest? (No) ----
table(dasc_out[[1]]$class, pdata$class)

# Correct the hidden batch effect
bc2_exprs3 <- sva::ComBat(bc_exprs3, batch = dasc_out[[1]]$class)

# PCA visualizations of hidden batch, class, and batch factors
autoplot(prcomp(t(bc2_exprs3), scale=F, center=T), data=as.data.frame(dasc_out[[1]]), color="class")
autoplot(prcomp(t(bc2_exprs3), scale=F, center=T), data=as.data.frame(pdata), color="class")
#autoplot(prcomp(t(bc2_exprs3), scale=F, center=T), data=as.data.frame(pdata), color="batch")


#### Level 5 of tree: Perform analysis on both to see if big changes are observed (Yes) ----

# t-test on the class factor and compare results
t.test.wrapper <- function(df, class_factor, output="features"){
  unique_classes <- unique(class_factor)
  class1 = df[,which(class_factor==unique_classes[1])]
  class2 = df[,which(class_factor==unique_classes[2])]
  result=c()
  for(i in 1:nrow(df))
  {
    if(sd(class1[i,])<1e-6 && sd(class2[i,])<1e-6)
    {
      class1[i,1]=jitter(class1[i,1])
    }
    a=t.test(class1[i,],class2[i,])
    result=rbind(result,a$p.value)
  }
  result=p.adjust(result,method = "bonferroni")
  result=as.data.frame(result)
  rownames(result)=rownames(class1)
  
  sig.feats = rownames(result)[which(result<0.05)]
  
  if (output == "features"){
    return(sig.feats)
  }
  if (output == "pvals"){
    return(result)
  }
}

# t-test on data without hidden batch effect
rownames(exprs) <- multipro$prot
bc_no_hidden <- sva::ComBat(na.omit(exprs), batch=batch_f) # DEG of dataset without hidden batches simulated, corrected for the main batch effect
correct_once_no_hidden_features <- t.test.wrapper(bc_no_hidden, class_factor = pdata$class, output = "features")

# t-test q-values (BH-adjusted p-values)
no_correct <- t.test.wrapper(exprs3, class_factor = pdata$class, output = "pvals")
correct_once <- t.test.wrapper(bc_exprs3, class_factor = pdata$class, output = "pvals")
correct_twice <- t.test.wrapper(bc2_exprs3, class_factor = pdata$class, output = "pvals")
correct_once_no_hidden <- t.test.wrapper(bc_no_hidden, class_factor = pdata$class, output = "pvals")

# Identify DEGs
no_correct_features <- t.test.wrapper(exprs3, class_factor = pdata$class, output = "features") # DEG of dataset without correction
correct_once_features <- t.test.wrapper(bc_exprs3, class_factor = pdata$class, output = "features") # DEG of dataset with correction for the main batch effect
correct_twice_features <- t.test.wrapper(bc2_exprs3, class_factor = pdata$class, output = "features") # DEG of dataset with correction for both main and hidden batch effects

length(intersect(correct_once_features, no_correct_features))
length(intersect(correct_once_features, correct_twice_features))
different_features <- setdiff(correct_twice_features, correct_once_features) # Which features are different?

# t-test on data without hidden batch effect
rownames(exprs) <- multipro$prot
bc_no_hidden <- sva::ComBat(na.omit(exprs), batch=batch_f) # DEG of dataset without hidden batches simulated, corrected for the main batch effect
correct_once_no_hidden_features <- t.test.wrapper(bc_no_hidden, class_factor = pdata$class, output = "features")

# Plot the barchart of total DEGs
diff_features <- c(`No correction`=length(no_correct_features),
                   `Corrected M`=length(correct_once_features),
                   `Corrected M+H`=length(correct_twice_features),
                   `Corrected M \n (no hidden)`= length(correct_once_no_hidden_features))
diff_features_mat <- data.frame(labs=factor(names(diff_features), levels=c("No correction",
                                                                           "Corrected M",
                                                                           "Corrected M+H",
                                                                           "Corrected M \n (no hidden)")), value=diff_features)

# Manuscript Figure 4D
ggplot(diff_features_mat, aes(x=labs, y=value, fill=labs)) +
  geom_bar(position="dodge", stat="identity", fill="#abc4ff", color="black") +
  geom_text(aes(y=1000), label=diff_features_mat$value, fontface="bold") +
  theme_bw() + labs(x="Data", y="DEG", title="Number of DEGs")


# Identify top 1000 DEGs
no_correct_top1000 <- no_correct %>% mutate(prots=rownames(.)) %>% .[order(.[,1], decreasing = F),] %>% rownames() %>% .[1:1000]
correct_once_top1000 <- correct_once %>% mutate(prots=rownames(.)) %>% .[order(.[,1], decreasing = F),] %>% rownames() %>% .[1:1000]
correct_twice_top1000 <- correct_twice %>% mutate(prots=rownames(.)) %>% .[order(.[,1], decreasing = F),] %>% rownames() %>% .[1:1000]
correct_once_no_hidden_top1000 <- correct_once_no_hidden %>% mutate(prots=rownames(.)) %>% .[order(.[,1], decreasing = F),] %>% rownames() %>% .[1:1000]

# Jaccard similarity function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# Determine similarity to the groundtruth using Jaccard index
jaccard_dat <- diff_features_mat
jaccard_dat$value <- c(
  jaccard(no_correct_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_once_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_twice_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_once_no_hidden_top1000, correct_once_no_hidden_top1000))

# Manuscript Figure 4E
ggplot(jaccard_dat, aes(x=labs, y=value, fill=labs)) +
  geom_bar(position="dodge", stat="identity", fill="#abc4ff", color="black") +
  geom_text(aes(y=0.1), label=round(jaccard_dat$value, 3), fontface="bold") +
  theme_bw() + labs(x="Data", y="Jaccard index", title="Top 1000 DEGs similarity to ground truth")

#### Level 6 of tree: Compare with other methods to assess stability (Results appear stable) ----
# Correct the hidden batch effect using Harman and SVA instead
# SVA correction function
SVAcorrect <- function(edata, pdata, n.batches=NULL, sva.batch=T){
  pheno = as.data.frame(pdata)
  
  # perform SVA
  mod = model.matrix(~as.factor(pheno$class), data=pheno)
  mod0 = model.matrix(~1,data=pheno)
  
  # filter and remove features that are completely 0
  exp_mean <- apply(edata,1,mean)
  filter_features <- which(exp_mean == 0)
  if(length(filter_features) > 0){
    f.edata = edata[-filter_features,]
  } else {
    f.edata = edata
  }
  if (sva.batch == TRUE){ # Select number of surrogate variables
    n.sv = num.sv(f.edata,mod,method="leek")
  } else {
    n.sv = n.batches
  }
  svobj = sva(f.edata,mod,mod0,n.sv=n.sv) # n.sv = of surrogate variables (1 batch)
  
  # Correct the data
  X = cbind(mod, svobj$sv)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(f.edata))
  rm(Hat)
  gc()
  P = ncol(mod)
  sva_edata = f.edata - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  return(sva_edata)
}

sva_exprs <- SVAcorrect(bc_exprs3, pdata[,c("class","hidden")], sva.batch=F, n.batches=length(unique(dasc_out[[1]]$class)))
correct_twice_features <- t.test.wrapper(sva_exprs, class_factor = pdata$class, output = "features") # DEG of dataset with correction for both main and hidden batch effects

# Plot the barchart of total DEGs
diff_features <- c(`No correction`=length(no_correct_features),
                   `Corrected M`=length(correct_once_features),
                   `Corrected M+H`=length(correct_twice_features),
                   `Corrected M \n (no hidden)`= length(correct_once_no_hidden_features))
diff_features_mat <- data.frame(labs=factor(names(diff_features), levels=c("No correction",
                                                                           "Corrected M",
                                                                           "Corrected M+H",
                                                                           "Corrected M \n (no hidden)")), value=diff_features)

ggplot(diff_features_mat, aes(x=labs, y=value, fill=labs)) +
  geom_bar(position="dodge", stat="identity", fill="#abc4ff", color="black") +
  geom_text(aes(y=1000), label=diff_features_mat$value, fontface="bold") +
  theme_bw() + labs(x="Data", y="DEG", title="Number of DEGs (SVA)")

# Identify top 1000 DEGs
correct_twice_top1000 <- correct_twice %>% mutate(prots=rownames(.)) %>% .[order(.[,1], decreasing = F),] %>% rownames() %>% .[1:1000]

# Determine similarity to the groundtruth using Jaccard index
jaccard_dat <- diff_features_mat
jaccard_dat$value <- c(
  jaccard(no_correct_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_once_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_twice_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_once_no_hidden_top1000, correct_once_no_hidden_top1000))

ggplot(jaccard_dat, aes(x=labs, y=value, fill=labs)) +
  geom_bar(position="dodge", stat="identity", fill="#abc4ff", color="black") +
  geom_text(aes(y=0.1), label=round(jaccard_dat$value, 3), fontface="bold") +
  theme_bw() + labs(x="Data", y="Jaccard index", title="Top 1000 DEGs similarity to ground truth (SVA)")

# Harman correction
harman_out <- harman(bc_exprs3, expt=class_f, batch=dasc_out[[1]]$class)
harman_exprs <- reconstructData(harman_out)
correct_twice_features <- t.test.wrapper(harman_exprs, class_factor = pdata$class, output = "features") # DEG of dataset with correction for both main and hidden batch effects

# Plot the barchart of total DEGs
diff_features <- c(`No correction`=length(no_correct_features),
                   `Corrected M`=length(correct_once_features),
                   `Corrected M+H`=length(correct_twice_features),
                   `Corrected M \n (no hidden)`= length(correct_once_no_hidden_features))
diff_features_mat <- data.frame(labs=factor(names(diff_features), levels=c("No correction",
                                                                           "Corrected M",
                                                                           "Corrected M+H",
                                                                           "Corrected M \n (no hidden)")), value=diff_features)

ggplot(diff_features_mat, aes(x=labs, y=value, fill=labs)) +
  geom_bar(position="dodge", stat="identity", fill="#abc4ff", color="black") +
  geom_text(aes(y=1000), label=diff_features_mat$value, fontface="bold") +
  theme_bw() + labs(x="Data", y="DEG", title="Number of DEGs (Harman)")

# Identify top 1000 DEGs
correct_twice_top1000 <- correct_twice %>% mutate(prots=rownames(.)) %>% .[order(.[,1], decreasing = F),] %>% rownames() %>% .[1:1000]

# Determine similarity to the groundtruth using Jaccard index
jaccard_dat <- diff_features_mat
jaccard_dat$value <- c(
  jaccard(no_correct_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_once_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_twice_top1000, correct_once_no_hidden_top1000),
  jaccard(correct_once_no_hidden_top1000, correct_once_no_hidden_top1000))

ggplot(jaccard_dat, aes(x=labs, y=value, fill=labs)) +
  geom_bar(position="dodge", stat="identity", fill="#abc4ff", color="black") +
  geom_text(aes(y=0.1), label=round(jaccard_dat$value, 3), fontface="bold") +
  theme_bw() + labs(x="Data", y="Jaccard index", title="Top 1000 DEGs similarity to ground truth (Harman)")

# Write session info
# writeLines(capture.output(sessionInfo()), "Thinking Points for Effective Batch Correction on Biomedical Data/sessionInfo.txt")
