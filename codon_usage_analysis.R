
#################################################
####    CYO PROJECT: CODON USAGE ANALYSIS    ####
#################################################

# Package installations
if(!require(DiagrammeR)) install.packages("DiagrammeR", repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(knitr)) install.packages("knitr", repos = "http://cran.us.r-project.org")
if(!require(gplots)) install.packages("gplots", repos = "http://cran.us.r-project.org")
if(!require(RColorBrewer)) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
if(!require(kernlab)) install.packages("kernlab", repos = "http://cran.us.r-project.org")
if(!require(randomForest)) install.packages("randomForest", repos = "http://cran.us.r-project.org")

# Load libraries
library(DiagrammeR)   # creation of flowcharts
library(tidyverse)
library(caret)        # training and tuning of models
library(knitr)        # table generation
library(gplots)       # heatmap
library(RColorBrewer) # color palettes

#  Global setting: number of significant digits
options(digits = 3)

# Downloading, unzipping and reading of the codon usage data set
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00577/codon_usage.csv.zip"
dl <- tempfile()
download.file(url, dl)
unzip(dl, "codon_usage.csv")
dat <- read.csv("codon_usage.csv")
  
# Dimension of the data set
dim(dat)

# Classes of the data set
sapply(dat,class)

# First 6 rows
head(dat)

# Structure of the data set
str(dat)

# Summary for each column
summary(dat)


#####################
# Genetic code table
#####################

Codon1 <- c("UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AUU", "AUC", "AUA", "AUG", "GUU", "GUC", "GUA", "GUG")
Aminoacid1 <- c("Phenylalanine", "Phenylalanine", "Leucine", "Leucine", "Leucine", "Leucine", "Leucine", "Leucine", "Isoleucine", "Isoleucine", "Isoleucine", "Methionine", "Valine", "Valine", "Valine", "Valine")
Codon2 <- c("UCU", "UCC", "UCA", "UCG", "CCU", "CCC", "CCA", "CCG", "ACU", "ACC", "ACA", "ACG", "GCU", "GCC", "GCA", "GCG")
Aminoacid2 <- c("Serine", "Serine", "Serine", "Serine", "Proline", "Proline", "Proline", "Proline", "Threonine", "Threonine", "Threonine", "Threonine", "Alanine", "Alanine", "Alanine", "Alanine")
Codon3 <- c("UAU", "UAC", "UAA", "UAG", "CAU", "CAC", "CAA", "CAG", "AAU", "AAC", "AAA", "AAG", "GAU", "GAC", "GAA", "GAG")
Aminoacid3 <- c("Tyrosine", "Tyrosine", "stop codon", "stop codon", "Histidine", "Histidine", "Glutamine", "Glutamine", "Asparagine", "Asparagine", "Lysine", "Lysine", "Aspartate", "Aspartate", "Glutamate", "Glutamate")
Codon4 <- c("UGU", "UGA", "UGA", "UGG", "CGU", "CGC", "CGA", "CGG", "AGU", "AGC", "AGA", "AGG", "GGU", "GGC", "GGA", "GGG")
Aminoacid4 <- c("Cysteine", "Cysteine", "stop codon", "Tryptophan", "Arginine", "Arginine", "Arginine", "Arginine", "Serine", "Serine", "Arginine", "Arginine", "Glycine", "Glycine", "Glycine", "Glycine")
codons <- data.frame(Codon1, Aminoacid1, Codon2, Aminoacid2, Codon3, Aminoacid3, Codon4, Aminoacid4)
codons %>% kable()


######################
# Tree of life figure
######################

# Creation of a flowchart using DiagrammeR package
grViz("
digraph boxes_and_circles {

  node [shape = box, fontsize = 12,
        fontname = Helvetica]
        
  Life Eukarya Fungi Plants[style = filled, fillcolor = LightSalmon]
  Animals Fishes Amphibians Reptiles Birds
  
  node [style = filled, fillcolor = LightSalmon]
        
  Archaea Bacteria Plasmids Viruses Bacteriophages Invertebrates
  Vertebrates Mammals Rodents Primates
  
  Life->{Bacteria Archaea Eukarya}
  Viruses->Bacteriophages
  Eukarya->{Fungi Plants Animals}
  Animals->{Invertebrates Vertebrates}
  Vertebrates->{Fishes Amphibians Reptiles Birds Mammals}
  Mammals->{Rodents Primates}
  
  subgraph {
  rank = same; Plasmids Eukarya Archaea Bacteria Viruses
  }
}
")


###############
#data cleaning
###############

# Re-formatting of columns and substitution of NAs by 0s
dat <- dat %>%
  mutate(Kingdom = as.factor(Kingdom), UUU = as.numeric(UUU), UUC = as.numeric(UUC)) %>%
  replace_na(list(UUU = 0, UUC = 0))

#check for remaining NAs
sum(is.na(dat))


############################################################
# Separation of a validation set as final hold-out test set
############################################################

# Validation set will be 10% of the data
set.seed(2005)
test_index <- createDataPartition(y = dat$Kingdom, times = 1, p = 0.1, list = FALSE)
train <- dat[-test_index,]
validation <- dat[test_index,]


################
# Data analysis
################

# Name vector for kingdom bar plot
Kingdoms <- c("Plasmid (plm)", "Archaea (arc)", "Primate (pri)", "Rodent (rod)",
              "Bacteriophage (phg)", "Mammal (mam)", "Invertebrate (inv)",
              "Vertebrate (vrt)", "Plant (pln)", "Virus (vrl)", "Bacteria (bct)")

# Bar plot: counts kingdom
train %>%
  group_by(Kingdom) %>%
  summarize(count = n()) %>%
  mutate(Kingdom = reorder(Kingdom, count)) %>%
  ggplot(aes(Kingdom, count)) +
  geom_bar(stat = "identity", fill = "steelblue3", width = 0.8) +
  coord_flip() +
  scale_x_discrete(labels = Kingdoms) +
  theme(axis.text.y = element_text(size = 10)) +
  xlab("") +
  ylab("Number of entries")

# Boxplot: Entries per SpeciesID
train %>%
  group_by(SpeciesID) %>%
  summarize(entries = n()) %>%
  group_by(entries) %>%
  summarize(count = n()) %>%
  ggplot(aes(entries, count)) +
  geom_bar(stat = "identity", fill = "steelblue3") +
  scale_y_continuous(trans = "log10") +
  xlab("Entries per SpeciesID") +
  ylab("Count")

# Table: counts per DNAtype
Name <- c("genomic", "mitochondrial", "chloroplast", "cyanelle", "plastid", "nucleomorph", "secondary_endosymbiont", "chromoplast", "NA", "apicoplast", "kinetoplast")
typetab <- train %>%
  group_by(DNAtype) %>%
  summarize(Count = n()) %>%
  mutate(Share = round(Count/sum(Count)*100, 2))
cbind(Name, typetab) %>% kable()
  
# Bar plot: DNAtype by Kingdom
train %>%
  group_by(Kingdom, DNAtype) %>%
  summarize(count = n()) %>%
  #mutate(Kingdom = reorder(Kingdom, count)) %>%
  ggplot(aes(Kingdom, count, fill = as.factor(DNAtype))) +
  geom_bar(position= position_fill(reverse = TRUE), stat = "identity") +
  scale_fill_brewer(name = "DNA type", palette = "Paired") +
  xlab("") +
  ylab("Proportion of DNA types") +
  theme(axis.text.x = element_text(size=11))

# Boxplot: Ncodons by Kingdom
train %>%
  ggplot(aes(Kingdom, Ncodons)) +
  geom_boxplot(fill = "steelblue3") +
  scale_y_continuous(trans = "log10") +
  xlab("") +
  ylab("Sum of codon numbers") +
  theme(axis.text.x = element_text(size=11))

# Histogram: sums of frequencies \
train %>%
  select(-Kingdom, -DNAtype, -SpeciesID, -Ncodons, -SpeciesName) %>%
  summarize(Row_sums = rowSums(.)) %>%
  ggplot(aes(Row_sums)) +
  geom_histogram(binwidth = 0.00002, col = "black", fill = "steelblue3")

# Boxplot: codon usage, sorted by mean frequency
train %>%
  select(-DNAtype, -SpeciesID, -Ncodons, -SpeciesName) %>%
  gather(Codon, Frequency, -Kingdom) %>%
  mutate(Codon = reorder(Codon, desc(Frequency))) %>%
  ggplot(aes(Codon, Frequency)) +
  geom_boxplot(outlier.size = 0.6, fill = "steelblue3") +
  scale_y_continuous(trans = "log10", limits = c(0.0001,0.3)) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1))
  
# Boxplot: faceting codon usage by vrl, rod and pri
train %>%
  select(-DNAtype, -SpeciesID, -Ncodons, -SpeciesName) %>%
  gather(Codon, Frequency, -Kingdom) %>%
  mutate(Codon = reorder(Codon, desc(Frequency))) %>%
  filter(Kingdom %in% c("vrl", "rod", "pri")) %>%
  ggplot(aes(Codon, Frequency)) +
  geom_boxplot(outlier.size = 0.6, fill = "steelblue3") +
  scale_y_continuous(trans = "log10", limits = c(0.0001,0.3)) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1)) +
  facet_grid(Kingdom ~ .)

# Bar plot: codons sorted by variance
train %>%
  select(-Kingdom, -DNAtype, -SpeciesID, -Ncodons, -SpeciesName) %>%
  gather(Codon, value) %>%
  group_by(Codon) %>%
  summarize(Codon = unique(Codon), Variance = var(value)) %>%
  arrange(desc(Variance)) %>%
  mutate(Codon = reorder(Codon, desc(Variance))) %>%
  ggplot(aes(Codon, Variance)) +
  geom_bar(stat = "identity", col = "black", fill = "steelblue3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1)) +
  xlab("")

# Transformation of predictors into matrix
usage <- train %>%
  select(-Kingdom, -DNAtype, -SpeciesID, -Ncodons, -SpeciesName) %>%
  as.matrix()

# Color palette with 11 colors from red to white to blue mixed to 100 transition colors
hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

# Heatmap of correlations
heatmap.2(cor(usage), trace = "none", col = hmcol, cexRow = 0.4, cexCol = 0.4)

# Correlation values for dark blue spots in the heatmap
cor(train$CUA, train$UGA)
cor(train$AUU, train$UUA)
cor(train$CCG, train$GCG)
cor(train$CCG, train$CGC)
cor(train$CGC, train$GCG)
cor(train$AAU, train$UAU)
cor(train$AUC, train$UUC)
cor(train$GGC, train$GCC)
cor(train$UCG, train$ACG)


##############################################
# Removing dispensable columns for prediction
##############################################

train <- train %>% select(-DNAtype, -SpeciesID, -Ncodons, -SpeciesName)
validation <- validation %>% select(-DNAtype, -SpeciesID, -Ncodons, -SpeciesName)



########################
####     Results    ####
########################


#################
# multinom model
#################

set.seed(2005)

# Control for the cross-validation method
control <- trainControl(method = "cv", number = 10, p = .9)

# Training and tuning the model
multimodel <- train(Kingdom ~ ., data = train, method = "multinom",
                    trControl = control,
                    trace = FALSE,    # suppress iterations output
                    tuneLength = 5)   # 5 different decay values
                    

# Plot of decay tuning
ggplot(multimodel, highlight = TRUE) + scale_x_continuous(trans = "log10")

# Best decay value
multimodel$bestTune

# Extraction of best accuracy
acc_multi <- max(multimodel$results$Accuracy)
acc_multi

# Accuracy table
acc_results <- tibble(Method = "Multinomial Logistic Regression Model",
                      Accuracy = acc_multi)
acc_results %>% kable()


############
# lda model
############

# Training the model
set.seed(2005)
control <- trainControl(method = "cv", number = 10, p = .9)
ldamodel <- train(Kingdom ~ ., data = train, method = "lda",
                  trControl = control)

# Extraction of accuracy
acc_lda <- ldamodel$results$Accuracy
acc_lda

# Accuracy table
acc_results <- bind_rows(acc_results,
                         tibble(Method="Linear Discriminant Analysis Model",
                                Accuracy = acc_lda))
acc_results %>% kable()


############
# knn model
############

# Training and tuning the model
set.seed(2005)
control <- trainControl(method = "cv", number = 10, p = .9)
knnmodel <- train(Kingdom ~ ., data = train, method = "knn",
                  trControl = control, 
                  tuneGrid = data.frame(k = seq(1, 7, 2)))

# Plot of k tuning
ggplot(knnmodel, highlight = TRUE)

# Best k value
knnmodel$bestTune

# Extraction of best accuracy
acc_knn <- max(knnmodel$results$Accuracy)
acc_knn

# Accuracy table
acc_results <- bind_rows(acc_results,
                         tibble(Method="k-Nearest Neighbors Model",
                                Accuracy = acc_knn))
acc_results %>% kable()


##################
# svmlinear model
##################

# Training and tuning the model
set.seed(2005)
control <- trainControl(method = "cv", number = 10, p = .9)
svmlinmodel <- train(Kingdom ~ ., data = train,  method = "svmLinear",
                     trControl = control,
                     tuneGrid = data.frame(C = seq(0.05, 0.3, 0.05)))

# Plot of C tuning
ggplot(svmlinmodel, highlight = TRUE)

# Best C value
svmlinmodel$bestTune

# Extraction of best accuracy
acc_svmlin <- max(svmlinmodel$results$Accuracy)
acc_svmlin

# Accuracy table
acc_results <- bind_rows(acc_results,
                         tibble(Method="Support Vector Machines Linear Model",
                                Accuracy = acc_svmlin))
acc_results %>% kable()


##################
# svmradial model
##################

# Training and tuning the model
set.seed(2005)
control <- trainControl(method = "cv", number = 5, p = .8)
svmradmodel <- train(Kingdom ~ ., data = train, method = "svmRadial",
                     trControl = control,
                     tuneGrid = expand.grid(sigma = seq(.018, 0.030, 0.004),
                                            C = seq(5, 9, 1)))

# Plot of sigma and C tuning
ggplot(svmradmodel, highlight = TRUE)

# Best sigma and C values
svmradmodel$bestTune

# Extraction of best accuracy
acc_svmrad <- max(svmradmodel$results$Accuracy)
acc_svmrad

# Accuracy table
acc_results <- bind_rows(acc_results,
                         tibble(Method="Support Vector Machines Radial Model",
                                Accuracy = acc_svmrad))
acc_results %>% kable()


#######################
# Random forests model
#######################

# Training the model and tuning for mtry parameter
set.seed(2005)
control <- trainControl(method = "cv", number = 5, p = .8)
mtrytune <- train(Kingdom ~ ., data = train, method = "rf",
                 trControl = control,
                 tuneGrid = data.frame(mtry = seq(6,20,2)))

# Plot of mtry tuning
ggplot(mtrytune, highlight = TRUE)

# Best mtry value
mtrytune$bestTune

# Training the model and tuning for ntree parameter
set.seed(2005)
control <- trainControl(method = "cv", number = 5, p = .8)
ntree <- c(500, 1000, 1500, 2000)
ntreetune <- lapply(ntree, function(nt){
  train(Kingdom ~ ., data = train, method = "rf",
        trControl = control,
        tuneGrid = data.frame(mtry = mtrytune$bestTune$mtry),
        ntree = nt)
})

# Extraction of accuracy results
rf_accuracy <- c(ntreetune[[1]]$results$Accuracy, ntreetune[[2]]$results$Accuracy, 
                 ntreetune[[3]]$results$Accuracy, ntreetune[[4]]$results$Accuracy)

# Plot of ntree tuning
qplot(ntree, rf_accuracy, geom=c("point", "line"), ylab = "Accuracy (Cross−Validation)")

# Best ntree value
ntree[which.max(rf_accuracy)]

# Best accuracy
acc_rf <- max(rf_accuracy)
acc_rf

# Accuracy table
acc_results <- bind_rows(acc_results,
                         tibble(Method="Random Forests Model",
                                Accuracy = acc_rf))
acc_results %>% kable()


##################################
# Analysis of variable importance
##################################

# Variable importance of the best rf model
imp <- varImp(ntreetune[[which.max(rf_accuracy)]])

# Conversion into a tibble for use in ggplot
varimp <- tibble(Codon = row.names(imp$importance),
                 Importance = imp$importance$Overall)

# Bar plot: codons sorted by variable importance
varimp %>% 
  mutate(Codon = reorder(Codon, desc(Importance))) %>%
  ggplot(aes(Codon, Importance)) +
  geom_bar(stat = "identity", col = "black", fill = "steelblue3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1)) +
  xlab("") +
  ylab("Variable importance")

# Extraction of the 6 codons with the highest varImp
topcodons <- varimp %>%
  arrange(desc(Importance)) %>% top_n(6) %>% .$Codon

# Boxplot of this 6 codons stratified by kingdom
train %>%
  select(Kingdom, all_of(topcodons)) %>%
  gather(Codon, Frequency, -Kingdom) %>%
  ggplot(aes(Kingdom, Frequency))+
  geom_boxplot(outlier.size = 0.5, fill = "steelblue3") +
  scale_y_continuous(trans = "log10", limits = c(0.0001,0.3)) +
  xlab("") +
  ylab("Frequency") +
  theme(axis.text.x = element_text(size=11)) +
  facet_wrap(. ~ Codon) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1))

# Boxplot: cysteine codon usage stratified by kingdom
train %>%
  select(Kingdom, UGU, UGC) %>%
  rowwise() %>% 
  mutate(Both = sum(c_across(UGU:UGC))) %>%
  gather(Codon, Frequency, -Kingdom) %>%
  ggplot(aes(Kingdom, Frequency, fill = Codon))+
  geom_boxplot(outlier.size = 0.6) +
  scale_y_continuous(trans = "log10", limits = c(0.0001,0.3)) +
  xlab("") +
  ylab("Frequency of cysteine codons") +
  theme(axis.text.x = element_text(size=11))


##############
# Final model
##############

# Training of the final model on the training set
finalmodel <- train(Kingdom ~ ., data=train, method="svmRadial",
                    trControl = trainControl(method="none"), # no cv
                    tuneGrid = data.frame(sigma = svmradmodel$bestTune$sigma,
                                          C = svmradmodel$bestTune$C)) # fixed parameters

# Prediction of taxa in the validation set
predicted_classes <- predict(finalmodel, validation)

# Accuracy of the final model
acc_final <- confusionMatrix(predicted_classes, validation$Kingdom)$overall["Accuracy"]

# Accuracy table
acc_results <- bind_rows(acc_results,
                         tibble(Method="Final Model on Validation Set",
                                Accuracy = acc_final))
acc_results %>% kable()

# Metrics
confusionMatrix(predicted_classes, validation$Kingdom)$byClass[,c("Sensitivity","Specificity", "Prevalence")] %>% kable()

# Confusion matrix
confusionMatrix(predicted_classes, validation$Kingdom)$table %>% kable()

# Boxplot: Comparison of codon usage between bacteria and plasmids
train %>%
  gather(Codon, Frequency, -Kingdom) %>%
  filter(Kingdom %in% c("bct", "plm")) %>%
  mutate(Codon = reorder(Codon, desc(Frequency))) %>%
  ggplot(aes(Codon, Frequency, fill = Kingdom)) +
  geom_boxplot(outlier.size = 0.4) +
  scale_y_continuous(trans = "log10", limits = c(0.0001,0.3)) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1), legend.position = "top") +
  scale_fill_discrete(name = "Taxon")
  
