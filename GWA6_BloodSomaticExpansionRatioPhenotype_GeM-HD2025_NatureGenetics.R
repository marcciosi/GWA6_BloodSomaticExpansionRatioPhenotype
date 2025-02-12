#This script processes the blood somatic expansion data from GeM-HD 2025 Nature Genetics
#The relevant data can be requested as "GWA123456 publication data" from CHDI
#the specific data file processed by the following script is 3999_individuals_GWA6_SER10_and_vars_associated_with_raw_phenotype.xlsx

#path to 3999_individuals_GWA6_SER10_and_vars_associated_with_raw_phenotype.xlsx
pathToData <- "PATH/TO/3999_individuals_GWA6_SER10_mapping_table_and_vars_associated_with_raw_phenotype.xlsx"
 

#load blood somatic expansion data for GWA6
GWA6_BloodSERdata <- readxl::read_excel(pathToData,
                                     sheet="Sheet1", na=c("")
                                                )
#turn variables alleleStructure1 and alleleStructureDesignation_in_paper into factors
GWA6_BloodSERdata$alleleStructure1 <- as.factor(GWA6_BloodSERdata$alleleStructure1)
GWA6_BloodSERdata$alleleStructureDesignation_in_paper <- as.factor(GWA6_BloodSERdata$alleleStructureDesignation_in_paper)

#fit Multiple linear regression "SERmodel1: ln(SER) ~ beta0 + beta1.CAG + beta2.age + beta
SERmodel1 <- lm(log(SER10) ~ centered.CAG1 + centered.Age,
                data = GWA6_BloodSERdata, na.action=na.exclude)

#fit Multiple linear regression "SERmodel2: ln(SER) ~ beta0 + beta1.CAG + beta2.age + beta3.(CAG x age) + beta"
SERmodel2 <- lm(log(SER10) ~ centered.CAG1*centered.Age,
                data = GWA6_BloodSERdata, na.action=na.exclude)

#fit Multiple linear regression "SERmodel3: ln(SER) ~ beta0 + beta1.CAG + beta2.age + beta3.(CAG x age) + beta4.CAG2 + beta5.age2 + beta"
SERmodel3 <- lm(log(SER10) ~ centered.CAG1*centered.Age + I(centered.CAG1^2) + I(centered.Age^2),
                data = GWA6_BloodSERdata, na.action=na.exclude)

#fit Multiple linear regression "SERmodel4: ln(SER) ~ beta0 + beta1.CAG + beta2.age + beta3.(CAG x age) + beta4.CAG2 + beta5.age2 + beta6.(CAG2 x age) + beta7.(CAG x age2) + beta"
SERmodel4 <- lm(log(SER10) ~ centered.CAG1*centered.Age + I(centered.CAG1^2) + I(centered.Age^2)
                + I(centered.CAG1^2):I(centered.Age) + centered.CAG1:I(centered.Age^2),
                data = GWA6_BloodSERdata, na.action=na.exclude)

#fit Multiple linear regression "SERmodel5: ln(SER) ~ beta0 + beta1.CAG + beta2.age + beta3.(CAG x age) + beta4.CAG2 + beta5.age2 + beta6.(CAG2 x age) + beta7.(CAG x age2) + beta8.(CAG2 x age2) + beta"
SERmodel5 <- lm(log(SER10) ~ centered.CAG1*centered.Age + I(centered.CAG1^2) + I(centered.Age^2)
                + I(centered.CAG1^2):I(centered.Age) + centered.CAG1:I(centered.Age^2)
                + I(centered.CAG1^2):I(centered.Age^2),
                data = GWA6_BloodSERdata, na.action=na.exclude)


#output AIC and adjusted R-squared for each of the 5 SERmodels above
AIC(SERmodel1)
AIC(SERmodel2)
AIC(SERmodel3)
AIC(SERmodel4)
AIC(SERmodel5)
summary(SERmodel1)$adj.r.squared
summary(SERmodel2)$adj.r.squared
summary(SERmodel3)$adj.r.squared
summary(SERmodel4)$adj.r.squared
summary(SERmodel5)$adj.r.squared


#fit Multiple linear regression "SERmodel4final: ln(SER) ~ beta0 + beta1.CAG + beta2.age + beta3.(CAG x age) + beta4.CAG2 + beta5.age2 + beta6.(CAG2 x age) + beta7.(CAG x age2)
#                                           + beta8.CCG + beta9.(CAG repeat-adjacent sequence) + beta10.PCRbatch + beta"
SERmodel4final <- lm(log(SER10) ~ centered.CAG1*centered.Age + I(centered.CAG1^2) + I(centered.Age^2)
                                            + I(centered.CAG1^2):I(centered.Age) + centered.CAG1:I(centered.Age^2)
                                            + CCG1 + (relevel(alleleStructure1, ref = "[CAG]n[CAA]1[CAG]1[CCG]1[CCA]1[CCG]n")) + as.factor(PCRbatch),
                                            data = GWA6_BloodSERdata, na.action=na.exclude)

#estimate effect of main "CAG repeat-adjacent sequences" on SER10 for Table1
summary(SERmodel4final)

#residual variation of SERmodel4final was used as the blood SER phenotype for the blood somatic expansion GWAS:
GWA6_BloodSERdata$bloodSERPhenotypeForTheSomaticExpansionGWAS <- residuals(SERmodel4final)

##################################################FIGURE 2A###################################
#make Figure2A

##### get predicted data for SER10 ######
#1) work out at what age the predicted values of SER10 reach an asymptote for any given value of CAG
# Given means
mean_CAG1 <- mean(GWA6_BloodSERdata$CAG1_MiSeqGla_Final)
mean_Age <- mean(GWA6_BloodSERdata$Age_at_sample_collection)

# Define the range of CAG values to evaluate
cagRange <- 40:50
# Create an empty list to store results
dataDummy <- list()

beta <- coef(SERmodel4final)
# Loop through each CAG value
for (actual_CAG1 in cagRange) {
  
  # Convert actual_CAG1 to centered value
  centered_CAG1 <- actual_CAG1 - mean_CAG1
  
  # Compute the centered Age asymptote
  denominator <- 2 * (beta["I(centered.Age^2)"] + beta["centered.CAG1:I(centered.Age^2)"] * centered_CAG1)
  
  if (denominator == 0) {
    warning(paste("Denominator is zero for CAG =", actual_CAG1, ". Check model fit or parameter values."))
    actual_Age_asymptote <- NA
  } else {
    centered_Age_asymptote <- - (beta["centered.Age"] + beta["centered.CAG1:centered.Age"] * centered_CAG1 + 
                                   beta["I(centered.CAG1^2):I(centered.Age)"] * centered_CAG1^2) / denominator
    actual_Age_asymptote <- centered_Age_asymptote + mean_Age
  }
  
  # Define age range dynamically
  # Define age range dynamically, capping at 92
  max_age <- ifelse(is.na(actual_Age_asymptote), 92, min(floor(actual_Age_asymptote), 92))# Use floor to ensure integer values
  ageRange <- 0:max_age

  # Store results for this CAG
  dataDummy[[as.character(actual_CAG1)]] <- data.frame(Age = ageRange, CAG = actual_CAG1)
}

dataDummy_df <- do.call(rbind, dataDummy) # Combine into a single dataframe
dataDummy_df <- data.frame("CAG1" = dataDummy_df$CAG, "Age_at_sample_collection" = dataDummy_df$Age) # rename columns
head(dataDummy_df,10)

#center CAG1 and Age_at_sampl_collection in dataDummy_df around the mean observed values for CAG1 and age_at_collection
dataDummy_df$centered.CAG1 <- dataDummy_df$CAG1 - mean(GWA6_BloodSERdata$CAG1_MiSeqGla_Final)
dataDummy_df$centered.Age <- dataDummy_df$Age_at_sample_collection - mean(GWA6_BloodSERdata$Age_at_sample_collection)
head(dataDummy_df,10)

# predict SM for the dataDummy_df values based on SERmodel4 from the real data
dataDummy_df$logSER_SERmodel4 <- predict(SERmodel4, dataDummy_df)
dataDummy_df$SER_SERmodel4 <- exp(dataDummy_df$logSER_SERmodel4)

dataDummy_df$CAG1.factor <- as.factor(dataDummy_df$CAG1)
GWA6_BloodSERdata$CAG1.factor <- as.factor(GWA6_BloodSERdata$CAG1_MiSeqGla_Final)


#for the 03/02/2025 final version of the paper
library(ggplot2)
Figure2A <-
  ggplot() + 
  geom_point(data=GWA6_BloodSERdata, mapping=aes(x=GWA6_BloodSERdata$Age_at_sample_collection, y=GWA6_BloodSERdata$SER10, color=GWA6_BloodSERdata$CAG1.factor), na.rm = TRUE)+
  geom_smooth(data=dataDummy_df,
              aes(x=dataDummy_df$Age_at_sample_collection, y=dataDummy_df$SER_SERmodel4, color = dataDummy_df$CAG1.factor), se = FALSE, na.rm = TRUE, method="loess") +
  xlab("age at sampling (years)") + 
  ylab("somatic expansion ratio") +
  labs(colour = "CAG", fill = "CAG") + 
  theme_bw(base_size = 14, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(.065,.735), #note: the above legend position is the one that works for 7x7inches pdf export
    axis.text = element_text(colour = "black", face = "bold"),
    panel.border = element_blank(), # to only have a border on the left and at the bottom, remove all borders first
    axis.line.x = element_line(colour = "black", size = 0.8),  # Add bottom border
    axis.line.y = element_line(colour = "black", size = 0.8))+   # Add left border  
  coord_cartesian(xlim =c(0, 80), ylim = c(0, 1.9))+ # use that instead of xlim and ylim to mask CAG50 dot at SER~2.5
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  guides(colour = guide_legend(reverse=T), fill = guide_legend(reverse=T)) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#B15928"))


#make Figure2B
#
SERmodel4final_withoutAllelestructures <- lm(log(SER10) ~ centered.CAG1*centered.Age + I(centered.CAG1^2) + I(centered.Age^2)
                                             + I(centered.CAG1^2):I(centered.Age) + centered.CAG1:I(centered.Age^2)
                                             + CCG1 + as.factor(PCRbatch),
                                             data = GWA6_BloodSERdata, na.action=na.exclude)
GWA6_BloodSERdata$AdjusterSERforFigure2B <- residuals(SERmodel4final_withoutAllelestructures)


#1) exclude 3 individuals with rare "CAG repeat-adjacent sequences"
GWA6_BloodSERdata_forFigure1B <- subset(GWA6_BloodSERdata, alleleStructureDesignation_in_paper != "other non-canonical variant mentioned in Table 1")
#drop level other non-canonical variant mentioned in Table 1" for variable alleleStructureDesignation_in_paper as that level is no longer present in GWA6_BloodSERdata_forFigure1B
GWA6_BloodSERdata_forFigure1B <- droplevels(GWA6_BloodSERdata_forFigure1B)
#check that the level "other non-canonical variant mentioned in Table 1" was dropped
levels(GWA6_BloodSERdata_forFigure1B$alleleStructureDesignation_in_paper)

# first reorder factor levels into order for figure2B into a new variable by putting canonical and CAACAG-dup first:
GWA6_BloodSERdata_forFigure1B$alleleStructureDesignation_in_paperFig2Border <- forcats::fct_relevel(GWA6_BloodSERdata_forFigure1B$alleleStructureDesignation_in_paper, "canonical", "CAACAG-dup", "CAA/CCA-loss")


#2) now make Figure2B
Figure2B <-
  ggplot(data = GWA6_BloodSERdata_forFigure1B,
         aes(alleleStructureDesignation_in_paperFig2Border, AdjusterSERforFigure2B, colour = alleleStructureDesignation_in_paperFig2Border)) +
  ggforce::geom_sina() +
  geom_boxplot(outlier.shape = NA, alpha = 0, colour = "black", lwd=0.7) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(), 
    legend.position="none",
    axis.text = element_text(colour = "black", face = "bold"),
    panel.border = element_blank(), # to only have a border on the left and at the bottom, remove all borders first
    axis.line.x = element_line(colour = "black", size = 0.8),  # Add bottom border
    axis.line.y = element_line(colour = "black", size = 0.8))+   # Add left border  
    labs(y = "Adjusted somatic expansion ratio", x = "")+
  scale_x_discrete(labels=c("Canonical", "CAACAG-dup","CAA/CCA-loss","CAA-loss","CCA-loss")) +
  scale_color_manual(values=c("#5c5c5c", "#5c5c5c", "#5c5c5c", "#5c5c5c", "#5c5c5c"))






