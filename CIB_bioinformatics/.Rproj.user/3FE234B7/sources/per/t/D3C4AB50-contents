############################################################################################
# This is a modified version of LOD-calculator.R (https://github.com/cmerkes/qPCR_LOD_Calc)
# Updated by Meghan M. Shea (April 2026) to accommodate ddPCR data as well as qPCR data
# Coding assistance provided by claude-4-5-sonnet accessed via Stanford AI Playground

# Other changes to the original script include:
# Allowing for user specified LOD threshold (instead of defaulting to 95%)
# Allowing for comparisons between different LOD/LOQ thresholds and ddPCR minimum positive droplets
# Saving all files + plots generated in a unique folder
############################################################################################


##READ BEFORE STARTING
############################################################################################
# Before running the code, you will need to modify the variables in the section below
# called "Items to Modify." This script includes interactive codes, and it is designed to be
# run using R studio, where you can click on "Source" in the upper right corner after modifying.

# This is a script to calculate Limit of Detection (LOQ) and Limit of Quantification (LOQ) for 
# eDNA assays using either qPCR or ddPCR. It requires you to assemble your data for many replicates
# of a dilution curve of low-copy standards into a single csv file with the following headings.

# For All Files:
# Target: The assay(s) being tested

# For qPCR:
# SQ (Starting Quantity): The expected values for your standards in copies/rxn
# Cq: qPCR quantification cycle

# For ddPCR: 
# SQ (Starting Quantity): The expected values for your standards in copies/uL (e.g. the concentration
# of gBlocks used in your serial dilution); this will be converted to copies/rxn by multiplying by
# the volume of template used. 
# Concentration: The measured concentration in copies/uL; this will be converted to copies/rxn by 
# multipying by the reaction volume. 
# TemplateVolume: volume of DNA template added to reaction
# ReactionVolume: volume of reaction
# PositiveDroplets: number of positive droplets

# Additional Notes: 
# - SQ values for NTCs should be NA, not 0.
# - Headings are not case-sensitive, but must be spelled the same.
# - Additional columns may be included for readability but are not used in calculations. 

# CONCEPTUAL NOTE: I have only added the template volume + reaction volume for ddPCR, so the qPCR 
# aspects of the script remain unmodified. But it probably makes sense to include for both.

##ITEMS TO MODIFY
############################################################################################

library(here)
## Load file
qPCR_AK_data <- readRDS(file=here("data","qPCR_results","data","qPCR_results","qPCR_AK_data.RDS"))

## There are 6 plates
qPCR_AK_data$Plate %>% unique()

## There are 4 targets (Chinook, Coho, Eulachon, cetacea) and there 2 positives (MOA and Haast)
qPCR_AK_data$`Target Name` %>% unique()

qPCR_AK_data %>% 
  filter(., str_detect(Task,"STANDARD")) %>% 
  dplyr::select(Plate, Target=`Target Name`, CT, Quantity) %>% 
  mutate(., Quantity=as.numeric(Quantity),
         CT=as.numeric(CT)) -> standards_all

library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(purrr)

plot_all_standards <- function(data) {
  
  # 1. Extract all unique combinations of Plate and Target present in the data
  combos <- data %>% 
    distinct(Plate, Target)
  
  # 2. Iterate over each unique combination to create a plot
  plot_list <- map2(combos$Plate, combos$Target, function(curr_plate, curr_target) {
    
    # Filter for the specific Plate and Target
    sub_df <- data %>% 
      filter(Plate == curr_plate, Target == curr_target)
    
    # Generate the ggplot
    p <- sub_df %>% 
      ggplot(aes(y = CT, x = Quantity)) + 
      geom_point() + 
      ggtitle(paste("Plate =", curr_plate, "Target =", curr_target)) +
      scale_x_log10(labels = label_log()) +
      geom_smooth(method = "lm", se = FALSE) +
      stat_cor(aes(label = after_stat(rr.label)), label.x = 1, label.y = 25)
    
    return(p)
  })
  
  # 3. Name each list element (e.g., "20m_1a_Chinook") for easy access later
  names(plot_list) <- paste(combos$Plate, combos$Target, sep = "_")
  
  return(plot_list)
}


standards_all %>% 
  filter(., Target %in% c("Chinook", "Coho","Eulachon","cetacea")) -> standards_targets

plot_all_standards(standards_targets) -> plot_all_standards_plots


## Define your CV threshold(s) for LoQ:
LOQ.Thresholds <- c(0.35)  # Vector of CV thresholds to compare

## Define your detection rate threshold(s) for LoD (e.g., 0.95 for 95% detection):
LOD.Thresholds <- c(0.95)  # Vector of detection thresholds to compare

## Define which logarithmic function to use for LoD model:
LOD.FCT <- "Best"
## Selecting "Best" will signal the code to automatically select the best fitting
##   model choice. Run the function getMeanFunctions() to print the list of all choices.
## Example: LOD.FCT <- W2.4()
##   This example will use the Weibull type II, 4 parameter function.

## Define which model to use for LoQ model:
LOQ.FCT <- "Best"
## Selecting "Best" will signal the code to automatically select the model with lowest
##   residual standard error. Change to "Decay" to use exponential decay model, "Linear"
##   to use linear model, "Pn" to use an nth-order polynomial model where n is numerical.
##   Example: "P2" will use a 2nd order polynomial model, "P3" will use 3rd order, and etc.
##   Selecting "Best" will test polynomial models up to 6th-order.

## Define your assay type:
ASSAY.TYPE <- "qPCR"  # Options: "qPCR" or "ddPCR"

## Define minimum positive droplets for ddPCR detection:
MIN.POSITIVE.DROPLETS.OPTIONS <- c(1, 2, 3)  # Vector of minimum droplet options to compare

############################################################################################


## Load packages:
library(ggplot2)
library(drc)
library(dplyr)

## Read and prepare data ONCE before loops:
DAT_RAW <- standards_targets
Targets <- unique(DAT_RAW$Target)
Plates <- unique(DAT_RAW$Plate)
i <-1

output_folder <- here("analysis","LoD")

## Check the data:
## For qPCR, expect: Plate, Target, Cq, SQ, Reaction Volume, Template Volume
## For ddPCR, expect: Plate, Target, Concentration, SQ (or Expected_Conc), Reaction Volume, Template Volume


if(sum(colnames(DAT_RAW)=="Target")!=1) { #Is there a "Target" column?
  A <- grep("target",colnames(DAT_RAW),ignore.case=TRUE)
  if(length(A)==1) { colnames(DAT_RAW)[A] <- "Target" } #Rename target column if it is mispelled but can be identified and there is only 1.
  if(length(A)!=1) { write("There is a problem with the 'Target' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) } #Add error message to analysis log.
  if(length(A)>1) { cat("ERROR: multiple 'Target' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
  if(length(A)==0) { print("ERROR: cannot detect 'Target' column.") }
}
DAT_RAW$Target <- as.factor(DAT_RAW$Target)


if(sum(colnames(DAT_RAW)=="SQ")!=1) { #Is there a "SQ" column?
  A <- grep("sq|copies|starting|quantity",colnames(DAT_RAW),ignore.case=TRUE)
  if(length(A)==1) { colnames(DAT_RAW)[A] <- "SQ" } #Rename SQ column if it is mispelled but can be identified and there is only 1.
  if(length(A)!=1) { write("There is a problem with the 'SQ' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) } #Add error message to analysis log.
  if(length(A)>1) { cat("ERROR: multiple 'SQ' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
  if(length(A)==0) { print("ERROR: cannot detect 'SQ' column.") }
}

DAT_RAW$SQ <- suppressWarnings(as.numeric(as.character(DAT_RAW$SQ))) #Non-numerical values (i.e. NTC) will be converted to NAs

if(ASSAY.TYPE == "qPCR") {
  if(sum(colnames(DAT_RAW)=="Cq")!=1) {
    A <- grep("cq|ct|cycle",colnames(DAT_RAW),ignore.case=TRUE)
    if(length(A)==1) { colnames(DAT_RAW)[A] <- "Cq" }
    if(length(A)!=1) { write("There is a problem with the 'Cq' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) }
    if(length(A)>1) { cat("ERROR: multiple 'Cq' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
    if(length(A)==0) { print("ERROR: cannot detect 'Cq' column.") }
  }
  DAT_RAW$Cq <- suppressWarnings(as.numeric(as.character(DAT_RAW$Cq)))
}

if(ASSAY.TYPE == "ddPCR") {
  if(sum(colnames(DAT_RAW)=="Concentration")!=1) {
    A <- grep("concentration|conc|copies",colnames(DAT_RAW),ignore.case=TRUE)
    if(length(A)==1) { colnames(DAT_RAW)[A] <- "Concentration" }
    if(length(A)!=1) { write("There is a problem with the 'Concentration' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) }
    if(length(A)>1) { cat("ERROR: multiple 'Concentration' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
    if(length(A)==0) { print("ERROR: cannot detect 'Concentration' column.") }
  }
  DAT_RAW$Concentration <- suppressWarnings(as.numeric(as.character(DAT_RAW$Concentration)))
  
  if(sum(colnames(DAT_RAW)=="TemplateVolume")!=1) {
    A <- grep("template|gBlock",colnames(DAT_RAW),ignore.case=TRUE)
    if(length(A)==1) { colnames(DAT_RAW)[A] <- "TemplateVolume" }
    if(length(A)!=1) { write("There is a problem with the 'TemplateVolume' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) }
    if(length(A)>1) { cat("ERROR: multiple 'TemplateVolume' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
    if(length(A)==0) { print("ERROR: cannot detect 'TemplateVolume' column.") }
  }
  DAT_RAW$TemplateVolume <- suppressWarnings(as.numeric(as.character(DAT_RAW$TemplateVolume)))
  
  if(sum(colnames(DAT_RAW)=="ReactionVolume")!=1) {
    A <- grep("reaction|rxn",colnames(DAT_RAW),ignore.case=TRUE)
    if(length(A)==1) { colnames(DAT_RAW)[A] <- "ReactionVolume" }
    if(length(A)!=1) { write("There is a problem with the 'ReactionVolume' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) }
    if(length(A)>1) { cat("ERROR: multiple 'ReactionVolume' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
    if(length(A)==0) { print("ERROR: cannot detect 'ReactionVolume' column.") }
  }
  DAT_RAW$ReactionVolume <- suppressWarnings(as.numeric(as.character(DAT_RAW$ReactionVolume)))
  
  if(sum(colnames(DAT_RAW)=="PositiveDroplets")!=1) {
    A <- grep("positive|pos.*drop",colnames(DAT_RAW),ignore.case=TRUE)
    if(length(A)==1) { colnames(DAT_RAW)[A] <- "PositiveDroplets" }
    if(length(A)!=1) { write("There is a problem with the 'PositiveDroplets' column.\n\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE) }
    if(length(A)>1) { cat("ERROR: multiple 'PositiveDroplets' columns detected.",colnames(DAT_RAW)[A],sep="\n") }
    if(length(A)==0) { print("ERROR: cannot detect 'PositiveDroplets' column (required when MIN.POSITIVE.DROPLETS > 1).") }
  }
  DAT_RAW$PositiveDroplets <- suppressWarnings(as.numeric(as.character(DAT_RAW$PositiveDroplets)))
  
  ##CHECK THAT THESE CALCULATIONS ARE RIGHT
  DAT_RAW$Concentration <- DAT_RAW$Concentration*DAT_RAW$ReactionVolume
  DAT_RAW$SQ <- DAT_RAW$SQ*DAT_RAW$TemplateVolume
} 


#if(sum(is.na(DAT_RAW$SQ))>0) {
#write(paste0("WARNING: ",sum(is.na(DAT_RAW$SQ))," DAT_RAWa points excluded without a valid starting quantity (SQ)!\nHere is a sample of the DAT_RAWa being excluded:\n"),
# file.path(output_folder, "Analysis Log.txt"),append=TRUE)
# suppressWarnings(write.table(head(DAT_RAW[is.na(DAT_RAW$SQ),]),file.path(output_folder, "Analysis Log.txt"),append=TRUE,
#                              sep="\t",eol="\n",row.names=FALSE,col.names=TRUE))
# write("\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE)
# print(paste0("WARNING: ",sum(is.na(DAT_RAW$SQ))," DAT_RAWa points excluded without a valid starting quantity (SQ)!"))
## print(head(DAT_RAW[is.na(DAT_RAW$SQ),]))
#}

## Check for wild outliers that the user should go back and review:


## Get matchups of all standards and markers used:
for (k in 1:length(Plates)) {
  for(i in 1:length(Targets)) {
    if(i==1) {
      Standards <- unique(DAT_RAW$SQ[DAT_RAW$Target==Targets[i]&!is.na(DAT_RAW$SQ)&DAT_RAW$Plate==Plates[k]])
      Target <- rep(as.character(Targets[i]),length(Standards))
      Plate <- rep(as.character(Plates[k]),length(Standards))
    }
    else {
      Standards <- c(Standards,unique(DAT_RAW$SQ[DAT_RAW$Target==Targets[i]&!is.na(DAT_RAW$SQ)&DAT_RAW$Plate==Plates[k]]))
      Target <- c(Target,rep(as.character(Targets[i]),
                             length(unique(DAT_RAW$SQ[DAT_RAW$Target==Targets[i]&!is.na(DAT_RAW$SQ)&DAT_RAW$Plate==Plates[k]]))))
      Plate <- c(Plate,rep(as.character(Plates[k]),
                            length(unique(DAT_RAW$SQ[DAT_RAW$Target==Targets[i]&!is.na(DAT_RAW$SQ)&DAT_RAW$Plate==Plates[k]]))))
    }
  }
}

OUTS <- data.frame(Plate, Target=Target,Standard=Standards,Outliers=NA)

## Identify any wells where the Cq value is more than 10% away from the median for
##   that standard.
for(i in 1:nrow(OUTS)) {
  MED <- median(DAT_RAW$Cq[DAT_RAW$SQ==OUTS$Standard[i]&DAT_RAW$Target==OUTS$Target[i]&DAT_RAW$Plate==OUTS$Plate[i]],na.rm=TRUE)
  A <- which(DAT_RAW$SQ==OUTS$Standard[i]&DAT_RAW$Target==OUTS$Target[i]&DAT_RAW$Plate==OUTS$Plate[i]&DAT_RAW$Cq<0.9*MED&!is.na(DAT_RAW$Cq))
  B <- which(DAT_RAW$SQ==OUTS$Standard[i]&DAT_RAW$Target==OUTS$Target[i]&DAT_RAW$Plate==OUTS$Plate[i]&DAT_RAW$Cq>1.1*MED&!is.na(DAT_RAW$Cq))
  if(length(c(A,B))>0) {
    OUTS$Outliers[i] <- paste(c(A,B),collapse=",")
  }
}

## If any outliers are detected, export the raw DAT_RAWa as csv and make a note in
##   the analysis log.
if(sum(!is.na(OUTS$Outliers))>0) {
  OUT.ROW <- paste(OUTS$Outliers[!is.na(OUTS$Outliers)],collapse=",")
  OUT.ROW2 <- unlist(strsplit(OUT.ROW,split=","))
  write.csv(DAT_RAW[OUT.ROW2,], file=file.path(output_folder, "Potential-Outliers.csv"), row.names=FALSE)
  write("Potential outliers have been detected. Please review the DAT_RAWa exported as
Potential-Outliers.csv, and determine if any DAT_RAWa points need to be excluded
or adjusted due to false positives or poorly normalized baselines.",
        file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
  write("\n",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
}


## Create master output folder with timestamp:
input_filename <- tools::file_path_sans_ext(basename("CIB_standards"))
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
master_output_folder <- here(output_folder,paste0(input_filename, "_comparison_", timestamp))
dir.create(master_output_folder, showWarnings = FALSE)

## Create comparison summary dataframe:
comparison_summary <- data.frame()

## Set MIN.POSITIVE.DROPLETS.OPTIONS based on assay type:
if(ASSAY.TYPE == "qPCR") {
  MIN.POSITIVE.DROPLETS.OPTIONS <- c(NA)  # Not applicable for qPCR
  cat("Note: Minimum positive droplets parameter is only applicable to ddPCR.\n\n")
} else {
  cat("Analyzing ddPCR data with minimum positive droplet thresholds:", 
      paste(MIN.POSITIVE.DROPLETS.OPTIONS, collapse=", "), "\n\n")
}

## Loop through all parameter combinations:
for(LOD.Threshold in LOD.Thresholds) {
  for(LOQ.Threshold in LOQ.Thresholds) {
    for(MIN.POSITIVE.DROPLETS in MIN.POSITIVE.DROPLETS.OPTIONS) {
      
      ## Create a working copy of the prepared data for this iteration:
      DAT <- DAT_RAW
      
      ## Create subfolder for this parameter combination:
      if(ASSAY.TYPE == "qPCR") {
        param_label <- paste0("LOD", LOD.Threshold*100, "_LOQ", LOQ.Threshold*100)
      } else {
        param_label <- paste0("LOD", LOD.Threshold*100, "_LOQ", LOQ.Threshold*100, "_MinPos", MIN.POSITIVE.DROPLETS)
      }
      output_folder <- file.path(master_output_folder, param_label)
      dir.create(output_folder, showWarnings = FALSE)
      
      ## Create an analysis log file:
      write(paste0("Analysis started: ", date(), "\n"), 
            file=file.path(output_folder, "Analysis Log.txt"))
      write(paste0("LOD Threshold: ", LOD.Threshold*100, "%\n"),
            file=file.path(output_folder, "Analysis Log.txt"), append=TRUE)
      write(paste0("LOQ CV Threshold: ", LOQ.Threshold, "\n"),
            file=file.path(output_folder, "Analysis Log.txt"), append=TRUE)
      if(ASSAY.TYPE == "ddPCR") {
        write(paste0("Minimum Positive Droplets: ", MIN.POSITIVE.DROPLETS, "\n\n"),
              file=file.path(output_folder, "Analysis Log.txt"), append=TRUE)
      } else {
        write("\n", file=file.path(output_folder, "Analysis Log.txt"), append=TRUE)
      }
      
      ## Apply minimum positive droplet threshold (ddPCR only, PARAMETER-DEPENDENT):
      if(ASSAY.TYPE == "ddPCR" && MIN.POSITIVE.DROPLETS > 1) {
        rows_filtered <- sum(!is.na(DAT$PositiveDroplets) & DAT$PositiveDroplets < MIN.POSITIVE.DROPLETS & !is.na(DAT$Concentration) & DAT$Concentration > 0)
        if(rows_filtered > 0) {
          DAT$Concentration[!is.na(DAT$PositiveDroplets) & DAT$PositiveDroplets < MIN.POSITIVE.DROPLETS] <- 0
          write(paste0("INFO: ",rows_filtered," data points set to 0 concentration due to having fewer than ",MIN.POSITIVE.DROPLETS," positive droplets.\n\n"),
                file.path(output_folder, "Analysis Log.txt"),append=TRUE)
          print(paste0("INFO: ",rows_filtered," data points filtered for insufficient positive droplets (threshold: ",MIN.POSITIVE.DROPLETS,")"))
        } else {
          write(paste0("INFO: No data points filtered. All detections had ",MIN.POSITIVE.DROPLETS," or more positive droplets.\n\n"),
                file.path(output_folder, "Analysis Log.txt"),append=TRUE)
        }
      } 
      
      ## Log NA values for this iteration:
      if(sum(is.na(DAT$SQ))>0) {
        write(paste0("WARNING: ",sum(is.na(DAT$SQ))," data points excluded without a valid starting quantity (SQ)!\nHere is a sample of the data being excluded:\n"),
              file.path(output_folder, "Analysis Log.txt"),append=TRUE)
        suppressWarnings(write.table(head(DAT[is.na(DAT$SQ),]),file.path(output_folder, "Analysis Log.txt"),append=TRUE,
                                     sep="\t",eol="\n",row.names=FALSE,col.names=TRUE))
        write("\n",file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      }
      
      ## Log outliers if they exist (from pre-loop analysis):
      if(ASSAY.TYPE == "qPCR" && exists("OUTS") && sum(!is.na(OUTS$Outliers))>0) {
        OUT.ROW <- paste(OUTS$Outliers[!is.na(OUTS$Outliers)],collapse=",")
        OUT.ROW2 <- unlist(strsplit(OUT.ROW,split=","))
        write.csv(DAT[OUT.ROW2,], file=file.path(output_folder, "Potential-Outliers.csv"), row.names=FALSE)
        write("Potential outliers have been detected. Please review the data exported as
Potential-Outliers.csv, and determine if any data points need to be excluded
or adjusted due to false positives or poorly normalized baselines.\n\n",
              file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      }
      
      
      
      ## Generate standard curves (qPCR) or linear models (ddPCR) using all data:
      curve.list <- ""
      DAT$Copy.Estimate <- rep(NA,nrow(DAT))
      DAT$Mod <- rep(0,nrow(DAT))
      
      if(ASSAY.TYPE == "qPCR") {
        for(k in 1:length(Plates)) {
          Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
        for(i in 1:length(Targets)) {
          STDS <- data.frame(S=unique(DAT$SQ[DAT$Target==Targets[i]&DAT$Plate==Plates[k]]),R=NA)
          ## Calculate detection rates for each standard:
          for(j in 1:nrow(STDS)) {
            STDS$R[j] <- sum(!is.na(DAT$Cq)&DAT$SQ==STDS$S[j]&DAT$Target==Targets[i]&DAT$Plate==Plates[k],na.rm=TRUE)/sum(DAT$SQ==STDS$S[j]&DAT$Target==Targets[i]&DAT$Plate==Plates[k],na.rm=TRUE)
          }
          ## Only use standards with 50% or greater detection rates for linear regression:
          if(sum(STDS$R>=0.5,na.rm=TRUE)>2) {
            STDS2 <- STDS$S[STDS$R>=0.5&!is.na(STDS$R)&!is.na(STDS$S)]
          }
          ## If there are not at least 3 standards with 50% or greater detection, use the top 3:
          if(sum(STDS$R>=0.5,na.rm=TRUE)<3) {
            STDS2 <- STDS$S[order(STDS$R,decreasing=TRUE)][1:3]
          }
          ## Identify the 2nd and 3rd quartiles of each used standard for inclusion in the
          ##   standard curve calculations
          for(j in 1:length(STDS2)) {
            D <- DAT$Cq[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&DAT$SQ==STDS2[j]]
            DAT$Mod[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&DAT$SQ==STDS2[j]&DAT$Cq>=quantile(D,na.rm=TRUE)[2]&DAT$Cq<=quantile(D,na.rm=TRUE)[4]&!is.na(DAT$SQ)] <- 1
          }
          if(length(unique(DAT$SQ[DAT$Target==Targets[i]&DAT$Plate==Plates[k]]))!=length(STDS2)) {
            ToWrite <- paste0("These standards not included in Plate:",Plates[k],"for Targets:",Targets[i],
                              " standard curve regression for copy estimate calculations, because they detected below 50%: ",
                              paste(setdiff(unique(DAT$SQ[DAT$Target==Targets[i]&DAT$Plate==Plates[k]]),STDS2),collapse=", "),"\n\n")
            write(ToWrite,file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
          }
          assign(paste0("curve_",k,"_",i),lm(Cq~log10(SQ),data=DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&DAT$Mod==1,]))
          curve.list <- c(curve.list,paste0("curve_",k,"_",i))
          Intercept <- coef(get(curve.list[i+1]))[1]
          Slope <- coef(get(curve.list[i+1]))[2]
          DAT$Copy.Estimate[DAT$Target==Targets[i]&DAT$Plate==Plates[k]] <- 10^((DAT$Cq[DAT$Target==Targets[i]&DAT$Plate==Plates[k]]-Intercept)/Slope)
        }
        
        }
      }
      
      if(ASSAY.TYPE == "ddPCR") {
        ## For ddPCR, concentration IS the copy estimate
        DAT$Copy.Estimate <- DAT$Concentration
        DAT$Mod <- rep(1, nrow(DAT))  # All data points are "modeled" in ddPCR
        
        ## Generate linear regression for expected vs measured concentrations
        for(i in 1:length(Plates)) {
          Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
        for(i in 1:length(Targets)) {
          # Filter for valid data points (non-NA and detected)
          valid_data <- DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k] & !is.na(DAT$SQ) & !is.na(DAT$Concentration) & DAT$Concentration > 0, ]
          
          if(nrow(valid_data) > 0) {
            assign(paste0("curve_",k,"_",i), lm(Concentration ~ SQ, data=valid_data))
            curve.list <- c(curve.list, paste0("curve_",k,"_",i))
          } else {
            # If no valid data, create a placeholder
            curve.list <- c(curve.list, NA)
            write(paste0("WARNING: No valid data points for linear regression of ", Targets[i], "\n\n"),
                  file=file.path(output_folder, "Analysis Log.txt"), append=TRUE)
          }
        }
        }
      }
      
      plater1 <- DAT %>% dplyr::select(Standards=SQ, Plate, Target) %>% distinct()
      DAT2_holder <- data.frame(Reps=NA,Detects=NA,
                         Copy.CV=NA)
      DAT2 <- cross_join(plater1,DAT2_holder)
      
      if(ASSAY.TYPE == "qPCR") {
        DAT2$Cq.mean <- NA
        DAT2$Cq.sd <- NA
        DAT2$Cq.CV <- NA
      }
      
      if(ASSAY.TYPE == "ddPCR") {
        DAT2$Conc.mean <- NA
        DAT2$Conc.sd <- NA
      }
      
      for(i in 1:nrow(DAT2)) {
        DAT2$Reps[i] <- sum(DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i],na.rm=TRUE)
        
        if(ASSAY.TYPE == "qPCR") {
          DAT2$Detects[i] <- sum(!is.na(DAT$Cq)&DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i],na.rm=TRUE)
          DAT2$Cq.mean[i] <- mean(DAT$Cq[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i]],na.rm=TRUE)
          DAT2$Cq.sd[i] <- sd(DAT$Cq[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i]],na.rm=TRUE)
          DAT2$Cq.CV[i] <- sqrt(2^(DAT2$Cq.sd[i]^2*log(2))-1)
        }
        
        if(ASSAY.TYPE == "ddPCR") {
          DAT2$Detects[i] <- sum(!is.na(DAT$Concentration)&DAT$Concentration>0&DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i],na.rm=TRUE)
          DAT2$Conc.mean[i] <- mean(DAT$Concentration[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i]],na.rm=TRUE)
          DAT2$Conc.sd[i] <- sd(DAT$Concentration[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i]],na.rm=TRUE)
        }
        
        DAT2$Copy.CV[i] <- sd(DAT$Copy.Estimate[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i]],na.rm=TRUE)/
          mean(DAT$Copy.Estimate[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]&DAT$Plate==DAT2$Plate[i]],na.rm=TRUE)
      }
      
      DAT2$Rate <- DAT2$Detects/DAT2$Reps
      
      
      if(ASSAY.TYPE == "qPCR") {
        DAT2$CV_metric <- DAT2$Cq.CV
        CV_label <- "Coefficient of variation for Cq-Values"
      } else {
        DAT2$CV_metric <- DAT2$Copy.CV
        CV_label <- "Coefficient of variation for Concentration Estimates"
      }
      
      
      write("Data Summary:",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      suppressWarnings(write.table(DAT2,file=file.path(output_folder, "Analysis Log.txt"),append=TRUE,sep="\t",eol="\n",
                                   row.names=FALSE,col.names=TRUE))
      
      ## Determine the lowest standard with user-defined threshold or greater detection:
      for(k in 1:length(Plates)) {
        Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
      for(i in 1:length(Targets)) {
        A <- min(DAT2$Standards[DAT2$Rate>=LOD.Threshold&DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])
        ToWrite <- paste0("For ",Plates[k]," and ",Targets[i],", the lowest standard with ",LOD.Threshold*100,"% or greater detection is: ",A," copies/reaction.")
        ToWrite2 <- ""
        if(length(which(DAT2$Rate<LOD.Threshold&DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]))>0) {
          B <- max(DAT2$Standards[DAT2$Rate<LOD.Threshold&DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])
          if(B>A) {
            ToWrite2 <- paste0("WARNING: For ",Plates[k]," and ",Targets[i],", ",B," copies/reaction standard detected at lower rate than ",A," copies/reaction standard.\nPlease retest.")
          }
        }
        if(length(which(DAT2$Rate<LOD.Threshold&DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]))==0) {
          ToWrite2 <- paste0("WARNING: LoD cannot be determined for ",Plates[k]," and ",Targets[i],", because it is lower than the lowest standard you tested.\nReport as <",A," copies/reaction, or retest with lower concentrations.")
        }
        write(paste0("\n\n",ToWrite,"\n"),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
        if(ToWrite2!="") { write(paste0(ToWrite2,"\n\n"),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE) }
        cat(ToWrite,ToWrite2,sep="\n")
      }
      }
      
      ## Determine LoD and LoQ by modeling, and summarize each assay:
      ## NOTE: LoD is now determined by dose-response modeling. Probit modeling code remains,
      ##         but has been converted to comments.
      
      if(ASSAY.TYPE == "qPCR") {
        DAT$Detect <- as.numeric(!is.na(DAT$Cq))
      }
      
      if(ASSAY.TYPE == "ddPCR") {
        DAT$Detect <- as.numeric(!is.na(DAT$Concentration) & DAT$Concentration > 0)
      }
      
      
      #LOD.list <- ""
      LOD.list2 <- ""
      LOD.list3 <- ""
      LOQ.list <- ""
      plater <- DAT2 %>% dplyr::select(Plate, Assay=Target) %>% distinct()
      dfer<-data.frame(R.squared=NA,Slope=NA,Intercept=NA,LOD.discrete=NA,
                 LOD=NA,LOQ=NA,rep2.LOD=NA,rep3.LOD=NA,rep4.LOD=NA,rep5.LOD=NA,rep8.LOD=NA,rep9.LOD=NA)
      DAT3 <- cross_join(plater, dfer)
      LOD.FCTS <- list(LL.2(),LL.3(),LL.3u(),LL.4(),LL.5(),W1.2(),W1.3(),W1.4(),W2.2(),W2.3(),
                       W2.4(),AR.2(),AR.3(),MM.2(),MM.3())
      for(k in 1:length(Plates)) {
        Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
      for(i in 1:length(Targets)) {
        ## Check input suitability for probit or dose-response modeling:
        if(sum(DAT2$Rate[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]!=1&DAT2$Rate[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]!=0)==0) {
          ToWrite <- paste0("WARNING: For ",Plates[k]," and ",Targets[i],", all standards detected fully or failed fully.  Therefore, the LoD model will not converge.")
          write(paste0(ToWrite,"\n\n"),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
          print(ToWrite)
        }
        if(sum(DAT2$Rate[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]!=1&DAT2$Rate[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]!=0)==1) {
          ToWrite <- paste0("WARNING: For ",Plates[k]," and ",Targets[i],", only 1 standard detected in the informative range (not 0% and not 100%).  Therefore, the LoD model results will be less reliable.")
          write(paste0(ToWrite,"\n\n"),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
          print(ToWrite)
        }
        ## Define probit model:
        #assign(paste0("LOD.mod",i),glm(Detect~SQ,data=DAT[DAT$Target==Targets[i],],
        #                               family=binomial(link="probit")))
        #LOD.list <- c(LOD.list,paste0("LOD.mod",i))
        ## Define LOQ model using lowest residual standard error selection:
        if(LOQ.FCT=="Best") {
          ## Remove previous marker LOQ models from environment if they exist:
          suppressWarnings(rm(LOQ1,LOQ2,LOQ3,LOQ4,LOQ5,LOQ6,LOQ7))
          tryCatch({ #skip if model cannot be determined.
            LOQ1 <- nls(CV_metric~SSasymp(log10(Standards),Asym,R0,lrc),
                        data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: decay LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          tryCatch({ #skip if model cannot be determined.
            LOQ2 <- lm(CV_metric~log10(Standards),data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: linear LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          tryCatch({ #skip if model cannot be determined.
            LOQ3 <- lm(CV_metric~poly(log10(Standards),2),data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: 2nd polynomial LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          tryCatch({ #skip if model cannot be determined.
            LOQ4 <- lm(CV_metric~poly(log10(Standards),3),data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: 3rd polynomial LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          tryCatch({ #skip if model cannot be determined.
            LOQ5 <- lm(CV_metric~poly(log10(Standards),4),data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: 4th polynomial LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          tryCatch({ #skip if model cannot be determined.
            LOQ6 <- lm(CV_metric~poly(log10(Standards),5),data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: 5th polynomial LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          tryCatch({ #skip if model cannot be determined.
            LOQ7 <- lm(CV_metric~poly(log10(Standards),6),data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],])
          }, error=function(e) {
            e
            cat("ERROR: 6th polynomial LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
          ## Determine which models were able to be determined:
          A <- sapply(c("LOQ1","LOQ2","LOQ3","LOQ4","LOQ5","LOQ6","LOQ7"),exists)
          B <- names(A)[A==TRUE]
          ## If at least 1 LOQ model was determined, select the one with the lowest
          ##   residual standard error:
          if(length(B)>0) {
            LOQ.res <- rep(NA,length(B))
            for(j in 1:length(B)) {
              LOQ.res[j] <- summary(get(B[j]))$sigma
            }
            C <- which(LOQ.res==min(LOQ.res,na.rm=TRUE))
            assign(paste0("LOQ.mod_",k,"_",i),get(B[C]))
            LOQ.list <- c(LOQ.list,paste0("LOQ.mod_",k,"_",i))
          }
        }
        ## Define LOQ model by exponential decay modeling:
        if(LOQ.FCT=="Decay") {
          tryCatch({ #skip if model cannot be determined.
            assign(paste0("LOQ.mod_",k,"_",i),nls(CV_metric~SSasymp(log10(Standards),Asym,R0,lrc),
                                           data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],]))
            LOQ.list <- c(LOQ.list,paste0("LOQ.mod_",k,"_",i))
          }, error=function(e) {
            e
            cat("ERROR: decay LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
        }
        ## Define LOQ model by linear modeling:
        if(LOQ.FCT=="Linear") {
          tryCatch({ #skip if model cannot be determined.
            assign(paste0("LOQ.mod_",i),lm(CV_metric~log10(Standards),
                                          data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],]))
            LOQ.list <- c(LOQ.list,paste0("LOQ.mod_",k,"_",i))
          }, error=function(e) {
            e
            cat("ERROR: linear LOQ model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
        }
        ## Define LOQ model by polynomial modeling:
        if(substr(LOQ.FCT,1,1)=="P") {
          Z <- as.numeric(substr(LOQ.FCT,2,nchar(LOQ.FCT)))
          tryCatch({ #skip if model cannot be determined.
            assign(paste0("LOQ.mod_",k,"_",i),lm(CV_metric~poly(log10(Standards),Z),
                                          data=DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],]))
            LOQ.list <- c(LOQ.list,paste0("LOQ.mod_",k,"_",i))
          }, error=function(e) {
            e
            cat("ERROR: ",Z,"-order polynomial LOQ model cannot be defined for ",as.character(Plates[k]),
                as.character(Targets[i]),sep="")
          })
        }
        ## Signal undetermined model with NA:
        if(length(LOQ.list)<i+1) {
          LOQ.list <- c(LOQ.list,NA)
        }
        ## Define the logarithmic model for LOD using user-selected function:
        if(is.list(LOD.FCT)==TRUE) {
          tryCatch({ #skip if model cannot be determined.
            assign(paste0("LOD.mod2_",k,"_",i),drm(Detect~SQ,data=DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k]],fct=LOD.FCT))
            LOD.list2 <- c(LOD.list2,paste0("LOD.mod2_",k,"_",i))
            LOD.list3 <- c(LOD.list3,LOD.FCT$name)
          }, error=function(e) {
            e
            cat("ERROR: LOD model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
          })
        }
        ## Define the logarithmic model with function automatically selected:
        if(is.character(LOD.FCT)) {
          if(LOD.FCT=="Best") {
            tryCatch({ #skip if model cannot be determined.
              ## Pull out data for specific assay:
              TEMP.DAT <- DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k],]
              ## Define a model to start with:
              LOD.mod <- drm(Detect~SQ,data=TEMP.DAT,fct=W2.4())
              ## Test all available models and select the best one:
              LOD.FCT2 <- row.names(mselect(LOD.mod,LOD.FCTS))[1]
              LOD.FCT3 <- getMeanFunctions(fname=LOD.FCT2)
              assign(paste0("LOD.mod2_",k,"_",i),drm(Detect~SQ,data=DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k],],fct=LOD.FCT3[[1]]))
              LOD.list2 <- c(LOD.list2,paste0("LOD.mod2_",k,"_",i))
              LOD.list3 <- c(LOD.list3,LOD.FCT2)
            }, error=function(e) {
              e
              cat("ERROR: LOD model cannot be defined for ",as.character(Plates[k]),as.character(Targets[i]),sep="")
            })
          }
        }
        ## Signal undetermined model with NA:
        if(length(LOD.list2)<i+1) {
          LOD.list2 <- c(LOD.list2,NA)
          LOD.list3 <- c(LOD.list3,NA)
        }
        ## Populate summary data:
        if(ASSAY.TYPE == "qPCR") {
          for(pp in 1:nrow(DAT3)) {
          DAT3$R.squared[pp] <- summary(get(curve.list[pp+1]))$r.squared
          DAT3$Slope[pp] <- coef(get(curve.list[pp+1]))[2]
          DAT3$Intercept[pp] <- coef(get(curve.list[pp+1]))[1]
          }
        }
        
        if(ASSAY.TYPE == "ddPCR") {
          if(!is.na(curve.list[i+1])) {
            for(pp in 1:nrow(DAT3)) {
            DAT3$R.squared[pp] <- summary(get(curve.list[pp+1]))$r.squared
            DAT3$Slope[pp] <- coef(get(curve.list[pp+1]))[2]
            DAT3$Intercept[pp] <- coef(get(curve.list[pp+1]))[1]
          }
            } else {
              for(pp in 1:nrow(DAT3)) {
            DAT3$R.squared[pp] <- NA
            DAT3$Slope[pp] <- NA
            DAT3$Intercept[pp] <- NA
              }
            }
        }
        DAT3$LOD.discrete[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- min(DAT2$Standards[DAT2$Rate>=LOD.Threshold&DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])
       
        
        ## Only get LOD values if the LOD model is defined:
        if(length(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]) !=0) {
          DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),LOD.Threshold,type="absolute")[1]
          DAT3$rep2.LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-sqrt(1-LOD.Threshold),type="absolute")[1]
          DAT3$rep3.LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^(1/3),type="absolute")[1]
          DAT3$rep4.LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^0.25,type="absolute")[1]
          DAT3$rep5.LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^0.2,type="absolute")[1]
          DAT3$rep8.LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^0.125,type="absolute")[1]
          DAT3$rep9.LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^(1/9),type="absolute")[1]
          ## Residual code using probit method:
          #DAT3$LOD[i] <- (qnorm(0.95)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
          #DAT3$rep2.LOD[i] <- (qnorm(0.50)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
          #DAT3$rep3.LOD[i] <- (qnorm(1/3)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
          #DAT3$rep4.LOD[i] <- (qnorm(0.25)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
          #DAT3$rep5.LOD[i] <- (qnorm(0.2)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
          #DAT3$rep8.LOD[i] <- (qnorm(0.125)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
        }
        ## Generate prediction data for LoQ:
        ## Only get LOQ if LOQ model is determined:
        if(length(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))]) !=0) {
          newData <- data.frame(Standards = seq(1, 10000))
          newData$CV_metric <- predict(get(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))]), newData)
          ## Determine what type of LOQ model is used and calculate LOQ accordingly:
          ## For exponential decay:
          if(as.character(get(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))])$call)[1]=="nls") {
            ## Look up lowest modeled standard below the CV threshold:
            DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- min(newData$Standards[newData$CV_metric<=LOQ.Threshold])
            ## Unless... If the background variation exceeds the CV threshold, adjust threshold:
            ## Determine the highest standard used:
            A <- max(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])
            if(min(newData$CV_metric[newData$Standards<=A])>LOQ.Threshold) {
              ## Set the adjusted threshold to 1.5x the lowest simulated Cq.CV 
              ##   within the range of data tested:
              B <- min(newData$CV_metric[newData$Standards<=A])
              DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- min(newData$Standards[newData$CV_metric<=B*1.5])
              ## Make a note of the adjusted threshold in the analysis log:
              ToWrite <- paste0("Note: All standards tested for ",Targets[i],
                                " yielded higher Cq.CV values than the user-defined CV threshold of ",
                                LOQ.Threshold,". The CV threshold has been adjusted to ",
                                B*1.5," for the LOQ of this marker.")
              write(paste0(ToWrite,"\n\n"),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
            }
          }
          if(as.character(get(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))])$call)[1]=="lm") {
            ## For polynomial:
            if(grepl("poly",as.character(get(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))])$call)[2])==TRUE) {
              ## Determine the highest standard used:
              A <- max(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])
              ## Adjust if the tested range does not cross below the CV threshold:
              if(min(DAT2$CV_metric[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],na.rm=TRUE)>LOQ.Threshold) {
                B <- min(DAT2$CV_metric[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],na.rm=TRUE)*1.5
                ## Make a note of the adjusted threshold in the analysis log:
                ToWrite <- paste0("Note: All standards tested for ",Plates[k]," and ",Targets[i],
                                  " yielded higher Cq.CV values than the user-defined CV threshold of ",
                                  LOQ.Threshold,". The CV threshold has been adjusted to ",
                                  B," for the LOQ of this marker.")
                write(paste0(ToWrite,"\n\n"),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
              }
              else {
                B <- LOQ.Threshold
              }
              ## Look up highest modeled standard below the CV threshold:
              C <- max(newData$Standards[newData$CV_metric<=B&newData$Standards<=A])
              ## Look up the highest modeled standard above the CV threshold...
              ##   and also below the highest standard below the CV threshold.
              ##   This captures the farthest right crossing point on a downward slope.
              D <- max(newData$Standards[newData$CV_metric>B&newData$Standards<C])
              ## LOQ is D + 1 to get back less than or equal to the CV threshold.
              DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- D+1
            }
            # For linear:
            else {
              ## Look up lowest modeled standard below the CV threshold:
              DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- min(newData$Standards[newData$CV_metric<=LOQ.Threshold])
            }
          }
          ## If modeled LOQ is calculated to be below the 95% LOD, set LOD as LOQ:
          if(!is.na(DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]])) {
            if(DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]<DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]) {
              DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]
            }
          }
          ## If modeled LOQ is calculated to be below the lowest standard tested,
          ##   set the lowest standard as the LOQ:
          if(DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]<min(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])) {
            DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]] <- min(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])
          }
        }
      }
      }
      write("Assay summary:",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("\nR.squared: The R-squared value of linear regression of all standards Cq-values vs log10 of the starting quantities.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("Slope: The slope of the linear regression.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("Intercept: The y-intercept of the linear regression.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write(paste0("\nLOD.discrete: The lowest standard with at least ",LOD.Threshold*100,"% positive detection."),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write(paste0("LOD: The ",LOD.Threshold*100,"% limit of detection as determined by dose-response modeling."),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write(paste0("LOQ: The limit of quantification as determined by decay modeling, using the user-selected CV threshold of: ",LOQ.Threshold),file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("\nrep2.LOD: The effective limit of detection if analyzing in 2 replicates.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("rep3.LOD: The effective limit of detection if analyzing in 3 replicates.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("rep4.LOD: The effective limit of detection if analyzing in 4 replicates.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("rep5.LOD: The effective limit of detection if analyzing in 5 replicates.",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("rep8.LOD: The effective limit of detection if analyzing in 8 replicates.\n\n",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write("rep9.LOD: The effective limit of detection if analyzing in 9 replicates.\n\n",file=file.path(output_folder, "Analysis Log.txt"),append=TRUE)
      write.csv(DAT3,file=file.path(output_folder, "Assay summary.csv"),row.names=FALSE)
      
      
      if(ASSAY.TYPE == "qPCR") {
        ## Plot Cq value vs Standard Concentration standard curves:
        DAT$Mod[DAT$Mod==0] <- "Excluded"
        DAT$Mod[DAT$Mod==1] <- "Modeled"
        for(k in 1:length(Plates)) {
          Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
        for(i in 1:length(Targets)) {
          ggOut <- ggplot(data=DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&is.na(DAT$SQ)==FALSE,],
                          aes(x=SQ,y=Cq,color=factor(Mod),shape=factor(Mod),size=factor(Mod))) + 
            geom_jitter(width=0.1,alpha=0.75) + 
            scale_shape_manual("",values=c(3,20),guide=FALSE) +
            scale_size_manual("",values=c(1,3)) +
            scale_x_log10() +
            scale_color_manual("",values=c("blue", "black")) +
            xlab("Standard Concentrations (Copies / Reaction)") +
            ylab("Cq-value") +
            geom_abline(intercept=coef(get(curve.list[which(sapply(curve.list, function(x) x == paste0("curve_",k,"_",i)))]))[1],
                        slope=coef(get(curve.list[which(sapply(curve.list, function(x) x == paste0("curve_",k,"_",i)))]))[2]) +
            geom_vline(xintercept=DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]],colour="red") +
            geom_vline(xintercept=DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]],linetype=2) +
            annotate("text",y=max(DAT$Cq[DAT$Target==Targets[i]],na.rm=TRUE)*0.99,color="red",
                     x=DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]*0.8,angle=90,label="LOD") +
            annotate("text",y=max(DAT$Cq[DAT$Target==Targets[i]],na.rm=TRUE)*0.94,
                     x=DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]*0.8,angle=90,label="LOQ") +
            theme_bw() + theme(legend.justification=c(1,1),legend.position=c(1,0.99)) +
            ggtitle(paste0("Standard curve for: ",Plates[k]," and ",Targets[i])) +
            theme(plot.title=element_text(hjust=0.5,size=20),
                  axis.title=element_text(size=16)) +
            theme(legend.title=element_blank(),
                  legend.text=element_text(size=11)) +
            annotate("text",y=min(DAT$Cq[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&is.na(DAT$SQ)==FALSE],na.rm=TRUE)*1.05,
                     x=min(DAT$SQ[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&is.na(DAT$SQ)==FALSE],na.rm=TRUE)*1.01,hjust=0,
                     label=(paste0("R-squared: ",DAT3$R.squared[i],"\ny = ",DAT3$Slope[i],"x + ",DAT3$Intercept[i])))
          ggsave(file.path(output_folder, paste0("StandardCurve_",Plates[k],"_and_",Targets[i], ".png")), 
                 plot=ggOut, width=8, height=6, dpi=300)
          print(ggOut)
          readline(prompt="Press [Enter] for next plot.")
          print("Calculating... Please wait.")
        }
        }
      }
      
      
      if(ASSAY.TYPE == "ddPCR") {
        ## Plot Concentration vs Expected Concentration:
        for(k in 1:length(Plates)) {
          Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
        for(i in 1:length(Targets)) {
          ggOut <- ggplot(data=DAT[DAT$Target==Targets[i]&DAT$Plate==Plates[k]&is.na(DAT$SQ)==FALSE,],
                          aes(x=SQ,y=Concentration)) + 
            geom_jitter(width=0.1,alpha=0.75,size=3) + 
            scale_x_log10() +
            scale_y_log10() +
            xlab("Expected Concentrations (Copies / Reaction)") +
            ylab("Measured Concentration (Copies / Reaction)") +
            geom_abline(intercept=0, slope=1, linetype="dashed", color="blue") +
            geom_vline(xintercept=DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]],colour="red") +
            geom_vline(xintercept=DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]],linetype=2) +
            annotate("text",y=max(DAT$Concentration[DAT$Target==Targets[i]&DAT$Plate==Plates[k]],na.rm=TRUE)*0.99,color="red",
                     x=DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]*0.8,angle=90,label="LOD") +
            annotate("text",y=max(DAT$Concentration[DAT$Target==Targets[i]&DAT$Plate==Plates[k]],na.rm=TRUE)*0.94,
                     x=DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]*0.8,angle=90,label="LOQ") +
            theme_bw() +
            ggtitle(paste0("ddPCR accuracy plot for: ",Plates[k]," and ",Targets[i])) +
            theme(plot.title=element_text(hjust=0.5,size=20),
                  axis.title=element_text(size=16))
          
          ggsave(file.path(output_folder, paste0("ddPCR_Accuracy_",Plates[k],"_and_", Targets[i], ".png")), 
                 plot=ggOut, width=8, height=6, dpi=300)
          print(ggOut)
          readline(prompt="Press [Enter] for next plot.")
          print("Calculating... Please wait.")
        }
        }
      }
      
      ## Plot the LOD models for each assay:
      for(k in 1:length(Plates)) {
        Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
      for(i in 1:length(Targets)) {
        if(length(is.na(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]))!=0) {
          png(file.path(output_folder, paste0("LOD_Plot_", Plates[k],"_and_",Targets[i], ".png")), 
              width=8, height=6, units="in", res=300)
          
          DAT4 <- rbind(ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),LOD.Threshold,interval="delta",type="absolute"),
                        ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-sqrt(1-LOD.Threshold),interval="delta",type="absolute"),
                        ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^(1/3),interval="delta",type="absolute"),
                        ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^0.25,interval="delta",type="absolute"),
                        ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^0.2,interval="delta",type="absolute"),
                        ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^0.125,interval="delta",type="absolute"),
                        ED(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),1-(1-LOD.Threshold)^(1/9),interval="delta",type="absolute"))
          if(substr(LOD.list3[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))],1,3)=="LL2") {
            DAT4 <- exp(DAT4)
          } 
    
          DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                        "5rep.LOD","8rep.LOD","9rep.LOD"),
                             Assay=rep(Targets[i],nrow(DAT4)))
          DAT4$Assay <- as.character(DAT4$Assay)
          if(i==1) {
            LOD.CI <- DAT4
          }
          if(i>1) {
            LOD.CI <- rbind(LOD.CI,DAT4)
          }
          if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
            DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
          }
          plot(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]),main=paste0("LoD Plot for: ",Plates[k]," and ",Targets[i]),
               ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
               xlim=c(min(DAT4[,1:4],na.rm=TRUE),max(DAT$SQ,na.rm=TRUE)))
          LODS <- sum(!is.na(DAT4[,1]))
          COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
                    rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
          PNTS <- c(15,16,17,18,25,3)
          YS <- c(LOD.Threshold,1-sqrt(1-LOD.Threshold),1-(1-LOD.Threshold)^(1/3),1-(1-LOD.Threshold)^0.25,1-(1-LOD.Threshold)^0.2,1-(1-LOD.Threshold)^0.125)
          LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
                     "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
          if(LODS<6) {
            LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
          }
          points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=1.2)
          for(j in 1:LODS) {
            lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=2)
            lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=2,col=COLS[j])
            lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=2,col=COLS[j])
          }
          legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS)
          Pval <- modelFit(get(LOD.list2[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))]))[[5]][2]
          mtext(paste0("FCT used: ",LOD.list3[which(sapply(LOD.list2, function(x) x == paste0("LOD.mod2_",k,"_",i)))],"    Lack of fit test: p = ",Pval),side=3)
          
          dev.off()
        }
        if(is.na(LOD.list2[i+1])) {
          png(file.path(output_folder, paste0("LOD_Plot_", Plates[k]," and ",Targets[i], "_unsolvable.png")), 
              width=8, height=6, units="in", res=300)
          
          plot(DAT2$Rate[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]~log10(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]),
               ylim=c(0,1),ylab="Detection Probability",
               xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
               main=paste0("LoD for: ",Plates[k]," and ",Targets[i]," unsolvable"))
          
          dev.off()
        }
        readline(prompt="Press [Enter] for next plot.")
        print("Calculating... Please wait.")
      }
      }
      
      LOD.CI <- LOD.CI[,c(6,5,1,3,4,2)]
      write.csv(LOD.CI,file=file.path(output_folder, "LOD_confint.csv"),row.names=FALSE)
      
      ## Plot the LoQ models of each assay:
     for(k in 1:length(Plates)) {
       Targets <- unique(DAT$Target[DAT$Plate==Plates[k]])
       for(i in 1:length(Targets)) {
        if(is.na(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))])==FALSE) {
          ## Re-generate prediction data for the model:
          newData <- data.frame(Standards = seq(1, 10000))
          newData$CV_metric <- predict(get(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))]), newData)
          
          ## Define LOQ polygon coordinates:
          PDAT <- data.frame(x=c(min(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],na.rm=TRUE),
                                 min(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],na.rm=TRUE),
                                 DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]],
                                 DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]),
                             y=c(min(c(DAT2$CV_metric[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],newData$CV_metric[newData$Standards<=max(DAT2$Standards[DAT2$Target==Targets[i]])&newData$Standards>=min(DAT2$Standards[DAT2$Target==Targets[i]])]),na.rm=TRUE)*0.9,
                                 newData$CV_metric[newData$Standards==DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]],
                                 newData$CV_metric[newData$Standards==DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]],
                                 min(c(DAT2$CV_metric[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],newData$CV_metric[newData$Standards<=max(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])&newData$Standards>=min(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])]),na.rm=TRUE)*0.9))
          if(DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]!=floor(DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]])) {
            PDAT$y[2:3] <- LOQ.Threshold
          }
        }
        
        Decay.Plot <- ggplot(DAT2[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k],], aes(x= Standards, y = CV_metric)) +
          geom_point(size=2) +
          scale_x_continuous(trans = 'log10') +
          ylab(CV_label) +
          xlab("Standard concentrations (Copies / Reaction)") +
          geom_vline(xintercept=DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]],color="red") +
          annotate("text",y=max(DAT2$CV_metric[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],na.rm=TRUE)*0.99,
                   x=DAT3$LOD[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]*0.8,angle=90,label="LOD",color="red") +
          theme(legend.position="none") +
          theme(plot.title=element_text(hjust=0.5))
        
        if(is.na(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))])==FALSE) {
          if(DAT3$LOQ[DAT3$Assay==Targets[i]&DAT3$Plate==Plates[k]]<=min(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]])) {
            PDAT$x[3:4] <- NA
            Decay.Plot <- Decay.Plot + 
              annotate("text",y=max(DAT2$CV_metric[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]],na.rm=TRUE)*0.99,
                       x=median(DAT2$Standards[DAT2$Target==Targets[i]&DAT2$Plate==Plates[k]]),
                       label="LOQ may be outside tested range.",hjust=0)
          }
          Decay.Plot <- Decay.Plot + geom_polygon(data=PDAT,aes(x=x,y=y,alpha=0.5))
          
          if(as.character(get(LOQ.list[which(sapply(LOQ.list, function(x) x == paste0("LOQ.mod_",k,"_",i)))])$call)[1]=="nls") {
            Decay.Plot <- Decay.Plot + 
              stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE) +
              ggtitle(paste0("Exponential Decay LOQ model for: ",Plates[k]," and ",Targets[i]))
          }
          
          if(as.character(get(LOQ.list[i+1])$call)[1]=="lm") {
            if(grepl("poly",as.character(get(LOQ.list[i+1])$call)[2])==TRUE) {
              B <- length(get(LOQ.list[i+1])$coefficients)-1
              Decay.Plot <- Decay.Plot +
                stat_smooth(method = "lm", formula = y ~ poly(x,B),se=FALSE) +
                ggtitle(paste0(B,"-order polynomial LOQ model for: ",Plates[k]," and ",Targets[i]))
            }
            else {
              Decay.Plot <- Decay.Plot +
                stat_smooth(method = "lm", formula = y ~ x,se=FALSE) +
                ggtitle(paste0("Linear LOQ model for: ",Plates[k]," and ",Targets[i]))
            }
          } 
        }
        
        if(is.na(LOQ.list[i+1])==TRUE) {
          Decay.Plot <- Decay.Plot + 
            ggtitle(paste0("LOQ model for: ",Plates[k]," and ",Targets[i]," not solvable."))
        }
        
        ggsave(file.path(output_folder, paste0("LOQ_Plot_", Plates[k]," and ",Targets[i], ".png")), 
               plot=Decay.Plot, width=8, height=6, dpi=300)
        
        print(Decay.Plot)
        readline(prompt="Press [Enter] for next plot.")
        print("Calculating... Please wait.")
      }
      } # End Targets loop for LOQ plots
      
      ## Add parameter info to DAT3 and store for comparison:
      DAT3$LOD_Threshold <- LOD.Threshold
      DAT3$LOQ_Threshold <- LOQ.Threshold
      if(ASSAY.TYPE == "ddPCR") {
        DAT3$Min_Pos_Droplets <- MIN.POSITIVE.DROPLETS
      }
      comparison_summary <- rbind(comparison_summary, DAT3)
      
    } # End MIN.POSITIVE.DROPLETS loop
  } # End LOQ.Threshold loop
} # End LOD.Threshold loop


## Save comparison summary:
write.csv(comparison_summary, 
          file=file.path(master_output_folder, "Parameter_Comparison_Summary.csv"),
          row.names=FALSE)

print("Creating comparison plots...")

comparison_summary %>% 
  as.tibble() %>% 
  mutate(., LOD_F= if_else(is.na(rep3.LOD),LOD.discrete,  rep3.LOD)) %>% 
  group_by(Assay) %>% 
  summarise(
    n = n(),
    min_LOD_F = min(LOD_F, na.rm = TRUE),
    max_LOD_F = max(LOD_F, na.rm = TRUE),
    median_LOD_F = median(LOD_F, na.rm = TRUE),
    mean_LOD_F = mean(LOD_F, na.rm = TRUE),
    sd_LOD_F = sd(LOD_F, na.rm = TRUE),
    
    # Calculate 95% CI bounds
    se = sd_LOD_F / sqrt(n),                     # Standard Error
    margin_error_LOD_F = qt(0.975, df = n - 1) * se, # Margin of Error (t-distribution)
    
    ci_lower_LOD_F = min_LOD_F - margin_error_LOD_F,
    ci_upper_LOD_F = max_LOD_F + margin_error_LOD_F,
    min_LOQ = min(LOQ, na.rm = TRUE),
    max_LOQ = max(LOQ, na.rm = TRUE),
    median_LOQ = median(LOQ, na.rm = TRUE),
    mean_LOQ = mean(LOQ, na.rm = TRUE),
    sd_LOQ = sd(LOQ, na.rm = TRUE),
    
    # Calculate 95% CI bounds
    se = sd_LOQ / sqrt(n),                     # Standard Error
    margin_error_LOQ = qt(0.975, df = n - 1) * se, # Margin of Error (t-distribution)
    
    ci_lower_LOQ = mean_LOQ - margin_error_LOQ,
    ci_upper_LOQ = mean_LOQ + margin_error_LOQ
  ) -> comparison_summary_assay

write.csv(comparison_summary_assay, 
          file=file.path(master_output_folder, "Parameter_Comparison_Summary_assay.csv"),
          row.names=FALSE)

## LOD comparison - faceted by LOD threshold:
if(ASSAY.TYPE == "qPCR") {
  lod_by_lod <- ggplot(comparison_summary, 
                       aes(x=Assay, y=LOD, 
                           fill=Assay)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOD_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOD Comparison by LOD Threshold",
         y="LOD (Copies/Reaction)",
         fill="LOQ CV Threshold") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
} else {
  lod_by_lod <- ggplot(comparison_summary, 
                       aes(x=Assay, y=LOD, 
                           fill=Assay)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOD_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOD Comparison by LOD Threshold",
         y="LOD (Copies/Reaction)",
         fill="LOQ_CV / MinPos") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
}

ggsave(file.path(master_output_folder, "LOD_by_LOD_Threshold.png"),
       plot=lod_by_lod, width=12, height=6, dpi=300)

## LOD comparison - faceted by LOQ threshold:
if(ASSAY.TYPE == "qPCR") {
  lod_by_loq <- ggplot(comparison_summary, 
                       aes(x=Assay, y=LOD, 
                           fill=Assay)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOQ_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOD Comparison by LOQ Threshold",
         y="LOD (Copies/Reaction)",
         fill="LOD Threshold") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
} else {
  lod_by_loq <- ggplot(comparison_summary, aes(x=Assay, y=LOD, fill=Assay)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOQ_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOD Comparison by LOQ Threshold",
         y="LOD (Copies/Reaction)",
         fill="LOD% / MinPos") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
}

ggsave(file.path(master_output_folder, "LOD_by_LOQ_Threshold.png"),
       plot=lod_by_loq, width=12, height=6, dpi=300)

## LOQ comparison - faceted by LOD threshold:
if(ASSAY.TYPE == "qPCR") {
  loq_by_lod <- ggplot(comparison_summary,
                       aes(x=Assay, y=LOQ,
                           fill=factor(LOQ_Threshold))) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOD_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOQ Comparison by LOD Threshold",
         y="LOQ (Copies/Reaction)",
         fill="LOQ CV Threshold") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
} else {
  loq_by_lod <- ggplot(comparison_summary,
                       aes(x=Assay, y=LOQ,
                           fill=interaction(factor(LOQ_Threshold),
                                            factor(Min_Pos_Droplets)))) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOD_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOQ Comparison by LOD Threshold",
         y="LOQ (Copies/Reaction)",
         fill="LOQ_CV / MinPos") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
}

ggsave(file.path(master_output_folder, "LOQ_by_LOD_Threshold.png"),
       plot=loq_by_lod, width=12, height=6, dpi=300)

## LOQ comparison - faceted by LOQ threshold:
if(ASSAY.TYPE == "qPCR") {
  loq_by_loq <- ggplot(comparison_summary,
                       aes(x=Assay, y=LOQ,
                           fill=factor(LOD_Threshold))) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOQ_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOQ Comparison by LOQ Threshold",
         y="LOQ (Copies/Reaction)",
         fill="LOD Threshold") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
} else {
  loq_by_loq <- ggplot(comparison_summary,
                       aes(x=Assay, y=LOQ,
                           fill=interaction(factor(LOD_Threshold),
                                            factor(Min_Pos_Droplets)))) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(Plate~LOQ_Threshold, labeller = label_both) +
    theme_bw() +
    labs(title="LOQ Comparison by LOQ Threshold",
         y="LOQ (Copies/Reaction)",
         fill="LOD% / MinPos") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(hjust=0.5, size=16),
          legend.position = "bottom")
}

ggsave(file.path(master_output_folder, "LOQ_by_LOQ_Threshold.png"),
       plot=loq_by_loq, width=12, height=6, dpi=300)


cat("\n\nComparison analysis complete! Results saved in:", master_output_folder, "\n")
cat("Total parameter combinations tested:", nrow(comparison_summary) / length(unique(comparison_summary$Assay)), "\n")
cat("\nGenerated plots:\n")
cat("  - LOD_by_LOD_Threshold.png\n")
cat("  - LOD_by_LOQ_Threshold.png\n")
if(ASSAY.TYPE == "ddPCR") {
  cat("  - LOD_by_MinPos.png\n")
}
cat("  - LOQ_by_LOD_Threshold.png\n")
cat("  - LOQ_by_LOQ_Threshold.png\n")
if(ASSAY.TYPE == "ddPCR") {
  cat("  - LOQ_by_MinPos.png\n")
}
if(length(unique(comparison_summary$Assay)) > 1) {
  cat("  - LOD_Grid.png (2D facet grid)\n")
  cat("  - LOQ_Grid.png (2D facet grid)\n")
}
cat("  - LOD_Sensitivity.png\n")


