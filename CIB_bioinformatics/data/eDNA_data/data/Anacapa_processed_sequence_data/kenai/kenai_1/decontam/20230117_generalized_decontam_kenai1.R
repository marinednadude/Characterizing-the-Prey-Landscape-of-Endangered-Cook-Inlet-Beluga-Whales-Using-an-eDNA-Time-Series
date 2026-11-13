# ---
# title: "Generalized Decontamination Script"
# author: "Zack Gold"
# date: "1/17/2023"
# output: html_document
# ---

# Load Libraries ----------------------------------------------------------------------------------------------------------------
library (tidyverse)
library (vegan)
library (proxy)
library(reshape2)
library(microDecon)
library(stringr)
library(knitr)
library(ranacapa)
library(dplyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(plotly)
library(optparse)
library(fitdistrplus)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
library(matrixStats)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source(here("data","eDNA_data","data","Anacapa_processed_sequence_data","twentymile","seqrun1","decontam","decontamination_utilities.R"))
here()
#Test inputs ----------------------------------------------------------------------------------------------------------------
working_directory <- here("data","eDNA_data","data","Anacapa_processed_sequence_data","kenai","kenai_1","decontam")
input_anacapa_path_1 <- here("data","eDNA_data","data","Anacapa_processed_sequence_data","kenai","kenai_1","c19_fishcard_taxonomy_tables","Summary_by_percent_confidence","60","c19_fishcard_ASV_raw_taxonomy_60.txt")

input_meta_path <- here("data","Metadata","metadata_kenai1_20220401.csv")
number_anacapa_tables <- "one"

read_type <- "merged_only"

# Step3 parameters
low_dist_probability_cutoff <- 0.001
minimum_read_cutoff <- 20000
high_dist_probability_cutoff <- 1

# Step4 parameters
swamping_cutoff <- 0.9

# Step6 parameters
som_level <- "double"
som_chains<-10
som_iterations<-1000
som_filter_choice <- "som_min_threshold_1"
som_filter_choice_2 <- "som_min_threshold_1_single"

# Step7 parameters
max_dist_probability_cutoff <- 0.96
#---

#---

#Load Data ----------------------------------------------------------------------------------------------------------------
setwd(working_directory)
dir.create("Output_csv")
dir.create("Output_R")
dir.create("Output_plots")

metadata <- read.table(input_meta_path, header = 1, sep = ",", stringsAsFactors = F)

anacapa_table_1 <- read.table(input_anacapa_path_1, header = 1, sep = "\t", stringsAsFactors = F)

#Format Anacapa Table 1
anacapa_table_1[1,1] %>% str_remove(.,"forward_") %>% str_sub(., end=-3) -> barcode_1
barcode_1_use <- paste0(barcode_1)
anacapa_table_1$Miseq_run <- barcode_1_use

###Fix Names
anacapa_table_1_names <- colnames(anacapa_table_1)
anacapa_table_1_names %>% as.data.frame() -> anacapa_table_1_names
anacapa_table_1_names[1,1] <- "seq_number"
colnames(anacapa_table_1_names) <- c("Seq_number")
anacapa_table_1_names %>% mutate(.,  Seq_number =str_remove(Seq_number,paste0(barcode_1,"_"))) -> anacapa_table_1_names

left_join(anacapa_table_1_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_1_names
colnames(anacapa_table_1) <- anacapa_table_1_names$final_names
anacapa_table_1 %>% 
  mutate(., seq_number = str_replace(seq_number, "c19_fishcard", barcode_1_use)) -> anacapa_table_1

#Generate Hash.key
anacapa_table_1 %>% 
  dplyr::select(seq_number, sum.taxonomy) %>% 
  filter(., str_detect(seq_number, "merged")) %>% 
  filter(., !str_detect(seq_number, "unmerged")) -> hash.key_1_c19

saveRDS(hash.key_1_c19,file= here(working_directory,"Output_R","hash.key.RDS"))

hash.key_updated <- readRDS( file=here(working_directory,"hash.key_updated.RDS"))

# Convert to Long Data

###Format for Long Data
anacapa_table_1$seq_number <- factor(anacapa_table_1$seq_number)
anacapa_table_1$Miseq_run <- factor(anacapa_table_1$Miseq_run)

columns <- colnames(anacapa_table_1)
remove <- c("seq_number","sum.taxonomy","Miseq_run")

gathercols <-  columns[! columns %in% remove] 

anacapa_table_1 %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="reads") -> anacapa_table_1_long

anacapa_table_1_long$reads <- as.numeric(anacapa_table_1_long$reads)
head(hash.key_updated)

anacapa_table_1_long %>% 
  dplyr::select(-sum.taxonomy) %>% 
  left_join(hash.key_updated) -> anacapa_table_1_long

#Merge All Tables
ASV.table_long <- anacapa_table_1_long

dim(ASV.table_long)
#---

#---

# Cleaning Process 1: Remove all Forward, Reverse, and Unmerged reads & Remove Singletons ----------------------------------------------------------------------------------------------------------------

#Filter Merged only reads

#Filter out forward, reverse, and unmerged reads
ASV.table_long %>% 
  filter(., str_detect(seq_number,"merged")) %>% 
  filter(., !str_detect(seq_number,"unmerged")) -> ASV.table_used

#Calculate % ASVs Kept
ASV.table_long %>%  dim() -> all_dim
ASV.table_used %>%  dim() -> used_only_dim

ASV.table_long %>% arrange(desc(reads)) %>% 
  filter(., !str_detect(seq_number,"merged")) %>% 
  dplyr::summarise(sum(reads))


paste0(round(used_only_dim[[1]]/all_dim[[1]]*100,2),"% ASVs retained")

saveRDS(ASV.table_long, file=here(working_directory,"ASV.table_long.RDS"))
saveRDS(ASV.table_used, file=here(working_directory,"ASV.table_used.RDS"))

#Total Initial Reads and ASVs
ASV.table_long %>% 
  group_by(Miseq_run) %>% 
  dplyr::summarise(n_distinct(seq_number), sum(reads))

#Total Merged Reads and ASVs
ASV.table_used %>% 
  group_by(Miseq_run) %>% 
  dplyr::summarise(n_distinct(seq_number), sum(reads))

# Technical Replicate Exploration -------
ASV.table_used <- readRDS(file=here(working_directory,"ASV.table_used.RDS"))

ASV.table_used %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(sample_read_depth= sum(reads)) -> sample_read_depth

#Mean Read Depth
mean(sample_read_depth$sample_read_depth)
median(sample_read_depth$sample_read_depth)

#SD Read Depth
sd(sample_read_depth$sample_read_depth)
#Min Read Depth
min(sample_read_depth$sample_read_depth)
#Max Read Depth
max(sample_read_depth$sample_read_depth)


#---

#---

ASV.table_used <- readRDS(file=here(working_directory,"ASV.table_used.RDS"))
hash.key_updated <- readRDS( file=here(working_directory,"hash.key_updated.RDS"))
metadata <- read.table(input_meta_path, header = 1, sep = ",", stringsAsFactors = F)



# Cleaning Process 2: Estimation of *Tag-jumping* or sample *Cross-talk* ----------------------------------------------------------------------------------------------------------------

## Step 1: Nest the dataset by origin of ASVs

###Identify Positives, Negatives, and Samples

###Create list of control samples
metadata %>% 
  filter(Sample_Control=="Control") %>% 
  dplyr::select(New_name) %>% unique() -> controls
controls <- controls$New_name

metadata %>% 
  filter(Control_Type=="Pos") %>% 
  dplyr::select(New_name) %>% unique()-> pos_controls
pos_controls <- pos_controls$New_name

metadata %>% 
  filter(Control_Type=="Neg") %>% 
  dplyr::select(New_name) %>% unique() -> neg_controls
neg_controls <- neg_controls$New_name

###New column that labels each ASV as from Positive (control) or Sample
ASV.table_used %>% 
  mutate(source = case_when(sample %in% pos_controls~"Positives",
                            sample %in% neg_controls~"Blanks",
                            TRUE ~"Samples")) -> ASV.table_used


###Convert to tibble
ASV.table_used <- as_tibble(ASV.table_used)

###Remove empty sequences
ASV.table_used %>% 
  filter(reads != 0)  -> ASV.table_used

###Rename Columns and remove seq_number
ASV.table_used %>%
  mutate(sample = as.character(sample),
         nReads = reads)  -> ASV.table_used

###ASVs in Positive Controls
ASV.table_used %>% 
  filter (source == "Positives") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.positives

hash.key_updated %>% 
  filter(seq_number %in% all.seqs.in.positives) %>% as_tibble() -> pos.contam.species

write.csv(pos.contam.species, file="pos.contam.species.csv")

###ASVs in Negative Controls

ASV.table_used %>% 
  filter (source == "Blanks") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.blanks

hash.key_updated %>% 
  filter(seq_number %in% all.seqs.in.blanks) %>% as_tibble() -> blank.contam.species

write.csv(blank.contam.species, file="blank.contam.species.csv")

### Visualize Read Counts Across Samples for Barcode_1
ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_1_use) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>% 
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() +ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_1

ggsave(plot=plot_1,here(working_directory,"Output_plots","Sample_Read_Depth_Barcode_1.png"), device = "png", width = 12, height = 8, units = "in")


###Nesting the dataset
ASV.table_used %>% 
  dplyr::group_by(Miseq_run, source) %>% 
  nest() %>% 
  pivot_wider(names_from=source, values_from=data) -> ASV.nested


####Summary.file.1
ASV.nested %>% 
  ungroup() %>% 
  dplyr::transmute(.,Miseq_run,Summary = purrr::map(Samples, ~ how.many(ASVtable = ., round = 0)))  -> ASV.summary

ASV.summary$Summary

## Step 2: Model the composition of the positive controls of each run 

###Jumping vector

ASV.nested %>% 
  mutate (contam.tibble = purrr::map(Positives, 
                                     function(.x){
                                       .x %>%
                                         ungroup() %>% 
                                         group_by(sample) %>%
                                         mutate (TotalReadsperSample = sum(nReads)) %>%
                                         mutate (proportion = nReads/TotalReadsperSample) %>%
                                         group_by(seq_number) %>%
                                         dplyr::summarise (vector_contamination = max(proportion))
                                     }) ) -> ASV.nested

###Vector Contamination in Barcode_1
ASV.nested$contam.tibble[[1]] %>% as.data.frame() %>% 
  ggplot(aes(x= vector_contamination))+
  geom_histogram() -> vc_plot_1 # Check how it looks
ggsave(plot=vc_plot_1,here(working_directory,"Output_plots","Vector_Contamination_Barcode_1.png"), device = "png", width = 12, height = 8, units = "in")


##Step 3: Substract the composition of the positive controls from the environment samples

ASV.nested %>% 
  ungroup() %>% 
  mutate(Step2.cleaned.tibble = map2(Samples, contam.tibble, function(.x,.y){ 
    .x %>%
      dplyr::group_by (sample) %>%
      mutate (TotalReadsperSample = sum (nReads)) %>%
      left_join(.y, by = "seq_number") %>%
      mutate (Updated_nReads = ifelse (!is.na(vector_contamination),  nReads - (ceiling(vector_contamination*TotalReadsperSample)), nReads)) %>%
      filter (Updated_nReads > 0) %>%
      ungroup() %>% 
      dplyr::select (sample, seq_number, nReads = Updated_nReads)
  })) -> ASV.nested

###Add this step to the summary table we were creating

####Summary.file.2
ASV.nested %>% 
  transmute(Miseq_run, Summary.2 = purrr::map(Step2.cleaned.tibble, ~ how.many(ASVtable = .,round = "2.Index_hopping"))) %>% 
  left_join(ASV.summary) %>% #use left join when there are many miseq runs to join
  mutate(Summary   = map2(Summary, Summary.2, bind_rows)) %>%
  dplyr::select(-Summary.2) -> ASV.summary 

ASV.summary$Summary

#---

#---

# Cleaning Process 3: **Remove Likely Lab Contaminants** ----------------------------------------------------------------------------------------------------------------

### Filter by Taxonomic Paths -----
# This list can be expanded upon additional work to identify known contaminants within OME PMEL samples
step3_species_to_remove <- c(
    "Eukaryota;Chordata;Aves;Accipitriformes;Accipitridae;Harpagornis;Harpagornis moorei",
    "Eukaryota;Chordata;Mammalia;Primates;Hominidae;Pan;Pan paniscus",
    "Eukaryota;Chordata;Mammalia;Primates;NA;NA;NA",
    "Eukaryota;Chordata;Mammalia;Artiodactyla;Suidae;Sus;Sus scrofa",
    "Eukaryota;Chordata;Mammalia;Artiodactyla;Bovidae;Bos;Bos taurus",
    "Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens",
    "Eukaryota;Chordata;Aves;Galliformes;Phasianidae;Gallus;Gallus gallus",
    "Eukaryota;Chordata;Mammalia;Carnivora;Canidae;Canis;Canis lupus",
    "Unknown;NA;NA;NA;NA;NA;NA",
    ""
  )

## Add to ASV.nested   ----
ASV.nested %>%
  mutate(
    Step3.remove.known.contaminants.tibble = purrr::map(
      Step2.cleaned.tibble,
      ~ left_join(., hash.key_updated) %>%
        filter(., !sum.taxonomy %in% step3_species_to_remove) %>% 
        filter(., str_detect(sum.taxonomy,"Chordata"))
    )
  ) -> ASV.nested


ASV.nested$Step3.remove.known.contaminants.tibble[[1]] %>% dim()
ASV.nested$Step2.cleaned.tibble[[1]] %>% dim()

## Summary.file.3 -----
ASV.nested %>%
  transmute(Miseq_run,
            Summary.1 = purrr::map(
              Step3.remove.known.contaminants.tibble,
              ~ how.many(ASVtable = ., round = "3.labcontaminants")
            )
  ) %>%
  left_join(ASV.summary) %>%
  mutate(Summary = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary


ASV.summary$Summary

# Cleaning Process 4: **Discarding PCR replicates with low number of reads** ----------------------------------------------------------------------------------------------------------------

###Pull out Sample Read Depth

ASV.nested %>% 
  dplyr::select(Miseq_run, Step3.remove.known.contaminants.tibble) %>% 
  unnest(Step3.remove.known.contaminants.tibble) %>% 
  group_by(Miseq_run,sample) %>%
  dplyr::summarise(tot = sum(nReads)) %>% 
  arrange(desc(tot))-> all.reps

all.reps %>%  
  filter(., Miseq_run==barcode_1_use) ->all.reps_1

all.reps_1 %>% 
  pull(tot) -> reads.per.sample_1

names(reads.per.sample_1) <- all.reps_1 %>% pull(sample)  

### Fit Normal Distribution

fit_1 <- fitdist(reads.per.sample_1, "gamma", lower=c(0,0), start=list(scale=1,shape=1))

all.reps_1 %>%  
  mutate(prob = pgamma(tot, shape = fit_1$estimate[[2]], scale = fit_1$estimate[[1]], lower.tail = TRUE,
                       log.p = FALSE)) -> all.reps_1

### Remove Outlier Samples
outliers_1 <- all.reps_1 %>% 
  filter(prob < low_dist_probability_cutoff  | tot < minimum_read_cutoff | prob > high_dist_probability_cutoff) # changed to 0.05 to save the two samples

outliers <- outliers_1
outliers %>% 
  unite(unique_sample,c("Miseq_run","sample"), remove=F) -> outliers

ASV.nested %>% 
  mutate (Step3.remove.known.contaminants.tibble.edited = purrr::map(Step3.remove.known.contaminants.tibble, 
                                                   function(.x){
                                                     .x %>%
                                                       mutate(., barcode_1_col= str_detect(seq_number,barcode_1_use)) %>%
                                                       mutate(Miseq_run = if_else(barcode_1_col=="FALSE",barcode_1_use,barcode_1_use)) %>% 
                                                       unite(unique_sample,c("Miseq_run","sample"), remove=F)
                                                   }) ) -> ASV.nested

ASV.nested %>% 
  mutate(Step.4.low.read.depth = purrr::map (Step3.remove.known.contaminants.tibble.edited, ~ filter(.,!unique_sample %in% outliers$unique_sample) %>% ungroup)) -> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step.4.low.read.depth, ~ how.many(ASVtable = .,round = "4.Low.Read.Depth"))) %>% 
  left_join(ASV.summary) %>% 
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary

#---

#---

# Cleaning Process 5: **Discarding PCR replicates with low complexity i.e. swamping ** ----------------------------------------------------------------------------------------------------------------

###Pull out Sample Read Depth

ASV.nested %>% 
  dplyr::select(Step.4.low.read.depth) %>% 
  unnest(Step.4.low.read.depth) %>% 
  group_by(unique_sample) %>%
  dplyr::summarise(tot = sum(nReads)) %>% 
  arrange(desc(tot))-> all.reps


ASV.nested %>% 
  dplyr::select(Step.4.low.read.depth) %>% 
  unnest(Step.4.low.read.depth) %>% 
  ungroup() %>% 
  left_join(hash.key_updated) %>% 
  group_by(unique_sample, sample, barcode_1_col, Miseq_run, sum.taxonomy) %>% 
  mutate(sum_reads = sum(nReads)) %>% 
  ungroup() %>% 
  left_join(all.reps) %>% 
  mutate(., prop_reads = sum_reads/tot) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>%  
  filter(., prop_reads > swamping_cutoff) %>%  
  dplyr::select(unique_sample, sum.taxonomy, prop_reads) %>% distinct() -> samples_2_discard

samples_2_discard

ASV.nested %>% 
  mutate(Step5.discard.swamped = purrr::map (Step.4.low.read.depth, ~ filter(.,!unique_sample %in% samples_2_discard$unique_sample) %>% ungroup)) -> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step5.discard.swamped, ~ how.many(ASVtable = .,round = "5.Swamped.Samples"))) %>% 
  left_join(ASV.summary) %>% 
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary

#---

#---

# Cleaning Process 6: **Site Occupancy Modeling** ----------------------------------------------------------------------------------------------------------------

##STAN Model

#Metadata test run issue: 2 sequencing runs used different names for the same site so metadata had repeated data. 
#Needed to make a metadata file that was distinct(), otherwise each sample gets duplicated in the spread
metadata %>%  dplyr::select(-Seq_number) %>%  distinct() -> metadata

#Hierarchical SOM 2 levels

####Set up data with different hierarchical levels and reduce down to unique patterns.
#Format 2 level hierarcjocal patterns
ASV.nested %>% 
  dplyr::select(Step5.discard.swamped) %>% 
  unnest(Step5.discard.swamped) %>% # This stage is a long table
  ungroup() %>%
  mutate(nReads = 1) %>% 
  left_join(metadata, by= c("sample"="New_name")) %>% 
  filter(.,!is.na(Bio_rep)) %>% 
  filter(., !is.na(Tech_rep)) %>% 
  filter(., Project =="Kenai") %>% 
  dplyr::select(seq_number, Site, Depth, Month,Bio_rep, Tech_rep, nReads) %>% 
  group_by(Site, Depth, Month) %>% 
  mutate(., K = n_distinct(Tech_rep)) -> data_som_int

data_som_int %>% 
  pivot_wider(values_from = nReads, names_from = Tech_rep, names_prefix = "Tech") %>% 
  ungroup() %>% 
  mutate(ndetections = rowSums(dplyr::select(., starts_with("Tech")),na.rm = TRUE)) %>% 
  dplyr::select(., !starts_with("Tech"))  %>% 
  pivot_wider(values_from = ndetections, names_from = Bio_rep, names_prefix = "Bio_rep") -> data_tech_summarized


#Now need to try to reduce the number of repeated tests. 
#Many species have identical patterns of presence (e.g. pattern of technical replicates a species was detected)
data_tech_summarized %>% 
  unite(repetition.level,starts_with("Bio_rep"), sep = ".", remove = FALSE) -> Pattern.of.presence_tech

#Select unique patterns
Pattern.of.presence_tech %>% 
  group_by(repetition.level) %>% 
  slice(1) -> unique_patterns

#Format for STAN Model
# K is number of technical replicates for each Biological replicate
# ndetections is the number of times a species was detected out of K trials

unique_patterns %>% 
  pivot_longer(., cols=starts_with("Bio_rep"), names_to="Bio_rep",values_to="ndetections") %>% 
  filter(., !is.na(ndetections)) %>% 
  mutate(Bio_rep =  str_replace_all(Bio_rep, "Bio_rep","")) -> unique_data

#Create unique identifier for species; for use in hierarchical modeling
#STAN model requires species names to be numbers.
SS_species <- unite(data = unique_data,
                    col = SS_species,
                    c("seq_number")
) %>% pull(SS_species)
unique_data$Species <- match(SS_species, unique(SS_species)) #index for unique site-species combinations

#Create unique identifier for combinations of site-biological_rep-species; for use in hierarchical modeling
SDS <- unite(data = unique_data,
             col = SDS,
             c("Site",   "Depth", "Month","Bio_rep", "seq_number")
) %>% pull(SDS)
unique_data$SiteDateSpecies <- match(SDS, unique(SDS)) #index for unique site-date-species combinations

#Create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = unique_data,
            col = SS,
            c("Site",   "Depth", "Month","seq_number")
) %>% pull(SS)
unique_data$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations

#Fix names to match STAN model format
unique_data %>%
  mutate(., N= ndetections) %>%
  dplyr::select(-ndetections)-> unique_data

unique_data$N <- as.numeric(unique_data$N)
##Stan Model
sink(here(working_directory,"Stan_SOM_hierarchical.stan"))
cat(
  "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    array[S] int<lower=1> Species;    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    array[S] int<lower=1> L;   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    array[S] int<lower=1> K;   // number of replicates per site (ncol)
    array[S] int<lower=0> N; // number of detections among these replicates
    array[S] int z;   // integer flag to help estimate psi parameter
}
parameters{/////////////////////////////////////////////////////////////////////
array[Nloc] real<lower=0,upper=1> psi;  //commonness parameter
array[Nspecies] real<lower=0,upper=1> p11; //true positive detection rate
array[Nspecies] real<lower=0,upper=1> p10; //false positive detection rate
}
transformed parameters{/////////////////////////////////////////////////////////////////////
}
model{/////////////////////////////////////////////////////////////////////
array[S] real p;

for (i in 1:S){
z[i] ~ bernoulli(psi[L[i]]);
p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
N[i] ~ binomial(K[i], p[i]);
}; 

//priors
psi ~ beta(2,2); 
p11 ~ beta(2,2); 
p10 ~ beta(1,10);
}
generated quantities{
array[S] real<lower=0,upper=1> Occupancy_prob;    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015

for (i in 1:S){
Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
/ ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
+ (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
);
}
}

",
fill=TRUE)
sink()

##Run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++

myHierarchicalModel <- stan(file = here(working_directory,"Stan_SOM_hierarchical.stan"), 
                            data = list(
                              S = nrow(unique_data),
                              Species = unique_data$Species,
                              Nspecies = length(unique(unique_data$Species)),
                              L = unique_data$SiteSpecies,
                              Nloc = length(unique(unique_data$SiteDateSpecies)),
                              K = unique_data$K,
                              N = unique_data$N,
                              z = ifelse(unique_data$N > 0, 1, 0)
                            ), 
                            chains = som_chains,   #number of chains
                            iter = som_iterations   #number of iterations per chain
)



myHierarchicalStanResults <- tidy(tibble(as.data.frame(myHierarchicalModel)))

#See histogram of occupancy probabilities for unique patterns of presence
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>% 
  ggplot() + geom_histogram( aes(mean)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "Unique Patterns of Presence") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_3

ggsave(plot=plot_3,here(working_directory,"Output_plots","Histogram_of_Mean_Occupancy_Probabilities_unique_patterns_of_presence.png"), device = "png", width = 12, height = 8, units = "in")

#Now Re-attach SOM to unique patterns
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>%
  separate(column, into=c("column","SiteDateSpecies"), sep="([\\[\\]])") %>% 
  pivot_wider(., names_from=column, values_from = c(mean, sd))-> myHierarchicalStanResults_wide
myHierarchicalStanResults_wide$SiteDateSpecies <- as.numeric(myHierarchicalStanResults_wide$SiteDateSpecies)

unique_data %>% 
  left_join(myHierarchicalStanResults_wide) %>% 
  group_by(seq_number,Site, Depth, Month) %>% 
  summarise(., max_Occupancy_prob = max(mean_Occupancy_prob)) -> unique_data_SOM

unique_patterns %>% 
  left_join(unique_data_SOM, by=c("seq_number","Site","Depth","Month")) %>% 
  dplyr::select(-seq_number,-Site,-Depth,-Month, -starts_with("Bio_rep"))-> unique_patterns_SOM_1

unique_patterns %>% 
  left_join(unique_data_SOM, by=c("seq_number","Site","Depth","Month")) %>% 
  dplyr::select(-seq_number,-Site,-Depth,-Month)-> unique_patterns_SOM


Pattern.of.presence_tech %>% 
  left_join(unique_patterns_SOM_1, by=c("K","repetition.level")) %>% 
  left_join(hash.key_updated)-> SOM_results

unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate(maxdetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., maxdetections == 1)  -> som_min_threshold_2_keep


#See histogram of occupancy probabilities for unique patterns of presence
SOM_results %>% 
  ggplot() + geom_histogram( aes(max_Occupancy_prob)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "All ASVs") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_4

ggsave(plot=plot_4,here(working_directory,"Output_plots","Histogram_of_Mean_Occupancy_Probabilities_all_ASVs.png"), device = "png", width = 12, height = 8, units = "in")


unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate (ndetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., ndetections ==1) %>% 
  summarise(max(max_Occupancy_prob)) -> som_min_threshold_1

unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate(maxdetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., maxdetections == 1) %>% 
  summarise(max(max_Occupancy_prob)) -> som_min_threshold_2

unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate(maxdetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., maxdetections == 3) %>% 
  summarise(min(max_Occupancy_prob)) -> som_min_threshold_3


list(som_min_threshold_1,som_min_threshold_2,som_min_threshold_3) -> som_min_thresh
names(som_min_thresh) <- c("som_min_threshold_1","som_min_threshold_2","som_min_threshold_3")

#Summarize
SOM_results %>% 
  mutate(., to_keep = case_when(max_Occupancy_prob < som_min_threshold_1$`max(max_Occupancy_prob)` ~ "Toss",
                                max_Occupancy_prob > som_min_threshold_2$`max(max_Occupancy_prob)` ~ "Keep",
                                TRUE ~ "Borderline")) %>%
  dplyr::count(to_keep) %>% 
  mutate(., total=sum(n)) %>% 
  mutate(per=paste0(round(100*n/total,2),'%')) %>% 
  dplyr::select(-total)

#Hierarchical SOM 1 level

#Format single level hierarcjocal patterns
ASV.nested %>% 
  dplyr::select(Step5.discard.swamped) %>% 
  unnest(Step5.discard.swamped) %>% # This stage is a long table
  ungroup() %>%
  mutate(nReads = 1) %>% 
  left_join(metadata, by= c("sample"="New_name")) %>% 
  filter(., is.na(Tech_rep)) %>% 
  dplyr::select(seq_number, Site, Depth, Month,Date_Sampled, Bio_rep, nReads) %>% 
  group_by(Site, Depth, Month,Date_Sampled) %>% 
  mutate(., K = n_distinct(Bio_rep)) -> data_som_int_single

data_som_int_single %>% 
  pivot_wider(values_from = nReads, names_from = Bio_rep, names_prefix = "Bio_rep") %>%
  ungroup() %>% 
  mutate (ndetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) -> data_tech_summarized_single


#Now need to try to reduce the number of repeated tests. 
#Many species have identical patterns of presence (e.g. pattern of technical replicates a species was detected)
data_tech_summarized_single %>% 
  unite(repetition.level,starts_with("Bio_rep"), sep = ".", remove = FALSE) -> Pattern.of.presence_tech_single

#Select unique patterns
Pattern.of.presence_tech_single %>% 
  group_by(repetition.level) %>% 
  slice(1) -> unique_patterns_single

#Format for STAN Model
# K is number of technical replicates for each Biological replicate
# ndetections is the number of times a species was detected out of K trials

unique_patterns_single %>% 
  filter(., !is.na(ndetections))  -> unique_data_single

#Create unique identifier for species; for use in hierarchical modeling
#STAN model requires species names to be numbers.
SS_species <- unite(data = unique_data_single,
                    col = SS_species,
                    c("seq_number")
) %>% pull(SS_species)
unique_data_single$Species <- match(SS_species, unique(SS_species)) #index for unique site-species combinations

#Create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = unique_data_single,
            col = SS,
            c("Site","Depth","Month", "Date_Sampled","seq_number")
) %>% pull(SS)
unique_data_single$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations

#Fix names to match STAN model format
unique_data_single %>%
  mutate(., N= ndetections) %>%
  dplyr::select(-ndetections)-> unique_data_single

unique_data_single$N <- as.numeric(unique_data_single$N)

sink("Stan_SOM_hierarchical_single.stan")
cat(
  "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    array[S] int<lower=1> Species;    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    array[S] int<lower=1> L;   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    array[S] int<lower=1> K;   // number of replicates per site (ncol)
    array[S] int<lower=0> N; // number of detections among these replicates
    array[S] int z;   // integer flag to help estimate psi parameter
}
parameters{/////////////////////////////////////////////////////////////////////
array[Nloc] real<lower=0,upper=1> psi;  //commonness parameter
array[Nspecies] real<lower=0,upper=1> p11; //true positive detection rate
array[Nspecies] real<lower=0,upper=1> p10; //false positive detection rate
}
transformed parameters{/////////////////////////////////////////////////////////////////////
}
model{/////////////////////////////////////////////////////////////////////
array[S] real p;

for (i in 1:S){
z[i] ~ bernoulli(psi[L[i]]);
p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
N[i] ~ binomial(K[i], p[i]);
}; 

//priors
psi ~ beta(2,2); 
p11 ~ beta(2,2); 
p10 ~ beta(1,10);
}
generated quantities{
array[S] real<lower=0,upper=1> Occupancy_prob;    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015

for (i in 1:S){
Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
/ ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
+ (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
);
}
}

",
fill=TRUE)
sink()

#run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++
myHierarchicalModel_single <- stan(file = "Stan_SOM_hierarchical_single.stan", 
                                   data = list(
                                     S = nrow(unique_data_single),
                                     Species = unique_data_single$Species,
                                     Nspecies = length(unique(unique_data_single$Species)),
                                     L = unique_data_single$SiteSpecies,
                                     Nloc = length(unique(unique_data_single$SiteSpecies)),
                                     K = unique_data_single$K,
                                     N = unique_data_single$N,
                                     z = ifelse(unique_data_single$N > 0, 1, 0)
                                   ), 
                                   chains = som_chains,   #number of chains
                                   iter = som_iterations   #number of iterations per chain
)



myHierarchicalStanResults_single <- tidy(tibble(as.data.frame(myHierarchicalModel_single)))

#See histogram of occupancy probabilities for unique patterns of presence
myHierarchicalStanResults_single %>% 
  filter(grepl("Occupancy_prob", column)) %>% 
  ggplot() + geom_histogram( aes(mean)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "Unique Patterns of Presence") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_3_single

ggsave(plot=plot_3_single,here(working_directory,"Output_plots","Histogram_of_Mean_Occupancy_Probabilities_unique_patterns_of_presence_single.png"), device = "png", width = 12, height = 8, units = "in")

#Now Re-attach SOM to unique patterns
myHierarchicalStanResults_single %>% 
  filter(grepl("Occupancy_prob", column)) %>%
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])") %>% 
  pivot_wider(., names_from=column, values_from = c(mean, sd))-> myHierarchicalStanResults_wide_single
myHierarchicalStanResults_wide_single$SiteSpecies <- as.numeric(myHierarchicalStanResults_wide_single$SiteSpecies)

unique_data_single %>% 
  left_join(myHierarchicalStanResults_wide_single) %>% 
  group_by(seq_number,Site,Depth,Month,Date_Sampled) %>% 
  summarise(., max_Occupancy_prob = max(mean_Occupancy_prob)) -> unique_data_SOM_single

unique_patterns_single %>% 
  left_join(unique_data_SOM_single, by=c("seq_number","Site","Depth","Month","Date_Sampled")) %>% 
  dplyr::select(-seq_number,-Site,-Depth,-Month,-Date_Sampled)  -> unique_patterns_SOM_single

unique_patterns_single %>% 
  left_join(unique_data_SOM_single, by=c("seq_number","Site","Depth","Month","Date_Sampled")) %>% 
  dplyr::select(-seq_number,-Site,-Depth,-Month,-Date_Sampled,-starts_with("Bio_rep"),-K)  -> unique_patterns_SOM_single_1

Pattern.of.presence_tech_single %>% 
  dplyr::select(., !starts_with("Bio_rep"), -ndetections) %>% 
  left_join(unique_patterns_SOM_single_1, by=c("repetition.level")) %>%
  left_join(hash.key_updated)-> SOM_results_single

#See histogram of occupancy probabilities for unique patterns of presence
SOM_results_single %>% 
  ggplot() + geom_histogram( aes(max_Occupancy_prob)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "All ASVs") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_4

ggsave(plot=plot_4,here(working_directory,"Output_plots","Histogram_of_Mean_Occupancy_Probabilities_all_ASVs_single.png"), device = "png", width = 12, height = 8, units = "in")


unique_patterns_SOM_single %>% 
  ungroup() %>% 
  mutate (ndetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., ndetections ==1) %>% 
  summarise(max(max_Occupancy_prob)) -> som_min_threshold_1_single

unique_patterns_SOM_single %>% 
  ungroup() %>% 
  mutate(maxdetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., ndetections == 2) %>% 
  summarise(min(max_Occupancy_prob)) -> som_min_threshold_2_single

unique_patterns_SOM_single %>% 
  ungroup() %>% 
  mutate(maxdetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>%
  filter(., ndetections == 3) %>% 
  summarise(min(max_Occupancy_prob)) -> som_min_threshold_3_single


list(som_min_threshold_1_single,som_min_threshold_2_single,som_min_threshold_3_single) -> som_min_thresh_single
names(som_min_thresh_single) <- c("som_min_threshold_1_single","som_min_threshold_2_single","som_min_threshold_3_single")

#Summarize
SOM_results_single %>% 
  mutate(., to_keep = case_when(max_Occupancy_prob < som_min_threshold_1_single$`max(max_Occupancy_prob)` ~ "Toss",
                                max_Occupancy_prob > som_min_threshold_2_single$`min(max_Occupancy_prob)` ~ "Keep",
                                TRUE ~ "Borderline")) %>% 
  dplyr::count(to_keep) %>% 
  mutate(., total=sum(n)) %>% 
  mutate(per=paste0(round(100*n/total,2),'%')) %>% 
  dplyr::select(-total)


#Save Data
SOM_results %>% 
  dplyr::select(-K, -repetition.level)-> SOM_results_to_keep_all

saveRDS(SOM_results_to_keep_all,"occupancy_results_all.RDS")

SOM_results_to_keep_all %>% 
  write.csv(paste0("Occupancy_all_",Sys.Date(),".csv"), row.names = F)

SOM_results_to_keep_all %>% 
  filter(., max_Occupancy_prob >= som_min_thresh[[som_filter_choice]][[1]]) -> SOM_results_to_keep_filtered

saveRDS(SOM_results_to_keep_filtered,"occupancy_results_filtered.RDS")

SOM_results_to_keep_filtered %>% 
  write.csv(paste0("Occupancy_all_filtered",Sys.Date(),".csv"), row.names = F)

SOM_results_single %>% 
  dplyr::select(-K, -repetition.level)-> SOM_results_to_keep_all_single

saveRDS(SOM_results_to_keep_all_single,"occupancy_results_all_single.RDS")

SOM_results_to_keep_all_single %>% 
  write.csv(paste0("Occupancy_all_single_",Sys.Date(),".csv"), row.names = F)

SOM_results_to_keep_all_single %>% 
  filter(., max_Occupancy_prob >= som_min_thresh_single[[som_filter_choice_2]][[1]]) -> SOM_results_to_keep_filtered_single

saveRDS(SOM_results_to_keep_filtered_single,"occupancy_results_filtered_single.RDS")

SOM_results_to_keep_filtered_single %>% 
  write.csv(paste0("Occupancy_all_single_filtered",Sys.Date(),".csv"), row.names = F)


SOM_results_to_keep_filtered %>% 
  mutate(Date_Sampled=NA) %>% 
  dplyr::select(seq_number,Site,Depth,Month,Date_Sampled,sum.taxonomy) -> double_hierarchy_keep

SOM_results_to_keep_filtered_single %>% 
  dplyr::select(seq_number,Site,Depth,Month,Date_Sampled,sum.taxonomy) -> single_hierarchy_keep

SOM_results_to_keep_filtered_combined <- rbind(double_hierarchy_keep,single_hierarchy_keep)

ASV.nested %>% 
  mutate (Step5.tibble.edited = purrr::map(Step5.discard.swamped, 
                                           function(.x){
                                             .x %>%
                                               left_join(hash.key_updated) %>% 
                                               dplyr::select(sample, seq_number, nReads, sum.taxonomy)
                                           }) ) -> ASV.nested
ASV.nested %>% 
  mutate(Step6.tibble = purrr::map (Step5.tibble.edited, ~ filter(.,sum.taxonomy %in% SOM_results_to_keep_filtered_combined$sum.taxonomy) %>% ungroup))-> ASV.nested


ASV.nested$Step5.tibble.edited[[1]] %>% 
  dplyr::select(seq_number) %>% distinct() -> asvs_pre_occ

ASV.nested$Step6.tibble[[1]] %>% 
  dplyr::select(seq_number) %>% distinct() -> asvs_post_occ

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step6.tibble, ~ how.many(ASVtable = .,round ="6.Occupancy"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary

ASV.summary$Summary

#---

#---

# Cleaning Process 7: **Dissimilarity between PCR (biological) replicates** ----------------------------------------------------------------------------------------------------------------

ASV.nested %>% 
  dplyr::select(Miseq_run,Step6.tibble) %>% 
  unnest(Step6.tibble) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.post_occ


ASV.nested %>% 
  dplyr::select(Miseq_run,Step5.tibble.edited) %>% 
  unnest(Step5.tibble.edited) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.pre_occ


## How many samples, how many ASVs
cleaned.tibble.post_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))

cleaned.tibble.pre_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))


## eDNA Index
cleaned.tibble.post_occ %>%
  dplyr::group_by(sample, Miseq_run) %>%
  mutate (Tot = sum(nReads),
          Normalized.reads = nReads / Tot) -> cleaned.tibble.post_occ #transforms raw number of reads to eDNA index

cleaned.tibble.pre_occ %>%
  dplyr::group_by(sample, Miseq_run) %>%
  mutate (Tot = sum(nReads),
          Normalized.reads = nReads / Tot) -> cleaned.tibble.pre_occ #transforms raw number of reads to eDNA index

#Calculate Vegan Distances
tibble_to_vegdist_all (cleaned.tibble.post_occ) -> all.distances.full.post
tibble_to_vegdist_all (cleaned.tibble.pre_occ) -> all.distances.full.pre

as.tibble(subset(melt(as.matrix(all.distances.full.post)))) -> all.distances.melted.post
as.tibble(subset(melt(as.matrix(all.distances.full.pre)))) -> all.distances.melted.pre

# Now, create a three variables for all distances, they could be PCR replicates, BIOL replicates, or from the same site
all.distances.melted.pre %>%
  separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Month = if_else(is.na(Month), Date_Sampled,Month)) %>% 
  mutate(., Site1=Site,
         Depth1=Depth,
         Month1=Month,
         Bio_rep1=Bio_rep,
         Tech_rep1=Tech_rep,
         Miseq_run1=Miseq_run,
         Project1=Project) %>% 
  dplyr::select(Var1,Miseq_run1,Site1,Depth1,Month1,Bio_rep1,Tech_rep1,Project1,Var2,value) %>% 
  separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Month = if_else(is.na(Month), Date_Sampled,Month)) %>% 
  mutate(., Site2=Site,
         Depth2=Depth,
         Month2=Month,
         Bio_rep2=Bio_rep,
         Tech_rep2=Tech_rep,
         Miseq_run2=Miseq_run,
         Project2=Project) %>% 
  dplyr::select(Var1,Var2,Site1,Depth1,Month1,Miseq_run1,Bio_rep1,Tech_rep1,Project1,Miseq_run2,Site2,Depth2,Month2,Bio_rep2,Tech_rep2,Project2,value) %>%
  unite( Site1,Depth1,Month1, Bio_rep1, col= "station1", remove=F) %>% 
  unite( Site2,Depth2,Month2, Bio_rep2, col= "station2", remove=F) %>% 
  unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
  unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
  unite( Miseq_run1, Site1,Depth1,Month1, col= "type_site1", remove=F) %>% 
  unite( Miseq_run2, Site2,Depth2,Month2, col= "type_site2", remove=F) %>% 
  ungroup() %>% 
  mutate(Distance.type = case_when(typers1 == typers2 ~ "PCR.replicates",
                                   type_site1 == type_site2 ~ "Biol.replicates",
                                   Project1==Project2  & Site1 == Site2 ~"Same Site",
                                   Project1==Project2  & Miseq_run1 == Miseq_run2 ~ "Same Project Different Site",
                                   TRUE ~ "Different Project")) %>% 
  dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
  filter (Sample1 != Sample2) -> all.distances.to.plot.pre



# Checking all went well
sapply(all.distances.to.plot.pre, function(x) summary(is.na(x)))

all.distances.to.plot.pre$Distance.type <- all.distances.to.plot.pre$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Site","Same Barcode Different Site","Different Project")

plot_5 <- ggplot (all.distances.to.plot.pre , aes (fill = Distance.type, x = value,after_stat(density))) +
  geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
  facet_wrap( ~ Distance.type) +
  labs (x = "Pairwise Dissimilarity", y = "Density" ,
        fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Pre Occupancy") +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                panel.background = element_blank())
ggsave(plot=plot_5,here(working_directory,"Output_plots","eDNA_Pairwise_Dissimilarity_Between_Samples_pre_occupancy.png"), device = "png", width = 12, height = 8, units = "in")



all.distances.melted.post %>%
  separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>%
  mutate(., Month = if_else(is.na(Month), Date_Sampled,Month)) %>% 
  mutate(., Site1=Site,
         Depth1=Depth,
         Month1=Month,
         Bio_rep1=Bio_rep,
         Tech_rep1=Tech_rep,
         Miseq_run1=Miseq_run,
         Project1=Project) %>% 
  dplyr::select(Var1,Miseq_run1,Site1,Depth1,Month1,Bio_rep1,Tech_rep1,Project1,Var2,value) %>% 
  separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Month = if_else(is.na(Month), Date_Sampled,Month)) %>% 
  mutate(., Site2=Site,
         Depth2=Depth,
         Month2=Month,
         Bio_rep2=Bio_rep,
         Tech_rep2=Tech_rep,
         Miseq_run2=Miseq_run,
         Project2=Project) %>% 
  dplyr::select(Var1,Var2,Site1,Depth1,Month1,Miseq_run1,Bio_rep1,Tech_rep1,Project1,Miseq_run2,Site2,Depth2,Month2,Bio_rep2,Tech_rep2,Project2,value) %>%
  unite( Site1,Depth1,Month1, Bio_rep1, col= "station1", remove=F) %>% 
  unite( Site2,Depth2,Month2, Bio_rep2, col= "station2", remove=F) %>% 
  unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
  unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
  unite( Miseq_run1, Site1,Depth1,Month1, col= "type_site1", remove=F) %>% 
  unite( Miseq_run2, Site2,Depth2,Month2, col= "type_site2", remove=F) %>% 
  ungroup() %>% 
  mutate(Distance.type = case_when(typers1 == typers2 ~ "PCR.replicates",
                                   type_site1 == type_site2 ~ "Biol.replicates",
                                   Project1==Project2  & Site1 == Site2 ~"Same Site",
                                   Project1==Project2  & Miseq_run1 == Miseq_run2 ~ "Same Project Different Site",
                                   TRUE ~ "Different Project"))  %>%
  dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
  filter (Sample1 != Sample2) -> all.distances.to.plot.post


# Checking all went well
sapply(all.distances.to.plot.post, function(x) summary(is.na(x)))

all.distances.to.plot.post$Distance.type <- all.distances.to.plot.post$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Site","Same Barcode Different Site", "Different Barcode")

plot_6 <- ggplot (all.distances.to.plot.post , aes (fill = Distance.type, x = value,after_stat(density))) +
  geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
  facet_wrap( ~ Distance.type) +
  labs (x = "Pairwise Dissimilarity", y = "Density" ,
        fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Post Site Occupancy") +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                      panel.background = element_blank())

ggsave(plot=plot_6,here(working_directory,"Output_plots","eDNA_Pairwise_Dissimilarity_Between_Samples_post_occupancy.png"), device = "png", width = 12, height = 8, units = "in")

# Instead of chosing based on the pairwise distances, we will use distance to centroid

# Double Hierarchy
# Need to Remove Single Samples

cleaned.tibble.post_occ  %>% 
  group_by(Miseq_run,Site,Depth,Month,Date_Sampled,Bio_rep) %>% 
  dplyr::summarise(cases = n_distinct(Tech_rep)) %>% 
  filter(., cases == 1) %>% 
  unite(col="location_id", Miseq_run,Site,Depth,Month,Date_Sampled,Bio_rep, remove = FALSE, sep=":") %>% ungroup() %>% 
  dplyr::select(location_id) -> singles
#Currently this method does not filter out these single tech reps. There are only 3

#Nest Data
cleaned.tibble.post_occ  %>% 
  unite(col="location_id", Miseq_run,Site,Depth,Month,Date_Sampled,Bio_rep, remove = FALSE, sep=":") %>% 
  filter(., !location_id %in% singles$location_id) %>% 
  group_by(location_id) %>% nest() -> nested.cleaning

nested.cleaning %>% 
  mutate(matrix = map(data, tibble_to_vegdist_bio_rep)) -> nested.cleaning

# Convert
nested.cleaning %>% mutate(ncomparisons = map(matrix, length)) -> nested.cleaning

#Calculate distances to centroid
nested.cleaning <- nested.cleaning %>% mutate(distances = map2(matrix, location_id, dist_to_centroid))
nested.cleaning %>% 
  separate(location_id, c("Miseq_run","Site","Depth","Month","Date_Sampled","Bio_rep"), sep=":") %>% 
  unnest_longer(distances) %>% 
  dplyr::select(Miseq_run,Site,Date_Sampled,Bio_rep, distances_id,distances) -> all_distances.groups #unnest data


#Calculate normal distribution of distances to centroid
normparams.step7 <- MASS::fitdistr(all_distances.groups$distances, "normal")$estimate

#Calculate Probability
probs.step7 <- pnorm(all_distances.groups$distances, normparams.step7[1], normparams.step7[2])

#Determine Outliers
outliers.step7 <- which(probs.step7>max_dist_probability_cutoff)

#Remove Outliers
discard.step7 <-names(all_distances.groups$distances[outliers.step7])

to_write_discarded.step7 <- tibble(distances_id = discard.step7,
                                   distance = all_distances.groups$distances[outliers.step7])

to_write_discarded.step7 <- to_write_discarded.step7 %>% left_join(tibble(distances_id = discard.step7,
                                                                          probs = probs.step7[outliers.step7]))
write_csv(to_write_discarded.step7 ,"step7.discared_samples.csv")


cleaned.tibble.post_occ %>% 
  unite(col="distance_id", Miseq_run,Site,Depth,Month,Date_Sampled,Bio_rep,Tech_rep, remove = FALSE, sep=":")  %>% 
  filter(., distance_id %in% discard.step7) %>% ungroup() %>% 
  dplyr::select(sample) %>% distinct() -> step7_samples_discard

# Single Hierarchy

nested.cleaning %>% 
  separate(location_id, c("Miseq_run","Site","Depth","Month","Date_Sampled","Bio_rep"), sep=":") %>% 
  unite(col="location_id", Miseq_run,Site,Depth,Month,Date_Sampled, remove = FALSE, sep=":") -> double_2_filter


#Need to Remove Single Samples
cleaned.tibble.post_occ  %>% 
  group_by(Miseq_run,Site,Depth,Date_Sampled,Month) %>% 
  dplyr::summarise(cases = n_distinct(Bio_rep)) %>% 
  ungroup() %>% 
  filter(., cases == 1) %>% 
  unite(col="location_id", Miseq_run,Site,Depth,Month,Date_Sampled,remove = FALSE, sep=":") %>% 
  ungroup() %>% 
  dplyr::select(location_id) %>% 
  filter(., !location_id %in% double_2_filter$location_id) -> singles_1
#Currently this method does not filter out these single tech reps. There are only 3

#Nest Data
cleaned.tibble.post_occ  %>% 
  unite(col="location_id", Miseq_run,Site,Depth,Month,Date_Sampled, remove = FALSE, sep=":") %>% 
  filter(., !location_id %in% singles_1$location_id) %>% 
  filter(., !location_id %in% double_2_filter$location_id)%>%
  group_by(location_id) %>% nest() -> nested.cleaning_1

nested.cleaning_1 %>% 
  mutate(matrix = map(data, tibble_to_vegdist_bio_rep_kenai)) -> nested.cleaning_1

# Convert
nested.cleaning_1 %>% mutate(ncomparisons = map(matrix, length)) -> nested.cleaning_1

#Calculate distances to centroid
nested.cleaning_1 <- nested.cleaning_1 %>% mutate(distances = map2(matrix, location_id, dist_to_centroid))

nested.cleaning_1 %>% 
  separate(location_id, c("Miseq_run","Site","Depth","Month","Date_Sampled"), sep=":") %>% 
  unnest_longer(distances) %>% 
  dplyr::select(Miseq_run,Site, Depth,Month,Date_Sampled,distances_id,distances) -> all_distances.groups_1 #unnest data


#Calculate normal distribution of distances to centroid
normparams.step7_1 <- MASS::fitdistr(all_distances.groups_1$distances, "normal")$estimate

#Calculate Probability
probs.step7_1 <- pnorm(all_distances.groups_1$distances, normparams.step7_1[1], normparams.step7_1[2])

#Determine Outliers
outliers.step7_1 <- which(probs.step7_1>max_dist_probability_cutoff)

#Remove Outliers
discard.step7_1 <-names(all_distances.groups_1$distances[outliers.step7_1])

to_write_discarded.step7_1 <- tibble(distances_id = discard.step7_1,
                                     distance = all_distances.groups_1$distances[outliers.step7_1])

to_write_discarded.step7_1 <- to_write_discarded.step7_1 %>% left_join(tibble(distances_id = discard.step7_1,
                                                                              probs = probs.step7_1[outliers.step7_1]))
write_csv(to_write_discarded.step7_1 ,"step7.discared_samples.csv")



cleaned.tibble.post_occ %>% 
  unite(col="distance_id", Miseq_run,Site,Depth,Month,Date_Sampled,Bio_rep, remove = FALSE, sep=":")  %>% 
  filter(., distance_id %in% discard.step7_1) %>% ungroup() %>% 
  dplyr::select(sample) %>% distinct() -> step7_1_samples_discard

rbind(step7_samples_discard,step7_1_samples_discard) ->discard.step7c

## Plot Final Replication Levels

cleaned.tibble.post_occ %>%
  filter(., !sample %in% discard.step7c$sample) %>% 
  dplyr::group_by(Miseq_run,Site,Depth,Month,Bio_rep,Date_Sampled) %>% 
  dplyr::summarise(Tech_reps = n_distinct(Tech_rep)) %>% 
  unite(Station,c("Site","Depth","Month","Date_Sampled"), remove=F) %>% 
  ggplot()+
  geom_raster(aes(x= Station, y = Bio_rep, fill = Tech_reps))+
  geom_text(aes(x= Station, y = Bio_rep, label = Tech_reps),color="white") +
  ggtitle("Sample Replication Level") +facet_wrap(~Miseq_run) +
  theme(axis.text.x = element_text(angle = 90))-> plot_7

ggsave(plot=plot_7,here(working_directory,"Output_plots","Final_sample_replication_level.png"), device = "png", width = 20, height = 8, units = "in")


ASV.nested %>% 
  mutate(Step7.tibble = purrr::map (Step6.tibble,
                                    ~ filter(.,!sample %in% discard.step7c$sample) %>% ungroup))-> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step7.tibble, ~ how.many(ASVtable = .,round ="7.Dissimilarity"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary

##Save Data
saveRDS(ASV.nested,file="ASV.nested_final.RDS")
saveRDS(ASV.summary,file="ASV.summary_final.RDS")

ASV.summary$Summary -> ASV.summary_table
write.csv(ASV.summary_table, file=here(working_directory,"ASV.summary_table.csv"))


#---

#---

#Generate Final Outputs ----------------------------------------------------------------------------------------------------------------

##Code for merging ASV tables

###Identify Unique Species for Sum.taxonomy merging

hash.key_updated %>% 
  distinct(.,sum.taxonomy) -> hashes_unique

hashes_unique$number <- row.names(hashes_unique)
hashes_unique$number <- paste0("taxon_",hashes_unique$number)
row.names(hashes_unique)<-hashes_unique$number

hash.key_updated %>% 
  left_join(hashes_unique, by="sum.taxonomy") -> hash.key.updated

#---

#---

###Post Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated$number %>% unique() -> total_taxa

ASV.nested$Step7.tibble[[1]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

total_kept_taxa <- barcode_1_pre_df_taxa %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_1_prop <- decostand(barcode_1_pre_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_1_prop_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)

####This gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)
####How both datasets are combined, on a common, comparable scale.

### Output Read Count Data
pre_results_reads = barcode_1_pre_df

hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

pre_results_reads$number <- rownames(pre_results_reads)

pre_results_reads %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> pre_results_reads

saveRDS(pre_results_reads,file=here(working_directory,"Output_R","post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS"))
write_csv(pre_results_reads ,here(working_directory,"Output_csv","post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv"))

### Output eDNA Index Data

combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file=here(working_directory,"Output_R","post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS"))
write_csv(combined_index ,here(working_directory,"Output_csv","post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv"))


### Twentymile Only - Post Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated$number %>% unique() -> total_taxa

ASV.nested$Step7.tibble[[1]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  filter(., str_detect(Sample,"TM")) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

total_kept_taxa <- barcode_1_pre_df_taxa %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_1_prop <- decostand(barcode_1_pre_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_1_prop_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)

####This gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)
####How both datasets are combined, on a common, comparable scale.

### Output Read Count Data
pre_results_reads = barcode_1_pre_df

hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

pre_results_reads$number <- rownames(pre_results_reads)

pre_results_reads %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> pre_results_reads

saveRDS(pre_results_reads,file=here(working_directory,"Output_R","tm_post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS"))
write_csv(pre_results_reads ,here(working_directory,"Output_csv","tm_post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv"))

### Output eDNA Index Data

combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file=here(working_directory,"Output_R","tm_post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS"))
write_csv(combined_index ,here(working_directory,"Output_csv","tm_post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv"))



### Kenai Only - Post Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated$number %>% unique() -> total_taxa

ASV.nested$Step7.tibble[[1]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  filter(., str_detect(Sample,"Bridge") | str_detect(Sample,"Dock")) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

total_kept_taxa <- barcode_1_pre_df_taxa %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_1_prop <- decostand(barcode_1_pre_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_1_prop_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)

####This gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)
####How both datasets are combined, on a common, comparable scale.

### Output Read Count Data
pre_results_reads = barcode_1_pre_df

hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

pre_results_reads$number <- rownames(pre_results_reads)

pre_results_reads %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> pre_results_reads

saveRDS(pre_results_reads,file=here(working_directory,"Output_R","kenai_post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS"))
write_csv(pre_results_reads ,here(working_directory,"Output_csv","kenai_post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv"))

### Output eDNA Index Data

combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file=here(working_directory,"Output_R","kenai_post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS"))
write_csv(combined_index ,here(working_directory,"Output_csv","kenai_post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv"))





