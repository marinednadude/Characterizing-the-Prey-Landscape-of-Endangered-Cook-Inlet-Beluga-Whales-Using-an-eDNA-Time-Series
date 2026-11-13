## NGN Run QM Model
# Author: Eily Allan - but really all from Ryan Kelly and Ole Shelton 
# Person running: Zack Gold
# Last modified: 8/24/24 by Zack - minor modifications to use Beluga data

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(rstan)
library(MCMCpack) #for rdirichelet function
library(here)
library(gridExtra)
library(unikn)
library(compositions)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here("analysis","Scripts","calibrate_metabarcoding.R"))

# read in taxa table 
taxa_table_2  <- readRDS( here("analysis","Output","combined_data_long_reads_c.RDS"))

mock_2 <- readRDS(here("data","Mock_community_data","data","mifish_mock_community_data.RDS"))

mock_2 %>%
  group_by(community, tech_rep) %>% 
  mutate(., totReads=sum(nReads),
         propReads = nReads/totReads,
         tot_start = sum(start_conc_ng),
         start_prop=start_conc_ng/tot_start) %>% 
  filter(., ID_mifish %in% c("Oncorhynchus gorbuscha" ,  "Oncorhynchus keta"   ,   "Oncorhynchus kisutch"   ,  "Oncorhynchus nerka" ,     
                             "Oncorhynchus tshawytscha")) %>% 
  ggplot(aes(x=tech_rep, y=propReads, fill = ID_mifish)) + geom_col() +facet_wrap(community~Cycles)

mock_2 %>%
  filter(., ID_mifish %in% c("Oncorhynchus gorbuscha" ,  "Oncorhynchus keta"   ,   "Oncorhynchus kisutch"   ,  "Oncorhynchus nerka" ,     
                             "Oncorhynchus tshawytscha")) %>% 
  filter(., nReads >0) %>% 
  group_by(community,Cycles, tech_rep) %>% 
  dplyr::summarise(n_distinct(ID_mifish))

enviro_2 <- taxa_table_2 %>% 
  filter(., River %in% c('Twentymile')) %>% 
  unite(col="station",c(Site,Depth), sep=".") %>% 
  dplyr::select(time =Sample_date,creek=River, station =station,biol=Bio_rep,tech=Tech_rep, species=Species, Nreads= nReads) %>% 
  dplyr::mutate(tech = replace_na(tech, 1))

# only focus on overlapping species

mock_2 <- mock_2 %>% 
  filter(., Cycles ==39) %>% 
  filter(., str_detect(community,"North")) %>% 
  mutate(., CommType = case_when(str_detect(community,"Even")~"Even",
                                 str_detect(community,"Skew_1")~"Skew_1",
                                 str_detect(community,"Skew_2")~"Skew_2")) %>%
  filter(., CommType !="Skew_1") %>% 
  filter(., ID_mifish %in% c("Oncorhynchus gorbuscha" ,  "Oncorhynchus keta"   ,   "Oncorhynchus kisutch"   ,  "Oncorhynchus nerka" ,     
                           "Oncorhynchus tshawytscha")) %>%  
  group_by(community, tech_rep) %>% 
  mutate(., totReads=sum(nReads),
         propReads = nReads/totReads,
         tot_start = sum(start_conc_ng),
         start_prop=start_conc_ng/tot_start) %>% 
  dplyr::select(Community=community, CommType,Tech_Rep=tech_rep,species=ID_mifish, nReads,totReads,propReads,start_prop, N_pcr_mock=Cycles) %>% 
  mutate(Tech_Rep = as.numeric(Tech_Rep),
         N_pcr_mock = as.numeric(N_pcr_mock)) 

# Prepare for stan model 
input_metabarcoding_RDS <- enviro_2
input_mock_comm_RDS <- mock_2
format_metabarcoding_data <- function(input_metabarcoding_RDS, input_mock_comm_RDS){
  require(tidyverse)
  
  Observation <- input_metabarcoding_RDS 
  Mock <- input_mock_comm_RDS %>% 
    mutate(time = 1,    #for formatting compatibility
           creek = 1,
           biol = 1) %>% 
    rename(Nreads = nReads) %>% 
    #unite(c(Genus, Species), col = "species", sep = " ") %>% 
    unite(c(Community, CommType), col = "Community", sep = "_") %>% 
    rename(station = Community) %>% 
    replace(is.na(.), 0) %>% 
    filter(start_prop > 0) %>% 
    filter(Nreads > 0) %>% #omit things that are absent from the mocks
    filter(! is.na(Nreads)) %>% 
    rename(tech = Tech_Rep,
           b_proportion = start_prop) %>% 
    #mutate(stationname = station) %>% 
    unite(c(time, creek, station,biol, tech), col = "Sample", sep = ".", remove = F) 
  #unite(c(time, creek, stationname, biol), col = "Sample", sep = ".", remove = F) %>%
  #rename(station = stationname) 
  
  # only keep species present in both mocks and observations
  keepSpecies <- intersect(Mock$species, Observation$species)
  Observation <- Observation %>% 
    filter(species %in% keepSpecies) 
  Mock <- Mock %>% filter(species %in% keepSpecies)

  # index species to a common standard 
  sp_list <- data.frame(
    species = c(Mock$species, Observation$species) %>% unique(),
    species_idx = NA)
  sp_list$species_idx <- match(sp_list$species, unique(sp_list$species)) 
  
  #reindex and renormalize to deal with omitted species
  Mock <- Mock %>% 
    left_join(sp_list) %>% 
    mutate(station = match(station, unique(station)),
           speciesname = species,
           species = species_idx) %>% 
    group_by(station, tech, time, creek, biol) %>% 
    mutate(b_proportion = b_proportion/sum(b_proportion)) %>% 
    ungroup()
  
  Observation <- Observation %>% 
    left_join(sp_list) %>% 
    mutate(stationname = station,
           speciesname = species,
           species = species_idx) %>% 
    mutate(station = case_when(stationname == "Knik Arm.Surface" ~ 1,
                               stationname == "Little Susitna.Surface" ~ 2,
                               stationname == "Susitna.Surface" ~ 3,
                               stationname == "Twentymile.Surface" ~ 4,
                               stationname == "Bridge.Surface" ~ 5,
                               stationname == "Dock.Surface" ~ 6,
                               stationname == "Bridge.5m" ~ 7,
                               stationname == "Dock.5m" ~ 8)) %>% 
    group_by(station, time, creek,biol) %>% 
    mutate(tech = match(tech, unique(tech)),
           tot_reads=sum(Nreads)) %>% 
    filter(., tot_reads >30000) %>% 
    ungroup() %>% 
    dplyr::select(-tot_reads)
  
  station_list <- data.frame(
    station = Observation$stationname %>% unique(),
    station_idx = Observation$station %>% unique())
  
  #list object containing named elements Observation, Mock, and Npcr
  
  return(
    metabarcoding_data <- list(
      Observation = Observation,
      Mock = Mock,
      N_pcr_mock = 39, 
      NSpecies = nrow(sp_list),
      station_list = station_list,
      sp_list = sp_list
    ))
}

qmdata <- format_metabarcoding_data(enviro_2, mock_2)
obs <- qmdata

obs$Observation %>% 
  ggplot(aes(x=time, y=log10(Nreads+1))) + geom_point() +facet_wrap(station~species)

obs$Mock %>%
  ggplot(aes(x=time, y=propReads, fill = speciesname)) + geom_col() +facet_wrap(Sample~.)


N_pcr_cycles <- 39
MOCK <- qmdata$Mock   

alrTransform <- function(MOCK){
  require(tidyverse)
  require(compositions)
  
  p_mock <- MOCK %>% 
    dplyr::select(species, station, tech, b_proportion) %>% 
    rename(tech_rep = tech) %>% 
    pivot_wider(names_from = species, values_from = b_proportion, values_fill = 1e-9) %>% 
    ungroup() 
  colnames(p_mock)[3:(length(unique(MOCK$species))+2)] <- paste0("alr_", 1:length(unique(MOCK$species)))
  
  p_mock <- alr(p_mock[,3:ncol(p_mock)]) %>% as.matrix() %>% as.data.frame()
  p_mock[,length(unique(MOCK$species))] <- 0  #add reference zero expressly
  names(p_mock)[length(unique(MOCK$species))] <- paste0("alr_", length(unique(MOCK$species)))
  
  p_mock <-  cbind(MOCK %>% dplyr::select(tech, station) %>% distinct(),
                   p_mock) %>% ungroup()
  names(p_mock)[1] <- "tech_rep"
  
  
  return(p_mock)
}


makeDesign <- function(obs, #obs is a named list with elements Observation, Mock, N_pcr_mock, sp_list
                       N_pcr_cycles){ #N_pcr_cycles is the number of PCR cycles in your experimental/enviro samples; currently a single value, could be made into a vector if this number varies
  #library(tidyverse)
  library(MCMCpack)
  library(compositions)
  library(rstan)
  library(dplyr)
  
  
  mock <- obs$Mock
  observed <- obs$Observation
  
  p_mock_all <- alrTransform(mock)
  
  mock_3 <- mock %>% 
    dplyr::select(species, station, tech, Nreads) %>% 
    ungroup() %>% 
    mutate(species = paste0("sp_", species)) %>% 
    pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
  N_pcr_mock <- rep(obs$N_pcr_mock, nrow(p_mock_all)) #assumes all have the same Npcr
  
  
  p_samp_all <- observed %>% 
    ungroup() %>% 
    unite(time, creek, station, biol, col = "station") %>% 
    dplyr::select(station, tech, species, Nreads) %>% 
    rename(tech_rep = tech) %>% 
    mutate(species = paste0("sp_", species)) %>% 
    arrange(species) %>% 
    group_by(station, tech_rep,species) %>% 
    dplyr::summarise(Nreads= sum(Nreads)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
  N_pcr_samp <- rep(N_pcr_cycles, nrow(p_samp_all))
  
  ########################################################################
  #### Create data frames that can be read into Stan model
  ########################################################################
  
  NOM <- as.name(colnames(p_mock_all)[1])
  formula_a <- eval(NOM) ~ N_pcr_mock -1
  model_frame <- model.frame(formula_a, p_mock_all)
  model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()
  N_pcr_mock_small <- cbind(N_pcr_mock, p_mock_all) %>%  filter(tech_rep == 1) %>% pull(N_pcr_mock)
  formula_b <- eval(NOM) ~ N_pcr_mock_small -1
  model_frame <- model.frame(formula_b, p_mock_all%>% filter(tech_rep==1))
  model_vector_a_mock_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
  
  N_obs_mock       <- nrow(p_mock_all)
  
  # unknown communities second
  # species compositions (betas)
  
  NOM <- as.name(colnames(p_samp_all)[1])    
  
  p_samp_all$station <- as.factor(p_samp_all$station)
  N_station = length(unique(p_samp_all$station))
  p_samp_all$tech_rep <- as.factor(p_samp_all$tech_rep)
  if(N_station == 1){
    formula_b <- eval(NOM) ~ 1  
  } else {
    formula_b <- eval(NOM) ~ station
  }
  
  model_frame <- model.frame(formula_b, p_samp_all)
  model_matrix_b_samp <- model.matrix(formula_b, model_frame)
  
  # choose a single representative for each station to make predictions to
  model_frame <- model.frame(formula_b, p_samp_all %>% filter(tech_rep==1))
  model_matrix_b_samp_small <- model.matrix(formula_b, model_frame)
  
  # efficiencies (alpha)
  formula_a <- eval(NOM) ~ N_pcr_samp -1
  model_frame <- model.frame(formula_a, p_samp_all)
  model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()
  N_pcr_samp_small <- cbind(N_pcr_samp, p_samp_all) %>% filter(tech_rep == 1) %>% pull(N_pcr_samp)
  formula_b <- eval(NOM) ~ N_pcr_samp_small -1
  
  model_frame <- model.frame(formula_b, p_samp_all %>% filter(tech_rep==1))
  model_vector_a_samp_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
  
  #counters 
  N_obs_samp_small <- nrow(model_matrix_b_samp_small)
  N_obs_samp <- nrow(p_samp_all)
  N_b_samp_col <- ncol(model_matrix_b_samp)
  
  
  #### Make Stan objects
  
  stan_data <- list(
    N_species = ncol(p_samp_all)-2,   # Number of species in data
    N_obs_samp = nrow(p_samp_all), # Number of observed community samples and tech replicates ; this will be Ncreek * Nt * Nbiol * Ntech * 2 [for upstream/downstream observations]
    N_obs_mock = nrow(p_mock_all), # Number of observed mock samples, including tech replicates
    N_obs_samp_small = nrow(p_samp_all %>% filter(tech_rep == 1)), # Number of unique observed community samples ; this will be Ncreek * Nt * Nbiol * 2 [for upstream/downstream observations]
    
    # Observed data of community matrices
    sample_data = p_samp_all %>% dplyr::select(contains("sp")),
    sample_vector = unique(p_samp_all$station),
    mock_data   = mock_3 %>% dplyr::select(contains("sp")),
    sp_list = obs$sp_list,
    
    # True proportions for mock community
    #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
    alr_mock_true_prop = p_mock_all %>% dplyr::select(contains("alr")),
    
    # vectors of PCR numbers
    N_pcr_samp = N_pcr_samp,
    N_pcr_mock = N_pcr_mock,
    
    # Design matrices: field samples
    N_b_samp_col = N_b_samp_col,
    model_matrix_b_samp = model_matrix_b_samp,
    model_matrix_b_samp_small = as.array(model_matrix_b_samp_small),
    model_vector_a_samp = model_vector_a_samp,
    model_vector_a_samp_small = as.array(model_vector_a_samp_small),
    
    # Design matrices: mock community samples
    model_vector_a_mock = as.array(model_vector_a_mock),
    
    # Priors
    alpha_prior = c(0,0.5),  # normal prior
    beta_prior = c(0,5),    # normal prior
    tau_prior = c(1,2)   # gamma prior
  )
  
  return(stan_data)
  
}


stan_metabarcoding_data <- makeDesign(qmdata, N_pcr_cycles = 39)
stan_metabarcoding_data$alr_mock_true_prop
stan_metabarcoding_data$mock_data
stan_metabarcoding_data$sample_vector

####################################################################
# Run ML QM model
####################################################################
stanmodelname <- here("analysis","scripts","quant_salmon_code","quant_metabar_rosetta_noSampleEta.stan" )
QM_likelihood <- function(stanmodelname, stan_metabarcoding_data){
  M <- stan_model(stanmodelname)
  
  stanOpt <- optimizing(M, data=stan_metabarcoding_data, iter=30000,draws=0,
                        verbose=T,  
                        tol_param=1e-40,
                        algorithm="LBFGS",
                        hessian = TRUE)
  
  MLest <- stanOpt$par[grep("int_samp_small", names(stanOpt$par))] %>%
    matrix(ncol = stan_metabarcoding_data$N_species) %>% 
    as.data.frame()
  names(MLest) <- stan_metabarcoding_data$sp_list$species
  rownames(MLest) <- stan_metabarcoding_data$sample_vector
  ML_a <- stanOpt$par[grep("alpha\\[", names(stanOpt$par))]  
  ML_a <- data.frame("alpha_est" = ML_a, 
                     "species" = stan_metabarcoding_data$sp_list$species)
  
  
  return(list(
    ML_modelfit = stanOpt,
    ML_estimates = MLest,
    ML_alpha_est = ML_a
  ))
  
}

ML_out <- QM_likelihood(stanmodelname, stan_metabarcoding_data)

#write_rds(ML_out, here("Output","metabarcoding","ML_out.RDS"))

ML_out$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE, sep="_") %>% 
  mutate(station = case_when(station ==1 ~"Knik Arm.Surface",
                             station ==2 ~"Little Susitna.Surface",
                             station ==3 ~"Susitna.Surface",
                             station ==4 ~"Twentymile.Surface",
                             station ==5 ~"Bridge.Surface",
                             station ==6 ~"Dock.Surface",
                             station ==7 ~"Bridge.5m",
                             station ==8 ~"Dock.5m")) %>% 
  filter(value > 0.001) %>%
  filter(., creek %in% c('Twentymile','Kenai')) %>% 
  ggplot(aes(x = biol, fill = species, y = value)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~creek ~station ~time) 

#ggsave(here("Figures","20230225_proportions_after_qm_ML.png"))

####################################################################
# Run Bayesian QM model
####################################################################

bayes_out <- QM_bayes(stanmodelname, stan_metabarcoding_data)
#save(bayes_out, file=here("analysis","Output","20250224_bayesout_salmon_only.RData")) # both Even and Skew 2
#load(here("analysis","Output","20250224_bayesout_salmon_only.RData"))

summaryout <- summary(bayes_out$Bayes_modelfit)$summary
write.csv(summaryout, here("analysis","Output","20261105_bayesout_salmon_eul_only_summary.csv"))

bayes_out$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") -> tibble_bayes

bayes_out$Bayes_25ci %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") -> tibble_bayes_25

bayes_out$Bayes_75ci %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") -> tibble_bayes_75

tibble_bayes$Bayes_25ci <- tibble_bayes_25$value
tibble_bayes$Bayes_75ci <- tibble_bayes_75$value

pars2 <- rstan::extract(bayes_out$Bayes_modelfit, par = "int_samp_small")

species_mapping <- stan_metabarcoding_data$sp_list
Stations <-  stan_metabarcoding_data$sample_vector %>% as.tibble()
Stations %>% 
  mutate(., station_idx =row_number())-> Stations

pars2$int_samp_small %>% 
    as_tibble() %>% 
    pivot_longer(., cols = `1.1`:`58.5`, names_to="station.species", values_to = "est") %>%
    separate(station.species, into = c("station_idx","species_idx"), remove = F, sep="\\.") %>%
    mutate(., station_idx = as.numeric(station_idx),
           species_idx = as.numeric(species_idx)) %>%
     left_join(species_mapping) %>%
       left_join(Stations) -> int_samp_small_tibble

saveRDS(int_samp_small_tibble, here("analysis","Output","int_samp_small_tibble.RDS"))

mean_est <- summary(bayes_out$Bayes_modelfit, par = "int_samp_small")$summary[,1] %>%
  matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
  as.data.frame()

names(mean_est) <- stan_metabarcoding_data$sp_list$species
rownames(mean_est) <- stan_metabarcoding_data$sample_vector


tibble_bayes %>% 
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE,sep = "_") %>% 
  mutate(station = case_when(station == "Knik Arm" ~ 1,
                             station == "Little Susitna" ~ 2,
                             station == "Susitna" ~ 3,
                             station == "Twentymile" ~ 4,
                             station == "Bridge" ~ 5,
                             station == "Dock" ~ 6)) %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) -> qm_data_plotting

saveRDS(qm_data_plotting, here("analysis","Output","qm_data_plotting.RDS"))

bayes_out$Bayes_alpha_est %>% as.tibble() %>%
  mutate(alpha_mult = exp(alpha_est*39)) -> alpha_est_tib


# 1. Divide observed proportions by alpha estimates

mock_2 %>% 
  left_join(alpha_est_tib) %>% 
  group_by(Community,CommType,Tech_Rep) %>%
  # Divide observed proportion by the exponentiated alpha
  mutate(unscaled_corrected = propReads / alpha_mult) %>%
  # Re-normalize the unscaled values so they sum to 1 for each sample
  mutate(bias_corrected_prop = unscaled_corrected / sum(unscaled_corrected, na.rm = TRUE)) %>% 
  mutate(., community = str_remove(CommType,"_2"))-> mock_bias_corrected

species_colors_salmon <- c("#efb28c","#fa8775","grey50","#5f865f","#a2d2df")

mock_bias_corrected %>% 
ggplot(aes(x=Tech_Rep, y=bias_corrected_prop, fill=species)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(.~community) +
  scale_fill_manual(values = species_colors_salmon, name = "Species") + 
  ylab("% Reads") + 
  xlab("Technical Replicate") +
  labs(title = "Bias Corrected Metabarcoding") +
  theme_bw(base_size = 32) +
  theme(plot.title = element_text(hjust = 0.5),strip.text = element_text(face = "bold"), strip.background = element_rect(fill = "white", color = "black")) +
  scale_y_continuous(labels = scales::percent) -> north_stacked_bar_plot_salmon_qm

mock_bias_corrected %>% 
  ggplot(aes(x=Tech_Rep, y=propReads, fill=species)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(.~community) +
  scale_fill_manual(values = species_colors_salmon, name = "Species") + 
  ylab("% Reads") + 
  xlab("Technical Replicate") +
  labs(title = "Metabarcoding") +
  theme_bw(base_size = 32) +
  theme(plot.title = element_text(hjust = 0.5),strip.text = element_text(face = "bold"), strip.background = element_rect(fill = "white", color = "black")) +
  scale_y_continuous(labels = scales::percent) -> north_stacked_bar_plot_salmon

# 1. Completely remove the legend from the first plot
p1_salmon <- north_stacked_bar_plot_salmon + 
  guides(fill = "none")

# 2. Apply the row/column formatting strictly to the second plot's legend
p2_salmon <- north_stacked_bar_plot_salmon_qm + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5))

# 3. Combine the plots and position the single collected legend at the bottom
p1_salmon + p2_salmon + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 20, b = 0, l = 0) 
  )

ggsave(file=here("analysis","Figures","Figure_S4.3_Stacked_bar_plot_salmon_qm.png"),
       width = 22,
       height = 14)

mock_bias_corrected %>% 
  mutate(., Raw = (propReads-start_prop)/start_prop,
         `Bias Corrected` = (bias_corrected_prop-start_prop)/start_prop) %>% 
  dplyr::select(CommType, Tech_Rep, species, Raw, `Bias Corrected`) %>% 
  pivot_longer(cols = `Raw`:`Bias Corrected`, values_to = "MRE", names_to = "Status" ) -> raw_int
  
raw_int$Status <- factor(raw_int$Status, levels = c("Raw","Bias Corrected"))

raw_int %>% 
  ggplot(aes(x= MRE, group=Status, color=Status, fill=Status)) +geom_density(alpha=0.4) +
  theme_pubclean() +
  xlab("Relative Error") +
  ylab("Density") +
  theme(legend.position ="bottom", 
        axis.text = element_text(size=24), 
        axis.title = element_text(size=32, face = "bold"),
        legend.text = element_text(size=24),
        legend.title = element_text(size=32, face = "bold")) +
  scale_color_manual(values=c("orange2","dodgerblue4")) +
  scale_fill_manual(values=c("orange2","dodgerblue4"))

ggsave(here("analysis","Figures","Figure_S4.4_Figure_bias_corrected_abundance.png"), width = 12, height = 8)

pseudocount <-0.001
# Pivot your data so each sample is a row and species are columns
matrix_true <- mock_bias_corrected %>%
  ungroup() %>% 
  mutate(., Sample_ID = str_c(CommType,":",Tech_Rep)) %>% dplyr::select(Sample_ID, species, start_prop) %>% 
  pivot_wider(names_from = species, values_from = start_prop) %>% 
  column_to_rownames("Sample_ID") + pseudocount

matrix_raw <- mock_bias_corrected %>%
  ungroup() %>% 
  mutate(., Sample_ID = str_c(CommType,":",Tech_Rep)) %>% dplyr::select(Sample_ID, species, propReads) %>% 
  pivot_wider(names_from = species, values_from = propReads) %>% 
  column_to_rownames("Sample_ID") + pseudocount

matrix_corr <- mock_bias_corrected %>%
  ungroup() %>% 
  mutate(., Sample_ID = str_c(CommType,":",Tech_Rep)) %>% dplyr::select(Sample_ID, species, bias_corrected_prop) %>% 
  pivot_wider(names_from = species, values_from = bias_corrected_prop) %>% 
  column_to_rownames("Sample_ID") + pseudocount

# 2. Calculate Aitchison Distances for each sample
# vegdist(..., method = "aitchison") calculates distances between rows
# We need the distance between a raw sample and its corresponding true sample
# A simple way to get pairwise row distances is looping or mapply:

calculate_row_dist <- function(mat1, mat2) {
  sapply(1:nrow(mat1), function(i) {
    # Combine the matched sample from both matrices into a 2-row matrix
    pair <- rbind(mat1[i, ], mat2[i, ])
    # Calculate Aitchison distance between those two rows
    vegdist(pair, method = "aitchison")[1]
  })
}

# Get the vectors of error (distance to truth)
error_raw <- calculate_row_dist(matrix_raw, matrix_true)
error_corr <- calculate_row_dist(matrix_corr, matrix_true)

# Create a dataframe to view the errors side-by-side
error_df <- data.frame(
  Sample_ID = rownames(matrix_true),
  Raw_Error = error_raw,
  Corr_Error = error_corr
)

# 3. Run the Paired Wilcoxon Signed-Rank Test
# alternative = "greater" tests if Raw Error is significantly GREATER than Corrected Error
test_result <- wilcox.test(
  x = error_df$Raw_Error, 
  y = error_df$Corr_Error, 
  paired = TRUE, 
  alternative = "greater" 
)

print(test_result)

error_df %>% 
  dplyr::summarise(mean(Raw_Error), sd(Raw_Error), mean(Corr_Error), sd(Corr_Error))
