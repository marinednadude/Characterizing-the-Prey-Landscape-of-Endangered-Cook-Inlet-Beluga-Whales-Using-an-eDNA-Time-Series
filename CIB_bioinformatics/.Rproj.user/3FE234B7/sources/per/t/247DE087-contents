## NGN Merging QM results with QPCR results
# Author: Eily Allan and Ryan Kelly
# Person running: Zack
# Last modified: 2/24/25 by Zack
# Date of run:  2/24/25 by Zack

# Set up  ---------------------------------------------------------------------------

## Load packages  ---------------------------------------------------------------------------
library(tidyverse)
library(here)
library(gtools)

## Pull posteriors of Bayesian model for QM  ---------------------------------------------------------------------------

#bayes_out <- readRDS(here("analysis","20231210_bayesout_salmon_only.RDS"))
load(file=here("analysis","Output","20250224_bayesout_salmon_only.RData"))

# My latest  ---------------------------------------------------------------------------

grabthese <- grep("int_samp_small", names(bayes_out$Bayes_modelfit@sim$samples[[1]]))
postSamples <- bayes_out$Bayes_modelfit@sim$samples[[1]][grabthese]

a <- as.data.frame(postSamples) 
mynames <- names(a)

b <- t(a) %>% 
  as_tibble() %>% 
  rownames_to_column("samp") %>% 
  mutate(samp = mynames) %>% 
  mutate(samp = str_remove_all(samp,  "int_samp_small\\.")) %>% 
  separate(samp, into = c("bottle", "species", "x"), "\\.") %>% 
  dplyr::select(-x)

b <- b %>% 
  group_by(bottle, species) %>% 
  nest()

b$bottle <- rep(row.names(bayes_out$Bayes_estimates),
                times = ncol(bayes_out$Bayes_estimates))
b$species <- rep(colnames(bayes_out$Bayes_estimates),
                 each = nrow(bayes_out$Bayes_estimates))

b %>% 
  filter(., str_detect(bottle,"Twentymile")) %>% 
  filter(., species =="Oncorhynchus kisutch") 


## Pull posteriors for qPCR data for cutthroat   ---------------------------------------------------------------------------


qMod_coho <- readRDS(file=here("data","qPCR_results","data","qPCR_results","qMod_coho.RDS"))
qMod_beluga <- readRDS(file=here("data","qPCR_results","data","qPCR_results","qMod_beluga.RDS"))
qMod_chinook <- readRDS(file=here("data","qPCR_results","data","qPCR_results","qMod_chinook.RDS"))
qMod_beluga_tibble <- readRDS(, file=here("data","qPCR_results","data","qPCR_results","qMod_beluga_tibble.RDS"))
qMod_coho_tibble <- readRDS(, file=here("data","qPCR_results","data","qPCR_results","qMod_coho_tibble.RDS"))
qMod_chinook_tibble <- readRDS(, file=here("data","qPCR_results","data","qPCR_results","qMod_chinook_tibble.RDS"))
qMod_eulachon_tibble <- readRDS(, file=here("data","qPCR_results","data","qPCR_results","qMod_eulachon_tibble.RDS"))
qMod_eulachon <- readRDS(file=here("data","qPCR_results","data","qPCR_results","qMod_eulachon.RDS"))

qPCR_AK_data <- readRDS(file=here("data","qPCR_results","data","qPCR_results","qPCR_AK_data.RDS"))

qPCR_AK_data %>% 
  mutate(., Sample_date=as.Date(Sample_date,"%m/%d/%y")) %>% 
  unite(col="station1",c(Site,Depth), sep=".") %>% 
  mutate(station = case_when(station1 == "Knik Arm.Surface" ~ 1,
                             station1 == "Little Susitna.Surface" ~ 2,
                             station1 == "Susitna.Surface" ~ 3,
                             station1 == "Twentymile.Surface" ~ 4,
                             station1 == "Bridge.Surface" ~ 5,
                             station1 == "Dock.Surface" ~ 6,
                             station1 == "Bridge.5m" ~ 7,
                             station1 == "Dock.5m" ~ 8)) %>% 
  unite(col="NTC_names", c(Plate,Sample.Name,Well.Position), remove=F) %>% 
  unite(col="Control_names", c(LabID,Sample.Name,Control_Type), remove=F) %>% 
  unite(col="Sample_names", c(Sample_date,River,station,Bio_rep),remove=F) %>% 
  mutate(., better_name = case_when(is.na(Sample_date) & str_detect(Sample.Name,"standard") ~ Sample.Name,
                                    Sample.Name =="NTC" ~ NTC_names,
                                    is.na(Sample_date) & str_detect(Sample_Control,"Control") ~Control_names ,
                                    TRUE ~ Sample_names)) %>% 
  dplyr::select(Sample.Name, better_name,Volume) %>% 
  distinct() -> name_holder

##### Coho  ---------------------------------------------------------------------------
 
qMod_coho_tibble %>% 
   left_join(name_holder) %>%  
   mutate(., dilution = 1) %>% 
   dplyr::select(c(Sample=better_name,dilution,Adj_Vol=Volume)) -> coho_names
 
 qgrabthese <- grep("envir_concentration", names(qMod_coho@sim$samples[[1]]))
 qpostSamples <- qMod_coho@sim$samples[[1]][qgrabthese]

qa <- as.data.frame(qpostSamples) 
 
 qb <- t(qa) %>% 
   as_tibble() %>% 
   add_column(coho_names) %>% 
   group_by(Sample, dilution, Adj_Vol) %>% 
   nest() %>% 
   rename(cohoqpcr = data) %>% 
   rename(bottle = Sample) 
 
 mergedb <- b %>% 
   left_join(qb, by="bottle") %>% 
   drop_na() %>% 
   #replace_na(list(cutqpcr = list(-0.30103), dilution=1, Adj_Vol = 2)) %>% ## back calculate -0.30103 to get final of 50 copies / L water for samples where there are reads but no qpcr data 
   mutate(coho_dnacopy_uL = map2(.x = cohoqpcr, .y = dilution, .f = function(x, y) y*(10^x))) %>% 
   mutate(coho_dnacopy_L = map2(.x = coho_dnacopy_uL, .y = Adj_Vol, .f = function(x,y) x*(200/y))) %>% 
   ungroup()
 
 mergedb <- mergedb %>%
   group_by(bottle, species) %>% 
   mutate(meanpropreads = mean(unlist(data))) %>%
   mutate(meancohoqpcr = mean(unlist(coho_dnacopy_uL))) %>% 
   mutate(meancohodna = mean(unlist(coho_dnacopy_L))) %>% 
   separate(bottle, into=c("time","creek","station", "bio"), remove=FALSE, sep = "_")
 
  #calculate mean total DNA, only using coho
 meantotdna <- mergedb %>% 
   filter(species == "Oncorhynchus kisutch") %>% 
   mutate(meantotdna = meancohodna/meanpropreads) %>% 
   dplyr::select(bottle, meantotdna)
 
 #going old-school to get this to work
 mergedc <- mergedb
 mergedc$meantotdna <- meantotdna$meantotdna[match(mergedc$bottle, meantotdna$bottle)]
 
 #check; are all values for total dna the same?
 mergedc %>% 
   dplyr::select(bottle, species, meanpropreads, meancohodna, meantotdna) %>% 
   filter(bottle == "2021-06-26_Kenai_6_1")
 
 mergedc <- mergedc %>% 
   ungroup() %>% 
   mutate(meandnaconc = meanpropreads*meantotdna)
 
 #check; do concentration values make sense?
 mergedc %>% 
   dplyr::select(bottle, species, meanpropreads, meancohodna, meantotdna, meandnaconc) %>% 
   slice(1:10)
 
 saveRDS(mergedc, file=here("analysis","Output","coho_anchor_predicted_DNA_conc.RDS"))
 
 quants_to_plot <- mergedc %>% 
   mutate(station = case_when(station ==  1~"Knik Arm.Surface" ,
                              station ==  2~"Little Susitna.Surface" ,
                              station ==  3~"Susitna.Surface" ,
                              station == 4~"Twentymile.Surface" ,
                              station == 5~"Bridge.Surface",
                              station == 6~"Dock.Surface" ,
                              station == 7~"Bridge.5m" ,
                              station == 8~"Dock.5m" )) %>%
   separate(time, into = c("month","year"), sep = 2) %>% 
   unite(newtime, c(year,month), sep="-") 
 
 ggplot(quants_to_plot, aes(x=newtime, y=log10(meandnaconc), color=species)) +
   geom_point() +
   facet_grid(~creek ~station ~species) +
   # scale_color_manual(values = pal_okabe_ito) +
   guides(color= 'none') +
   labs(x="Date (YY-MM)", y= "Log10 copies/L water") +
   theme_bw() +
   scale_x_discrete(guide = guide_axis(angle = -45))
 
 
 ggsave(here("analysis","Figures","quant_ts_after_qm_coho_anchor_20230227.png"), units="in", width=12, height=8)
 
 ggplot(quants_to_plot, aes(x=newtime, y=log10(meantotdna))) +
   geom_point() +
   facet_grid(~creek ~station) +
   labs(x="Date (YY-MM)", y= "Log10 copies total  DNA/L water") +
   theme_bw() +
   scale_x_discrete(guide = guide_axis(angle = -45))
 
 ggsave(here("analysis","Figures","totaldna_ts_after_qm_coho_anchor_20230227.png"))
 
 
 simple <- quants_to_plot %>%
   dplyr::select(c(bottle, newtime, creek, station, bio, species, meanpropreads, meancohodna, meantotdna, meandnaconc))
 
saveRDS(simple, here("analysis","Output","coho_anchor_abs_abundance_notflowcorrected.RDS"))
 

##### Chinook  ---------------------------------------------------------------------------

qMod_chinook_tibble %>% 
  left_join(name_holder) %>%  
  mutate(., dilution =1) %>% 
  dplyr::select(c(Sample=better_name,dilution,Adj_Vol=Volume)) -> chinook_names

qgrabthese <- grep("envir_concentration", names(qMod_chinook@sim$samples[[1]]))
qpostSamples <- qMod_chinook@sim$samples[[1]][qgrabthese]

qa <- as.data.frame(qpostSamples) 

qb <- t(qa) %>% 
  as_tibble() %>% 
  add_column(chinook_names) %>% 
  group_by(Sample, dilution, Adj_Vol) %>% 
  nest() %>% 
  rename(chinookqpcr = data) %>% 
  rename(bottle = Sample) 


mergedb <- b %>% 
  left_join(qb, by="bottle") %>% 
  drop_na() %>% 
  #replace_na(list(cutqpcr = list(-0.30103), dilution=1, Adj_Vol = 2)) %>% ## back calculate -0.30103 to get final of 50 copies / L water for samples where there are reads but no qpcr data 
  mutate(chinook_dnacopy_uL = map2(.x = chinookqpcr, .y = dilution, .f = function(x, y) y*(10^x))) %>% 
  mutate(chinook_dnacopy_L = map2(.x = chinook_dnacopy_uL, .y = Adj_Vol, .f = function(x,y) x*(200/y))) %>% 
  ungroup()

mergedb <- mergedb %>%
  group_by(bottle, species) %>% 
  mutate(meanpropreads = mean(unlist(data))) %>%
  mutate(meanchinookqpcr = mean(unlist(chinook_dnacopy_uL))) %>% 
  mutate(meanchinookdna = mean(unlist(chinook_dnacopy_L))) %>% 
  separate(bottle, into=c("time","creek","station", "bio"), remove=FALSE, sep = "_")

mergedb %>% 
  dplyr::select(bottle, species, meanpropreads, meanchinookqpcr, meanchinookdna)

#calculate mean total DNA, only using chinook
meantotdna <- mergedb %>% 
  filter(species == "Oncorhynchus tshawytscha") %>% 
  mutate(meantotdna = meanchinookdna/meanpropreads) %>% 
  dplyr::select(bottle, meantotdna)

#going old-school to get this to work
mergedc <- mergedb
mergedc$meantotdna <- meantotdna$meantotdna[match(mergedc$bottle, meantotdna$bottle)]

#check; are all values for total dna the same?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meanchinookdna, meantotdna)

mergedc <- mergedc %>% 
  ungroup() %>% 
  mutate(meandnaconc = meanpropreads*meantotdna)

#check; do concentration values make sense?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meanchinookdna, meantotdna, meandnaconc) %>% 
  filter(bottle == "2020-11-25_Kenai_5_3")

saveRDS(mergedc, file=here("analysis","Output","chinook_anchor_predicted_DNA_conc.RDS"))

quants_to_plot <- mergedc %>% 
  mutate(station = case_when(station ==  1~"Knik Arm.Surface" ,
                             station ==  2~"Little Susitna.Surface" ,
                             station ==  3~"Susitna.Surface" ,
                             station == 4~"Twentymile.Surface" ,
                             station == 5~"Bridge.Surface",
                             station == 6~"Dock.Surface" ,
                             station == 7~"Bridge.5m" ,
                             station == 8~"Dock.5m" )) %>%
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") 

ggplot(quants_to_plot, aes(x=newtime, y=log10(meandnaconc), color=species)) +
  geom_point() +
  facet_grid(~creek ~station ~species) +
  # scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') +
  labs(x="Date (YY-MM)", y= "Log10 copies/L water") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45))


ggsave(here("analysis","Figures","quant_ts_after_qm_chinook_anchor_20230227.png"), units="in", width=12, height=8)

ggplot(quants_to_plot, aes(x=newtime, y=log10(meantotdna))) +
  geom_point() +
  facet_grid(~creek ~station) +
  labs(x="Date (YY-MM)", y= "Log10 copies total  DNA/L water") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45))

ggsave(here("analysis","Figures","totaldna_ts_after_qm_chinook_anchor_20230227.png"))


simple <- quants_to_plot %>%
  dplyr::select(c(bottle, newtime, creek, station, bio, species, meanpropreads, meanchinookdna, meantotdna, meandnaconc))

saveRDS(simple, here("analysis","Output","chinook_anchor_abs_abundance_notflowcorrected.RDS"))
simple$species %>%  unique()


### Eulachon ----------

qMod_eulachon_tibble %>%
  left_join(name_holder) %>%
  mutate(., dilution =1) %>%
dplyr::select(c(Sample=better_name,dilution,Adj_Vol=Volume)) -> eulachon_names

qgrabthese <- grep("envir_concentration", names(qMod_eulachon@sim$samples[[1]]))
qpostSamples <- qMod_eulachon@sim$samples[[1]][qgrabthese]

qa <- as.data.frame(qpostSamples)

qb <- t(qa) %>%
  as_tibble() %>%
  add_column(eulachon_names) %>%
  group_by(Sample, dilution, Adj_Vol) %>%
  nest() %>%
  rename(eulqpcr = data) %>%
  rename(bottle = Sample)
qb$eulqpcr[[1]] %>% dim()

eul_est_qPCR_values <- qb %>% 
  #filter(species == "Thaleichthys pacificus") %>%
  drop_na() %>%
  #replace_na(list(cutqpcr = list(-0.30103), dilution=1, Adj_Vol = 2)) %>% ## back calculate -0.30103 to get final of 50 copies / L water for samples where there are reads but no qpcr data
  mutate(eul_dnacopy_uL = map2(.x = eulqpcr, .y = dilution, .f = function(x, y) y*(10^x))) %>%
  mutate(eul_dnacopy_L = map2(.x = eul_dnacopy_uL, .y = Adj_Vol, .f = function(x,y) x*(200/y))) %>%
  ungroup()

mergedb <- eul_est_qPCR_values %>%
  mutate(., species = "Thaleichthys pacificus") %>% 
  group_by(bottle, species) %>%
  mutate(meaneulqpcr = mean(unlist(eul_dnacopy_uL))) %>%
  mutate(meaneuldna = mean(unlist(eul_dnacopy_L))) %>%
  separate(bottle, into=c("time","creek","station", "bio"), remove=FALSE, sep = "_")

mergedb %>%
  filter(species == "Thaleichthys pacificus") %>%
  filter(., str_detect(bottle,"Twentymile")) %>%
  ggplot(aes(x = as.Date(time), fill = species)) +
  geom_point(aes(y = meaneuldna)) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(species~creek, scales = "free") +
  ylab("Total DNA after correction") +
  labs(x="Date", y="Total DNA after correction", fill="Species", color="") +
  theme_bw()

saveRDS(mergedb, file=here("analysis","Output","eul_anchor_predicted_DNA_conc.RDS"))

quants_to_plot <- mergedb %>% 
  filter(!time %in% c("52385", "NA","Black Tuna")) %>%
  mutate(station = case_when(station ==  1~"Knik Arm.Surface" ,
                             station ==  2~"Little Susitna.Surface" ,
                             station ==  3~"Susitna.Surface" ,
                             station == 4~"Twentymile.Surface" ,
                             station == 5~"Bridge.Surface",
                             station == 6~"Dock.Surface" ,
                             station == 7~"Bridge.5m" ,
                             station == 8~"Dock.5m" )) %>%
  separate(time, into = c("month","year"), sep = 2) %>%
  unite(newtime, c(year,month), sep="-") 


ggplot(quants_to_plot, aes(x=newtime, y=log10(meaneuldna), color=species)) +
  geom_point() +
  facet_grid(~creek) +
  # scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') +
  labs(x="Date (YY-MM)", y= "Log10 copies/L water") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45))


ggsave(here("analysis","Figures","quant_ts_after_qm_eul_anchor_20230227.png"), units="in", width=12, height=8)

simple <- quants_to_plot %>%
  dplyr::select(c(bottle, newtime, creek, station, bio, species, meaneuldna))

saveRDS(simple, here("analysis","Output","eul_qPCR_abs_abundance_notflowcorrected.RDS"))



### Beluga ----------

qMod_beluga_tibble %>%
  left_join(name_holder) %>%
  mutate(., dilution =1) %>%
  dplyr::select(c(Sample=better_name,dilution,Adj_Vol=Volume)) -> beluga_names

qgrabthese <- grep("envir_concentration", names(qMod_beluga@sim$samples[[1]]))
qpostSamples <- qMod_beluga@sim$samples[[1]][qgrabthese]

qa <- as.data.frame(qpostSamples)

qb <- t(qa) %>%
  as_tibble() %>%
  add_column(beluga_names) %>%
  group_by(Sample, dilution, Adj_Vol) %>%
  nest() %>%
  rename(belugaqpcr = data) %>%
  rename(bottle = Sample)
qb$belugaqpcr[[1]] %>% dim()

beluga_est_qPCR_values <- qb %>% 
  drop_na() %>%
  mutate(beluga_dnacopy_uL = map2(.x = belugaqpcr, .y = dilution, .f = function(x, y) y*(10^x))) %>%
  mutate(beluga_dnacopy_L = map2(.x = beluga_dnacopy_uL, .y = Adj_Vol, .f = function(x,y) x*(200/y))) %>%
  ungroup()

mergedb <- beluga_est_qPCR_values %>%
  mutate(., species = "Delphinapterus leucas") %>% 
  group_by(bottle, species) %>%
  mutate(meanbelugaqpcr = mean(unlist(beluga_dnacopy_uL))) %>%
  mutate(meanbelugadna = mean(unlist(beluga_dnacopy_L))) %>%
  separate(bottle, into=c("time","creek","station", "bio"), remove=FALSE, sep = "_")

mergedb %>%
  filter(species == "Delphinapterus leucas") %>%
  filter(., str_detect(bottle,"Twentymile")) %>%
  ggplot(aes(x = as.Date(time), fill = species)) +
  geom_point(aes(y = meanbelugadna)) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(species~creek, scales = "free") +
  ylab("Total DNA after correction") +
  labs(x="Date", y="Total DNA after correction", fill="Species", color="") +
  theme_bw()

saveRDS(mergedb, file=here("analysis","Output","beluga_anchor_predicted_DNA_conc.RDS"))

quants_to_plot <- mergedb %>% 
  filter(!time %in% c("52385", "NA","Black Tuna")) %>%
  mutate(station = case_when(station ==  1~"Knik Arm.Surface" ,
                             station ==  2~"Little Susitna.Surface" ,
                             station ==  3~"Susitna.Surface" ,
                             station == 4~"Twentymile.Surface" ,
                             station == 5~"Bridge.Surface",
                             station == 6~"Dock.Surface" ,
                             station == 7~"Bridge.5m" ,
                             station == 8~"Dock.5m" )) %>%
  separate(time, into = c("month","year"), sep = 2) %>%
  unite(newtime, c(year,month), sep="-") 
  

ggplot(quants_to_plot, aes(x=newtime, y=log10(meanbelugadna), color=species)) +
  geom_point() +
  facet_grid(~creek) +
  # scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') +
  labs(x="Date (YY-MM)", y= "Log10 copies/L water") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45))


ggsave(here("analysis","Figures","quant_ts_after_qm_beluga_anchor_20230227.png"), units="in", width=12, height=8)

simple <- quants_to_plot %>%
  dplyr::select(c(bottle, newtime, creek, station, bio, species, meanbelugadna))

saveRDS(simple, here("analysis","Output","beluga_qPCR_abs_abundance_notflowcorrected.RDS"))


