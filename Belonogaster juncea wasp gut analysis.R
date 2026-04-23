library(tidyverse)
library(vegan)
library(viridis)
library(gt)
library(bipartite)
#read all the files into R and combine into one dataset
plate1<-read_csv("Plate 1 ASV_filtered_read_count.csv")
plate2<-read_csv("plate 2 ASV_filtered_read_count.csv")
plate3<-read_csv("plate 3 ASV_filtered_read_count.csv")
plate4<-read_csv("plate 4 ASV_filtered_read_count.csv")
plate5<-read_csv("plate 5 ASV_filtered_read_count.csv")
plate6<-read_csv("plate 6 ASV_filtered_read_count.csv")
plate7<-read_csv("plate 7 ASV_filtered_read_count.csv")
plate8<-read_csv("plate 8 ASV_filtered_read_count (1).csv")
plate9<-read_csv("plate 9 ASV_filtered_read_count.csv")
plates<-rbind(plate1,plate2,plate3,plate4,plate5,plate6,plate7,plate8,plate9)
#separate village and sample info into separate columns
plates<- plates %>%
  separate(sample_name, into = c("plate", "village", "num1", "num2"),
           sep = "_") %>%
  mutate(sample = paste(num1, num2, sep = "_")) %>%
  select(-num1, -num2)
#filtering
filtered<- plates%>% filter(!sample == "NA_NA")
filtered<-filtered%>% filter(!village == "EXT")
filtered<-filtered%>% filter(!village == "BLANK")
filtered<-filtered%>% filter(!village == "POS")
filtered<-filtered%>% filter(!species == "Belonogaster juncea")
filtered<-filtered%>% filter(phylum == "Arthropoda")
#quick counting of how many reads were attributed to hymenoptera
wasps<-africanspecies%>% filter(order == "Hymenoptera")
nrow(wasps)
#importing GBIF dataset of african arthropod occurance records 
africa<-`0060548.260226173443078`
#creating a list of all arthropods identified in africa
africanarthropods<-unique(africa$species)
#creating a list of all genuses identified in africa
arthgenus<-unique(africa$genus)
#making those lists into data frames
africanarthropods<-as.data.frame(africanarthropods)
arthgenus<-as.data.frame(arthgenus)
#renaming columns to allow joining
africanarthropods<- africanarthropods %>% rename(species= africanarthropods)
arthgenus<- arthgenus %>% rename(genus=arthgenus)
#joining species lists to filtered dataset, those who are present in africa 
#will be maintained in  new dataset, those not present in africa will be 
#filtered out
africanspecs<-right_join(filtered, africanarthropods, by="species")
africanspecies<-right_join(filtered, arthgenus, by="genus")
#trimming dataset of the arthropods that we did not find in Belonogaster diet
africanspecs<- africanspecs%>%filter(!is.na(genus))
africanspecies<- africanspecies%>%filter(!is.na(phylum))
#checking output
nrow(table(africanspecies$genus))
view(africanspecies)
unique(africanspecs$species)
view(africanspecs)
nrow(table(filtered$species))
nrow(table(africanspecs$species))
nrow(table(africanspecies$species))

########### primer comparison ##################################################
#filtering out hymenopteran species (one primer was a wasp exclusion primer)
nowasps<-plates%>% filter(!order == "Hymenoptera")
#grouping reads by primer
nowaspGEN<-nowasps%>%filter(plate == c("plate01","plate02","plate03","plate04"))
nowaspWEX<-nowasps%>%filter(plate == c("plate05","plate06","plate07","plate08"))
#count number of attributed non-wasp reads
primerscomp<- c((nrow(nowaspGEN)),(nrow(nowaspWEX)))
#combining counts of attributed reads for each primer
primerscomp <- data.frame(
  name=c("General Primers","Wasp Exclusion Primers") ,  
  value=c(nrow(nowaspGEN),nrow(nowaspWEX)))
#plotting a comparison 
ggplot(primerscomp, aes(x = name, y = value, fill = name)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_d()

########### species accumulation curves ########################################
#making a dataset with only VILLAGE, SAMPLE AND SPECIES IN IT
SACdata<- plates%>% select(village, sample, species)
#Create species presence/absence dataset
species_pa <- SACdata %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = species, 
              values_from = Presence, 
              values_fill = 0) %>% 
  select(-village, -sample)
species_pa
#calculating overall SAC
accum <- specaccum(species_pa, method = "random")
#plotting overall SAC
plot(accum, ci.type = "polygon", ci.col = "lightblue",
     main = "Species Accumulation Curve",
     xlab = "Quadrats", ylab = "Species Richness")
#making presence/absence datasets per village
village_list <- SACdata %>% group_split(village)
for(data in village_list){
  v <- unique(data$village)
  species_pa <- data %>%
    distinct() %>%
    mutate(Presence = 1) %>%
    pivot_wider(names_from = species,
                values_from = Presence,
                values_fill = 0) %>%
    select(-village, -sample)
#converting to matrix
species_pa <- as.matrix(species_pa)
#skiping villages with insufficient data
if(nrow(species_pa) < 2 | ncol(species_pa) < 2){
  message(paste("Skipping village:", v, "- not enough samples or species"))
  next}
#species accumulation calculation, sampled randomly 
accum <- specaccum(species_pa, method = "random")
#plotting the SACs
plot(accum,
       ci.type = "polygon",
       ci.col = "lightblue",
       main = paste("Species Accumulation Curve -", v),
       xlab = "Quadrats",
       ylab = "Species Richness")}

######### species compositions analysis ########################################
#separating out dataset to add individual nest info to each sample
comp<- filtered %>%
  separate(sample, into = c("nest", "sample"), sep = "_")
#counting the number of attributed reads to each species
speccount<-table(filtered$species,filtered$class, filtered$order)
speccount<-as.data.frame(speccount)
#plotting table of the most common species
TOPJOB<- gt(head(speccount[order(-speccount$Freq),], 50)) %>% tab_header(
  title = md("**50 Most common species in Belonogaster diet**")) %>%
  data_color(Var3, method= "factor", apply_to = "fill", palette= "viridis")
#saving the table
gtsave(TOPJOB, "50commonspecies2.docx", path = NULL)
#plotting table in R
TOPJOB
#filtering
speccy<-speccount%>%filter(Freq>5)
#pie chart of all species attribution composition
ggplot(speccy, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()
#making an INSECT ONLY order chart
#making subset of data for arthropods
Arty<-filtered%>% filter(phylum == "Arthropoda")
#removing sample species from count
Arty<-Arty%>% filter(!species == "Belonogaster juncea")
#counting the number of attributed reads to each insect order
artyspeccount<-table(africanspecies$order)
artyspeccount<- as.data.frame(artyspeccount)
#plotting the composition of order in B. juncea diet
ggplot(artyspeccount, aes(x="", y=Freq, fill= Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + scale_color_viridis() +
  theme_void()
#counting the number of attributed reads to each family (not used in project)
famcount<-table(comp$family)
famcount<-as.data.frame(famcount)
#plotting
ggplot(famcount, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()
#counting the number of attributed reads to each class
classcount<-table(comp$class)
classcount<-as.data.frame(classcount)
#plotting
ggplot(classcount, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + scale_color_viridis()

############# counting pests in diet ###########################################
#looking for pests
#counting instances of pest species in the diet in filtered dataset
pests2<-c(length(which(africanspecs$species=="Helicoverpa armigera")),  
length(which(africanspecs$species=="Zonocerus variegatus")),
length(which(africanspecs$species=="Spodoptera frugiperda")),
length(which(africanspecs$species=="Maruca vitrata")),
length(which(africanspecs$species=="Spodoptera eridania")),
length(which(africanspecs$species=="Lagria hirta")) ,
length(which(africanspecs$species=="Taphronota calliparea")),
length(which(africanspecs$species=="Chrysodeixis acuta")),
length(which(africanspecs$species=="Spoladea recurvalis")), 
length(which(africanspecs$species=="Anomis involuta")), 
length(which(africanspecs$species=="Spodoptera littoralis")),
length(which(africanspecs$species=="Macrotermes herus")),
length(which(africanspecs$species=="Hymenia perspectalis")),
length(which(africanspecs$species=="Acraea eponina")),
length(which(africanspecs$species=="Paysandisia archon")))

#counting instances of pest species in the diet in raw dataset (not used in project)
pests1<-c(length(which(plates$species=="Helicoverpa armigera")),  
          length(which(plates$species=="Zonocerus variegatus")),
          length(which(plates$species=="Spodoptera frugiperda")),
          length(which(plates$species=="Maruca vitrata")),
          length(which(plates$species=="Spodoptera eridania")),
          length(which(plates$species=="Lagria hirta")) ,
          length(which(plates$species=="Taphronota calliparea")),
          length(which(plates$species=="Chrysodeixis acuta")),
          length(which(plates$species=="Spoladea recurvalis")), 
          length(which(plates$species=="Anomis involuta")), 
          length(which(plates$species=="Spodoptera littoralis")),
          length(which(plates$species=="Macrotermes herus")),
          length(which(plates$species=="Hymenia perspectalis")),
          length(which(plates$species=="Acraea eponina")),
          length(which(plates$species=="Paysandisia archon")))

#putting the names of the pests into an object
pests<-c("Helicoverpa armigera",
  "Zonocerus variegatus",
  "Spodoptera frugiperda",
  "Maruca vitrata",
  "Spodoptera eridania",
  "Lagria hirta",
  "Taphronota calliparea", 
  "Chrysodeixis acuta", 
  "Spoladea recurvalis", 
  "Anomis involuta", 
  "Spodoptera littoralis",
  "Macrotermes herus",
  "Hymenia perspectalis",
  "Acraea eponina", 
  "Paysandisia archon")

#combining the raw dataset pest species counts with names
Pests2<-cbind(pests,pests2)
Pests2<-as.data.frame(Pests2)
Pests2<-Pests2 %>% rename(Species = `pests`, Count = `pests2`)
#combining the filtered dataset pest species counts with names
Meep2<-cbind(africanspecs$order, africanspecs$species) 
Meep2<-as.data.frame(Meep2)
Meep2<-Meep2%>%rename(Species = V2, Order = V1)
#making a unified data frame with filterd
Pests2<-Pests2%>%left_join(Meep2, by = "Species")
#making table
gt_tbl2 <- gt(unique(Pests2)) %>% tab_header(
  title = md("**Pest species in Belonogaster diet**"))
#plot table
gt_tbl2

############# Effect of density on diet diversity ##############################
#comparing feeding diversity of nests in competition
#creating a per nest diversity of diet
#creating new data set combining village and nest numbers
compet<- plates %>%
  separate(ASV, into = c("plate", "village", "num1", "num2"), sep = "_") %>%
  mutate(village_nest = paste(village, num1, sep = "_")) %>%
  select(-village, -num1)
#calculating simpson's diversity(paper said it was good for comparing sites)
compet%>% 
  group_by(village_nest) %>%
  table(compet$species)%>%
  summarise(diversity=diversity(index= "simpson")) %>%
              ungroup()
#finding a per-nest diversity
#making lists of all the names of villages
villages <- unique(comp$village)
nests <- unique(comp$nest)
#making an object for the results to go into
results <- data.frame(village = character(),
                      nest = character(),
                      simpson = numeric(),
                      stringsAsFactors = FALSE)
#actual calculation of diversity of diet in villages
for (v in villages) {
  for (n in nests) {
    temp <- comp %>% dplyr::filter(village == v, nest == n)
    f (nrow(temp) > 0) {
      simpson_val <- diversity(table(temp$species), index = "simpson")
      results <- rbind(results,
                       data.frame(village = v,
                                  nest = n,
                                  simpson = simpson_val))}}}
#finding an average per-nest diversity for villages
village_diversity<-results%>% group_by(village) %>%
  summarise(mean_diversity = mean(simpson)) %>%
  ungroup()
#finding number of nests per village
number_nests<-comp%>% group_by(village) %>%
    summarise(nests = n_distinct(nest))%>% 
  ungroup()
#combining datasets
nest.diversity<-right_join(number_nests,village_diversity)
#preliminary checking of output
boxplot(nest.diversity$mean_diversity ~ nest.diversity$nests)
#final plot to show results
nest.diversity %>%
  ggplot(aes(x=factor(nests), y=mean_diversity, fill=factor(nests))) +
  geom_boxplot() + theme_classic() + 
  scale_fill_manual(values=c("darkseagreen","firebrick1", "plum1", "black",
                              "black", "black", "black")) +
  theme(legend.position="none",plot.title = element_text(size=11)) +
  xlab("Number of Nests") +
  ylab("Average Simpson Diversity")
#finding average per-nest diversity for number of nests in a village
nestno.diversity<-nest.diversity%>% group_by(nests) %>%
  summarise(mean = mean(mean_diversity)) %>% 
  ungroup()
view(nestno.diversity)
boxplot(nestno.diversity)
#statistics
nestno.divsersity_aov <- aov(nest.diversity$mean_diversity 
                             ~ nest.diversity$nests)
summary(nestno.divsersity_aov)
fitted.model <- glm(mean_diversity ~ nests, 
                    family= gaussian ,data= nest.diversity)
summary(fitted.model)
#checking assumptions
plot(nestno.divsersity_aov, pch = 16, col = "blue")
kruskal.test(nest.diversity$mean_diversity ~ nest.diversity$nests)
######
#comparison on villge diversity 
results <- comp %>%
  group_by(village, nest) %>%
  summarise( simpson = diversity(table(species), index = "simpson", .groups = "drop" )
#plotting graph
results %>%
  ggplot( aes(x=village, y=simpson, fill=village)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("village")
#ANOVA to check for statistical difference
vilcomp<-aov(results$simpson~ results$village)
summary(vilcomp)
############# NMDS analysis ####################################################
#making individual nest NMDS analysis
#creating a dataset formatted for NMDS
#creating species counts for each nest
species.n.villages3<-compet%>% group_by(village_nest) %>%
  count(species)%>%
  ungroup()
#pivot to create the wide version of the data
NMDSdat3<-species.n.villages3%>%
  pivot_wider( names_from = species,values_from = n,values_fill = 0)
#get rid of character variables
NMDSdat3 <- column_to_rownames(NMDSdat3, var = "village_nest")
#NMDS calculation
NMDS3=metaMDS(NMDSdat3,k=2)
#plotting NMDS
#creating the treatment categories for grouping on plot
temp2<-as.data.frame(NMDSdat3)
treatments<-temp2 %>%
  rownames_to_column("id") %>%
  separate(id, into = c("village", "nest"), sep = "_")
treat=treatments$village
#plotting NMDS with grouping village
ordiplot(NMDS3,type="n")
ordihull(NMDS3,groups=treat,draw="polygon",col="grey90",label=F)
#######
#plotting individual sample NMDS
#separating so that there will be individual sample code
indivs<- plates %>%
  separate(ASV, into = c("plate", "village", "num1", "num2"), sep = "_") %>%
  mutate(village_nest_indiv = paste(village, num1, num2, sep = "_")) %>%
  mutate(village = paste(village)) %>%
  select(-num1, -num2)
#filtering
indivs<-indivs%>% filter(!sample == "NA_NA")
indivs<-indivs%>% filter(!village == "EXT")
indivs<-indivs%>% filter(!village == "BLANK")
#make the individual species detection list
species.n.villages4<-indivs%>% group_by(village_nest_indiv) %>%
  count(species)%>%
  ungroup()
view(species.n.villages4)
#pivotting wider
NMDSdat4<-species.n.villages4%>%
  pivot_wider( names_from = species,values_from = n,values_fill = 0)
view(NMDSdat4)
#get rid of character variables
NMDSdat4 <- column_to_rownames(NMDSdat4, var = "village_nest_indiv")
#NMDS calculation
NMDS4=metaMDS(NMDSdat4,k=2)
#creating the treatment categories for grouping on plot
temp3<-as.data.frame(NMDSdat4)
treatments<-temp3 %>%
  rownames_to_column("id") %>%
  separate(id, into = c("village", "nest"), sep = "_")
treat=treatments$village
view(treatments)
#plotting the NMDS and combining, plotting polygons to group villages
plot(NMDS4)
ordihull(NMDS4,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(NMDS4,display="sites",col="red",air=0.01)
orditorp(NMDS4,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)
#regional compaersion NMDS
#creating a new dataset separating samples by region
Yaounde <-c("CM2","CM3","CM4","CM6","CM7","CM9","CM11")
Bafousam <- c("CM13", "CM14", "CM15","CM16","CM18","CM19")
Yaounde<- as.data.frame(Yaounde)
Yaounde<- Yaounde %>% mutate(region= "Yaounde")
Yaounde<- Yaounde %>% rename(village = Yaounde)
Bafousam<- as.data.frame(Bafousam)
Bafousam<- Bafousam %>% mutate(region= "Bafousam")
Bafousam<- Bafousam %>% rename(village = Bafousam)
regions<-rbind(Yaounde,Bafousam)
indivvies<- filtered %>%
  separate(ASV, into = c("plate", "village", "num1", "num2"), sep = "_") %>%
  mutate(village_nest_indiv = paste(village, num1, num2, sep = "_")) %>%
  mutate(village = paste(village)) %>%
  select(-num1, -num2)
indivvies<-left_join(indivvies, regions)
indivvies<- indivvies%>% mutate(region_village_nest_indiv= paste(region,village_nest_indiv, sep = "_"))
#creating a dataset that has the species compositions of each sample
species_mat <- indivvies %>%
  group_by(village_nest_indiv, species) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = species,
              values_from = n,
              values_fill = 0)
#remove character variables
NMDSdat4 <- column_to_rownames(species_mat, "village_nest_indiv")
#create the treatment groups
treat <- indivvies %>%
  distinct(village_nest_indiv, region) %>%
  filter(village_nest_indiv %in% rownames(NMDSdat4)) %>%
  arrange(match(village_nest_indiv, rownames(NMDSdat4))) %>%
  pull(region)
#NMDS calculation
NMDS4 <- metaMDS(NMDSdat4, distance = "bray", k = 2)
#PERMANOVA statistic calculation
adonis2(
  NMDSdat4 ~ treat,
  method = "bray",
  permutations = 999
)
#PERMDISP statistic calculation
bd <- betadisper(vegdist(NMDSdat4, method = "bray"), treat)
anova(bd)
############# distribution of pest predation ###################################
#quick-checking if pest species are found across all the villages. 
view(filter((plates %>% group_by(village)%>% 
  count(species) %>% 
  ungroup()), species == "Spodoptera eridania"))
#making a data set of just read counts of pest species
pesty<-plates%>% filter(species == "Helicoverpa armigera" |species == 
                       "Zonocerus variegatus"|species == 
                        "Spodoptera frugiperda"|species == 
                        "Maruca vitrata"|species == 
                        "Spodoptera eridania"|species == 
                        "Lagria hirta"|species == 
                        "Taphronota calliparea"|species == 
                        "Chrysodeixis acuta"|species == 
                        "Spoladea recurvalis"|species == 
                        "Anomis involuta"|species == 
                        "Spodoptera littoralis"|species == 
                        "Macrotermes herus"|species == 
                        "Hymenia perspectalis"|species == 
                        "Acraea eponina"|species == 
                        "Paysandisia archon")
#filtering
pesty<-pesty%>% filter(!sample == "NA_NA")
pesty<-pesty%>% filter(!village == "EXT")
pesty<-pesty%>% filter(!village == "BLANK")
pesty<-pesty%>% filter(!village == "POS")
#counting number of pest read counts per village
pestvillages<-table(pesty$village)
pestvillages<-as.data.frame(pestvillages)
view(pestvillages)
#plotting the number of pest read identifications in each village
ggplot(pestvillages, aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  scale_fill_viridis_d() +
  guides(fill = FALSE) +
  labs(x = "Village",y = "Count")
#repeating coding method but calculating a per sample value for pest reads
#essentially a percentage of total attributed reads that are pests
pestease<-plates%>% filter(!sample == "NA_NA")
pestease<-pestease%>% filter(!village == "EXT")
pestease<-pestease%>% filter(!village == "BLANK")
pestease<-pestease%>% filter(!village == "POS")
pestvillages<-pestvillages%>% rename(village = Var1)
pestease2<- pestease%>% group_by(village) %>%
  tally() %>%
  ungroup()
pestease2<-pestease2%>% right_join(pestvillages)
pestease2<-as.data.frame(pestease2)
#this is the step that creates the proportion value
pestease2<- pestease2%>%mutate(Freq/n)
#plotting
ggplot(pestease2, aes(x=village, y=Freq/n, fill=village)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  scale_fill_viridis_d() +
  guides(fill = FALSE) +
  labs(x = "Village",y = "Count")
######
#pest predation across region comparison
pestyys<-indivvies%>% filter(species == "Helicoverpa armigera" |species == 
                          "Zonocerus variegatus"|species == 
                          "Spodoptera frugiperda"|species == 
                          "Maruca vitrata"|species == 
                          "Spodoptera eridania"|species == 
                          "Lagria hirta"|species == 
                          "Taphronota calliparea"|species == 
                          "Chrysodeixis acuta"|species == 
                          "Spoladea recurvalis"|species == 
                          "Anomis involuta"|species == 
                          "Spodoptera littoralis"|species == 
                          "Macrotermes herus"|species == 
                          "Hymenia perspectalis"|species == 
                          "Acraea eponina"|species == 
                          "Paysandisia archon")

#creating regional tags for the samples
woahhh<-pestyys%>%filter(region == "Bafousam")
woahhhh<-as.data.frame(table(woahhh$species))
table(woahhh$species)
nrow(woahhhh)

hahah<-pestyys%>%filter(region == "Yaounde")
hahaha<-as.data.frame(table(hahah$species))
table(hahah$species)
nrow(hahaha)

woahhhh<-woahhhh%>%mutate(region= "Bafoussam")
hahaha<-hahaha%>%mutate(region= "Yaoundé")
boofers<-rbind(woahhhh,hahaha)
view(boofers)
#plotting graph
ggplot(boofers, aes(fill=region, y=Freq, x=Var1)) + 
  geom_bar(position='dodge', stat='identity') +
  scale_color_viridis() + 
  theme(text = element_text(size =15),
  axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x ="Pest species", y= "Count")
############## bipartite network ###############################################
#checking number of reads attributed to each  order
table(africanspecs$order)
#creating the data frame for the bipartite analysis 
#adding predator information (all = Belonogaster)
africanspecs<-africanspecs%>%mutate(predator = "Belonogaster juncea")
network<-cbind(africanspecs$species, africanspecs$predator)
network<-as.data.frame(network)
network<-network%>% rename(species = V1,predator= V2)
#counting number of attributed reads to each species 
specieses<-table(network$species)
specieses<-as.data.frame(specieses)
specieses<-specieses%>% rename(species = Var1)
head(specieses)
#combining the counts to the predator and species list
network<-network%>%left_join(specieses, by = "species")
#Pivoting wider to give summations of each species
networks<-network%>%pivot_wider(names_from = predator, values_from = Freq,
                                values_fn = dplyr::first)
#trimming unnecessary columns
networks <- networks%>%column_to_rownames("species")
#transform to format the bipartite package can understand
web <- as.matrix(networks)
#plotting the bipartite network
plotweb(web, text_size = 0.6, spacing= 0.5)
