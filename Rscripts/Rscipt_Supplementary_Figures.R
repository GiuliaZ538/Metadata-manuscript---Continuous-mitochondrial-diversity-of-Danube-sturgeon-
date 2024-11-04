#packages
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(readr)
library(scales)
library(data.table)
library(forcats)

#1.Mapping statistics

##Figure S1.1a,b
###load datasets
raw_reads <- read.csv("/path/to/Table_rawreads.csv", sep = ",")
data_duplicates <- read.csv("/path/to/Duplicates.csv", sep = ",")

###plot
rawreads <- raw_reads$total_raw_reads_pair1
adaptremoved <-raw_reads$After_PolyT_Remnants
sample <- raw_reads$SampleID

raw2 <- data.frame(rawreads, adaptremoved, sample)
raw3 <- melt(raw2, id.vars='sample')
head(raw3)

a1.1 <- raw3 %>%
  mutate(sample = fct_relevel(sample,
                              "PA-01","PA-04","PA-06","PA-11","PA-16","LV-19","LV-24","LV-28","LV-29","KU-32","KU-34","KU-35","DB-38","DB-39","VM-43","VM-44","DK-47","DK-49","DP-52","DP-53","ST-54","ST-55","ST-56"))%>%
  ggplot(aes(x=sample, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge')+
  scale_fill_manual(values = c("rawreads" = "#104E8B", "adaptremoved" ="#CD6600"), 
    name = "Reads",labels = c("rawreads" = "Total","adaptremoved" = "After adapter removal"))+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
a1.1

b1.1 <- data_duplicates %>%
  mutate(SampleID = fct_relevel(SampleID,
                                "PA-01","PA-04","PA-06","PA-11","PA-16","LV-17", "LV-19", "LV-22", "LV-24", "LV-28", "LV-29", "KU-32", "KU-34", "KU-35", "DB-38", "DB-39", "VM-43", "VM-44", "DK-47", "DK-49","DP-52","DP-53", "ST-54","ST-55","ST-56"))%>%
  ggplot(aes(y=Percentage_dupl, x=SampleID)) +
  geom_point(color = "blue", size =3) + 
  geom_segment(aes(x=SampleID, xend=SampleID, y=0, yend=Percentage_dupl))+
  geom_text(aes(label = Percentage_dupl), vjust = -0.7, color = "black")+
  labs(x = "Sample ID", y = "Percentage of duplicates (%)")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
b1.1

pdf(file="FigureS1.1.pdf", width = 12, height = 8)
grid.arrange(a1.1, b1.1,
             ncol = 2, nrow = 1)
dev.off()

##Figure S1.2a,b,c
Coverage_data <- read.csv("/path/to/Coverage_depth_Bwa_Bowtie.csv", sep = ",")

Coverage_data <- Coverage_data %>%
  mutate(Type = factor(Type, 
                       levels = c("Bowtie_pad", "Bwa_pad"), 
                       labels = c("Bowtie", "BWA"))) 
###Depth of coverage
a1.2 <- Coverage_data %>%
  mutate(SampleID = fct_relevel(SampleID,
                                "PA-01","PA-04","PA-06","PA-11","PA-16","LV-17", "LV-19", "LV-22", "LV-24", "LV-28", "LV-29", "KU-32", "KU-34", "KU-35", "DB-38", "DB-39", "VM-43", "VM-44", "DK-47", "DK-49","DP-52","DP-53", "ST-54","ST-55","ST-56"))%>%
  ggplot(aes(y=SampleID, x=Coverage_Depth, fill=Type)) +
  geom_bar(stat='identity', position='dodge')+
  labs(x = "Depth of coverage (%)", y = "Sample ID")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")
a1.2

###Breadth of coverage
b1.2 <- Coverage_data %>%
  mutate(SampleID = fct_relevel(SampleID,
                                "PA-01","PA-04","PA-06","PA-11","PA-16","LV-17", "LV-19", "LV-22", "LV-24", "LV-28", "LV-29", "KU-32", "KU-34", "KU-35", "DB-38", "DB-39", "VM-43", "VM-44", "DK-47", "DK-49","DP-52","DP-53", "ST-54","ST-55","ST-56"))%>%
  ggplot(aes(y=SampleID, x=Breadth_Coverage, fill=Type)) +
  geom_bar(stat='identity', position='dodge')+
  labs(x = "Breadth of coverage (%)", y = "Sample ID")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")
b1.2

###Covered bases
c1.2 <- Coverage_data %>%
  mutate(SampleID = fct_relevel(SampleID,
                                "PA-01","PA-04","PA-06","PA-11","PA-16","LV-17", "LV-19", "LV-22", "LV-24", "LV-28", "LV-29", "KU-32", "KU-34", "KU-35", "DB-38", "DB-39", "VM-43", "VM-44", "DK-47", "DK-49","DP-52","DP-53", "ST-54","ST-55","ST-56"))%>%
  ggplot(aes(y=SampleID, x=Cov_bases, fill=Type)) +
  geom_bar(stat='identity', position='dodge')+
  labs(x = "Covered bases (bp)", y = "Sample ID")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c1.2

pdf(file="FigureS1.2.pdf", width = 12, height = 12)
grid.arrange(a1.2, b1.2, c1.2, 
             ncol = 2, nrow = 2)
dev.off()

##Figure S1.3 a,b
###Merge context_data with dataframe(adding species information)
Coverage_data$new <- context_data$Species[match(Coverage_data$SampleID, context_data$SampleID)]
names(Coverage_data)[names(Coverage_data) == 'new'] <- 'Species'

a1.3 <- Coverage_data %>%
  mutate(SampleID = fct_relevel(SampleID,
                                "PA-01","PA-04","PA-06","PA-11","PA-16","LV-17", "LV-19", "LV-22", "LV-24", "LV-28", "LV-29", "KU-32", "KU-34", "KU-35", "DB-38", "DB-39", "VM-43", "VM-44", "DK-47", "DK-49","DP-52","DP-53", "ST-54","ST-55","ST-56"))%>%
  ggplot(aes(y=Breadth_Coverage, x=SampleID)) +
  facet_wrap(~Species, scales = "free", ncol = 2) +
  geom_bar(stat="identity", position = "stack", width=0.5)+
  labs(x = "Sample ID", y = "Breadth of Coverage (%)")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
a1.3 

b1.3 <- Coverage_data %>%
  mutate(SampleID = fct_relevel(SampleID,
                                "PA-01","PA-04","PA-06","PA-11","PA-16","LV-17", "LV-19", "LV-22", "LV-24", "LV-28", "LV-29", "KU-32", "KU-34", "KU-35", "DB-38", "DB-39", "VM-43", "VM-44", "DK-47", "DK-49","DP-52","DP-53", "ST-54","ST-55","ST-56"))%>%
  ggplot(aes(y=Coverage_Depth, x=SampleID)) +
  facet_wrap(~Species, scales = "free", ncol = 2) +
  geom_bar(stat="identity", position = "stack", width=0.5)+
  labs(x = "Sample ID", y = "Depth of Coverage (%)")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
b1.3 

pdf(file="FigureS1.3.pdf", width = 12, height = 12)
grid.arrange(a1.3,b1.3, 
             ncol = 2, nrow = 1)
dev.off()

#2.Primary statistics from metaDMG analysis
##Figure S2.1 Damage and mean fragment length against period

###Loading files
output_metadmg <- read_csv("/path/to/final_output_metadmg.csv")
df1 <- output_metadmg

context_metadata <- read.csv ("/path/to/metadata.csv", sep = ",")

###Merge context_data and depth_data with dataframe(adding new column)
df1$new <-context_metadata$Site[match(df1$sample, context_metadata$SampleID)]
names(df1)[names(df1) == 'new'] <- 'Site'

df1$new <- context_metadata$Period[match(df1$sample, context_metadata$SampleID)]
names(df1)[names(df1) == 'new'] <- 'Period'

df1$new <- context_metadata$Species[match(df1$sample, context_metadata$SampleID)]
names(df1)[names(df1) == 'new'] <- 'Species'

###Rename long names
df1$sample[df1$sample == "no_dupl_sort_DB-39-stellatus_merged"] <- "DB-39"
df1$sample[df1$sample == "no_dupl_sort_DB-38-Lib-BLEDD-gueld_S17_merged"] <- "DB-38"
df1$sample[df1$sample == "no_dupl_sort_DK-49-Lib-DD-Acipens_S24_merged"] <- "DK-49"
df1$sample[df1$sample == "no_dupl_sort_DP-52-Lib-BLEDD-Acipenser_S25_merged"] <- "DP-52"
df1$sample[df1$sample == "no_dupl_sort_DP-53-Lib-BLEDD-gueld_S26_merged"] <- "DP-53"
df1$sample[df1$sample == "no_dupl_sort_KU-32-Lib-DD-Huso_S14_merged"] <- "KU-32"
df1$sample[df1$sample == "no_dupl_sort_KU-34-Lib-BLEDD-stellatus_S15_merged"] <- "KU-34"
df1$sample[df1$sample == "no_dupl_sort_KU-35-Lib-DD-Acipens_S16_merged"] <- "KU-35"
df1$sample[df1$sample == "no_dupl_sort_LV-19-Lib-DD-Acipens_S9_merged"] <- "LV-19"
df1$sample[df1$sample == "no_dupl_sort_LV-24-Lib-DD-Huso_S11_merged"] <- "LV-24"
df1$sample[df1$sample == "no_dupl_sort_LV-28-Lib-DD-Huso_S12_merged"] <- "LV-28"
df1$sample[df1$sample == "no_dupl_sort_LV-29-Lib-DD-Acipenser_S13_merged"] <- "LV-29"
df1$sample[df1$sample == "no_dupl_sort_PA-1-Lib-DD-huso_S1_merged"] <- "PA-01"
df1$sample[df1$sample == "no_dupl_sort_PA-11-Lib-DD-Acipens_S4_merged"] <- "PA-11"
df1$sample[df1$sample == "no_dupl_sort_PA-16-Lib-DD-gueld_S6_merged"] <- "PA-16"
df1$sample[df1$sample == "no_dupl_sort_PA-4-Lib-DD-gueld_S2_merged"] <- "PA-04"
df1$sample[df1$sample == "no_dupl_sort_PA-6-Lib-DD-stellatus_S3_merged"] <- "PA-06"
df1$sample[df1$sample == "no_dupl_sort_ST-54-Lib-BLEDD-gueld_S27_merged"] <- "ST-54"
df1$sample[df1$sample == "no_dupl_sort_ST-55-Lib-BLEDD-Huso_S28_merged"] <- "ST-55"
df1$sample[df1$sample == "no_dupl_sort_ST-56-Lib-DD-Huso_S29_merged"] <- "ST-56"
df1$sample[df1$sample == "no_dupl_sort_VM-43-Lib-BLEDD-Huso_S21_merged"] <- "VM-43"
df1$sample[df1$sample == "no_dupl_sort_VM-44-Lib-BLEDD-Huso_S22_merged"] <- "VM-44"
df1$sample[df1$sample == "no_dupl_sort_DK-47-Lib-BLEDD-Huso_S23_merged"] <- "DK-47"
       
output_metadmg_short <- df1

###plot
pd <- position_dodge(0.5)
a2.1 <- output_metadmg_short %>%
  mutate(Period = fct_relevel(Period,
                              "Mesolithic","Mesolithic-Neolithic Transition","Neolithic","Roman","Late Roman","High Medieval","Late Medieval")) %>%
  ggplot(aes(x=D_max, y=Period, colour=Site)) + geom_errorbar(aes(xmin=D_max-D_max_std, xmax=D_max+D_max_std, colour=Site, group=Site), position=pd) +
  geom_point(size = 2.2, position=pd)+
  scale_colour_manual(values =c("#CD5B45", "#EEB422", "#458B00", "#53868B", "#7AC5CD", "#104E8B", "#68228B", "#EEA2AD")) + 
  theme_minimal()+
  theme(legend.position = "none")
a2.1

b2.1 <- output_metadmg_short %>%
  mutate(Period = fct_relevel(Period,
                              "Mesolithic","Mesolithic-Neolithic Transition","Neolithic","Roman","Late Roman","High Medieval","Late Medieval")) %>%
  ggplot(aes(x=mean_L, y=Period, colour=Site)) + geom_errorbar(aes(xmin=mean_L-std_L, xmax=mean_L+std_L, colour=Site, group=Site), position=pd) +
  geom_point(size = 2.2, position=pd)+ 
  scale_colour_manual(values =c("#CD5B45", "#EEB422", "#458B00", "#53868B", "#7AC5CD", "#104E8B", "#68228B", "#EEA2AD")) + theme_minimal()
b2.1

pdf(file="FigureS2.1.pdf", width = 12, height = 8)
grid.arrange(a2.1, b2.1,
             ncol = 2, nrow = 1)
dev.off()

##Figure S2.2 Damage and mean fragment length (only Mesolithic and Neolithic sites)
###filter dataset
output_metadmg_filtered <- output_metadmg_short %>% filter(N_reads > 100, grepl("Mesolithic", Period) | grepl("Neolithic", Period))

###plot
a2.2 <- ggplot(output_metadmg_filtered, aes(x=D_max, y=sample, colour=Site))+ 
geom_errorbar(aes(xmin=D_max-D_max_std, xmax=D_max+D_max_std, colour=Site, group=Site), position=pd)+
geom_point(size = 2.2, position=pd)+
scale_colour_manual(values =c("#CD5B45", "#EEB422", "#458B00", "#53868B", "#7AC5CD", "#104E8B", "#68228B", "#EEA2AD"))+ 
theme_minimal()+
theme(legend.position = "none")
a2.2

b2.2 <- ggplot(output_metadmg_filtered, aes(x=mean_L, y=sample, colour=Site))+ 
geom_errorbar(aes(xmin=mean_L-std_L, xmax=mean_L+std_L, colour=Site, group=Site), position=pd)+
geom_point(size = 2.2, position=pd)+ 
scale_colour_manual(values =c("#CD5B45", "#EEB422", "#458B00", "#53868B", "#7AC5CD", "#104E8B", "#68228B", "#EEA2AD")) +
theme_minimal()
b2.2

pdf(file="FigureS2.2.pdf", width = 12, height = 8)
grid.arrange(a2.2, b2.2,
             ncol = 2, nrow = 1)
dev.off()



