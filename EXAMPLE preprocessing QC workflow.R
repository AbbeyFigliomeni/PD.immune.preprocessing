# FLOW CYTOMETRY POOR QUALITY SAMPLE IDENTIFICATION AND EXCLUSION ----
# ABBEY FIGLIOMENI
# 2022-11-25

# Libraries ----
library(CytoML)
library(dplyr)
library(flowWorkspace)
library(flowCore)
library(stringr)
library(ggcyto)
library(ggplot2)

setwd("[working directory]")

# Reading FULL Panel 1 workspace file (.wsp) ----
ws <- open_flowjo_xml("20221123 gDT mono Treg preprocessing.wsp")
ws  # shows groups in workspace, group ID on left

# Extracting sample info from WS ----
sample_info <- fj_ws_get_samples(ws, group_id = 8) # can view samples contained within certain group
sample_info

# import entire workspace GROUP with FCS files ----
# NOTE - SAMPLES AND WSP MUST BE IN SAME FILE DIRECTORY AND THERE CAN BE NO FCS FILES WITH SAME NAME
gs <- flowjo_to_gatingset(ws, name = "full stains and FMOs", subset = grep("full stain", sample_info$name))
sampleNames(gs)
length(sampleNames(gs))

# change sample names (derived from $FIL) to sample ID
sampleNames(gs) <- as.character(str_extract(sampleNames(gs), "[0-9][0-9][0-9][0-9]")) 
sampleNames(gs)
pData(gs)$ID <- sampleNames(gs)
pData(gs)

# visualize gating nodes
gs_get_pop_paths(gs, path="auto")
plot(gs)

# CRRITERIA 1 & 2: samples with viability <25%, CD3/CD14 viability <90% ----
exclusion_df <- pData(gs)
sample_list <- row.names(exclusion_df)  
sample_list

pop_list <- c("Live", "liveCD3pos", "liveCD14pos")

for (s in sample_list){
  for (p in pop_list){
    stat <- gh_pop_get_stats(gs[[s]], nodes=p, type="percent")
    exclusion_df[s, p] <- stat$percent*100
  }
}


# CRITERIA 3: visualise and identify  abnormal cell morphology ----
autoplot(gs[1:12], "PanScatter") + facet_wrap(~ID) # 2023, 3137, 2070 look abnormal
autoplot(gs[13:24], "PanScatter") + facet_wrap(~ID) # 2008 abnormal
autoplot(gs[25:36], "PanScatter") + facet_wrap(~ID)
autoplot(gs[37:48], "PanScatter") + facet_wrap(~ID)
autoplot(gs[49:60], "PanScatter") + facet_wrap(~ID) 

# subset dataframe based on criteria 1/2/3
exclusion_df$Morphology <- "ok"
poor_morph <- c("2070", "3137", "2008", "2023", "2041") 
exclusion_df[exclusion_df$ID %in% poor_morph,]$Morphology <- "bad"

poor_samples <- subset(exclusion_df, Live <25 | liveCD3pos <90 | liveCD14pos <90 | Morphology == "bad")
poor_samples


# CONCLUSIONS ----

## CHECK IN FLOWJO - SEE PDF FOR STAINING ##
poor_samples$Staining <- "ok"
poor_staining <- c("3137", "2070", "2008", "2023")
poor_samples[poor_samples$ID %in% poor_staining, ]$Staining <- "bad"


poor_samples$Exclude <- "N"
exclude_samples <- c("2070", "3137", "2008", "2023")
poor_samples[poor_samples$ID %in% exclude_samples, ]$Exclude <- "Y"  # all excluded met 2 or more of the criteria. Keep watch on quality of 2041, 2039, 3162, 3139
excluded <- subset(poor_samples, Exclude == "Y")

# Save file for later reference ---- 
setwd("[working directory]")

write.csv(poor_samples, "20221125 gDT Panel 1 Poor Quality Exclusion Decision FINAL.csv")
write.csv(excluded, "20230110 gDT Panel 1 Panel 1 Excluded Samples.csv")


