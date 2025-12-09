### The purpose of this script is use SnapMine to check PSIs of CEs in public RNA-Seq datasets
### SnapMine package is still in development but is functional. It can be downloaded from GitHub.

wdGit <- "/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/NMDi_TDP43kd_repo/"
wdPersonal <- "/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Projects/NMD/NotGitHub/"

### Install Packages
update.packages("BiocManager")
library("BiocManager")
packages <- c("tidyverse","ggplot2","data.table","biomaRt", "pheatmap","viridis", "progress")

for(i in packages){
  if(!(i %in% rownames(installed.packages()))){
    BiocManager::install(pkgs = i, ask = F)
  } else{
    BiocManager::install(pkgs = i, ask = F, update=T)
  }
  library(package = i, character.only = T, quietly = T)
  rm(i)
}

rm(packages)

## Install SnapMine package from GitHub
remotes::install_github("irikas/SnapMine_Paper", subdir = "snapmine")
library(snapmine)

## List Package Versions
versions = installed.packages()[names(sessionInfo()$otherPkgs), "Version"]
data.frame(versions)
rm(versions)

## versions
# clusterProfiler     4.16.0
# ggstats             0.11.0
# ggpubr               0.6.2
# snapmine        0.0.0.9000
# progress             1.2.3
# viridis              0.6.5
# viridisLite          0.4.2
# pheatmap            1.0.13
# biomaRt             2.64.0
# data.table          1.17.8
# lubridate            1.9.4
# forcats              1.0.1
# stringr              1.6.0
# dplyr                1.1.4
# purrr                1.2.0
# readr                2.1.6
# tidyr                1.3.1
# tibble               3.3.0
# ggplot2              4.0.1
# tidyverse            2.0.0
# BiocManager        1.30.27


### Load files and run Snaptron query - Copied from introduction vignette!

## Load df
df <- fread(paste0(wdGit,"/data/251126_CEdf.csv"))

## Original df did not have strand info so need to use biomaRt
## Update df
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
info <- getBM(attributes = c("ensembl_gene_id", "strand"), values = df$`Ensembl ID`, mart = mart)
df <- df %>%
  left_join(info, by = join_by(`Ensembl ID` == ensembl_gene_id)) %>%
  mutate(strand = ifelse(strand == 1, "+", "-"))
rm(mart, info)

## Create info_df file
info_df <- data.frame(
  novel_junc_id = paste(df$`Gene Symbol`, df$`ID #`, sep = "_"),
  compilation = "sra_human",
  strand = df$strand,
  novel_junc_left_coord = df$LeftJunction,
  novel_junc_right_coord = df$RightJunction,
  canon_junc_coord = df$`Canonical Splice Junction`
)

saveRDS(info_df, "/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/NMDi_TDP43kd_repo/data/251206_SnapMine_info_df.RDS")

## Load data.frames with SnapMine results - as generated through intro vignette
#result_df_long <- readRDS("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Projects/NMD/NotGitHub/output_tables/251204_snapmine_results.RDS")
result_df_wide <- readRDS("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Projects/NMD/NotGitHub/output_tables/251204_snapmine_result_df_wide.RDS")

### Visualize CE inclusion across human SRA samples
## Subset to samples with avgPSI > 15 for at least one CE
## Subset to avgPSI > 5 for at least 10% tested CEs
result_df_wide_subset <- result_df_wide %>%
  mutate(maxPSI = apply(result_df_wide[, -1],
                        MARGIN = 1, function(x) max(x, na.rm = T))) %>%
  filter(maxPSI > 5) %>% dplyr::select(-maxPSI)


%>%
  mutate(testFilt = apply(result_df_wide_subset[, -1],
                    MARGIN = 1,
                    function(x) length(which(x > 5)))) %>%
  filter(testFilt > round(0.10 * nrow(info_df), 0)) %>%
  dplyr::select(-testFilt) %>%
  column_to_rownames(var = "sampleID")


result_df_wide_subset_overall <- result_df_wide_subset

m.results <- as.matrix(result_df_wide_subset_overall)
pheatmap(m.results, show_rownames = F, fontsize = 10,
         color = viridis(100,option="B"))

## Remove CEs with no values higher than 5%
colKeep <- apply(X = result_df_wide_subset_overall, MARGIN = 2, FUN = function(x) max(x) > 15)
result_df_wide_subset_overall_CEsubset <- result_df_wide_subset_overall[,colKeep]

pheatmap(as.matrix(result_df_wide_subset_overall_CEsubset), show_rownames = F, fontsize = 8,
         color = viridis(100,option="B"), legend_breaks = c(0,20,40,60,80,100))

# CE with medium-high inclusion
## Subset to avgPSI > 15 for at least 20
result_df_TDPneg <- result_df_wide_subset_overall%>%
  mutate(n5 = apply(result_df_wide_subset_overall,
                    MARGIN = 1,
                    function(x) length(which(x > 15)))) %>%
  filter(n5 > 20) %>%
  dplyr::select(-n5)

colKeep <- apply(X = result_df_TDPneg, MARGIN = 2, FUN = function(x) max(x) > 50)
result_df_TDPneg <- result_df_TDPneg[,colKeep]
pheatmap(as.matrix(result_df_TDPneg),
         #show_rownames = F,
         fontsize = 8,
         color = viridis(100,option="B"), legend_breaks = c(0,20,40,60,80,100))


# Samples with high CE inclusion
## Subset to avgPSI > 20 for at least 42 targets
## Subset to CEs within those samples with a max inclusion greater than 50
result_df_TDPneg <- result_df_wide_subset_overall%>%
  mutate(n5 = apply(result_df_wide_subset_overall,
                    MARGIN = 1,
                    function(x) length(which(x > 15)))) %>%
  filter(n5 > 42) %>%
  dplyr::select(-n5)

colKeep <- apply(X = result_df_TDPneg, MARGIN = 2, FUN = function(x) max(x) > 50)
result_df_TDPneg <- result_df_TDPneg[,colKeep]

samples <- data.frame(sampleID = rownames(result_df_TDPneg))
sampleMeta <- snapmine::add_SnaptronMeta(samples, compilation = "sra_human") %>%
  mutate(rail_id = as.character(rail_id))
samples <- samples %>% full_join(sampleMeta, by = join_by("sampleID" == "rail_id"))

#manual
rowAnn = data.frame(sampleID = samples$sampleID,
                    Condition = c(rep("HeLa SMG6/7KD UPF1 CLIP",2),
                                  rep("MDA-MB231 shTDP43/shSRSF3",2),
                                  "HeLa shUPF1",
                                  "hMN siTDP-43",
                                  rep("MRC-5 MERS-CoV",5))) %>%
  column_to_rownames("sampleID")


pheatmap(as.matrix(result_df_TDPneg),
         show_rownames = T,
         fontsize = 12,
         color = viridis(100,option="B"), legend_breaks = c(0,20,40,60,80,100),
         border_color = "white",
         annotation_row = rowAnn,
         annotation_colors = list(Condition =
                                    c(`HeLa SMG6/7KD UPF1 CLIP` = "#788AA3",
                                      `MDA-MB231 shTDP43/shSRSF3` = "#C83E4D",
                                      `HeLa shUPF1` = "#FFCF99",
                                      `hMN siTDP-43` = "#C83E4D",
                                      `MRC-5 MERS-CoV` = "#302F4D")))

# Samples with high CE inclusion
## Subset to avgPSI > 5 for at least 25% targets
## Subset to CEs within those samples with a max inclusion greater than 15 for at least 25% of samples
result_df_TDPneg <- result_df_wide_subset_overall%>%
  mutate(n5 = apply(result_df_wide_subset_overall,
                    MARGIN = 1,
                    function(x) length(which(x > 5)))) %>%
  filter(n5 > (0.25*423)) %>%
  dplyr::select(-n5)

colKeep <- (apply(result_df_TDPneg,
                  MARGIN = 2,
                  function(x) length(which(x > 15)))>(0.25*nrow(result_df_TDPneg)))

result_df_TDPneg <- result_df_TDPneg[,colKeep]

## Add CE values from i3N
ces_kept <- str_split_i(colnames(result_df_TDPneg),"_",2)
ces_kept <- data.frame(ID = ces_kept) %>%
  left_join(df %>% dplyr::select(`ID #`,avgPSI_sgTARDBP, `avgPSI_sgTARDBP_shXRN1/shUPF1`),
            by = join_by(ID == `ID #`))

result_df_TDPneg[8,] = ces_kept$avgPSI_sgTARDBP
rownames(result_df_TDPneg)[8] = "sgTDP-43"

result_df_TDPneg[9,] = ces_kept$`avgPSI_sgTARDBP_shXRN1/shUPF1`
rownames(result_df_TDPneg)[9] = "sgTDP-43 + shUPF1/shXRN1"

samples <- data.frame(sampleID = rownames(result_df_TDPneg)[-c(8:9)])
sampleMeta <- snapmine::add_SnaptronMeta(samples, compilation = "sra_human") %>%
  mutate(rail_id = as.character(rail_id))
samples <- samples %>% full_join(sampleMeta, by = join_by("sampleID" == "rail_id"))

#manual
rowAnn = data.frame(sampleID = c(samples$sampleID, "sgTDP-43","sgTDP-43 + shUPF1/shXRN1"),
                    Condition = c(rep("MDA-MB231 shTDP43",2),
                                  rep("hMN siTDP43",2),
                                  rep("MERS-CoV",3),
                                  rep("i3N",2))) %>%
  column_to_rownames("sampleID")

pheatmap(as.matrix(result_df_TDPneg),
         show_rownames = T,
         fontsize = 12,
         color = viridis(100,option="B"), legend_breaks = c(0,20,40,60,80,100),
         border_color = "black",
         annotation_row = rowAnn,
         annotation_colors = list(Condition =
                                    c(`MERS-CoV` = "#788AA3",
                                      `MDA-MB231 shTDP43` = "#FFCF99",
                                      `hMN siTDP43` = "#C83E4D",
                                      `i3N` = "#302F4D")))

