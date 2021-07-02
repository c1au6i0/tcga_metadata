# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Script to assamble RNA table for Luigi ----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# @@@@@@@@@@@@@@@@@@@@@@@@
# Load TCGA Metadata -----
# @@@@@@@@@@@@@@@@@@@@@@@@

library(here)
source(here('code', "libraries.R"))

load(file = here("objs", "annotations_methylation.RDA"))


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Load WGS and methylation -----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# tcga tables
load( file = here("objs", "methylation_wgs_tabs.RDA"))


# @@@@@@@@@@@@@@@@@@@
# RNA summaries -----
# @@@@@@@@@@@@@@@@@@@

# this was taken from the gdc portal
rna_set <- read_tsv(here("data", "rna_seq", "case_set_rna_seq.2021-03-03.tsv"))


# table of type of unique samples
biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  filter(case_id %in% !!rna_set$id) %>%
  inner_join(biospecimen_cl$sample[, c("sample_submitter_id", "sample_type")]) %>%
  group_by(case_submitter_id, sample_type) %>%
  summarize(n = n()) %>%
  group_by(sample_type) %>%
  summarize(n = n())

names(biospecimen_cl$analyte)

# samples per donor
biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  filter(case_id %in% !!rna_set$id) %>%
  inner_join(biospecimen_cl$sample[, c("sample_submitter_id", "sample_type")]) %>%
  group_by(case_submitter_id, sample_type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = sample_type, values_from = n) %>%
  View()

biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  filter(case_id %in% !!rna_set$id) %>%
  inner_join(biospecimen_cl$aliquot[, c("aliquot_submitter_id", "analyte_submitter_id")]) %>%
  mutate(center = str_extract(aliquot_submitter_id, ".{2}$")) %>%
  distinct(center)

biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  filter(case_id %in% !!rna_set$id) %>%
  inner_join(biospecimen_cl$sample[, c("sample_submitter_id", "sample_type")]) %>%
  group_by(sample_submitter_id) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# are there replicates?
# yes there are 4 replicates
rna_replicates <- biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  filter(case_id %in% !!rna_set$id) %>%
  mutate(sample_id = gsub(".$", "", sample_submitter_id)) %>%
  group_by(sample_id) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  filter(n > 1)

reg_expr <- paste0(rna_replicates$sample_id, collapse = '|')

# we removed the second
rna_replicates_to_remove <- biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA")  %>%
  filter(case_id %in% !!rna_set$id) %>%
  filter(str_detect(sample_submitter_id, reg_expr)) %>%
  filter(str_detect(sample_submitter_id, "B$")) %>%
  pull(sample_submitter_id)


# @@@@@@@@@@
# Tab -----
# @@@@@@@@@

View(tab_dan_wgs_rna_PPCG)


# we need to use the same PPCG ids used for methylation
tab_join <- tab_dan_wgs_rna_PPCG %>%
  mutate(case_submitter_id = str_sub(Meth_Local_Sample_ID, 1, 12)) %>%
  select(PPCG_Donor_ID, case_submitter_id, WGS_Local_Sample_ID, wgs_dan_PPCG_Sample_ID)


ppcg_donor_ids <-
   tab_join %>%
   select(PPCG_Donor_ID, case_submitter_id) %>%
   distinct(case_submitter_id, .keep_all = TRUE)




ppcg_wgs_ids <-
  tab_join %>%
  select(WGS_Local_Sample_ID, wgs_dan_PPCG_Sample_ID) %>%
  distinct(WGS_Local_Sample_ID, .keep_all = TRUE) %>%
  filter(!is.na(WGS_Local_Sample_ID)) %>%
  mutate(wgs_donor = str_sub(WGS_Local_Sample_ID, 1, 12))


rna_tab <- biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  filter(case_id %in% !!rna_set$id) %>%

  # we take only the first replicate
  filter(!sample_submitter_id %in% rna_replicates_to_remove) %>%

  inner_join(biospecimen_cl$sample[, c("sample_submitter_id", "sample_type")]) %>%
  select(case_submitter_id, sample_submitter_id, analyte_submitter_id, sample_type) %>%
  # first we check the donors
  left_join(ppcg_donor_ids, by = "case_submitter_id") %>%

  # then we join by sample with the ppcg id of wgs
  left_join(ppcg_wgs_ids[, 1:2], by = c( "sample_submitter_id" = "WGS_Local_Sample_ID")) %>%
  rename(Matching_WGS_Sample = wgs_dan_PPCG_Sample_ID) %>%
  mutate(WGS_Local_Sample_ID = if_else(!is.na(Matching_WGS_Sample), sample_submitter_id, NA_character_)) %>%
  mutate(RNA_Local_Sample_ID = sample_submitter_id) %>%
  mutate(`RNA_Sample_Type` = recode(sample_type,
    "Solid Tissue Normal" = "BEN",
    "Primary Tumor" = "TUM",
    "Metastatic" = "MET"
  )) %>%
  mutate(RNA_Distance_To_WGS_Sample = if_else(!is.na(Matching_WGS_Sample), "MATCH",  NA_character_)) %>%

  # now we check if the sample is among the methilation samples
  mutate(Meth_RNAseq_annot = if_else(
    str_sub(sample_submitter_id, 1, 15) %in%
    str_sub(!!tab_dan_wgs_rna_PPCG$Meth_Local_Sample_ID, 1, 15),
    "MATCH",
    NA_character_
    )) %>%
  mutate(center = "USA/WCM",
         PPCG_Meth_Sample_ID = NA,
         Platform = "RNA_seq",
         MetastaticSite = NA) %>%
  select(center,
         PPCG_Donor_ID,
         Matching_WGS_Sample,
         WGS_Local_Sample_ID,
         PPCG_Meth_Sample_ID,
         Platform,
         RNA_Local_Sample_ID,
         RNA_Sample_Type,
         MetastaticSite,
         RNA_Distance_To_WGS_Sample,
         Meth_RNAseq_annot
         )


rna_tab %>%
  write.csv(here("reports", "rna_table.csv"), row.names = FALSE)


rna_tab %>%
  group_by(PPCG_Donor_ID, RNA_Sample_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = RNA_Sample_Type, values_from = n) %>%
  View()
# so 52 patients with tum and ben
# one with tum and met


col_names_rna <- c("Center",
                   "PPCG_Donor_ID",
                   "Matching_WGS_Sample[ID/ADD/NA]",
                   "WGS_Local_Sample_ID[ID/ADD/NA]",
                   "PPCG_Meth_Sample_ID",
                   "Platform[RNAseq; Affymetrix; Beadchip; ...]",
                   "RNA_Local_Sample_ID",
                   "RNA_Sample_Type[NOR/BEN/TUM/MET/NA]",
                   "MetastaticSite[LN/BONE/LIVER/BRAIN/LUNG/BLADDER/ADRENAL/other]",
                   "RNA_Distance_To_WGS_Sample[MATCHED/NEAR-MM/DISTANT/UNMATCH/NA]",
                   "Meth_RNAseq_annot[MATCHED/NEAR-MM/DISTANT/NO_RNA/NA]",
                   "fastq_filename",
                   "bam_filename",
                   "Date_of_scan_or_Sequencing[Y/M/D]",
                   "Comments")

names(rna_tab)

unique(rna_tab$PPCG_Donor_ID)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Read Methylation IDs of Daniel ----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


ppcg_new_donor_ids <-
  readxl::read_xlsx(path = here("data", "ppcg_meth_ids.xlsx")) %>%
    dplyr::rename(ppcg_donor_id = PPCG_Donor_Id) %>%
    mutate(tcga_donor_id = str_sub(Meth_Local_Sample_ID, 1, 12)) %>%
    select(ppcg_donor_id, tcga_donor_id) %>%
    distinct(ppcg_donor_id, .keep_all = TRUE)




write.csv(ppcg_new_donor_ids, file = here("data", "ppcg_tcga_matching_donors_ids.csv"))


ppcg_new_donor_ids %>%
  group_by(tcga_donor_id) %>%
  summarize(n = n()) %>%
  filter(n > 1)




rna_tab_with_ids <-
  rna_tab %>%
  mutate(rna_donor_id = str_sub(RNA_Local_Sample_ID, 1, 12)) %>%
  left_join(ppcg_new_donor_ids, by = c("rna_donor_id" = "tcga_donor_id")) %>%
  mutate(PPCG_Donor_ID = if_else(!is.na(ppcg_donor_id), ppcg_donor_id, PPCG_Donor_ID)) %>%
  select(-ppcg_donor_id)


rna_tab %>%
  group_by(PPCG_Donor_ID, RNA_Sample_Type) %>%
  summarise(n = n()) %>%
  filter(n > 1)


write.csv(rna_tab_with_ids, here("reports", "rna_table.csv"), row.names = FALSE)


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Centers that analyze the samples ----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ls_c <- biospecimen_cl$analyte  %>%
  filter(analyte_type == "RNA") %>%
  pull(analyte_id)

biospecimen_cl$aliquot %>%
  filter(analyte_id %in% ls_c) %>%
  distinct(source_center)

# Source Center 22 & 23

biospecimen_cl$portion %>%
  select(is_ffpe) %>%
  table()



# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes






