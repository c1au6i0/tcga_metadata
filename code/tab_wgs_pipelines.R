# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TAB WGS PIPELINES
# create the tab neeeded for the Sanger Pipeline
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(here)
source(here('code', "libraries.R"))


# This are the columns that we need for the Sanger Pipeline

# column 1: Donor_ID
# column 2: Tissue_ID
# column 3: is_normal
# column 4: is_normal_for_donor
# column 5: Sample_ID
# column 6:


library(readxl)
ppcg_wgs_main <- read_excel(here("reports", "ppcg_wgs_main_final.xlsx"))

ppcg_wgs_us <-  ppcg_wgs_main %>%
    filter(Country == "USA")


names(ppcg_wgs_us)

# With selected col

# Local_Sample_Id is used to join with dataframe with names files.

ppcg_wgs_us_sl <- ppcg_wgs_us %>%
  select(PPCG_Donor_ID, Sample_Type, PPCG_Sample_ID, Is_Normal, Is_Normal_for_Donor, Local_Sample_ID)

  # file names as reported in bionimbus, that is the portal where we downloaded the data.
  # so we
  library(jsonlite)
bionimbus_wgs <- fromJSON(here("data", "wgs", "metadata.json"))

sanger_tab <- bionimbus_wgs %>%
    select(file_name, sample_barcode) %>%
    rename(Local_Sample_ID = sample_barcode,
           relative_file_path = file_name
           ) %>%

   inner_join(ppcg_wgs_us_sl) %>%
   relocate("relative_file_path", .after = Is_Normal_for_Donor) %>%
   mutate(relative_file_path = paste0("/omics/data/ppcg_wgs/", relative_file_path)) %>%
   mutate(relative_file_path = gsub(".bam", "_unmapped.bam", relative_file_path)) %>%
   select(-Local_Sample_ID)


sanger_tab <- sanger_tab %>%
  rename(Donor_ID = PPCG_Donor_ID,
         Tissue_ID = Sample_Type,
         is_normal = Is_Normal,
         is_nornal_for_donor = Is_Normal_for_Donor,
         Sample_ID = PPCG_Sample_ID
         ) %>%
  relocate(Sample_ID, .before = relative_file_path)
write.table(sanger_tab,
            file = here("reports", "sanger_wgs.tsv"),
            sep='\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)
