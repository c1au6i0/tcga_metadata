# @@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Get Metadata From TCGA -----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@
# the following code was used to import data.

# ++++++++++++++++++++
# Dataframes and Lists
# ++++++++++++++++++++
# all_ann:
#    dataframe with all metadata taken from GDC
# biospecimen_cl:
#    list of dataframes regariding biospecimens
# clinical:
#    list of tabs regarding clinical data

library(here)
source(here("code", "libraries.R"))


# @@@@@@@@@@@@@@@@@@@@@@
# Clinical
# @@@@@@@@@@@@@@@@@@@@@@

# Download Biospecimen and Clinical data (in `csv`) from the
# [TCGA-PRAD website](https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-PRAD%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22DNA%20Methylation%22%5D%7D%7D%5D%7D&searchTableTab=cases)
# in `data/methylation` folder.
#
untar(here("data", "methylation","clinical.cases_selection.2021-01-25.tar.gz"), exdir = here("data", "methylation","clinical"))

to_imp <- list.files(here("data", "methylation","clinical"), full.names = TRUE)

clinical <- lapply(to_imp, read.delim)
names(clinical) <- sub("\\.tsv", "", list.files(here("data", "methylation","clinical")))

clean_data <- function(dat){
  dat[dat == "Not Reported"| dat == "'--"] <- NA
  # # dat[dat == "'--"] <- NA
  dat %>%
    select(where(~!all(is.na(.))))
}
clinical <- lapply(clinical, clean_data)


# @@@@@@@@@@@@@@@@
# Biospecimen Info
# @@@@@@@@@@@@@@@@

untar(here("data", "methylation", "biospecimen.cases_selection.2021-01-25.tar.gz"), exdir = here("data", "methylation","biospecimen"))

to_imp <- list.files(here("data", "methylation","biospecimen"), full.names = TRUE)
biospecimen <- lapply(to_imp, read.delim)
names(biospecimen) <- sub("\\.tsv", "", list.files(here("data", "methylation","biospecimen")))

biospecimen_cl <- lapply(biospecimen, clean_data)
biospecimen_cl$aliquot$analyte_type <- NULL

unique(biospecimen_cl$aliquot$analyte_type)

clean_data(biospecimen$portion)


# @@@@@@@@@@@@@@@
# All annotations
# @@@@@@@@@@@@@@@


lapply(biospecimen_cl, names)

ids_merge1 <- grep("id", intersect(
  names(biospecimen_cl$aliquot),
  names(biospecimen_cl$analyte)), value = TRUE
  )

analyte <- full_join( biospecimen_cl$analyte,
                      biospecimen_cl$aliquot,
                      by = ids_merge1,
                      suffix = c("_analyte", "_aliquot")
                      )

# there are portions without id


ids_merge2 <- grep("id", intersect(
  names(analyte),
  names(biospecimen_cl$portion)),
  value = TRUE
)

portion <- full_join( analyte,
                      biospecimen_cl$portion,
                      by = ids_merge2[1:6],
                      suffix = c("_analyte", "_portion")
)

length(unique(analyte$portion_id))
length(unique(biospecimen_cl$portion$portion_id))

nrow(analyte)
nrow(portion)

# IDS shared among datasets
# ids <- Reduce(intersect, lapply(biospecimen_cl, names))
#
# lapply(biospecimen_cl, nrow)

# all_bios <- plyr::join_all(biospecimen_cl, by = ids, type = "full")
# ids_2 <- Reduce(intersect, lapply(list(clinical$clinical, all_bios), names))
#
# all_ann <- full_join(clinical$clinical, all_bios, by = ids_2)
#
# View(all_ann)
#
# all_ann <- all_ann %>%
#   mutate(analyte_type_id = recode(analyte_type_id,
#                                   "W" = "Whole Genome Amplification",
#                                   "D" = "DNA",
#                                   "H" = "mirVana RNA",
#                                   "R" = "RNA",
#                                   "T" = "total RNA",
#                                   "x" = "Whole Genome Amplification second reaction"
#   ))




save(biospecimen_cl, clinical,
     file = here("objs", "annotations_methylation.RDA"))
