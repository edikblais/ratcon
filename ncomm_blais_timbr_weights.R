
options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)
source("ncomm_helper.R")

library(limma)
library(Biobase)
library(reshape2)
library(dplyr)
library(ggplot2)
library(xlsx)

path.matlab.directory = ""
path.rcode.directory = ""
path.manuscript.directory = ""
path.expression.directory = ""

gpr.color = c(`non-enzymatic` = "#7F7F7F", shared = "#674EA7", 
              `rat-specific` = "#C0504D", `human-specific` = "#4F81BD",
              rno = "#C0504D", hsa = "#4F81BD", 
              similar = "#674EA7", opposite = "#F79646")

rxn.info.load = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table3.xlsx"),
                           sheetIndex = 1, startRow = 2) %>% as.tbl

rxn.info = rxn.info.load %>% 
  mutate(enabled = ifelse(!is.na(enabled), enabled == "true", F)) %>% 
  mutate(enabled_hsa = ifelse(!is.na(enabled_hsa), enabled_hsa == "true", F)) %>% 
  mutate(enabled_rno = ifelse(!is.na(enabled_rno), enabled_rno == "true", F))%>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))

rxn.gene = rxn.info %>% filter(enabled) %>% 
  select(rxn_id, hsa = gpr_hsa, rno = gpr_rno) %>%
  melt(c("rxn_id")) %>% ef_df %>% 
  # Excel annoyingly transforms GPR rules into date/time values
  mutate(value = gsub("193.5625","4645:30",value,fixed = T)) %>%
  mutate(value = gsub("[\\(\\)\\;\\:]+",";",value)) %>%
  ef_split_value(";") %>% ef_df %>% 
  rename(organism_id = variable, gene_id = value) %>% 
  filter(nchar(gene_id) > 0) %>% distinct 

# all gene_id values should be integers
rxn.gene %>% count(grepl("^[0-9]+$",gene_id))


differential.expression.info = data_frame(
  efit_root = c(#"hsa_hep_t2_one_caffeine_d1_gene_efit_robust",
    "hsa_hep_t2_one_caffeine_d2_gene_efit_robust",
    "hsa_hep_t2_one_caffeine_d3_gene_efit_robust",
    "rno_hep_t2_one_caffeine_d1_gene_efit_robust",
    "rno_hep_t2_one_caffeine_d2_gene_efit_robust",
    "rno_hep_t2_one_caffeine_d3_gene_efit_robust")) %>% 
  mutate(efit_file = paste0(path.expression.directory,efit_root,".txt.gz")) %>%
  mutate(organism_id = gsub("_.*","",efit_root),
         dose_id = gsub(".*_","",gsub("_gene_.*","",efit_root)),
         drug_id = gsub("_.*","",gsub(".*_one_","",efit_root)),
         time_id = gsub("_.*","",gsub(".*_hep_","",efit_root)))

differential.expression.load = differential.expression.info %>% 
  with(setNames(efit_file, efit_root)) %>% 
  lapply(read.table,sep = "\t",quote = "", header = T,comment.char = "",
         stringsAsFactors = F, check.names = F) %>% 
  bind_rows(.id = "efit_root") %>% as.tbl %>% 
  mutate(gene_id = as.character(gene_id), 
         feature_id = as.character(feature_id)) %>%
  left_join(differential.expression.info)

# 2175 human genes and 1927 rat genes 
# map to these two expression datasets
differential.expression.load %>% 
  select(organism_id, gene_id) %>% distinct %>%
  count(organism_id, gene_id %in% c(rxn.gene[["gene_id"]]))

fdr.cutoff = 0.1

metabolic.differential.expression = differential.expression.load %>%
  semi_join(rxn.gene %>% select(gene_id, organism_id) %>% distinct) %>% 
  mutate(significant = fdr < fdr.cutoff,
         direction = ifelse(significant, sign(logfc),0)) %>%
  group_by(efit_id, efit_root, organism_id, drug_id, time_id, dose_id) %>%
  mutate(metabolic_gene_count = n(),
         n_up = sum(direction > 0), 
         n_dn = sum(direction < 0), 
         n_significant = sum(direction != 0)) %>% ungroup %>%
  mutate(pct_significant = 100 * n_significant / metabolic_gene_count) %>%
  mutate(efit_ok = pct_significant > 1)

# Doses are not necessarily comparable between species, 
# so we manually selected treatments that induced similar 
# numbers of differentially expressed genes to make 
# TIMBR predictions more comparable. For caffeine, we selected 
# dose 2 (d2) for rats and dose 3 (d3) for humans.
metabolic.differential.expression %>% 
  select(efit_id, n_up, n_dn, n_significant, pct_significant, efit_ok) %>% 
  distinct %>% data.frame


rxn.gene %>% select(rxn_id, rxn_enzymatic) %>% distinct

rxn.pubmed = rxn.info %>% filter(enabled) %>% select(rxn_id, pubmed_id) %>% 
  melt(c("rxn_id")) %>% ef_df %>% 
  mutate(value = gsub("\\-","",value)) %>% 
  mutate(value = gsub("PMID[\\:\\-]*",";PMID:",value)) %>% 
  ef_split_value(";") %>% ef_df %>% 
  filter(grepl("PMID|DOI|UNIPROT", value)) %>%
  filter(grepl("[0-9]+",value)) %>% distinct

timbr.weights.default = rxn.info %>% filter(enabled) %>% 
  select(rxn_id,rxn_class) %>% distinct %>% 
  left_join(rxn.pubmed %>% count(rxn_id) %>% ungroup %>% rename(pubmed_count = n)) %>% 
  mutate(pubmed_count = ifelse(!is.na(pubmed_count), pubmed_count, 0)) %>% 
  left_join(bind_rows(list(
    rxn.gene %>%
      mutate(variable = paste0(organism_id, "_count_all")) %>% 
      group_by(rxn_id, variable) %>% 
      summarize(value = length(unique(setdiff(gene_id,"0")))) %>% ungroup,
    rxn.gene %>%
      semi_join(metabolic.differential.expression %>% select(organism_id, gene_id) %>% distinct) %>% 
      mutate(variable = paste0(organism_id, "_count")) %>% 
      group_by(rxn_id, variable) %>% 
      summarize(value = length(unique(setdiff(gene_id,"0")))) %>% ungroup)) %>%
      dcast(rxn_id ~ variable, value.var = "value", fill = 0) %>% ef_df) %>% 
  mutate(rno_count = ifelse(!is.na(rno_count), rno_count, 0),
         rno_count_all = ifelse(!is.na(rno_count_all), rno_count_all, 0),
         hsa_count = ifelse(!is.na(hsa_count), hsa_count, 0),
         hsa_count_all = ifelse(!is.na(hsa_count_all), hsa_count_all, 0)) %>%
  mutate(rxn_enzymatic = hsa_count > 0 | rno_count > 0,
         rxn_referenced = pubmed_count > 0 ) %>%
  mutate(timbr_weight_enzymatic = ifelse(rxn_enzymatic,1,2),
         timbr_weight_referenced = ifelse(rxn_referenced,1,2),
         timbr_weight_class = ifelse(rxn_class == "boundary",2, ifelse(rxn_class == "transport",2,1)),
         timbr_weight_default = timbr_weight_enzymatic * timbr_weight_referenced * timbr_weight_class)

timbr.rxn.setup = timbr.weights.default %>% select(rxn_id, timbr_weight_default) %>%
  left_join(rxn.info %>% select(rxn_id, rno = gpr_rno, hsa = gpr_hsa)) %>%
  melt(c("rxn_id","timbr_weight_default")) %>% ef_df %>%
  mutate(value = gsub("193.5625","4645:30",value, fixed = T))

timbr.expression.setup = metabolic.differential.expression %>% 
  mutate(limma_id = paste0(organism_id, "_", drug_id, "_", time_id, "_", dose_id),
         limma_ok = efit_ok) %>% filter(limma_ok) %>% 
  ef_df_slice("limma_id")

timbr.weights.list = timbr.expression.setup %>% lapply(ef_timbr_weights,timbr.rxn.setup)
timbr.weights = timbr.weights.list %>% bind_rows

timbr.weights %>% with(qplot(timbr_weight_ctl, timbr_weight_trt))

rxn.irreversible = bind_rows(list(
  timbr.weights %>% select(rxn_id) %>% distinct %>% mutate(irxn_id = paste0(rxn_id,"_f")),
  timbr.weights %>% select(rxn_id) %>% distinct %>% mutate(irxn_id = paste0(rxn_id,"_r"))))

timbr.weights.irreversible = bind_rows(list(
  timbr.weights %>% mutate(rxn_irreversible = paste0(rxn_id,"_f")),
  timbr.weights %>% mutate(rxn_irreversible = paste0(rxn_id,"_r"))))

rno.timbr.weights = bind_rows(list(
  timbr.weights.irreversible  %>%
    filter(organism_id == "rno") %>% 
    mutate(timbr_id = paste0(limma_id, "_ctl")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_ctl),
  timbr.weights.irreversible  %>%
    filter(organism_id == "rno") %>% 
    mutate(timbr_id = paste0(limma_id, "_trt")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_trt))) %>%
  dcast(organism_id + rxn_id + rxn_irreversible ~ timbr_id, value.var = "rxn_weight") %>% ef_df

hsa.timbr.weights = bind_rows(list(
  timbr.weights.irreversible  %>%
    filter(organism_id == "hsa") %>% 
    mutate(timbr_id = paste0(limma_id, "_ctl")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_ctl),
  timbr.weights.irreversible  %>%
    filter(organism_id == "hsa") %>% 
    mutate(timbr_id = paste0(limma_id, "_trt")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_trt))) %>%
  dcast(organism_id + rxn_id + rxn_irreversible ~ timbr_id, value.var = "rxn_weight") %>% ef_df

write.table(rno.timbr.weights, file = paste0(path.matlab.directory,"ncomm_blais_timbr_weights_rno.txt"), 
            sep = "\t", quote = F, row.names = F)
write.table(hsa.timbr.weights, file = paste0(path.matlab.directory,"ncomm_blais_timbr_weights_hsa.txt"), 
            sep = "\t", quote = F, row.names = F)
