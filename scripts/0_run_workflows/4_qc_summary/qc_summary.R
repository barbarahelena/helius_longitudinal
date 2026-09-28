## QC summary: sequencing depth (16S, shotgun) and MAG quality
## Revision: R3 #13, R2 minor #7 (issue #16)
## Reports QC as applied in the existing analysis — no re-filtering, no new thresholds.
## Numbers that cannot be derived from local files are flagged, not estimated.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(openxlsx)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

#### Output folder ####
out_dir <- "results/0_data_cleaning/qc_summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants (as applied in the analysis — do not change) ####
RARE_DEPTH         <- 13000  # --rarelevel, scripts/0_run_workflows/1_run_16s_pipeline/run_vsearch.sh
FILTER_MAXEE       <- 1      # filter_maxee, 16S pipeline params
MAG_COMPL_MIN      <- 70     # completeness > 70%, alistipes_bins_annotation/filter_samplesheets_by_quality.py
MAG_CONT_MAX       <- 10     # contamination < 10%, idem
GTDBTK_COMPL_MIN   <- 50     # nf-core/mag GTDB-Tk input filter (bins below are not classified)
GTDBTK_CONT_MAX    <- 10
PRESENCE_THRESHOLD <- 1      # re-mapping depth for MAG presence, 4_alistipes_anno/utils.R
MIN_CLADE_N        <- 5      # MIN_CLADE_SIZE, 4_alistipes_anno/3_draw_tree.R

tp_levels <- c("baseline", "follow-up")
tp_from_id <- function(x) if_else(str_starts(x, "HELIBA"), "baseline", "follow-up")
study_sample <- "^HELI(BA|FU)_[0-9]+$"

# Summary statistics of reads per group
depth_stats <- function(df) {
    df %>%
        filter(!is.na(reads)) %>%
        summarise(n = n(),
                  median = median(reads), q1 = quantile(reads, 0.25), q3 = quantile(reads, 0.75),
                  min = min(reads), max = max(reads), .groups = "drop") %>%
        mutate(across(c(median, q1, q3, min, max), round),
               `Median (IQR)` = sprintf("%s (%s–%s)", format(median, big.mark = ","),
                                        format(q1, big.mark = ","), format(q3, big.mark = ",")),
               Range = sprintf("%s–%s", format(min, big.mark = ","), format(max, big.mark = ",")))
}

med_iqr <- function(x, digits = 1) {
    sprintf("%.*f (%.*f–%.*f)", digits, median(x), digits, quantile(x, 0.25), digits, quantile(x, 0.75))
}

#### Clinical data / analysed sets ####
clin <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
eth <- clin %>% distinct(ID, EthnicityTot)

# 16S analysed set: participants with paired 16S + clinical data (datacleaning.R)
ids_16s_paired <- read.csv2("data/16s/ids_16s_paired.csv")$x

#### 16S: per-sample reads at each QC stage ####
print('16S read counts per stage..')
res16s <- "data/16s/results_helius_paired"

raw16s <- read_tsv(file.path(res16s, "fastqc_analysis/fastqc_summary.tsv"), show_col_types = FALSE) %>%
    mutate(sampleID = str_remove(sample, "_[12]$")) %>%
    group_by(sampleID) %>%
    summarise(raw = max(total_sequences), .groups = "drop")

filter_files <- list.files(file.path(res16s, "vsearch"), pattern = "\\.filter_stats\\.txt$", full.names = TRUE)
filt16s <- map_dfr(filter_files, function(f) {
    line <- str_subset(readLines(f), "sequences kept")
    tibble(sampleID = str_remove(basename(f), "\\.filter_stats\\.txt$"),
           kept = as.numeric(str_extract(line, "^[0-9]+")),
           discarded = as.numeric(str_match(line, "([0-9]+) sequences discarded")[, 2]))
}) %>%
    mutate(merged = kept + discarded)

mapped16s <- read_tsv(file.path(res16s, "vsearch/mapping_rate_summary.tsv"), show_col_types = FALSE) %>%
    dplyr::select(sampleID = sample, mapped = mapped_reads)

ps_decontam <- readRDS(file.path(res16s, "phyloseq/decontam/phyloseq_decontam.RDS"))
decontam16s <- tibble(sampleID = sample_names(ps_decontam), decontam = sample_sums(ps_decontam))

ps_rare <- readRDS("data/16s/phyloseq_rarefied.RDS")
rare16s <- tibble(sampleID = sample_names(ps_rare), rarefied = sample_sums(ps_rare))

reads16s <- raw16s %>%
    left_join(filt16s %>% dplyr::select(sampleID, merged, kept), by = "sampleID") %>%
    left_join(mapped16s, by = "sampleID") %>%
    left_join(decontam16s, by = "sampleID") %>%
    left_join(rare16s, by = "sampleID") %>%
    filter(str_detect(sampleID, study_sample)) %>%
    mutate(timepoint = tp_from_id(sampleID),
           ID = str_c("S", str_remove(sampleID, "^HELI.._")),
           passed_qc = !is.na(rarefied),
           analysed = ID %in% ids_16s_paired & passed_qc)

stages16s <- c(raw = "1. Raw read pairs",
               merged = "2. Merged pairs",
               kept = sprintf("3. Quality-filtered (maxEE %s)", FILTER_MAXEE),
               mapped = "4. Mapped to ASVs",
               decontam = "5. After decontam",
               rarefied = "6. Rarefied")

#### Shotgun: per-sample read pairs at each QC stage ####
print('Shotgun read counts per stage..')
# nf-core/mag preprocessing: fastp -> Bowtie2 host removal -> Bowtie2 PhiX removal
read_sg_batch <- function(b) {
    d <- sprintf("data/shotgun/multiqc_data_%d", b)
    gs <- yaml::read_yaml(file.path(d, "multiqc_general_stats.yaml"))
    raw <- map_dfr(str_subset(names(gs), "_raw_1$"), ~ tibble(
        sampleID = str_remove(.x, "_run[0-9]+_raw_1$"),
        raw = round(gs[[.x]][["fastqc_raw_reads-total_sequences"]] * 1e6)))
    # bowtie2-1 = host removal (input = fastp-passed pairs); bowtie2 = PhiX removal
    host <- yaml::read_yaml(file.path(d, "multiqc_bowtie2_bowtie2-1.yaml"))
    phix <- yaml::read_yaml(file.path(d, "multiqc_bowtie2.yaml"))
    bt <- function(y, field) map_dfr(str_subset(names(y), "^HELI"), ~ tibble(
        sampleID = str_remove(.x, "_run[0-9]+$"), value = y[[.x]][[field]]))
    raw %>%
        left_join(bt(host, "total_reads") %>% rename(trimmed = value), by = "sampleID") %>%
        left_join(bt(phix, "paired_aligned_none") %>% rename(host_removed = value), by = "sampleID") %>%
        mutate(batch = b)
}
readssg <- map_dfr(1:3, read_sg_batch)

# Sanity check of stage order: host-removal input must equal fastp-passed reads / 2
fastp_pass <- map_dfr(1:3, function(b) {
    gs <- yaml::read_yaml(sprintf("data/shotgun/multiqc_data_%d/multiqc_general_stats.yaml", b))
    map_dfr(str_subset(names(gs), "^HELI(BA|FU)_[0-9]+_run[0-9]+$"), ~ tibble(
        sampleID = str_remove(.x, "_run[0-9]+$"),
        fastp_pairs = gs[[.x]][["fastp-filtering_result_passed_filter_reads"]] * 1e6 / 2))
})
stopifnot(all(abs(inner_join(readssg, fastp_pass, by = "sampleID") %>%
                      with(trimmed - fastp_pairs)) < 1))

sg_abundance <- rownames(readRDS("data/shotgun/shotgun_abundance.RDS"))
readssg <- readssg %>%
    filter(str_detect(sampleID, study_sample)) %>%
    mutate(timepoint = tp_from_id(sampleID),
           ID = str_c("S", str_remove(sampleID, "^HELI.._")),
           passed_qc = sampleID %in% sg_abundance)
sg_paired_ids <- readssg %>%
    filter(passed_qc, sampleID %in% clin$sampleID) %>%
    count(ID) %>% filter(n == 2) %>% pull(ID)
readssg <- readssg %>% mutate(analysed = passed_qc & ID %in% sg_paired_ids)

stagessg <- c(raw = "1. Raw read pairs",
              trimmed = "2. After trimming (fastp)",
              host_removed = "3. After host + PhiX removal")

# The MetaPhlAn profiles come from an earlier preprocessing run (2024) of the same raw
# files; compare raw read pairs where that run's FastQC output is available locally.
old_raw <- read_tsv("data/shotgun/multiqc/multiqc_data/multiqc_fastqc.txt", show_col_types = FALSE) %>%
    filter(str_ends(Sample, "_1")) %>%
    transmute(sampleID = str_remove(Sample, "_1$"), raw_old = `Total Sequences`)
raw_check <- inner_join(old_raw, readssg %>% dplyr::select(sampleID, raw), by = "sampleID")
n_raw_identical <- sum(abs(raw_check$raw_old - raw_check$raw) < 1)

#### Depth tables ####
depth_table <- function(reads, stages) {
    long <- reads %>%
        pivot_longer(all_of(names(stages)), names_to = "stage", values_to = "reads") %>%
        mutate(stage = factor(stages[stage], levels = stages),
               timepoint = factor(timepoint, levels = tp_levels))
    bind_rows(
        long %>% mutate(set = "All study samples sequenced"),
        long %>% filter(analysed) %>% mutate(set = "Analysed (paired, after QC)")
    ) %>%
        group_by(set, timepoint, stage) %>%
        depth_stats() %>%
        arrange(set, stage, timepoint) %>%
        dplyr::select(Set = set, Timepoint = timepoint, Stage = stage, n,
                      Median = median, Q1 = q1, Q3 = q3, Min = min, Max = max,
                      `Median (IQR)`, Range)
}
depth16s <- depth_table(reads16s, stages16s)
depthsg  <- depth_table(readssg, stagessg)

#### Descriptive checks: depth by timepoint and by ethnic group ####
# Final pre-normalisation QC stage: 16S after decontam; shotgun after host removal
depth_checks <- function(reads, col, label) {
    d <- reads %>% filter(analysed) %>%
        dplyr::select(ID, timepoint, reads = all_of(col)) %>%
        left_join(eth, by = "ID")
    w <- d %>% dplyr::select(ID, timepoint, reads) %>%
        pivot_wider(names_from = timepoint, values_from = reads) %>%
        drop_na()
    wt <- wilcox.test(w$baseline, w$`follow-up`, paired = TRUE, exact = FALSE)
    paired_row <- tibble(
        Test = "Wilcoxon signed-rank, baseline vs follow-up (paired)",
        Timepoint = "both", Stage = label, n = nrow(w),
        Detail = sprintf("median baseline %s; follow-up %s; median within-person difference (FU - BA) %s",
                         format(round(median(w$baseline)), big.mark = ","),
                         format(round(median(w$`follow-up`)), big.mark = ","),
                         format(round(median(w$`follow-up` - w$baseline)), big.mark = ",")),
        Statistic = unname(wt$statistic), p = wt$p.value)
    eth_rows <- map_dfr(tp_levels, function(tp) {
        dd <- d %>% filter(timepoint == tp, !is.na(EthnicityTot), EthnicityTot != "Other") %>% droplevels()
        kt <- kruskal.test(reads ~ EthnicityTot, data = dd)
        meds <- dd %>% group_by(EthnicityTot) %>%
            summarise(m = format(round(median(reads)), big.mark = ","), .groups = "drop")
        tibble(Test = "Kruskal-Wallis, by ethnic group (excl. Other)",
               Timepoint = tp, Stage = label, n = nrow(dd),
               Detail = str_c(meds$EthnicityTot, ": ", meds$m, collapse = "; "),
               Statistic = unname(kt$statistic), p = kt$p.value)
    })
    bind_rows(paired_row, eth_rows) %>%
        mutate(Statistic = round(Statistic, 1), p = signif(p, 3))
}
checks16s <- depth_checks(reads16s, "decontam", stages16s[["decontam"]])
checkssg  <- depth_checks(readssg, "host_removed", stagessg[["host_removed"]])

#### Exclusions and paired participants ####
exclusions <- bind_rows(
    reads16s %>%
        group_by(Modality = "16S", Timepoint = timepoint) %>%
        summarise(`n sequenced` = n(),
                  `Low-depth threshold` = sprintf("< %s reads after decontam (rarefaction depth)",
                                                  format(RARE_DEPTH, big.mark = ",")),
                  `n excluded: low depth` = sum(!is.na(decontam) & decontam < RARE_DEPTH),
                  `n excluded: other QC` = sum(is.na(decontam)),
                  `Other QC reason` = if_else(`n excluded: other QC` > 0,
                                              "absent after DADA2/VSEARCH processing", "–"),
                  `n retained after QC` = sum(passed_qc),
                  `n in analysed paired set` = sum(analysed), .groups = "drop"),
    readssg %>%
        group_by(Modality = "Shotgun", Timepoint = timepoint) %>%
        summarise(`n sequenced` = n(),
                  `Low-depth threshold` = "None applied",
                  `n excluded: low depth` = 0L,
                  `n excluded: other QC` = sum(!passed_qc),
                  `Other QC reason` = if_else(`n excluded: other QC` > 0,
                                              str_c("MetaPhlAn profile all NA (datacleaning.R): ",
                                                    str_c(sampleID[!passed_qc], collapse = ", ")), "–"),
                  `n retained after QC` = sum(passed_qc),
                  `n in analysed paired set` = sum(analysed), .groups = "drop")
)
# 16S check: low-depth exclusions must be exactly the samples missing from the rarefied object
stopifnot(all(with(reads16s, (!is.na(decontam) & decontam < RARE_DEPTH) == (!is.na(decontam) & !passed_qc))))

paired_after_qc <- function(reads) reads %>% filter(passed_qc) %>% count(ID) %>% filter(n == 2) %>% nrow()
paired <- tibble(
    Modality = c("16S", "Shotgun"),
    `Participants sequenced` = c(n_distinct(reads16s$ID), n_distinct(readssg$ID)),
    `Paired after QC` = c(paired_after_qc(reads16s), paired_after_qc(readssg)),
    `Paired in analysed set (clinical data, no antibiotics)` =
        c(n_distinct(reads16s$ID[reads16s$analysed]), n_distinct(readssg$ID[readssg$analysed]))
)

#### MAGs ####
print('MAG quality..')
tier_std <- function(compl, cont) {
    case_when(compl >= 90 & cont < 5 ~ "High (>=90%, <5%)",
              compl >= 50 & cont < 10 ~ "Medium (>=50%, <10%)",
              TRUE ~ "Low")
}
tier_levels <- c("High (>=90%, <5%)", "Medium (>=50%, <10%)", "Low")

checkm2 <- map_dfr(1:3, ~ read_tsv(sprintf("data/shotgun/summaries/batch%d/checkm2_summary.tsv", .x),
                                   show_col_types = FALSE) %>% mutate(batch = .x)) %>%
    mutate(tier = factor(tier_std(Completeness, Contamination), levels = tier_levels),
           pass_analysis = Completeness > MAG_COMPL_MIN & Contamination < MAG_CONT_MAX)
gtdbtk <- map_dfr(1:3, ~ read_tsv(sprintf("data/shotgun/summaries/batch%d/gtdbtk_summary.tsv", .x),
                                  show_col_types = FALSE, col_types = cols(.default = "c"))) %>%
    transmute(Name = str_remove(user_genome, "\\.fa$"), classification)
checkm2 <- checkm2 %>% left_join(gtdbtk, by = "Name")

# Analysed Alistipes putredinis MAGs (post-filter batch tables used in 4_alistipes_anno)
batch_files <- sprintf("data/shotgun/alistipes_annotation/bins_alistipes_batch%d.csv", 1:3)
mag_batch <- map_dfr(batch_files, ~ read.csv(.x, check.names = FALSE))
mags <- mag_batch %>%
    transmute(bin_name = sub("\\.fa$", "", bin), Completeness, Contamination) %>%
    mutate(tier = factor(tier_std(Completeness, Contamination), levels = tier_levels))
stopifnot(all(mags$Completeness > MAG_COMPL_MIN & mags$Contamination < MAG_CONT_MAX))

trans <- read.delim("data/shotgun/alistipes_annotation/bin_translation_table.tsv") %>%
    mutate(subject_id = as.character(subject_id))
depth_own <- mag_batch %>%
    transmute(bin_name = sub("\\.fa$", "", bin), across(starts_with("Depth "))) %>%
    pivot_longer(starts_with("Depth "), names_to = "sampleID", values_to = "depth") %>%
    mutate(sampleID = sub("^Depth ", "", sampleID), depth = replace_na(depth, 0)) %>%
    inner_join(trans %>% dplyr::select(bin_name, subject_id), by = "bin_name") %>%
    filter(sampleID == str_c("HELIBA_", subject_id) | sampleID == str_c("HELIFU_", subject_id)) %>%
    mutate(timepoint = tp_from_id(sampleID)) %>%
    dplyr::select(bin_name, subject_id, timepoint, depth) %>%
    pivot_wider(names_from = timepoint, values_from = depth, values_fill = 0)

clades <- readRDS("results/3_species_change/4_alistipes_anno/bin_quality_clade.RDS") %>%
    dplyr::select(bin_name, clade)
mags <- mags %>%
    left_join(depth_own, by = "bin_name") %>%
    left_join(clades, by = "bin_name") %>%
    mutate(clade = replace_na(clade, "Not placed in phylogeny"),
           present_ba = baseline >= PRESENCE_THRESHOLD,
           present_fu = `follow-up` >= PRESENCE_THRESHOLD)
stopifnot(!any(is.na(mags$subject_id)))

n_aput_classified <- sum(str_detect(checkm2$classification, "s__Alistipes putredinis$"), na.rm = TRUE)
n_unclassified    <- sum(is.na(checkm2$classification))
mag_counts <- tibble(
    Item = c("Total bins (CheckM2, all taxa, 3 batches)",
             str_c("  ", tier_levels),
             sprintf("  Passing analysis filter (>%s%% completeness, <%s%% contamination)", MAG_COMPL_MIN, MAG_CONT_MAX),
             sprintf("  Classified by GTDB-Tk (input filter >=%s%%, <=%s%%)", GTDBTK_COMPL_MIN, GTDBTK_CONT_MAX),
             "  Not classified by GTDB-Tk (below input filter)",
             "Alistipes putredinis bins (GTDB-Tk)",
             "  Passing analysis filter = analysed MAGs",
             str_c("  Analysed MAGs: ", tier_levels),
             "  Analysed MAGs placed in phylogeny (clade assigned)",
             "  Analysed MAGs not placed in phylogeny",
             "A. putredinis bins among bins not classified by GTDB-Tk"),
    n = c(nrow(checkm2),
          as.integer(table(checkm2$tier)[tier_levels]),
          sum(checkm2$pass_analysis),
          sum(!is.na(checkm2$classification)),
          n_unclassified,
          n_aput_classified,
          nrow(mags),
          as.integer(table(mags$tier)[tier_levels]),
          sum(mags$clade != "Not placed in phylogeny"),
          sum(mags$clade == "Not placed in phylogeny"),
          NA)
) %>%
    mutate(`%` = case_when(str_detect(Item, "^  Analysed MAGs: ") ~ round(100 * n / nrow(mags), 1),
                           str_starts(Item, "  ") & row_number() <= 7 ~ round(100 * n / nrow(checkm2), 1),
                           TRUE ~ NA_real_),
           Note = case_when(
               str_starts(Item, "A. putredinis bins among") ~
                   sprintf("NOT DERIVABLE: %s bins below the GTDB-Tk input filter have no taxonomy", format(n_unclassified, big.mark = ",")),
               Item == "  Passing analysis filter = analysed MAGs" ~
                   "All classified A. putredinis bins passed the analysis filter",
               TRUE ~ ""))

quality_summary <- function(df, set) {
    bind_rows(df %>% mutate(tier = "All"), df %>% mutate(tier = as.character(tier))) %>%
        mutate(tier = factor(tier, levels = c("All", tier_levels))) %>%
        group_by(Tier = tier, .drop = TRUE) %>%
        summarise(n = n(),
                  `Completeness median (IQR)` = med_iqr(Completeness),
                  `Contamination median (IQR)` = med_iqr(Contamination, 2),
                  .groups = "drop") %>%
        mutate(Set = set, .before = 1)
}
mag_quality <- bind_rows(
    quality_summary(mags, "Analysed A. putredinis MAGs"),
    quality_summary(checkm2, "All bins (all taxa)")
)

clade_order <- c(sort(str_subset(unique(mags$clade), "^Clade")),
                 sort(str_subset(unique(mags$clade), "^Unclassified")),
                 "Not placed in phylogeny")
per_clade_rows <- function(df) {
    df %>%
        summarise(`n MAGs` = n(),
                  `n participants` = n_distinct(subject_id),
                  `MAG present at baseline` = n_distinct(subject_id[present_ba]),
                  `MAG present at follow-up` = n_distinct(subject_id[present_fu]),
                  `MAG present at both` = n_distinct(subject_id[present_ba & present_fu]),
                  `Completeness median (IQR)` = med_iqr(Completeness),
                  `Contamination median (IQR)` = med_iqr(Contamination, 2),
                  .groups = "drop")
}
mag_clade <- bind_rows(
    mags %>% group_by(Clade = clade) %>% per_clade_rows(),
    mags %>% per_clade_rows() %>% mutate(Clade = "Total")
) %>%
    mutate(Clade = factor(Clade, levels = c(clade_order, "Total")),
           `In clade analyses` = case_when(Clade == "Total" ~ "",
                                           str_starts(Clade, "Clade") & `n MAGs` >= MIN_CLADE_N ~ "Yes",
                                           TRUE ~ "No"),
           .after = Clade) %>%
    arrange(Clade)
analysed_clades <- mag_clade %>% filter(`In clade analyses` == "Yes")

#### Excel ####
print('Writing Excel..')
st_head <- createStyle(fgFill = "#1F4E79", fontColour = "#FFFFFF", textDecoration = "bold",
                       border = "TopBottomLeftRight", borderColour = "#BFBFBF", wrapText = TRUE, valign = "center")
st_sub  <- createStyle(fgFill = "#D6E4F0", textDecoration = "bold", fontColour = "#1F4E79")
st_band <- createStyle(fgFill = "#F2F7FB")
st_body <- createStyle(border = "TopBottomLeftRight", borderColour = "#D9D9D9")
st_num  <- createStyle(numFmt = "#,##0")

# Write a sheet as a stack of blocks, each with a subheader title row
write_blocks <- function(wb, sheet, blocks) {
    addWorksheet(wb, sheet)
    row <- 1
    ncol_max <- max(map_int(blocks, ~ ncol(.x$df)))
    for (b in blocks) {
        writeData(wb, sheet, b$title, startRow = row, startCol = 1)
        mergeCells(wb, sheet, cols = 1:ncol_max, rows = row)
        addStyle(wb, sheet, st_sub, rows = row, cols = 1:ncol_max, gridExpand = TRUE)
        row <- row + 1
        df <- b$df
        writeData(wb, sheet, df, startRow = row, startCol = 1, headerStyle = st_head)
        body_rows <- seq_len(nrow(df)) + row
        addStyle(wb, sheet, st_body, rows = body_rows, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
        band <- body_rows[seq_along(body_rows) %% 2 == 0]
        if (length(band)) addStyle(wb, sheet, st_band, rows = band, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
        int_cols <- which(map_lgl(df, ~ is.numeric(.x) && all(.x == round(.x), na.rm = TRUE)))
        if (length(int_cols)) addStyle(wb, sheet, st_num, rows = body_rows, cols = int_cols, gridExpand = TRUE, stack = TRUE)
        row <- row + nrow(df) + 2
    }
    setColWidths(wb, sheet, cols = 1:ncol_max, widths = "auto")
}
notes <- function(...) list(title = "Notes", df = tibble(Note = c(...)))

wb <- createWorkbook()
write_blocks(wb, "Depth_16S", list(
    list(title = "16S rRNA gene sequencing: reads per sample at each QC stage", df = depth16s),
    list(title = "Descriptive checks (analysed paired set; not main results)", df = checks16s),
    notes("Stage 1: read pairs per sample (FastQC, R1). Stage 2: merged pairs entering the quality filter. Stage 3: merged reads passing VSEARCH maxEE filter.",
          "Stage 4: reads mapped to ASVs (usearch_global, 97% identity). Stage 5: ASV counts after removal of 83 decontam contaminant ASVs.",
          sprintf("Stage 6: rarefied to %s reads; samples below this depth were excluded.", format(RARE_DEPTH, big.mark = ",")),
          "Analysed set: participants with paired 16S data and clinical data, excluding antibiotic use (data/16s/ids_16s_paired.csv).",
          "Descriptive checks use stage 5 (after decontam, before rarefaction).")
))
write_blocks(wb, "Depth_shotgun", list(
    list(title = "Shotgun metagenomics: read pairs per sample at each QC stage", df = depthsg),
    list(title = "Descriptive checks (analysed paired set; not main results)", df = checkssg),
    notes("Read counts from the nf-core/mag preprocessing (MultiQC, 3 batches): FastQC raw -> fastp -> Bowtie2 host removal -> Bowtie2 PhiX removal.",
          sprintf("The MetaPhlAn profiles were generated in an earlier preprocessing run of the same raw files; raw read pairs were identical for %d of %d samples for which that run's FastQC output is available locally.", n_raw_identical, nrow(raw_check)),
          "Post-trimming and post-host-removal counts of that earlier run are not available locally for all samples; the counts shown are from the nf-core/mag run.",
          "Analysed set: participants with shotgun data at both timepoints and clinical data (antibiotic users excluded).",
          "Descriptive checks use stage 3 (after host + PhiX removal).")
))
write_blocks(wb, "Exclusions", list(
    list(title = "Samples excluded during QC, per modality and timepoint", df = exclusions),
    list(title = "Participants with paired baseline and follow-up data", df = paired),
    notes("'Paired after QC': both timepoints passed sequencing QC. 'Paired in analysed set' additionally requires clinical data and no antibiotic use (datacleaning.R).",
          "No sequencing-depth threshold was applied to shotgun samples.")
))
write_blocks(wb, "MAG_summary", list(
    list(title = "Bin and MAG counts by quality tier", df = mag_counts),
    list(title = "Completeness and contamination (CheckM2), median (IQR)", df = mag_quality),
    notes(sprintf("Analysis filter applied to MAGs: completeness >%s%% and contamination <%s%% (filter_samplesheets_by_quality.py). High/medium tiers are reported for reference only.", MAG_COMPL_MIN, MAG_CONT_MAX),
          "Low tier: bins not meeting the medium criteria.",
          "Not placed in phylogeny: analysed MAGs without a clade assignment in 3_draw_tree.R.")
))
write_blocks(wb, "MAG_per_clade", list(
    list(title = "Alistipes putredinis MAGs per clade", df = mag_clade),
    notes(sprintf("MAG present at a timepoint: re-mapping depth of the participant's own MAG >= %sx in that timepoint's sample (PRESENCE_THRESHOLD, utils.R).", PRESENCE_THRESHOLD),
          sprintf("In clade analyses: named clades with >= %s MAGs (MIN_CLADE_SIZE).", MIN_CLADE_N),
          "Each participant contributed one MAG.")
))
saveWorkbook(wb, file.path(out_dir, "qc_summary.xlsx"), overwrite = TRUE)

#### Supplementary figure ####
print('Plotting..')
tier_cols <- c("High (>=90%, <5%)" = "#1F78B4", "Medium (>=50%, <10%)" = "#FDBF6F", "Low" = "grey60")
pl_mag <- ggplot(mags, aes(x = Completeness, y = Contamination, colour = tier)) +
    geom_vline(xintercept = c(MAG_COMPL_MIN, 90), linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = c(5, MAG_CONT_MAX), linetype = "dashed", colour = "grey60") +
    geom_point(size = 2, alpha = 0.8) +
    scale_colour_manual(values = tier_cols, drop = TRUE) +
    scale_x_continuous(limits = c(MAG_COMPL_MIN - 2, 100)) +
    scale_y_continuous(limits = c(0, MAG_CONT_MAX)) +
    labs(x = "Completeness (%)", y = "Contamination (%)", colour = "",
         title = sprintf("A. putredinis MAGs (n = %d)", nrow(mags))) +
    theme_Publication()

depth_plot <- function(reads, col, title, ylab, hline = NULL) {
    d <- reads %>% filter(analysed) %>%
        mutate(reads = .data[[col]], timepoint = factor(timepoint, levels = tp_levels,
                                                        labels = c("Baseline", "Follow-up")))
    p <- ggplot(d, aes(x = timepoint, y = reads)) +
        geom_line(aes(group = ID), colour = "grey50", alpha = 0.05) +
        geom_violin(aes(fill = timepoint), alpha = 0.6, colour = NA) +
        geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
        scale_y_log10(labels = scales::label_comma()) +
        scale_fill_manual(values = c("Baseline" = "#4E79A7", "Follow-up" = "#F28E2B"), guide = "none") +
        labs(x = "", y = ylab, title = title) +
        theme_Publication()
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "firebrick")
    p
}
pl_16s <- depth_plot(reads16s, "decontam",
                     sprintf("16S (n = %d pairs)", paired$`Paired in analysed set (clinical data, no antibiotics)`[1]),
                     "Reads after decontam (log10)", hline = RARE_DEPTH)
pl_sg  <- depth_plot(readssg, "host_removed",
                     sprintf("Shotgun (n = %d pairs)", paired$`Paired in analysed set (clinical data, no antibiotics)`[2]),
                     "Read pairs after host removal (log10)")

(pl_qc <- ggarrange(pl_mag, ggarrange(pl_16s, pl_sg, ncol = 2, labels = c("B", "C")),
                    nrow = 2, labels = c("A", ""), heights = c(1, 1)))
ggsave(file.path(out_dir, "suppl_fig_qc.pdf"), pl_qc, width = 10, height = 11)

#### Key numbers for the response letter ####
fmt <- function(x) format(round(x), big.mark = ",")
med_line <- function(tab, stage) {
    tab %>% filter(Set == "Analysed (paired, after QC)", Stage == stage) %>%
        mutate(txt = sprintf("  %s: median %s (IQR %s–%s), n = %d", Timepoint, fmt(Median), fmt(Q1), fmt(Q3), n)) %>%
        pull(txt)
}
n_high <- sum(mags$tier == tier_levels[1])
key <- c(
    "QC summary — key numbers (generated by scripts/0_run_workflows/4_qc_summary/qc_summary.R)",
    "",
    "MAGs",
    sprintf("  Total bins (all taxa): %s; A. putredinis MAGs analysed: %d (all >%s%% completeness, <%s%% contamination)",
            fmt(nrow(checkm2)), nrow(mags), MAG_COMPL_MIN, MAG_CONT_MAX),
    sprintf("  High quality (>=90%%, <5%%): %d/%d (%.1f%%) of analysed MAGs; %.1f%% of all bins",
            n_high, nrow(mags), 100 * n_high / nrow(mags), 100 * mean(checkm2$tier == tier_levels[1])),
    sprintf("  Completeness median (IQR): %s%%; contamination: %s%%", med_iqr(mags$Completeness), med_iqr(mags$Contamination, 2)),
    sprintf("  Minimum per analysed clade (n >= %d MAGs): %d participants (%s); present at both timepoints: min %d (%s)",
            MIN_CLADE_N, min(analysed_clades$`n participants`),
            analysed_clades$Clade[which.min(analysed_clades$`n participants`)],
            min(analysed_clades$`MAG present at both`),
            analysed_clades$Clade[which.min(analysed_clades$`MAG present at both`)]),
    sprintf("  FLAG: A. putredinis bins that failed QC cannot be counted — %s bins below the GTDB-Tk input filter were not classified.",
            fmt(n_unclassified)),
    "",
    "Sequencing depth, analysed paired set",
    sprintf("  16S reads after decontam (rarefied to %s):", fmt(RARE_DEPTH)),
    med_line(depth16s, stages16s[["decontam"]]),
    "  Shotgun read pairs after host removal:",
    med_line(depthsg, stagessg[["host_removed"]]),
    "",
    "Paired participants",
    sprintf("  16S: %s paired after QC; %s in analysed set", fmt(paired$`Paired after QC`[1]), fmt(paired$`Paired in analysed set (clinical data, no antibiotics)`[1])),
    sprintf("  Shotgun: %s paired after QC; %s in analysed set", fmt(paired$`Paired after QC`[2]), fmt(paired$`Paired in analysed set (clinical data, no antibiotics)`[2])),
    "",
    "Exclusions",
    sprintf("  16S low depth (<%s reads): %s", fmt(RARE_DEPTH),
            str_c(exclusions$Timepoint[exclusions$Modality == "16S"], " ",
                  exclusions$`n excluded: low depth`[exclusions$Modality == "16S"], collapse = ", ")),
    sprintf("  Shotgun: no depth threshold; %d sample(s) excluded for an all-NA MetaPhlAn profile", sum(exclusions$`n excluded: other QC`[exclusions$Modality == "Shotgun"])),
    "",
    "Depth checks (descriptive)",
    sprintf("  16S %s: %s p = %s", checks16s$Timepoint, checks16s$Test, format(checks16s$p)),
    sprintf("  Shotgun %s: %s p = %s", checkssg$Timepoint, checkssg$Test, format(checkssg$p))
)
writeLines(key, file.path(out_dir, "qc_key_numbers.txt"))
cat(key, sep = "\n")
