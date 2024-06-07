# EdgeR and LIMMA
library(edgeR)
matbase <- matbase[,colSums(matbase) != 0]
clindf <- readRDS("data/clinicaldata_wide.RDS")
clindf <- clindf %>% filter(ID %in% colnames(matbase)) %>% droplevels(.)
dge <- DGEList(matbase)
str(dge)
dge$genes <- rownames(dge)
all(colnames(dge) == clindf$ID) # TRUE
dge$samples$ethnicity <- clindf$Ethnicity

lcpm <- cpm(dge, log=TRUE)
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)

lcpm <- cpm(dge, log=TRUE)
par(mfrow=c(1,2))
col.group <- clindf$Ethnicity
levels(col.group) <-  c("firebrick", "royalblue4")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=clindf$Ethnicity, col=col.group)
title(main="Ethnicity")

ethnicity <- clindf$Ethnicity
design <- model.matrix(~0+ethnicity)
colnames(design) <- gsub("Ethnicity", "", colnames(design))
design

contr.matrix <- makeContrasts(
    DutchvsSur = Dutch-Surinamese, 
    levels = c("Dutch", "Surinamese"))
contr.matrix

v <- voom(dge, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))
dutchsur <- topTreat(efit, coef=1, n=Inf)
head(dutchsur)
library(Glimma)
glMDPlot(efit, coef=1, status=dt, main=colnames(efit)[1],
         side.main="ENTREZID", counts=lcpm, groups=ethnicity, launch=TRUE)

## Make baseline + follow-up set for Ethnicity analyses
rownames(df)
which(duplicated(colnames(df)))
sums <- colSums(df[1:ncol(df)])
length(sums)
abundant <- names(sums[which(sums > 300)])
prevalent <- names(df[,colSums(df > 0) > 100])
abprev <- abundant[abundant %in% prevalent]
length(abprev)
dftotal <- df %>% dplyr::select(all_of(abprev)) %>% 
    mutate_all(as.integer)
dim(dftotal) #  genes left
mattot <- t(as.matrix(dftotal))
rownames(mattot)
mattot <- mattot[,colSums(mattot) != 0]
clindf <- readRDS("data/clinicaldata_long.RDS")
clindf <- clindf %>% filter(sampleID %in% colnames(mattot)) %>% droplevels(.)
clindf <- clindf[order(clindf$sampleID, colnames(mattot)),]
all(clindf$sampleID == colnames(mattot)) # TRUE

timepoint <- fct_recode(clindf$timepoint, "followup" = "follow-up")
ethnicity <- clindf$Ethnicity
ethtime <- as.factor(str_c(ethnicity, ".", timepoint))
dge <- DGEList(counts = mattot, group = timepoint)
str(dge)
dge$genes <- rownames(dge)
all(colnames(dge) == clindf$sampleID) # TRUE
dge$samples$ethnicity <- ethnicity
dge$samples$timepoint <- timepoint
dge$samples$eth.time <- ethtime

design <- model.matrix(~0+ethtime)
fit <- voomLmFit(dge, design, plot=TRUE)
fit <- eBayes(fit)
cont.dif <- makeContrasts(BaseFU = (ethtimeSurinamese.baseline-ethtimeSurinamese.followup)-(ethtimeDutch.baseline-ethtimeDutch.followup), levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
tab <- topTable(fit2, adjust = "BH", number = 50) %>% filter(P.Value < 0.05)
diffids <- tab$ID
diffmat <- mattot[rownames(mattot) %in% tab$ID,]
diffdf <- as.data.frame(t(diffmat)) %>% 
    mutate(sampleID = rownames(.)) %>% 
    left_join(., clindf, by = "sampleID")

head(diffdf)[1:5,1:5]

diffdf <- diffdf %>% mutate(across(1:8, ~log10(.x+1)))
diffmeans <- diffdf %>% pivot_longer(., 1:8, names_to = "ARG", values_to = "count")  %>% 
    mutate(ARG = as.factor(ARG)) %>% 
    group_by(ARG, timepoint, Ethnicity) %>% 
    summarise(meancount = mean(count), 
              sdcount = sd(count),
              ncount = length(count))
gene <- levels(diffmeans$ARG)

library(lme4)
library(afex)
linearmixed_gene <- function(data, var){
    statres <- c()
    data1 <- data %>% mutate(var = .data[[var]])
    model <- lmer(var ~ Ethnicity*timepoint + (1|ID), data = data1)
    res <- summary(model)
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres <- cbind(gene = var, group1 = "baseline", group2 = "follow-up", pval)
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    return(statres)
}

res <- c()
for (a in gene) {
    print(a)
    res <- rbind(res, linearmixed_gene(diffdf, var = a))
}

genesig <- res %>% filter(p_signif != "") %>% dplyr::select(gene) %>% mutate_all(as.factor) %>% droplevels(.)

for(a in 1:nlevels(genesig$gene)){
    var <- paste0(genesig$gene[[a]])
    plota <- ggplot(data = diffdf, aes(x = timepoint, y = .data[[var]], fill = Ethnicity, 
                                       group = interaction(timepoint,Ethnicity))) +
        geom_violin(alpha = 0.2) +
        geom_boxplot(fill = "white",width = 0.1) +
        # stat_pvalue_manual(res, y.position = 150, label = "p_signif", 
        #                    remove.bracket = TRUE, bracket.size = 0) +
        scale_fill_jco() + 
        # scale_y_continuous(limits = c(90,170), breaks = seq(from = 90, to = 170, by = 10)) +
        theme_Publication() +
        labs(x = "Timepoint", y = "Log10(counts+1)", title = var, 
             color = "Ethnicity")
    print(plota)
}

## Make baseline + follow-up set for Diabetes
clindf2 <- clindf %>% filter(!is.na(DM))
dm <- fct_recode(clindf2$DM, "Diabetes" = "Yes", "NoDiabetes" = "No")
timepoint <- fct_recode(clindf2$timepoint, "followup" = "follow-up")
dmtime <- as.factor(str_c(dm, ".", timepoint))
mattot2 <- mattot[,colnames(mattot) %in% clindf2$sampleID]
dge <- DGEList(counts = mattot2)
all(colnames(dge) == clindf2$sampleID) # TRUE
dge$samples$diabetes <- dm
dge$samples$timepoint <- timepoint
dge$samples$dm.time <- dmtime

design <- model.matrix(~0+dmtime)
fit <- voomLmFit(dge, design, plot=TRUE)
fit <- eBayes(fit)
cont.dif <- makeContrasts(BaseFU = (dmtimeDiabetes.followup-dmtimeDiabetes.baseline)-(dmtimeNoDiabetes.baseline-dmtimeNoDiabetes.followup), levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTable(fit2, adjust = "BH", number = 50)


pl <- plotMDS(mattot)
x <- pl$x
y <- pl$y
mdsdf <- data.frame(dim1 = pl$x, dim2 = pl$y, sampleID = rownames(pl$distance.matrix.squared))
mdsdf <- left_join(mdsdf, clindf, by = "sampleID") %>% 
    filter(!is.na(DM)) %>% 
    mutate(DM = fct_recode(DM, "Diabetes" = "Yes", "NoDiabetes" = "No"))
ev <- pl$var.explained[1:2]*100

ggplot(mdsdf, aes(x = dim1, y = dim2, color = interaction(Ethnicity, timepoint), shape = timepoint)) +
    geom_point(alpha = 0.5) +
    xlab(str_c("LogFC dim 1 ", round(ev[1], 1), "%")) +
    ylab(str_c("LogFC dim 2 ", round(ev[2], 1), "%")) +
    ggsci::scale_color_jco() +
    ggsci::scale_fill_jco() +
    stat_ellipse(geom = "polygon", aes(color = interaction(Ethnicity, timepoint), fill = interaction(Ethnicity, timepoint)), 
                 alpha = 0.1, type = "norm") +
    theme_Publication() +
    labs(shape = "", color = "", fill = "")

ggplot(mdsdf, aes(x = dim1, y = dim2, color = interaction(DM, timepoint), shape = timepoint)) +
    geom_point(alpha = 0.5) +
    xlab(str_c("LogFC dim 1 ", round(ev[1], 1), "%")) +
    ylab(str_c("LogFC dim 2 ", round(ev[2], 1), "%")) +
    ggsci::scale_color_jco() +
    ggsci::scale_fill_jco() +
    stat_ellipse(geom = "polygon", aes(color = interaction(DM, timepoint), 
                                       fill = interaction(DM, timepoint)), 
                 alpha = 0.2, type = "norm") + 
    theme_Publication() +
    labs(shape = "", color = "", fill = "")

ggplot(mdsdf, aes(x = dim1, y = dim2, color = timepoint)) +
    geom_point(alpha = 0.5) +
    xlab(str_c("LogFC dim 1 ", round(ev[1], 1), "%")) +
    ylab(str_c("LogFC dim 2 ", round(ev[2], 1), "%")) +
    ggsci::scale_color_jco() +
    ggsci::scale_fill_jco() +
    stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), 
                 alpha = 0.1, type = "norm") + 
    theme_Publication()
