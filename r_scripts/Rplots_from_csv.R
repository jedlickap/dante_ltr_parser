#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Rplots_from_csv
csv_path <- args[1]
longest_seq_length <- as.numeric(args[2])

directory_path <- dirname(csv_path)
# csv_path <- "Saaz_LTR_retrotransposons_annotation_TE_characteristics.csv"
t <- read.csv(csv_path, header = T)

# load used libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(webshot2)

# Output directory
out_folder <- paste0(directory_path,"/tabs_and_plots/")

# Check if the directory "tabs_and_plots/" exists
if (!dir.exists(out_folder)) {
  # If the directory does not exist, create it
  dir.create(out_folder)
  
  cat("Directory 'tabs_and_plots/' created.\n")
} else {
  cat("Directory 'tabs_and_plots/' already exists.\n")
}

# Assuming 't' is your data frame
summary_table <- t %>%
  group_by(te_sfam, te_fam) %>%  # Group by te_sfam and te_fam
  summarize(
    Count = n(),                                 # Count the number of entries in each group
    TE_length = mean(te_length, na.rm = TRUE),
    LTR_length = mean(ltr_avg_len, na.rm = TRUE),    # Average ltr_avg_len for each group
    LTR_identity = mean(ltr_identity, na.rm = TRUE),  # Average ltr_identity for each group
    LTR_K80 = mean(K80, na.rm = TRUE),               # Average K80 value for each group
    TE_age_MYA = mean(MYA, na.rm = TRUE),               # Average MYA for each group
    TSD_PBS = sum(tsd == "True" & pbs == "True"),
    TSD_only = sum(tsd == "True" & pbs == "False"),
    PBS_only = sum(tsd == "False" & pbs == "True"),
    NO_TSD_PBS = sum(tsd == "False" & pbs == "False"),
    Autonomous = sum(autonomy_stat == "autonomous"),     # Count of autonomous TEs
    Non_autonomous = sum(autonomy_stat == "non_autonomous") # Count of non-autonomous TEs
  )

sumCopia <-subset(summary_table,summary_table$te_sfam=="Ty1/copia")
copiaOrder <- sumCopia[order(sumCopia$Count,decreasing = TRUE),]$te_fam
sumGypsy <-subset(summary_table,summary_table$te_sfam=="Ty3/gypsy")
gypsyOrder <- sumGypsy[order(sumGypsy$Count,decreasing = TRUE),]$te_fam

desired_order <- c(copiaOrder,gypsyOrder)

# Convert te_fam to a factor with the specified order
summary_table <- summary_table %>%
  mutate(te_fam = factor(te_fam, levels = desired_order)) %>%
  arrange(te_fam)


# Summarize the required columns: sums and means
summary_row <- summary_table %>%
  summarise(
    te_sfam = unique(te_sfam), 
    te_fam = "Summary",  # to make sure this is compatible with character columns
    Count = sum(Count, na.rm = TRUE),
    TE_length = mean(TE_length, na.rm = TRUE),
    LTR_length = mean(LTR_length, na.rm = TRUE),
    LTR_identity = mean(LTR_identity, na.rm = TRUE),
    LTR_K80 = mean(LTR_K80, na.rm = TRUE),
    TE_age_MYA = mean(TE_age_MYA, na.rm = TRUE),
    TSD_PBS = sum(TSD_PBS, na.rm = TRUE),
    TSD_only = sum(TSD_only, na.rm = TRUE),
    PBS_only = sum(PBS_only, na.rm = TRUE),
    NO_TSD_PBS = sum(NO_TSD_PBS, na.rm = TRUE),
    Autonomous = sum(Autonomous, na.rm = TRUE),
    Non_autonomous = sum(Non_autonomous, na.rm = TRUE)
  )

# Bind the summary row to the original dataframe
final_table <- bind_rows(summary_table, summary_row)
# Round floats in table
final_table <- final_table %>% mutate(across(c('TE_length', 'LTR_length','LTR_identity','LTR_K80','TE_age_MYA'), round, 2))
# Save summary table into csv
write.csv(final_table,paste0(out_folder,"summary_table.csv"), row.names = FALSE)

# table format to html and png
colNr <- length(colnames(final_table))
table_md <- final_table %>%
  kbl() %>%
  kable_styling(full_width = FALSE, fixed_thead = TRUE) %>%
  column_spec(c(1:colNr), width = "3cm") %>%   # Adjust column width, e.g., for the first column
  row_spec(nrow(final_table)-1, bold = TRUE, background = "lightgray") %>%
  row_spec(nrow(final_table), bold = TRUE, background = "lightgray")

# Save the table as an HTML file (intermediate step)
save_kable(table_md, paste0(out_folder,"kable_summary_table.html"))
# Convert the HTML table to a PNG image with a larger viewport for wider columns
webshot(paste0(out_folder,"kable_summary_table.html"), file = paste0(out_folder,"kable_summary_table.png"), vwidth = 1600, vheight = 800, cliprect = "viewport")
# Display the reordered table

# tsd_pbs_stat column
t <- t %>%
  mutate(tsd_pbs_stat = case_when(
    tsd == "True" & pbs == "True" ~ "TSD_PBS",
    tsd == "True" & pbs == "False" ~ "TSD_only",
    tsd == "False" & pbs == "True" ~ "PBS_only",
    tsd == "False" & pbs == "False" ~ "NONE"
  ))
# filter superfamilies
t_copia <- subset(t,t$te_sfam=="Ty1/copia")
t_copia$te_fam <- factor(t_copia$te_fam,levels = copiaOrder)
t_gypsy <- subset(t,t$te_sfam=="Ty3/gypsy")
t_gypsy$te_fam <- factor(t_gypsy$te_fam,levels = gypsyOrder)

# Family specific (non)autonomous TEs
copiaPl <- ggplot(t_copia,aes(x=te_fam,fill=autonomy_stat)) +
  geom_bar() +
  scale_fill_manual(values = c("red","blue"),
                    name="Prot. domains status:") +
  ggtitle("Ty1/copia") +
  ylab("# of TEs") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="None")

gypsyPl <- ggplot(t_gypsy,aes(x=te_fam,fill=autonomy_stat)) +
  geom_bar() +
  scale_fill_manual(values = c("red","blue"),
                    name="Prot. domains status:") +
  ggtitle("Ty3/gypsy") +
  ylab("# of TEs") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank())

png(paste0(out_folder,"autnom_nonauton_summary.png"),res=200,units = 'in',width = 9,height = 5)
copiaPl + gypsyPl
dev.off()

t_copia$tsd_pbs_stat <- factor(t_copia$tsd_pbs_stat, levels = c("TSD_PBS","TSD_only","PBS_only","NONE"))
t_gypsy$tsd_pbs_stat <- factor(t_gypsy$tsd_pbs_stat, levels = c("TSD_PBS","TSD_only","PBS_only","NONE"))

copiaPl2 <- ggplot(t_copia,aes(x=te_fam,fill=tsd_pbs_stat)) +
  geom_bar() +
  scale_fill_manual(values = c("#f4895f","#f8e16f","#95cf92","#369acc"),
                    name="TSD & PBS status:") +
  ggtitle("Ty1/copia") +
  ylab("# of TEs") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="None")

gypsyPl2 <- ggplot(t_gypsy,aes(x=te_fam,fill=tsd_pbs_stat)) +
  geom_bar() +
  scale_fill_manual(values = c("#f4895f","#f8e16f","#95cf92","#369acc"),
                    name="TSD & PBS status:") +
  ggtitle("Ty3/gypsy") +
  ylab("# of TEs") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank())

png(paste0(out_folder,"tsd_pbs_summary.png"),res=200,units = 'in',width = 9,height = 5)
copiaPl2 + gypsyPl2
dev.off()

# LTR identity
copiaFamCnt <- length(copiaOrder)
copiaPl3 <- ggplot(t_copia, aes(x=te_fam,y=ltr_identity,fill=te_fam)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')[1:copiaFamCnt]) +
  ggtitle("Ty1/copia") +
  ylab("LTR identity [%]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position="None")

gypsyFamCnt <- length(gypsyOrder)
gypsyPl3 <- ggplot(t_gypsy, aes(x=te_fam,y=ltr_identity,fill=te_fam)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')[1:gypsyFamCnt]) +
  ggtitle("Ty3/gypsy") +
  ylab("LTR identity [%]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position="None")

# TE age MYA
copiaPl4 <- ggplot(t_copia, aes(x=te_fam,y=MYA,fill=te_fam)) +
  geom_violin() +
  scale_fill_manual(values = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')[1:copiaFamCnt]) +
  scale_y_reverse() +
  #ggtitle("Ty1/copia") +
  ylab("TE age [MYA]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="None")

gypsyPl4 <- ggplot(t_gypsy, aes(x=te_fam,y=MYA,fill=te_fam)) +
  geom_violin() +
  scale_fill_manual(values = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')[1:gypsyFamCnt]) +
  scale_y_reverse() +
  #ggtitle("Ty3/gypsy") +
  ylab("TE age [MYA]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="None")

png(paste0(out_folder,"ltr_ident_mya.png"),res=200,units = 'in',width = 12,height = 7)
# (p1 | p2)/(p3 | p4)
(copiaPl3 | gypsyPl3) /(copiaPl4 | gypsyPl4)
dev.off()

# TE length
copiaPl5 <- ggplot(t_copia, aes(x=te_fam,y=te_length,fill=te_fam)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')[1:copiaFamCnt]) +
  ggtitle("Ty1/copia") +
  ylab("TE length [bp]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position="None")

gypsyPl5 <- ggplot(t_gypsy, aes(x=te_fam,y=te_length,fill=te_fam)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')[1:gypsyFamCnt]) +
  ggtitle("Ty3/gypsy") +
  ylab("TE length [bp]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position="None")

# LTR length
copiaPl6 <- ggplot(t_copia, aes(x=te_fam,y=ltr_avg_len,fill=te_fam)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')[1:copiaFamCnt]) +
  #ggtitle("Ty1/copia") +
  ylab("LTR length [bp]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="None")

gypsyPl6 <- ggplot(t_gypsy, aes(x=te_fam,y=ltr_avg_len,fill=te_fam)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')[1:gypsyFamCnt]) +
  #ggtitle("Ty3/gypsy") +
  ylab("LTR length [bp]") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="None")

png(paste0(out_folder,"te_ltr_lengths.png"),res=200,units = 'in',width = 12,height = 7)
# (p1 | p2)/(p3 | p4)
(copiaPl5 | gypsyPl5) /(copiaPl6 | gypsyPl6)
dev.off()

# Chromosome specific abundance
## pick-up bins
# Determine column names based on the value of longest_seq_length
# less than 80Mbp -> 100kbp and more than 80Mbp 1Mbp
if (longest_seq_length <= 80) {
  xc_col <- t_copia$bin_100kbp
  yc_col <- t_copia$inters_100kbp
  sec <- '100kbp'
} else {
  xc_col <- t_copia$bin_1Mbp
  yc_col <- t_copia$inters_1Mbp
  sec <- '1Mbp'
}

if (longest_seq_length <= 80) {
  xg_col <- t_gypsy$bin_100kbp
  yg_col <- t_gypsy$inters_100kbp
} else {
  xg_col <- t_gypsy$bin_1Mbp
  yg_col <- t_gypsy$inters_1Mbp
}


twentyCols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#f5fffa')
gypsyCols <- twentyCols[1:gypsyFamCnt]
cs <- gypsyFamCnt + 1
ce <- (gypsyFamCnt+copiaFamCnt)+1
copiaCols <-twentyCols[cs:ce]

copiaPl7 <- ggplot(t_copia,aes(x=xc_col,y=yc_col,fill=te_fam)) +
  geom_bar(stat='identity', position = 'stack') +
  ggtitle("Ty1/copia") +
  scale_fill_manual(values = copiaCols,
                    name="TE family:") +
  ylab(paste0("TE density [bp per ",sec,"]")) + 
  xlab(paste0("Chromosome sections [",sec,"]")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0)) +
  facet_grid(rows = vars(chromosome))

gypsyPl7 <- ggplot(t_gypsy,aes(x=xg_col,y=yg_col,fill=te_fam)) +
  geom_bar(stat='identity', position = 'stack') +
  ggtitle("Ty3/gypsy") +
  scale_fill_manual(values = gypsyCols,
                    name="TE family:") +
  ylab(paste0("TE density [bp per ",sec,"]")) + 
  xlab(paste0("Chromosome sections [",sec,"]")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0)) +
  facet_grid(rows = vars(chromosome))

png(paste0(out_folder,"te_cnt_chr_specific_density.png"),res=200,units = 'in',width = 17,height = 18)
copiaPl7 / gypsyPl7
dev.off()

# Plot age category
ageCatCols <- rev(c('#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026'))
copiaPl8 <- ggplot(t_copia,aes(x=xc_col,y=yc_col,fill=age_cat)) +
  geom_bar(stat='identity', position = 'stack') +
  ggtitle("Ty1/copia") +
  scale_fill_discrete(name = "Insertion time [MYA]:",
                      type = ageCatCols,
                      labels = c("0.0-1.0","1.0-2.0","2.0-3.0","3.0-4.0","4.0-5.0","5.0-10.0",">10.0")) +
  ylab(paste0("TE density [bp per ",sec,"]")) + 
  xlab(paste0("Chromosome sections [",sec,"]")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0)) +
  facet_grid(rows = vars(chromosome))


gypsyPl8 <- ggplot(t_gypsy,aes(x=xg_col,y=yg_col,fill=age_cat)) +
  geom_bar(stat='identity', position = 'stack') +
  ggtitle("Ty3/gypsy") +
  scale_fill_discrete(name = "Insertion time [MYA]:",
                      type = ageCatCols,
                      labels = c("0.0-1.0","1.0-2.0","2.0-3.0","3.0-4.0","4.0-5.0","5.0-10.0",">10.0")) +
  ylab(paste0("TE density [bp per ",sec,"]")) + 
  xlab(paste0("Chromosome sections [",sec,"]")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0),
        legend.title = element_blank(),
        legend.position="None") +
  facet_grid(rows = vars(chromosome))

png(paste0(out_folder,"te_age_chr_specific_density.png"),res=200,units = 'in',width = 17,height = 18)
copiaPl8 / gypsyPl8
dev.off()

# TE age chr and fam specific

copiaPl9 <- ggplot(t_copia,aes(x=xc_col,y=yc_col,fill=age_cat)) +
  geom_bar(stat='identity', position = 'stack') +
  ggtitle("Ty1/copia") +
  scale_fill_discrete(name = "Insertion time [MYA]:",
                      type = ageCatCols,
                      labels = c("0.0-1.0","1.0-2.0","2.0-3.0","3.0-4.0","4.0-5.0","5.0-10.0",">10.0")) +
  ylab(paste0("TE density [bp per ",sec,"]")) + 
  xlab(paste0("Chromosome sections [",sec,"]")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0)) +
  facet_grid(rows = vars(te_fam),
             cols = vars(chromosome))


gypsyPl9 <- ggplot(t_gypsy,aes(x=xg_col,y=yg_col,fill=age_cat)) +
  geom_bar(stat='identity', position = 'stack') +
  ggtitle("Ty3/gypsy") +
  scale_fill_discrete(name = "Insertion time [MYA]:",
                      type = ageCatCols,
                      labels = c("0.0-1.0","1.0-2.0","2.0-3.0","3.0-4.0","4.0-5.0","5.0-10.0",">10.0")) +
  ylab(paste0("TE density [bp per ",sec,"]")) + 
  xlab(paste0("Chromosome sections [",sec,"]")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),  axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0),
        legend.title = element_blank(),
        legend.position="None") +
  facet_grid(rows = vars(te_fam),
             cols = vars(chromosome))

png(paste0(out_folder,"te_age_chr_fam_spec_density.png"),res=200,units = 'in',width = 25,height = 25)
copiaPl9 / gypsyPl9
dev.off()