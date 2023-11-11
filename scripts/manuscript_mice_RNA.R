# Libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)
library(fgsea)
library(gplots)
library(cowplot)

# Load data & functions
load("data/ATR_AZD_data.RData")
load("scripts/functions/manuscript_functions.RData")

# Volcano
############

# Label same genes as in Nat Comm paper
genes_to_label<- sort(c("BRCA1", "BRCA2", "E2F1", "E2F2", "E2F8", "TOP2A", "MCM4", "MCM5", "MCM6", "MCM10", "MKI67", "CCNB1", "CCNB2", "CCNA2", "FOXM1", "BRIP1", "BARD1", "FANCD2", "FANCI", "FANCB", "CDKN1A","FAS", "TP53INP1", "GDF15", "NOTCH1"))
genes_to_label_MGI<- MGI_to_HGNC[MGI_to_HGNC$HGNC.symbol%in%genes_to_label,"MGI.symbol"]
# genes_to_label<- sort(c("BRCA1", "BRCA2", "E2F1", "E2F2", "E2F8", "TOP2A", "MCM4", "MCM5", "MCM6", "MCM10", "MKI67", "CCNB1", "CCNB2", "CCNA2", "FOXM1", "BRIP1", "BARD1", "FANCD2", "FANCI", "FANCB", "CDKN1A","FAS", "TP53INP1", "GDF15", "NOTCH1"))

# Plot volcano
genes_DE<- list()
res_proc<- as.data.frame(res_diff_expr_mice$AZD)
p_tmp<- plot_volcano(DE_results = res_proc, gene = genes_to_label_MGI, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5, plotTH=F)
p_tmp<- p_tmp + 
  # ggtitle(label = "Alk-F1178S;Th-MYCN") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7)
  ) +
  scale_x_continuous(name = "log2(Fold Change)", limits = c(-10,10)) +
  scale_y_continuous(name = "-log10(Padj)")
p_volc<- p_tmp
p_volc

# n DE?
res_proc<- na.omit(res_proc)
genes_DE<- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$MGI, th_logFC = 1.5, th_logP = 2, curve = 0.5)

sapply(genes_DE, "length")
# up down   DE 
# 1252  487 1739 

# fGSEA + RS plots main pathways
##################################

pws<- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_P53_PATHWAY"
)

res_proc<- as.data.frame(res_diff_expr_mice$AZD)
stat<- get_GSEA_stat(res_proc)
GSEA_res<- do_fGSEA(Ha_MGI_ls, stat)
GSEA_res$isLabel<- F
# GSEA_res$isLabel[GSEA_res$padj< 1e-40]<- T
GSEA_res$isLabel[1]<- T # Label top one

p<- ggplot(GSEA_res, aes(x = NES, y = -log10(padj), key=pathway, colour=isLabel)) +
  # ggtitle(label = "Alk-F1178S;Th-MYCN") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Normalized Enrichment Score") + 
  ylab("-log10(Padj)") + 
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
  ) +
  geom_point(size=2)

# Label
p<- p + 
  geom_text_repel(data = subset(GSEA_res, isLabel),
                  aes(label = pathway, fontface=3),
                  size = 2,
                  colour="black",
                  box.padding   = 0.5, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_color_manual(values=c("grey", "black"))

p_GSEA<- p
p_GSEA

# RS
for(j in 1:length(pws)){
  p_pw<- format(signif(GSEA_res[GSEA_res$pathway==pws[j],"padj"],3),scientific=T)
  p<- plotEnrichment(Ha_MGI_ls[[pws[j]]],stat,ticksSize=.1,ticksLength = 0.2) +
    # ggtitle(paste0(pws[j],"\n(Padj=",p_pw,")")) + 
    geom_line(size=0.5, col="green") +
    theme(
      plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7),
    ) +
    scale_x_continuous(name = "Rank") +
    scale_y_continuous(name = "Enrichment Score", limits = c(-1,1))
  assign(paste0("p_RS_",j), p)
  }

# Transcription factors?
##################################

GSEA_res<- do_fGSEA(TFT_MGI_ls, stat)
GSEA_res$isLabel<- F
# GSEA_res$isLabel[GSEA_res$padj< 1e-40]<- T
GSEA_res$isLabel[1:5]<- T # Label top one

p<- ggplot(GSEA_res, aes(x = NES, y = -log10(padj), key=pathway, colour=isLabel)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Normalized Enrichment Score") + 
  ylab("-log10(Padj)") + 
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
  ) +
  geom_point(size=2)

# Label
p<- p + 
  geom_text_repel(data = subset(GSEA_res, isLabel),
                  aes(label = gsub("_TARGET_GENES","",pathway), fontface=3),
                  size = 2,
                  colour="black",
                  box.padding   = 0.5, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_color_manual(values=c("grey", "black"))

p_GSEA_TFT<- p
p_GSEA_TFT

# RS
TFs<- c(
  "ING1_TARGET_GENES"
)

for(j in 1:length(TFs)){
  p_pw<- format(signif(GSEA_res[GSEA_res$pathway==TFs[j],"padj"],3),scientific=T)
  p<- plotEnrichment(TFT_MGI_ls[[TFs[j]]],stat,ticksSize=.1,ticksLength = 0.2) +
    # ggtitle(paste0(pws[j],"\n(Padj=",p_pw,")")) + 
    geom_line(size=0.5, col="green") +
    theme(
      plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7),
    ) +
    scale_x_continuous(name = "Rank") +
    scale_y_continuous(name = "Enrichment Score", limits = c(-1,1))
  assign(paste0("p_RS_FTF_",j), p)
}


# Previous chemical and genetic pertubations?
################################################

GSEA_res<- do_fGSEA(CGP_MGI_ls, stat)
GSEA_res$isLabel<- F
GSEA_res$isLabel[1]<- T # Label top one

p<- ggplot(GSEA_res, aes(x = NES, y = -log10(padj), key=pathway, colour=isLabel)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Normalized Enrichment Score") + 
  ylab("-log10(Padj)") + 
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
  ) +
  geom_point(size=2)

# Label
p<- p + 
  geom_text_repel(data = subset(GSEA_res, isLabel),
                  aes(label = pathway, fontface=3),
                  size = 2,
                  colour="black",
                  box.padding   = 0.5, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_color_manual(values=c("grey", "black"))

p_GSEA_CGP<- p
p_GSEA_CGP

# RS
CGPs<- c(
  "LAZARO_GENETIC_MOUSE_MODEL_HIGH_GRADE_LARGE_CELL_NEUROENDOCRINE_LUNG_CARCINOMA_UP"
)

for(j in 1:length(CGPs)){
  p_pw<- format(signif(GSEA_res[GSEA_res$pathway==CGPs[j],"padj"],3),scientific=T)
  p<- plotEnrichment(CGP_MGI_ls[[CGPs[j]]],stat,ticksSize=.1,ticksLength = 0.2) +
    # ggtitle(paste0(pws[j],"\n(Padj=",p_pw,")")) + 
    geom_line(size=0.5, col="green") +
    theme(
      plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7),
    ) +
    scale_x_continuous(name = "Rank") +
    scale_y_continuous(name = "Enrichment Score", limits = c(-1,1))
  assign(paste0("p_RS_CGP_",j), p)
}

# Create final figure
######################
p_GSEA_panels<- plot_grid(
  p_GSEA, p_RS_1 + scale_y_continuous(name = "Enrichment Score", limits = c(-1,0.1)),
  ncol = 2,
  labels = c("H","I"))

p_RNA_panels<- plot_grid(
  p_volc,
  p_GSEA_panels,
  ncol = 1,
  rel_heights = c(2,1),
  labels = c("G",NA)
  )

# p_TFT_panels<- plot_grid(
#   p_GSEA_TFT, p_RS_FTF_1 + scale_y_continuous(name = "Enrichment Score", limits = c(-1,0.1)),NA,
#   ncol = 1,
#   labels = c("K","L"))
# 
p<- plot_grid(
  p_RNA_panels, NA, NA,
  rel_widths = c(2,1,1),
  ncol = 3
)

ggsave("results/figs/manuscript_fig5GI_mice_RNA.pdf", p, device = cairo_pdf, width = 178, height = 2/5 * 265, units = "mm")

# Add results to suppl. table
##################################

res_BAY<- as.data.frame(res_diff_expr_mice$BAY)[,c(1,2,5,6)]
res_AZD<- as.data.frame(res_diff_expr_mice$AZD)[,c(1,2,5,6)]

# Save
WriteXLS::WriteXLS(c("res_BAY","res_AZD"),"results/tables/manuscript_summary_table_RNA_mice.xlsx",row.names = T, SheetNames = c("elimusertib","ceralasertib"))


