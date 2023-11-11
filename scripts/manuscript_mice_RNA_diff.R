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

# GSEA: differentiation?
#########################

# MSigDB: Chemical & Genetic perturbations? Cellular signature?
# Panglao: Cells?

p_ls<- list()
GSEA_ls<- list()
gsets<- c("CGP_MGI_ls", "CSign_MGI_ls", "PDB_MGI_ls", "TFT_MGI_ls")
gsets_names<- c("chemical and genetic perturbations", "cell type signature genesets", "PanglaoDB", "TFT")
for(i in 1:length(gsets)){
  gs<- gsets[i]
  gs_name<- gsets_names[i]
  for(drug in c("BAY","AZD")){
    res_proc<- as.data.frame(res_diff_expr_mice[[drug]])
    stat<- get_GSEA_stat(res_proc)
    GSEA_res<- do_fGSEA(get(gs), stat)
    GSEA_res$isLabel<- F
    if(gs=="TFT_MGI_ls") GSEA_res$isLabel[GSEA_res$padj< 1e-20]<- T
    GSEA_res$isLabel[grepl("NEURON",toupper(GSEA_res$pathway))]<- T
    GSEA_res$isLabel[grepl("SCHWANN",toupper(GSEA_res$pathway))]<- T
    GSEA_res$isLabel[grepl("NEUROENDOCRINE",toupper(GSEA_res$pathway))]<- T
    GSEA_ls[[drug]][[gs]]<- GSEA_res
    
    p<- ggplot(GSEA_res, aes(x = NES, y = -log10(padj), key=pathway, colour=isLabel)) +
      ggtitle(label = gs_name) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("normalized enrichment score") + 
      ylab("-log10(Padj)") + 
      theme(
        legend.position="none",
        plot.title = element_text(hjust = 0.5, size=8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size=6),  
        axis.title = element_text(size=7),
      ) +
      geom_point(size=2) +
      geom_point(data = subset(GSEA_res, subset = isLabel==T), mapping = aes(x = NES, y = -log10(padj), colour=isLabel), size=2)
    
    # Label
    GSEA_subset<- subset(GSEA_res, isLabel)
    GSEA_subset$pathway<- tolower(GSEA_subset$pathway)
    GSEA_subset$pathway<- gsub("descartes_organogenesis","DO",GSEA_subset$pathway)
    GSEA_subset$pathway<- gsub("_target_genes","",GSEA_subset$pathway)
    p<- p + 
      geom_text_repel(data = GSEA_subset,
                      aes(label = pathway, fontface=3),
                      size = 2,
                      colour="black",
                      box.padding   = 0.5, 
                      point.padding = 0.5,
                      segment.color = 'grey50',
                      max.overlaps = Inf) +
      scale_color_manual(values=c("grey", "black"))
    p_ls[[drug]][[gs]]<- p
  }

}

# p_ls$BAY$CGP_MGI_ls
# p_ls$BAY$CSign_MGI_ls
# p_ls$BAY$PDB_MGI_ls
# p_ls$BAY$TFT_MGI_ls
# 
# p_ls$AZD$CGP_MGI_ls
# p_ls$AZD$CSign_MGI_ls # Less neurons
# p_ls$AZD$PDB_MGI_ls

# View(GSEA_ls$BAY$CSign_MGI_ls)
# View(GSEA_ls$BAY$PDB_MGI_ls)

# RS plots main pathways
##################################

GSEA_dbs<- c("CSign_MGI_ls", "CSign_MGI_ls", "CGP_MGI_ls", "CGP_MGI_ls", "PDB_MGI_ls","PDB_MGI_ls","TFT_MGI_ls")
gs_names<- c(
  "DESCARTES_ORGANOGENESIS_SCHWANN_CELL_PRECURSOR", 
  "DESCARTES_ORGANOGENESIS_EXCITATORY_NEURONS",
  "LAZARO_GENETIC_MOUSE_MODEL_HIGH_GRADE_LARGE_CELL_NEUROENDOCRINE_LUNG_CARCINOMA_UP",
  "LEIN_NEURON_MARKERS",
  "Neurons",
  "Schwann_cells",
  "ING1_TARGET_GENES"
)

p_ls2<- list()
for(i in 1:length(GSEA_dbs)){
  for(drug in c("AZD","BAY")){
    gs_name<- gs_names[i]
    GSEA_db<- GSEA_dbs[i]
    gs_sel<- get(GSEA_db)[[gs_name]]
    res_proc<- as.data.frame(res_diff_expr_mice[[drug]])
    stat<- get_GSEA_stat(res_proc)
    
    p_pw<- as.numeric(GSEA_ls[[drug]][[GSEA_db]][GSEA_ls[[drug]][[GSEA_db]]$pathway==gs_name,"padj"])
    ES_sign<- sign(as.numeric(GSEA_ls[[drug]][[GSEA_db]][GSEA_ls[[drug]][[GSEA_db]]$pathway==gs_name,"NES"]))
    ymax<- round(max(abs(as.numeric(GSEA_ls[[drug]][[GSEA_db]][GSEA_ls[[drug]][[GSEA_db]]$pathway==gs_name,"ES"]))),digits = 1)+0.1
    gs_name<- tolower(gs_name)
    gs_name<- gsub("descartes_organogenesis","DO",gs_name)
    p<- plotEnrichment(gs_sel,stat,ticksSize=.1,ticksLength = 0.2) +
      # ggtitle(paste0(drug, "\n",gs_name,"\n(Padj=",p_pw,")")) +
      ggtitle(gs_name) +
      geom_line(size=0.5, col="green") +
      theme(
        plot.title = element_text(hjust = 0.5, size=6, face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size=6),
        axis.title = element_blank()
        # axis.title = element_text(size=7),
      ) +
      # scale_x_continuous(name = "rank") +
      if(ES_sign==1) scale_y_continuous(limits = c(-0.1,ymax))
      if(ES_sign== -1) scale_y_continuous(limits = c(-ymax,0.1))
    p_ls2[[drug]][[gs_name]]<- p
  }
}

# Merge & save plots
####################

# Main figure: cellular signatures
p_GSEA_main<- plot_grid(
  p_ls$BAY$PDB_MGI_ls,
  p_ls$BAY$CSign_MGI_ls,
  ncol = 2,
  labels = c("C", "D")    
)

p_rs_main<- plot_grid(
  p_ls2$BAY$neurons,
  p_ls2$BAY$schwann_cells,  
  p_ls2$BAY$DO_schwann_cell_precursor,
  p_ls2$BAY$DO_excitatory_neurons,  
  ncol = 4
)

p<- plot_grid(
  p_GSEA_main, p_rs_main,
  rel_heights = c(2,1.2),
  ncol = 1
)

ggsave("results/figs/manuscript_fig6DE_mice_RNA_diff.pdf", p, device = cairo_pdf, width = 2/3*178, height = 265/3, units = "mm")

# TFT: not in manuscript
p_RS_TFT<- p_ls2$BAY$ing1_target_genes + scale_y_continuous(name = "Enrichment Score", limits = c(-1,0.1))


p_TFT<- plot_grid(
  p_ls$BAY$TFT_MGI_ls, p_RS_TFT,NA,
  ncol = 1,
  labels = c("A","B")
)

p_TFT<- plot_grid(
  p_TFT, NA,
  ncol = 2,
  rel_widths = c(1,3)
)

ggsave("results/figs/manuscript_mice_RNA_TFT.pdf", p_TFT, device = cairo_pdf, width = 178, height = 265/2.5, units = "mm")

# Expression?
sapply(res_diff_expr_mice, function(x) x["Ing1",c("log2FoldChange","pvalue")])
# $BAY
# log2FoldChange      pvalue
# -0.669872 5.44769e-05
# 
# $AZD
# log2FoldChange      pvalue
# -0.779928 1.84151e-05

# Suppl figure: CGP
p_rs_suppl<- plot_grid(
  p_ls2$BAY$lazaro_genetic_mouse_model_high_grade_large_cell_neuroendocrine_lung_carcinoma_up + scale_y_continuous(name = "Enrichment Score", limits = c(-1,0.1)),
  p_ls2$BAY$lein_neuron_markers + scale_y_continuous(name = "Enrichment Score", limits = c(-0.1,1)),  
  ncol = 2
)

p_suppl<- plot_grid(
  p_ls$BAY$CGP_MGI_ls, p_rs_suppl,NA,
  ncol = 1,
  rel_heights = c(2,1,3)
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Fig S5: ATR#2, 2023",
    fontface = 'italic',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p_suppl<- plot_grid(
  p_suppl, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

ggsave("results/figs/manuscript_figS5_mice_RNA_CGP.pdf", p_suppl, device = cairo_pdf, width = 178, height = 265, units = "mm")


