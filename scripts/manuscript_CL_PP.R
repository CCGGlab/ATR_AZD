# Load libraries, functions & data
###################################
library(ggplot2)
library(cowplot)
library(DEGreport)
library(ggrepel)
library(ggvenn)
library(ggpubr)
load("data/ATR_AZD_data.RData")
load("scripts/functions/manuscript_functions.RData")

# General description
#####################
prot_sites<- rownames(normalized_counts_CL_PP)
prots<- unique(gsub("_.*","",prot_sites))
sites<- gsub(".*_","",prot_sites)
aa<- substr(sites,1,1)

length(prot_sites) # 11026
length(prots) # 3244
table(aa)
# S    T    Y 
# 9655 1351   20 

# Plot Volcano
##############

p_volc_ls<- list()
sites_DP_ls<- list()
for(cond in c("BAY1895344_50nM", "AZD6738_50nM", "AZD6738_1uM")){
    
  # Get data
  res_proc<- res_diff_expr_CL_PP[[cond]]
  
  # Get DE
  sites_DP<- get_DE(res_proc$logFC,res_proc$FDR,rownames(res_proc),th_logFC=0.3,th_logP= -log10(0.05),curve=0.1)
  cat(cond, "DP:", length(sites_DP$DE), "Up:", length(sites_DP$up), "Down:", length(sites_DP$down),"Prots DP:",length(unique(gsub("_.*","",sites_DP$DE))),"Prots up:",length(unique(gsub("_.*","",sites_DP$up))),"Prots down:",length(unique(gsub("_.*","",sites_DP$down))),"\n")
  sites_DP_ls[[cond]]<- sites_DP 
    
  # Use protein names for labels (not gene names)
  rownames(res_proc)<- gsub("PRKDC_","DNAPK_",rownames(res_proc))
  rownames(res_proc)<- gsub("CHEK1_","CHK1_",rownames(res_proc))
  rownames(res_proc)<- gsub("CHEK2_","CHK2_",rownames(res_proc))
  
  # Sites in discussion and/or labelled
  prot_sel<- rownames(res_proc)[grep("ATR_|ATM_|CHK1_|CHK2_|DNAPK_",rownames(res_proc))]

  res_proc$dummy<- NA # Format for plot_volcano function,
  p_volc<- plot_volcano(DE_results = res_proc, gene=prot_sel, isProt=T, useHyperbolicTH = T, p_cu = 0.05, logFC_cu = 0.3, curve = 0.1, labelAll = T, labelCol="black", plotTH = F, plot_nominal_p = T) # Plot unadjusted P for visualization purposes in volcano
  p_volc<- p_volc +
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
    scale_x_continuous(name = "log2(Fold Change)", limits = c(-2.0,2.8)) +
    scale_y_continuous(name = "-log10(P)")
  
  # Label ATM/ATR in red
  idx_sel<- grep("ATR_|ATM_|DNAPK_",rownames(res_proc))
  df_sel<- cbind(res_proc[idx_sel,],gene_id=rownames(res_proc[idx_sel,]))
  colnames(df_sel)[c(1,3)]<- c("log2FoldChange","padj")
  p_volc<- p_volc + geom_point(data=df_sel, colour="red", size=2) # this adds a red point
  
  # add title
  p_volc<- p_volc + ggtitle(paste0(gsub("_"," (",cond),")"))
  
  # Same scales
  p_volc_ls[[cond]]<- p_volc + 
    ylim(0,17)
}

# BAY1895344_50nM DP: 618 Up: 521 Down: 97 Prots DP: 368 Prots up: 301 Prots down: 79 
# AZD6738_50nM DP: 19 Up: 8 Down: 11 Prots DP: 16 Prots up: 6 Prots down: 10 
# AZD6738_1uM DP: 276 Up: 242 Down: 34 Prots DP: 173 Prots up: 146 Prots down: 29 

sapply(res_diff_expr_CL_PP, function(x) x["ATR_T1989",])
# AZD6738_1uM  AZD6738_50nM BAY1895344_50nM
# logFC  -0.416       -0.0228      -1.74          
# pvalue 8.383041e-06 0.7311025    6.668585e-15   
# FDR    0.194        0.949        2.25e-14   


# PK predictions 
#################

# Load PK substrate data
PP_Atlas <- as.data.frame(readxl::read_excel("downloads/PP_Atlas/41586_2022_5575_MOESM5_ESM.xlsx", 
                                             sheet = "Supplementary Table 3"))
PP_Atlas$site<- paste0(PP_Atlas$Gene,"_",PP_Atlas$Phosphosite)

PK_res_ls<- list()
for(cond in names(res_diff_expr_CL_PP)){
  cat(cond,"\n")
  # Get result
  res_PP<- res_diff_expr_CL_PP[[cond]]
  # PK prediction
  PK_pred<- do_PK_enrich_volcano(PP_atlas = PP_atlas, res_PP = res_PP, cu_pct = 99, plot_labels = c("ATM","ATR","DNAPK","SMG1","P70S6K","P70S6KB","RSK4"))
  # add to results
  PK_res_ls[[cond]]<- PK_pred
}

# PK_res_ls$BAY1895344_50nM$p
# PK_res_ls$AZD6738_1uM$p
# View(PK_res_ls$BAY1895344_50nM$res)
# View(PK_res_ls$AZD6738_1uM$res)

# Show PK sites
################
p_PK_sites_ls<- list()
for(PK in c("ATM","ATR","DNAPK","RPS6K")){
  # Get sites
  PK_sites<- PP_Atlas$site[which(PP_Atlas[[paste0(PK,"_percentile")]] > 99)]
  if(PK=="RPS6K") PK_sites<- PP_Atlas$site[unique(c(
    which(PP_Atlas[[paste0("RSK4","_percentile")]] > 99),
    which(PP_Atlas[[paste0("P70S6K","_percentile")]] > 99),
    which(PP_Atlas[[paste0("P70S6KB","_percentile")]] > 99)))]
  
  # Get PP data
  for(cond in names(res_diff_expr_CL_PP)){
    res_proc<- res_diff_expr_CL_PP[[cond]]
    res_proc$site<- rownames(res_proc)
    res_proc$isLabel<- res_proc$site%in%PK_sites
    
    # Plot
    p_volc<- ggplot() +
      geom_point(data = subset(res_proc,isLabel==F), mapping = aes(x = logFC, y = -log10(pvalue)), colour="#d0d3d4") +
      geom_point(data = subset(res_proc,isLabel==T), mapping = aes(x = logFC, y = -log10(pvalue)), colour="#3498db") +
      # ggtitle(paste0(cond, ": ",PK, " sites")) +
      scale_color_manual(values = c("#d0d3d4", "#3498db")) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size=6),  
        axis.title = element_text(size=7)
      ) +
      scale_x_continuous(name = "log2(Fold Change)", limits = c(-2.0,3.0)) +
      scale_y_continuous(name = "-log10(P)", limits = c(0,17)) 
    # +
      # geom_text_repel(data = subset(res_proc, res_proc$isLabel),
      #                 aes(x = logFC, y = -log10(pvalue), label = site),
      #                 size = 2,
      #                 colour="black",
      #                 box.padding   = 0.5,
      #                 point.padding = 0.5,
      #                 max.overlaps = 25,
      #                 segment.color = 'grey50')
    p_PK_sites_ls[[paste0(cond,"_",PK)]]<- p_volc
    res_proc_sel<- res_proc[res_proc$isLabel,]
    assign(paste0("PK_sites_",cond,"_",PK), as.data.frame(res_proc_sel[order(res_proc_sel$logFC),]))
  }
}

# Correlation between data?
############################

# Data frame
common_sites<- intersect(rownames(res_diff_expr_CL_PP$BAY1895344_50nM), rownames(res_diff_expr_CL_PP$AZD6738_1uM))
cor_df<- data.frame(
  site = common_sites,
  BAY_50nM = res_diff_expr_CL_PP$BAY1895344_50nM[common_sites,"logFC"],
  AZD_50nM = res_diff_expr_CL_PP$AZD6738_50nM[common_sites,"logFC"],
  AZD_1µM = res_diff_expr_CL_PP$AZD6738_1uM[common_sites,"logFC"]
)

# Use protein names for labels (not gene names)
cor_df$site<- gsub("PRKDC_","DNAPK_",cor_df$site)
cor_df$site<- gsub("CHEK1_","CHK1_",cor_df$site)
cor_df$site<- gsub("CHEK2_","CHK2_",cor_df$site)

# Label
cor_df$isLabel<- grepl("ATR_|ATM_|CHK1_|CHK2_|DNAPK_",cor_df$site)

# Plot
p_BAY_AZD50 <- ggplot(cor_df, aes(BAY_50nM,AZD_50nM, label=site)) +
  geom_point() + 
  geom_smooth(method = lm) + # Add regression line
  xlab(paste("Log2FC elimusertib (50 nM)")) +
  ylab(paste("Log2FC ceralasertib (50 nM)")) +
  # geom_cor(method = "pearson",cex=3) +
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
  geom_point(data=subset(cor_df, isLabel), colour="red", size=2) + # this adds a red point
  geom_label_repel(data=subset(cor_df, isLabel), aes(label=site), colour="red",label.size=NA,fontface = "italic",size=2,box.padding = 0.05,label.padding = 0.05) +
  stat_cor(method = "pearson",r.digits=2, p.digits = 2,size=2)

p_BAY_AZD1 <- ggplot(cor_df, aes(BAY_50nM,AZD_1µM, label=site)) +
  geom_point() + 
  geom_smooth(method = lm) + # Add regression line
  xlab(paste("Log2FC elimusertib (50 nM)")) +
  ylab(paste("Log2FC ceralasertib (1 µM)")) +
  # geom_cor(method = "pearson",cex=3) +
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
  geom_point(data=subset(cor_df, isLabel), colour="red", size=2) + # this adds a red point
  geom_label_repel(data=subset(cor_df, isLabel), aes(label=site), colour="red",label.size=NA,fontface = "italic",size=2,box.padding = 0.05,label.padding = 0.05) +
  stat_cor(method = "pearson",r.digits=2, p.digits = 2,size=2) 


# Create figure
###############

# Volcano panel
p_volc<- plot_grid(
  p_volc_ls$BAY1895344_50nM + ggtitle("elimusertib (50 nM)"), 
  p_volc_ls$AZD6738_50nM + ggtitle("ceralasertib (50 nM)"), 
  p_volc_ls$AZD6738_1uM + ggtitle("ceralasertib (1 µM)"),
  ncol=3
)
# 
# # Cor panel
# p_cor<- plot_grid(
#   p_BAY_AZD50, p_BAY_AZD1,NA,
#   ncol=3,
#   rel_widths = c(1,1,1.5)
# )

# PK panel
# equal y scale
p_PK_BAY<- PK_res_ls$BAY1895344_50nM$p + ggtitle("elimusertib (50nM)") + theme(strip.background = element_blank(),strip.text.x = element_blank())
p_PK_AZD<- PK_res_ls$AZD6738_1uM$p + ggtitle("ceralasertib (1µM)") + theme(strip.background = element_blank(),strip.text.x = element_blank())
ymax<- ceiling(max(c(layer_scales(p_PK_BAY)$y$range$range,layer_scales(p_PK_AZD)$y$range$range)))
p_PK_BAY<- p_PK_BAY + ylim(0,ymax)
p_PK_AZD<- p_PK_AZD + ylim(0,ymax)

p_PK<- plot_grid(
  p_PK_BAY, 
  p_PK_AZD,
  NA,
  ncol=3,
  rel_widths = c(1,1,1.5)
)

# Venn diagram
names(sites_DP_ls)<- c("elimusertib (50 nM)", "ceralasertib 50 nM", "ceralasertib 1 µM")
p_venn<- ggvenn(sapply(sites_DP_ls,function(x) x$DE), show_percentage = F, stroke_size = 0, set_name_size = 2, text_size = 2)

# Merge venn with corr plot
p_venn_cor<- plot_grid(
  p_venn, p_BAY_AZD1,NA,
  ncol=3,
  rel_widths = c(1,1,1.5),
  labels=c("B","C")
)

# Merge with volc
p_all<- plot_grid(
  p_volc, p_venn_cor, p_PK,
  ncol=1,
  rel_heights = c(1,0.8,0.8),
  labels = c("A",NA,"D")
)

# save
ggsave(paste0("results/figs/manuscript_fig4AD_CL_PP.pdf"),  p_all, device = cairo_pdf, width = 178, height = 2/3*265, units = "mm")

# Supplement with PK sites
###########################
p_PK_sites<- plot_grid(
  p_PK_sites_ls$BAY1895344_50nM_ATM, p_PK_sites_ls$AZD6738_50nM_ATM, p_PK_sites_ls$AZD6738_1uM_ATM,
  p_PK_sites_ls$BAY1895344_50nM_ATR, p_PK_sites_ls$AZD6738_50nM_ATR, p_PK_sites_ls$AZD6738_1uM_ATR,
  p_PK_sites_ls$BAY1895344_50nM_DNAPK, p_PK_sites_ls$AZD6738_50nM_DNAPK, p_PK_sites_ls$AZD6738_1uM_DNAPK,
  p_PK_sites_ls$BAY1895344_50nM_RPS6K, p_PK_sites_ls$AZD6738_50nM_RPS6K, p_PK_sites_ls$AZD6738_1uM_RPS6K,
  ncol=3
)
ggsave(paste0("results/figs/manuscript_fig4_S_PKsites.pdf"),  p_PK_sites, device = cairo_pdf, width = 178, height = 2/3*265, units = "mm")
ggsave(paste0("results/figs/manuscript_fig4_S_PKsites.png"),  p_PK_sites, width = 178, height = 2/3*265, units = "mm")

# Save list to pick
all_df<- paste0("PK_sites_",rep(c("BAY1895344_50nM", "AZD6738_50nM", "AZD6738_1uM"),4),"_",rep(c("ATM","ATR","DNAPK","RPS6K"),each=3))
WriteXLS::WriteXLS(all_df,"temp/table_sites.xlsx",row.names = T,SheetNames = gsub("PK_sites_","",all_df))


# Summary table
#################
BAY<- as.data.frame(res_diff_expr_CL_PP$BAY1895344_50nM)
AZD50<- as.data.frame(res_diff_expr_CL_PP$AZD6738_50nM)
AZD1000<- as.data.frame(res_diff_expr_CL_PP$AZD6738_1uM)
WriteXLS::WriteXLS(c("BAY","AZD50","AZD1000"),"results/tables/manuscript_tableS1_PP.xlsx",row.names = T, SheetNames = c("elimusertib 50nM","ceralasertib 50nM","ceralasertib 1µM"))

