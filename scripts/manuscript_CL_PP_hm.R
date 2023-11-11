library(ggplot2)
library(cowplot)
library(DEGreport)
library(ggrepel)
library(pheatmap)
library(dendextend)
load("data/ATR_AZD_data.RData")
load("scripts/functions/manuscript_functions.RData")

# Select sites
################

# Heatmap focussing on all genes active in ATM/ATR pathways 
ATM_genes<- CP_ls$PID_ATM_PATHWAY
ATR_genes<- CP_ls$PID_ATR_PATHWAY
DNAPK_genes<- CP_ls$PID_DNA_PK_PATHWAY
MTOR_genes<- CP_ls$PID_MTOR_4PATHWAY
ATMATR_genes<- unique(c(ATR_genes,ATM_genes,DNAPK_genes,MTOR_genes))
res_proc<- normalized_counts_CL_PP[gsub("_.*","",rownames(normalized_counts_CL_PP))%in%ATMATR_genes,]

# Plot heatmap?
###############

# Use pheatmap, see tutorial at https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/#comment-8232

data<- na.omit(res_proc)
data<- log2(data/rowMeans(data[,grep("DMSO",colnames(data))],na.rm=T))

sample_names<- gsub("THYM_REL_","",colnames(data))
sample_names<- substr(sample_names,1,nchar(sample_names)-2)

# Column colors
sample_col<- data.frame(
  row.names = colnames(data),
  treatment = sample_names
)

# Get main clusters
my_hclust_gene <- hclust(dist(data), method = "ward.D2")
gene_clusters <- cutree(tree = as.dendrogram(my_hclust_gene), k = 3)

# Row colors
gene_col<- data.frame(
  row.names = rownames(data),
  ATR = c("O","ATR")[as.numeric(gsub("_.*","",rownames(data))%in%ATR_genes)+1],
  ATM = c("O","ATM")[as.numeric(gsub("_.*","",rownames(data))%in%ATM_genes)+1],
  MTOR = c("O","MTOR")[as.numeric(gsub("_.*","",rownames(data))%in%MTOR_genes)+1],
  DNAPK = c("O","DNAPK")[as.numeric(gsub("_.*","",rownames(data))%in%DNAPK_genes)+1],
  # SQ = c("O","SQ")[as.numeric(CL_PP_SQ_matrix[rownames(data),"isSQ"])+1],
  cluster = paste0("cluster","_",gene_clusters[rownames(data)])
)

my_colour = list(
  ATR = c(ATR = "black", O = "white"),
  ATM = c(ATM = "black", O = "white"),
  MTOR = c(MTOR = "black", O = "white"),
  DNAPK = c(DNAPK = "black", O = "white")
  # SQ = c(SQ = "black", O = "white")
)

# Use protein names for labels (not gene names)
rownames(data)<- gsub("PRKDC_","DNAPK_",rownames(data))
rownames(data)<- gsub("CHEK1_","CHK1_",rownames(data))
rownames(data)<- gsub("CHEK2_","CHK2_",rownames(data))

# Determine which sites to label
idx_genes_sel<- grepl("ATR_|ATM_|CHK1_|CHK2_|\\bRRM2|RPTOR_|ATRX_|AEBP2_|CASC3_|DCK_|E2F3_|EXO1_|FANCD2_|FANCI_|NBN_|RPA2_|SLBP_|TOPBP1_|UTP14A_|WRN_|DNAPK",rownames(data))
genes_sel<- rownames(data)
genes_sel[!idx_genes_sel]<- ""

rg<- 2 # Maximal color scale
p<- pheatmap(data, cutree_rows = 3, cutree_cols = 4, annotation_row = gene_col,annotation_col = sample_col, annotation_colors = my_colour, clustering_method="ward.D2", color = colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-rg, rg, length.out = 50),show_colnames = F, show_rownames = F)
p_withSites<- pheatmap(data, cutree_rows = 3, cutree_cols = 4, annotation_row = gene_col,annotation_col = sample_col, annotation_colors = my_colour, clustering_method="ward.D2", color = colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-rg, rg, length.out = 50),show_colnames = F, show_rownames = T, labels_row = genes_sel, fontsize_row = 6)

# Clustering significant?
#########################
clust_enrich_df<- data.frame(
  clust=factor(as.character(rep(1:3,4)),levels=rev(c("3","1","2"))),
  PK=factor(rep(c("ATM","ATR","DNAPK","MTOR"),each=3),levels=c("DNAPK","MTOR","ATM","ATR")),
  p=NA,
  OR=NA
)

for(i in 1:nrow(clust_enrich_df)){
  ft<- fisher.test(table(gene_clusters==clust_enrich_df$clust[i], gsub("_.*","",names(gene_clusters))%in%get(paste0(clust_enrich_df$PK[i],"_genes"))))
  clust_enrich_df[i,c("p","OR")]<- c(ft$p.value,ft$estimate)
}
clust_enrich_df$sign_label<- ""
clust_enrich_df$sign_label[clust_enrich_df$p>0.05]<- "X"
clust_enrich_df

# clust    PK            p        OR
# 1      1   ATM 5.687817e-01 0.8668715
# 2      2   ATM 1.246176e-06 0.2077150
# 3      3   ATM 3.183606e-08 5.2594439
# 4      1   ATR 1.441257e-03 0.3621472
# 5      2   ATR 3.351720e-03 2.7659325
# 6      3   ATR 3.345745e-01 1.4101749
# 7      1 DNAPK 2.268373e-01 0.4421026
# 8      2 DNAPK 1.000000e+00 0.8068475
# 9      3 DNAPK 5.827679e-02 3.1749039
# 10     1  MTOR 2.956090e-03 2.0146423
# 11     2  MTOR 3.666076e-02 1.7916422
# 12     3  MTOR 6.823130e-10 0.1058650

p_enrich<- ggplot(clust_enrich_df, aes(x=paste0(PK," pw"),y=clust,fill=log2(OR),label=sign_label)) +
  geom_point(shape = 21,size=3*abs(log10(clust_enrich_df$p)),stroke=0,colour="white") +
  scale_fill_gradient2(low = "blue", high = "red", 
                              mid = "white", midpoint = 0, space = "Lab") +
  geom_text(size=2) + 
  xlab("") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylab("cluster") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text.x = element_text(size=7),  
    axis.text.y = element_text(size=7),
    axis.title = element_text(size=7),
    legend.text = element_text(size=6),
    legend.title = element_text(size=7),
    legend.key.size = unit(3, 'mm')
  )

# Barplot mean ATR/ATM per cluster
####################################
conds<- names(res_diff_expr_CL_PP)
clust<- 1:3
logFC_df<- data.frame(
  cond = rep(conds, each=length(clust)),
  clust= rep(rep(clust,length(conds))),
  avg=NA,
  ci_l=NA,
  ci_u=NA,
  p=NA
)

for(c in conds){
    tt<- tapply(res_diff_expr_CL_PP[[c]]$logFC, gene_clusters[rownames(res_diff_expr_CL_PP[[c]])],"t.test")
    logFC_df$avg[logFC_df$cond==c]<- sapply(tt,function(x) x$estimate)
    logFC_df$ci_l[logFC_df$cond==c]<- sapply(tt,function(x) x$conf.int)[1,]
    logFC_df$ci_u[logFC_df$cond==c]<- sapply(tt,function(x) x$conf.int)[2,]
    logFC_df$p[logFC_df$cond==c]<- sapply(tt,function(x) x$p.value)
}

logFC_df$clust<- paste0("Clust ",logFC_df$clust) 
logFC_df$clust<- factor(logFC_df$clust, levels=c("Clust 1","Clust 3","Clust 2","Clust 4"))
logFC_df$cond_short<- gsub("AZD6738_","AZD (",logFC_df$cond)
logFC_df$cond_short<- gsub("BAY1895344_","BAY (",logFC_df$cond_short)
logFC_df$cond_short<- paste0(logFC_df$cond_short,")")

p_clust<- ggplot(logFC_df, aes(x=cond_short,y=avg)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_grid(clust~., scales="free", space="free") +
  geom_errorbar(aes(x=cond_short, ymin=ci_l, ymax=ci_u), width=0.4) +
  xlab("") +
  ylab("logFC") +
  theme(
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    strip.text = element_text(size = 8),
  )
p_clust

# Some numbers
###############

table(gene_clusters)
# gene_clusters
# 1   2   3 
# 201  65  63 

logFC_df
# cond   clust          avg        ci_l        ci_u            p cond_short
# 1     AZD6738_1uM Clust 1  0.066460587  0.04674164  0.08617953 2.796750e-10  AZD (1uM)
# 2     AZD6738_1uM Clust 2 -0.377830769 -0.42855481 -0.32710673 1.871016e-22  AZD (1uM)
# 3     AZD6738_1uM Clust 3  0.694936508  0.59334801  0.79652500 2.220999e-20  AZD (1uM)
# 4    AZD6738_50nM Clust 1  0.051994229  0.03370963  0.07027883 6.767439e-08 AZD (50nM)
# 5    AZD6738_50nM Clust 2 -0.206658462 -0.26025453 -0.15306239 1.076847e-10 AZD (50nM)
# 6    AZD6738_50nM Clust 3  0.009324635 -0.03772002  0.05636929 6.933086e-01 AZD (50nM)
# 7 BAY1895344_50nM Clust 1  0.042189453  0.02149858  0.06288032 8.216403e-05 BAY (50nM)
# 8 BAY1895344_50nM Clust 2 -0.366349231 -0.44171888 -0.29097958 3.314598e-14 BAY (50nM)
# 9 BAY1895344_50nM Clust 3  0.967571429  0.83708566  1.09805720 4.745367e-22 BAY (50nM)

# Save plots
##########

pdf("results/figs/manuscript_fig4C_CL_PP_hm_sites.pdf")
grid::grid.newpage()
grid::grid.draw(p_withSites$gtable)
dev.off()

ggsave("results/figs/manuscript_fig4C_CL_PP_hm_clust.pdf", p_clust, device = cairo_pdf, width = 178/4, height = 265/4, units = "mm")

# Enrichment plot
ggsave("results/figs/manuscript_fig4D_CL_PP_hm_clust_enrich.pdf", p_enrich, device = cairo_pdf, width = 178/3, height = 265/6, units = "mm")
