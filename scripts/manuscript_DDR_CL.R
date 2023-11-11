# Explore DDR in RNA-Seq  data from different NB cell lines treated with BAY or lorlatinib

library(ggplot2)
library(gplots)

# Load data
load("data/ATR_AZD_data.RData")

# How many DDR genes?
length(DDR_ls$DNA_Damage_Response) # 276

# Compare logFC to ctrl
########################
CLs<- sort(unique(sample_info_CL$CL))
logFC_matrix<- NULL
for(CL in CLs){
  res_tmp<- res_diff_expr_CL[[CL]]
  for(tr in c("Lorlatinib","BAY1895344")){
    if(!tr%in%names(res_tmp)) next
    logFC_tmp<- sapply(res_tmp[[tr]], function(x) x$log2FoldChange)
    colnames(logFC_tmp)<- paste0(CL,"_",tr,"_",colnames(logFC_tmp)) 
    rownames(logFC_tmp)<- res_tmp[[tr]][[1]]$HGNC
    if(!is.null(logFC_matrix)&&nrow(logFC_tmp)!=nrow(logFC_matrix)) next
    logFC_matrix<- cbind(logFC_matrix, logFC_tmp)
  }
}

cond<- colnames(logFC_matrix)

logFC_df<- data.frame(
  cond=cond,
  pw=rep("DNA_Damage_Response",each=length(cond)),
  avg=NA,
  ci_l=NA,
  ci_u=NA,
  p=NA
)
logFC_df$CL<- gsub("_.*","",logFC_df$cond)
logFC_df$time<- as.numeric(gsub(".*_","",logFC_df$cond))
logFC_df$treatment<- "lorlatinib"
logFC_df$treatment[grep("BAY",logFC_df$cond)]<- "elimusertib"
logFC_df$names<- paste0(logFC_df$CL," (",logFC_df$time,"h)")
logFC_df$names<- gsub("CLB","CLB-",logFC_df$names)

logFC_df$names<- factor(logFC_df$names, levels = rev(unique(logFC_df$names)))

for(c in cond){
  pw<- "DNA_Damage_Response"
  logFC_df$avg[logFC_df$cond==c&logFC_df$pw==pw]<- mean(logFC_matrix[rownames(logFC_matrix)%in%DDR_ls[[pw]],c],na.rm=T)
  tt<- t.test(logFC_matrix[rownames(logFC_matrix)%in%DDR_ls[[pw]],c])
  logFC_df$ci_l[logFC_df$cond==c&logFC_df$pw==pw]<- tt$conf.int[1]
  logFC_df$ci_u[logFC_df$cond==c&logFC_df$pw==pw]<- tt$conf.int[2]
  logFC_df$p[logFC_df$cond==c&logFC_df$pw==pw]<- tt$p.value
}

logFC_df
#           cond                  pw           avg        ci_l         ci_u            p     CL time   treatment         names
# 1   CLBBAR_Lorlatinib_1 DNA_Damage_Response  2.783231e-03 -0.03297107  0.038537529 8.783037e-01 CLBBAR    1  lorlatinib  CLB-BAR (1h)
# 2  CLBBAR_Lorlatinib_24 DNA_Damage_Response -9.428589e-02 -0.17066568 -0.017906110 1.573926e-02 CLBBAR   24  lorlatinib CLB-BAR (24h)
# 3  CLBBAR_BAY1895344_24 DNA_Damage_Response -4.572690e-02 -0.08728652 -0.004167276 3.117455e-02 CLBBAR   24 elimusertib CLB-BAR (24h)
# 4  CLBBAR_BAY1895344_48 DNA_Damage_Response -1.905815e-01 -0.24498277 -0.136180258 3.943451e-11 CLBBAR   48 elimusertib CLB-BAR (48h)
# 5    CLBGE_Lorlatinib_1 DNA_Damage_Response  8.022605e-03 -0.02943210  0.045477312 6.735605e-01  CLBGE    1  lorlatinib   CLB-GE (1h)
# 6   CLBGE_Lorlatinib_24 DNA_Damage_Response -1.511631e-01 -0.19581565 -0.106510629 1.515231e-10  CLBGE   24  lorlatinib  CLB-GE (24h)
# 7   CLBGE_BAY1895344_24 DNA_Damage_Response -2.409988e-01 -0.31151592 -0.170481682 1.086367e-10  CLBGE   24 elimusertib  CLB-GE (24h)
# 8   CLBGE_BAY1895344_48 DNA_Damage_Response -5.083479e-01 -0.61034173 -0.406354143 1.549634e-19  CLBGE   48 elimusertib  CLB-GE (48h)
# 9    IMR32_Lorlatinib_1 DNA_Damage_Response  7.320931e-03 -0.01829574  0.032937601 5.741027e-01  IMR32    1  lorlatinib    IMR32 (1h)
# 10   IMR32_Lorlatinib_6 DNA_Damage_Response -3.215852e-02 -0.07058893  0.006271898 1.006131e-01  IMR32    6  lorlatinib    IMR32 (6h)
# 11  IMR32_Lorlatinib_24 DNA_Damage_Response -7.439986e-05 -0.02127035  0.021121547 9.944907e-01  IMR32   24  lorlatinib   IMR32 (24h)
# 12     NB1_Lorlatinib_1 DNA_Damage_Response  9.241242e-02  0.05009828  0.134726560 2.408290e-05    NB1    1  lorlatinib      NB1 (1h)
# 13     NB1_Lorlatinib_6 DNA_Damage_Response -5.072981e-01 -0.60176087 -0.412835244 5.360856e-22    NB1    6  lorlatinib      NB1 (6h)
# 14    NB1_Lorlatinib_24 DNA_Damage_Response -6.134614e-01 -0.72638873 -0.500534154 2.147418e-22    NB1   24  lorlatinib     NB1 (24h)
# 15   SKNAS_Lorlatinib_1 DNA_Damage_Response -1.852418e-02 -0.04833487  0.011286521 2.222328e-01  SKNAS    1  lorlatinib    SKNAS (1h)
# 16   SKNAS_Lorlatinib_6 DNA_Damage_Response -2.904908e-02 -0.09066421  0.032566056 3.540945e-01  SKNAS    6  lorlatinib    SKNAS (6h)
# 17  SKNAS_Lorlatinib_24 DNA_Damage_Response -2.790232e-02 -0.06157379  0.005769145 1.039525e-01  SKNAS   24  lorlatinib   SKNAS (24h)

p<- ggplot(logFC_df, aes(x=names,y=avg)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_grid(treatment~., scales="free", space="free") +
  geom_errorbar(aes(x=names, ymin=ci_l, ymax=ci_u), width=0.4) +
  xlab("") +
  ylab("log2FC 276 DDR genes") +
  theme(
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    strip.text = element_text(size = 8)
  )
p
ggsave("results/figs/manuscript_fig2E_DDR_CL.pdf",p, width = 178/2, height = 265/3, unit="mm")  

