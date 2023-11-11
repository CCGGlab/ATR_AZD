# Explore DDR in RNA-Seq data from different NB cell lines

library(ggplot2)
library(gplots)
library(cowplot)

# Load data
load("data/ATR_AZD_data.RData")
load("CLEAN/data/ALK_all_studies.RData") # Data available in CLEAN

# Merge genesets to focus on
GS_ls<- list(
  DNA_Damage_Response=DDR_ls$DNA_Damage_Response,
  HALLMARK_G2M_CHECKPOINT=Ha_ls$HALLMARK_G2M_CHECKPOINT,
  HALLMARK_E2F_TARGETS=Ha_ls$HALLMARK_E2F_TARGETS
)
pws<- names(GS_ls)

# Compare logFC to ctrl
########################
CLs<- sort(unique(sample_info$RNA$CL))
logFC_matrix<- NULL
for(CL in CLs){
  res_tmp<- res_diff_expr$RNA[[CL]]
  for(tr in names(res_tmp)){
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
  pw=rep(pws,each=length(cond)),
  avg=NA,
  ci_l=NA,
  ci_u=NA,
  p=NA
)
logFC_df$CL<- gsub("_.*","",logFC_df$cond)
logFC_df$time<- as.numeric(gsub(".*_","",logFC_df$cond))
logFC_df$names<- sapply(logFC_df$cond, function(x) paste(unlist(strsplit(x,"_"))[-1],collapse = "_"))
logFC_df$names<- sapply(logFC_df$names, function(x) paste(rev(rev(unlist(strsplit(x,"_")))[-1]),collapse = "_"))
logFC_df$names[!is.na(logFC_df$time)]<- paste0(logFC_df$names[!is.na(logFC_df$time)]," (",logFC_df$time[!is.na(logFC_df$time)],"h)")
logFC_df$names<- factor(logFC_df$names, levels = rev(unique(logFC_df$names)))

# Get study reference
logFC_df$study<- as.character(sapply(logFC_df$cond, function(x) unique(sample_info$RNA$study[grep(x,sample_info$RNA$sampleID)])))

for(c in cond){
  for(pw in pws){
    logFC_df$avg[logFC_df$cond==c&logFC_df$pw==pw]<- mean(logFC_matrix[rownames(logFC_matrix)%in%GS_ls[[pw]],c],na.rm=T)
    tt<- t.test(logFC_matrix[rownames(logFC_matrix)%in%GS_ls[[pw]],c])
    logFC_df$ci_l[logFC_df$cond==c&logFC_df$pw==pw]<- tt$conf.int[1]
    logFC_df$ci_u[logFC_df$cond==c&logFC_df$pw==pw]<- tt$conf.int[2]
    logFC_df$p[logFC_df$cond==c&logFC_df$pw==pw]<- tt$p.value
  }
}

p<- ggplot(logFC_df, aes(x=names,y=avg)) +
  geom_bar(stat="identity", ) +
  coord_flip() +
  facet_grid(CL~pw, scales="free", space="free") +
  geom_errorbar(aes(x=names, ymin=ci_l, ymax=ci_u), width=0.4) +
  xlab("") +
  ylab("logFC") +
  theme(
    axis.text = element_text(size=4),  
    axis.title = element_text(size=7),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 6, angle=0)
  )
p

# Footer
footer<- ggdraw() + 
  draw_label(
    "Fig Sx: ATR#2, 2023",
    fontface = 'italic',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

ggsave("results/figs/manuscript_figSx_DDR_CL_suppl.pdf",p, width = 178, height = 265, unit="mm")  

