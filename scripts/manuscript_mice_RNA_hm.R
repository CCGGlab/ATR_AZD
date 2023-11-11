################################
# manuscript_mice_RNA_hm.R
################################

# Visualize responses in DDR genes upon ATRi treatment in mice

# Libraries
library(ggplot2)
library(gplots)

# Load data
load("data/ATR_AZD_data.RData")

# Create heatmap
#################

pdf("results/figs/manuscript_fig5J_mice_RNA_hm.pdf")
genes_pw<- DDR_MGI_ls$DNA_Damage_Response # 273

# Take log(x+1) of count
data<- data.matrix(na.omit(log2(normalized_counts_mice[genes_pw,]+1)))

# Rank genes based on DE results
res_proc<- res_diff_expr_mice$BAY[genes_pw,]
data<- data[order(res_proc[rownames(data),"stat"]),]

# Define Colors for sidebars
cond<- unique(sample_info_mice$condition)
cond_cols<- rainbow(length(cond))
names(cond_cols)<- cond

# Plot
heatmap.2(main = "DDR", data,Rowv = F, dendrogram = "col", col=bluered(75),ColSideColors=  cond_cols[sample_info_mice[colnames(data),"condition"]],symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="row",density.info='none', cexRow=0.5, cexCol=0.8,key.title = "z score")
legend(y=1.1, x=.75, xpd=TRUE,     
       # legend = names(cond_cols),
       legend = c("Ctrl","BAY1895344","AZD"),
       col = cond_cols[c("ALK","ALK_ATRi","ALK_ATRi_AZD")], 
       lty= 1,             
       lwd = 5,           
       cex=.5
)
dev.off()

