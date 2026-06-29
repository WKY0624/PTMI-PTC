
library(limma)
library(GSEABase)
library(GSVA)
library(reshape2)
library(ggplot2)

expFile="05.TPM100.txt"    
riskFile="totalRisk.txt"
gmtFile="16.h.all.v2023.1.Hs.symbols.gmt"  

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
param=gsvaParam(
  data,
  geneSets,
  assay = NA_character_,
  annotation = NULL,
  minSize = 10,
  maxSize = 500,
  kcdf = c("Gaussian"),
  tau = 1,
  maxDiff = TRUE,
  absRanking = FALSE
)
gsvaResult=gsva(param,verbose = TRUE)

data=t(gsvaResult)
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)


cli=read.table("06.time_DSS.txt", header=T, sep="\t", check.names=F)
cox_result <- data.frame()

for (pathway in colnames(data)) {
  dat <- data.frame(
    score = data[, pathway],
    DSS = cli$DSS,
    Event = cli$Event
  )
  fit <- coxph(Surv(DSS, Event) ~ score, data = dat)
  tmp <- summary(fit)
  cox_result <- rbind(cox_result, data.frame(
    pathway = pathway,
    coef = tmp$coefficients[1],
    lower = tmp$conf.int[,"lower .95"],
    upper = tmp$conf.int[,"upper .95"]
  ))
}


pathway_group_info <- data.frame(
  pathway = c(
    "HALLMARK_GLYCOLYSIS", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_XENOBIOTIC_METABOLISM",
    "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_HEME_METABOLISM",
    
    "HALLMARK_PEROXISOME","HALLMARK_APICAL_JUNCTION","HALLMARK_APICAL_SURFACE", 
    
    "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT","HALLMARK_P53_PATHWAY",
    
    "HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_PROTEIN_SECRETION",
    
    "HALLMARK_DNA_REPAIR","HALLMARK_UV_RESPONSE_UP", "HALLMARK_UV_RESPONSE_DN",
    
    "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_ALPHA_RESPONSE", 
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_COAGULATION",
    
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_NOTCH_SIGNALING","HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_ANDROGEN_RESPONSE",  "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MTORC1_SIGNALING",
    
    
    
    "HALLMARK_SPERMATOGENESIS","HALLMARK_ADIPOGENESIS","HALLMARK_MYOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_PANCREAS_BETA_CELLS",'HALLMARK_ANGIOGENESIS'
    
  ),
  group = c(
    rep("Metabolism", 7),
    rep('Cellular component',3),
    rep("Proliferation", 6),
    rep("Pathway", 5),
    rep('DNA damage',3),
    rep("Immune", 8),
    rep("Signaling", 12),
    rep("Development",6)
  ))

cox_result$group <- pathway_group_info$group[match(cox_result$pathway, pathway_group_info$pathway)]


ggplot(cox_result, aes(x = coef, y = reorder(pathway, coef), color = group)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15, size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
  scale_color_manual(values = mycol) +
  xlim(-10, 40) +
  theme_bw() +
  labs(
    x = "Cox Coefficient",
    y = "Hallmarks",
    color = "Group"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10,  color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )
ggsave("GSVAcox.pdf", height = 9, width = 8)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F] 

outTab=data.frame()
for(Geneset in colnames(data)){
  for(gene in colnames(risk)){
    x=as.numeric(data[,Geneset])
    y=as.numeric(risk[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Geneset=Geneset, cor, text, pvalue))
  }
}


outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
ysort = unique(outTab$Geneset)
xface = c('plain','plain','plain','plain','plain','plain','plain','plain','plain','plain','plain','plain',"bold")


ggplot(outTab, aes(Gene, Geneset)) + 
  geom_tile(aes(fill = cor), colour = "white", size = 0.5)+ 
  scale_fill_gradient2(low = "#8FB4DC", mid = "white", high = "#e26844",  
                       breaks=c(-0.5,0,0.5),labels=c(-1,0,1),   
                       limits=c(-0.75,0.68)) +    
  geom_text(aes(label=text),col ="black",size = 3) +  
  scale_y_discrete(limits=factor(sort(ysort, decreasing = T)),position = 'left') +
  scale_x_discrete(limits=c('UBE2C','NOS3','SOD3','LPCAT1','ELOVL6','H4C9','GRIN2A','LRAT','TYMS','MGLL','GRIN2D','ATAD2','RiskScore')) +
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.title.x=element_blank(), 
        axis.ticks = element_line(linewidth = 0.3),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   colour = 'black',
                                   face = xface,
                                   size = 8),    
        axis.text.y = element_text(size = 8,
                                   colour = 'black')) +  
  ggtitle("Hallmarks GSVA correlation", ) +
  theme(plot.title = element_text(size =10, hjust = 0.5)) + 
  guides(fill = guide_colorbar(title = paste0("*** P<0.001","\n", " ** P<0.01","\n", "  * P<0.05","\n","\n","Correlation\nCoefficient"),
                               title.position = 'top',
                               title.theme = element_text(size = 8,face = "plain",colour = "black"),
                               label = T,   
                               label.theme = element_text(size = 8,face = "plain",colour = "black"),
                               raster = T,
                               frame.colour = NULL,    
                               barwidth = unit(3,"mm"),
                               barheight = unit(18,"mm"),
                               nbin = 50,   
                               ticks = T,   
                               draw.ulim = T,   
                               draw.llim = T,  
  )
  )

ggsave("GSVAcor.pdf", height = 9, width = 5.8)

