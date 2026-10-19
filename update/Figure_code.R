## ==== Fig 1B、C、D、E ======
if(F){
  load( "./12_immune_cell/P11s_all_newcelltype_seu.obj.Rdata")
  
  seu.obj<-P11s_all_seu.obj
  DimPlot(seu.obj,
          reduction = "tsne",
          label = FALSE,
          label.size = 4,
          group.by = "Celltype2",
          split.by = "group_tissue",
          raster=FALSE)+
    scale_color_manual(values = pal_main_celltype)
  ggsave(filename = paste0("./14_immune_Mye/02_all_cell_splitby_tissue_tsne_",t,".pdf"),width = 16,height = 5.5)
  
  
  ## 分组织来源以及病人来源进行tsne的绘制
  seu.obj<-P11s_all_seu.obj
  DimPlot(seu.obj,
          reduction = "tsne",
          label = FALSE,
          label.size = 4,
          group.by = "group_tissue",
          raster=FALSE)+
    scale_color_manual(values = mypal.tissue)
  ggsave(filename = paste0("./14_immune_Mye/02_all_cell_tissue_tsne_",t,".pdf"),width = 7,height = 5.5)
  

  DimPlot(seu.obj,
          reduction = "tsne",
          label = FALSE,
          label.size = 4,
          group.by = "Celltype2",
          raster=FALSE)+
    scale_color_manual(values = pal_main_celltype)
  ggsave(filename = paste0("./14_immune_Mye/02_all_cell_tsne_",t,".pdf"),width = 7,height = 5.5)
  
  DimPlot(seu.obj,
          reduction = "tsne",
          label = FALSE,
          label.size = 4,
          group.by = "group_patient",
          raster=FALSE)+
    scale_color_manual(values = mypal.patient)
  ggsave(filename = paste0("./14_immune_Mye/02_all_cell_patient_tsne_",t,".pdf"),width = 7,height = 5.5)
  
}

## ===== Fig 1F、G =====
if(F){
  DefaultAssay(seu.obj) <- "RNA"
  feature<-list(
    "CTC"=c("PTPRC","CD9","PPBP"),
    "Epithelial"=c("EPCAM","KRT8","CDH1"),
    "Fibroblast"=c("COL1A1","DCN"),
    "Endothelial"=c("CDH5","VWF"),
    "T cell"=c("CD3E","CD3D"),
    "NK"=c("GNLY","KLRD1"),
    "B cell"=c("CD79B","MS4A1"),
    
    "Monocyte"=c("CD14","S100A12"),
    "Macrophage"=c("FCGR3A","CD68","CD163"),
    "Neutrophil"=c("FCGR3B","FPR1"),
    "DC"=c("CD1C","FCER1A"), ## CD1C是DC细胞的marker
    "Mast"=c("TPSAB1","CPA3")
  )
  
  genes<-unlist(feature)       
  names(genes)<-gsub("[0-9]$","",names(genes))
  names(genes)<-gsub("[0-9]$","",names(genes))
  genes[!genes %in% rownames(seu.obj@assays$RNA@data)]
  
  ### 1 小提琴图
  cellorder <- levels(seu.obj$Celltype2)
  ## data to use
  count <- as.matrix(GetAssayData(object = seu.obj, slot = "data")[genes, ])
  count <- as.data.frame(t(count))
  # reshape
  count$cluster <- seu.obj$Celltype2
  df <- reshape2::melt(count, id = c("cluster"), 
                       variable.name = "gene", 
                       value.name = "expression") #将宽数据转换成长数据,id中是要保留不整合的列.这里选择cluster.
  
  width <- 0.7 * length(cellorder) 
  height <-  0.4* length(genes) 
  
  table(df$cluster)
  ## plot
  pdf(file = paste0("./14_immune_Mye/03_all_cell_marker_gene_violin_",t,".pdf"), width = 7, height = 12)
  ggplot(df, aes(x = cluster, y = expression)) +
    geom_violin(aes(fill = cluster), scale = "width") + 
    theme_bw() + 
    xlab("")+
    ylab("Expression")+
    theme(panel.grid = element_blank(), 
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), 
          axis.title = element_text(color = "black", size = 10), 
          axis.text.x = element_text(color = "black", size = 10,angle = 45, hjust = 1, vjust = 1), 
          # axis.text.y = element_text(color = "black", size = 10), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          strip.background = element_blank(), 
          strip.text.y = element_text(size = 10, angle = 0)) + 
    facet_grid(rows = vars(gene), scales = "free_y") + 
    scale_fill_manual(values = pal_main_celltype)+ # 配色
    guides(fill="none")
  dev.off()
  
  ## marker gene 的featureplot
  feature<-list(
    "CTC"=c("PPBP","PF4"),
    "Epithelial"=c("EPCAM","KRT8","CDH1","MUC1"),
    "Fibroblast"=c("FAP","COL1A1","DCN"),
    "Endothelial"=c("ENG","CDH5","VWF","CLDN5"),
    "T cell"=c("CD3E","CD3D","CD3G"),
    "NK"=c("GNLY","KLRD1","NCAM1"),
    "B cell"=c("CD79A","CD79B","MS4A1"),
    "Monocyte"=c("CD14","S100A12"),
    "Macrophage"=c("FCGR3A","CD68","CD163"),
    "Neutrophil"=c("FCGR3B","FPR1"),
    "DC"=c("CD1C","FCER1A","HLA-DQA1"), ## CD1C是DC细胞的marker
    "Mast"=c("TPSAB1","CPA3","TPSB2")
  )
  for(celltype in names(feature)){
    seu.obj<-P11s_all_seu.obj
    genes<-feature[[celltype]]
    markers<-paste0('(',paste(genes,collapse = ","),")")
    # 3. 提取表达矩阵 (使用 GetAssayData) 并计算每个细胞 (col) 的均值
    # 注：Matrix::colMeans 效率高，适合稀疏矩阵
    mean_expr <- Matrix::colMeans(GetAssayData(seu.obj)[genes, ])
    
    # 4. 将计算好的均值存入 meta.data 中
    seu.obj$GeneSet_Mean <- mean_expr
    print(FeaturePlot(seu.obj, features = "GeneSet_Mean", reduction = "tsne", raster = F))+
      ggtitle(paste0(celltype, "\n", markers))
    ggsave(filename = paste0("./14_immune_Mye/03_allcell_",celltype,"_marker_gene_featureplot_",t,".pdf"), width = 4.6, height = 4.3)
  }
}

## ===== Fig 1H、I、J、K =====
if(F){
  seu.obj<-P11s_Ep_seu.obj
  
  p<-DimPlot(seu.obj,
              reduction = "tsne",
              label = FALSE,
              label.size = 4,
              group.by = "Celltype",
              raster=FALSE)+
    scale_color_manual(values = pal_main_celltype)
  
  p<-DimPlot(seu.obj,
              reduction = "tsne",
              label = FALSE,
              label.size = 4,
              group.by = "group_tissue",
              raster=FALSE)+
    scale_color_manual(values = mypal.tissue)
  
  p<-DimPlot(seu.obj,
              reduction = "tsne",
              label = FALSE,
              label.size = 4,
              group.by = "label",
              raster=FALSE)+
    scale_color_manual(values = pal_label)
  
  p<-DimPlot(seu.obj,
              reduction = "tsne",
              label = FALSE,
              label.size = 4,
              group.by = "group_patient",
              raster=FALSE)+
    scale_color_manual(values = mypal.patient)
  
  CTC_seu.obj<-subset(seu.obj,subset=label=="CTC")
  p<-DimPlot(CTC_seu.obj,
              reduction = "tsne",
              label = FALSE,
              label.size = 4,
              group.by = "group_patient",
              raster=FALSE)+
    scale_color_manual(values = mypal.patient)
}

## ===== Fig 2A、B =====
if(F){
  seu.obj<-P11s_CTC_tumor_cells_seu.obj
  Rank<-"Rank0.05"
  cells_AUC<-cells_AUC.list[[Rank]]
  AUCell_res<-as.data.frame(getAUC(cells_AUC))
  
  cluster<-"tissue"
  df <- data.frame(group= seu.obj@meta.data[,cluster],
                   cell = colnames(seu.obj))
  ha.t <- HeatmapAnnotation(Group = seu.obj@meta.data[,cluster], 
                            Patient=seu.obj$group_patient,
                            col = list(Group=pal_CTC,
                                       Patient=mypal.patient))
  data<-AUCell_res[,df$cell]
  
  pathways<-c("Ferroptosis",
              "Apoptosis",
              "Mitotic_CD_RCD",
              "Autosis_RCD",
              "Immunogenic_cell_death",
              "Anoikis",
              "Cuproptosis",
              "REACTOME_AUTOPHAGY",
              "Necroptosis_RCD",
              "Paraptosis",
              "Pyroptosis",
              "Efferocytosis_RCD",
              "Entosis",
              "Disulfidptosis_RCD",
              "Lysosome-dependentcelldeath",
              # "Alkaliptosis",
              "Alkaliptosis_RCD",
              # "Oxeiptosis",
              "Oxeiptosis_RCD",
              "NETosis"
  )
  
  df<-geneSets.df[geneSets.df$PCD %in% pathways,]
  df$PCD<-gsub("_RCD|HALLMARK_|REACTOME_","",df$PCD)
  df$PCD<-gsub("_"," ",df$PCD)
  df$PCD<-gsub("dependentcelldeath","dependent cell death",df$PCD)
  df$PCD<-str_to_title(tolower(df$PCD))
  colnames(df)[1]<-"symbol"
  write.xlsx(df,file = "./12_Patients11/PCD_genesets.xlsx")
  
  pathways<-sigPathways[rownames(sigPathways) %in% pathways,]
  pathways<-pathways[order(pathways$logFC),]
  
  df.plot <- data[rownames(pathways),,drop=F ]
  row_group<-sigPathways[rownames(df.plot),]
  
  rownames(df.plot)<-gsub("_RCD|HALLMARK_|REACTOME_","",rownames(df.plot))
  rownames(df.plot)<-gsub("_"," ",rownames(df.plot))
  rownames(df.plot)<-gsub("dependentcelldeath","dependent cell death",rownames(df.plot))
  rownames(df.plot)<-str_to_title(tolower(rownames(df.plot)))
  
  ## 绘制PDC通路的barplot
  fig<-list()
  for(gene in rownames(df.plot)){
    plot_data<-df.plot[gene,,drop=F]
    plot_data<-as.data.frame(t(plot_data))
    colnames(plot_data)[1]<-"Expr"
    plot_data<-plot_data[rownames(seu.obj@meta.data),,drop=F]
    plot_data$group<-seu.obj$tissue
    lable.y<-plot_data %>% 
      dplyr::group_by(group) %>%
      dplyr::summarise(mean = mean(Expr,na.rm = T),)
    
    fig[[gene]]<-ggbarplot(plot_data, x = "group", y="Expr", add="mean_se",fill="group",  
                           palette = pal_CTC,
                           position=position_dodge(0.8)
                           
    ) +
      stat_compare_means(
        ref.group = "CTC",
        label = "p.signif",
        label.y = max(lable.y$mean)
      )+
      # facet_wrap(~Group)+
      xlab("") +
      ylab(gene) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title  = element_text(size = 16)) +
      guides(fill = F)
  }
  fig<-fig[!sapply(fig, is.null)]
  # w<-ceiling(sqrt(length(fig)-1))
  w<-5
  p<-do.call(grid.arrange,c(fig,ncol=w))
  ggsave(p, filename = paste0("./12_Patients11/05_scRNA_CTC_select_PCDpathways_barplot_",t,".pdf"),height = length(fig)/w*4.5,width=w*3.5)
  
  ## 绘制PDC通路的热图
  df.plot <- t(scale(t(df.plot)))
  cols = colorRamp2(seq(from = -2, to = 2, length.out = 11), 
                    c("#053061", "#2166AC", "#4393C3",
                      "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                      "#F4A582", "#D6604D","#B2182B", "#67001F"))
  
  mat_plot <- as.matrix(df.plot)
  
  x.length <- unlist(lapply(strsplit(rownames(df.plot), split = ""), length))
  width = (max(x.length) * 0.12) + 8
  height = (0.13 * length(rownames(df.plot)) + 3)
  
  right_anno = rowAnnotation(
    "log2FC" = anno_barplot(
      pathways$logFC, 
      bar_width = 0.8, 
      gp = gpar(fill = "steelblue"), # 设置柱子颜色
      border = T                # 是否显示边框
    ),
    show_annotation_name = TRUE     # 是否显示注释名称
  )
  
  
  pdf(file = paste0("./13_20260612/06_scRNA_CTC_AUCell_select_PCDpathways_heatmap_raster8_",t,".pdf"), height = height, width = width)
  print(
    Heatmap(mat_plot, 
            cluster_rows = F, 
            cluster_columns = F, 
            show_column_dend = F, 
            show_column_names = F, 
            top_annotation = ha.t,
            right_annotation = right_anno,
            column_split = group$group, 
            cluster_column_slices = F, 
            row_split = row_group$type,
            col = cols, 
            column_title_gp = gpar(fontsize = 10), 
            row_names_gp = gpar(fontsize = 10), 
            # 1. 强制开启栅格化
            use_raster = TRUE,
            
            # 2. 使用 CairoPNG 配合高 DPI（需安装 Cairo 包，无图形界面服务器极佳）
            raster_device = "CairoPNG",
            raster_device_param = list(res = 300), 
            
            # 3. 提高清晰度倍数
            raster_quality = 8,
            name = "Z score", 
            row_names_max_width = max_text_width(rownames(df.plot), gp = gpar(fontsize = 10)))
  )
  dev.off()
}

## ===== Fig 2C =====
if(F){
  seu.obj<-P11s_CTC_tumor_cells_seu.obj
  
  score2 <- fread(file ="./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/C2_CP_KEGG_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score3 <- fread(file ="./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/C2_CP_REACTOME_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score4 <- fread(file = "./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/C5_BP_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score1 <- fread(file = "./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/H_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score5 <- fread(file = "./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/PCD_18_pathways_02_Pt11s_Ep_RNA_data.tsv",header = T, data.table = F)
  
  gsvascore <- cbind(score1,score2, score3, score4,score5)
  # gsvascore<-as.data.frame(gsvascore)
  rownames(gsvascore)<-gsvascore$id
  gsvascore<-gsvascore[,-1]
  gsvascore<-gsvascore[colnames(seu.obj),]
  
  pathways<-data.frame(pathways=colnames(gsvascore))
  rownames(pathways)<-pathways$pathways
  
  pathways<-pathways[c("GO_SEQUESTERING_OF_METAL_ION",
                       "GO_CELLULAR_TRANSITION_METAL_ION_HOMEOSTASIS",
                       "GO_TRANSITION_METAL_ION_HOMEOSTASIS",
                       "GO_CELLULAR_ION_HOMEOSTASIS",
                       "GO_IRON_ION_HOMEOSTASIS",
                       "GO_IRON_ION_TRANSPORT",
                       
                       "KEGG_FATTY_ACID_METABOLISM",
                       "GO_FATTY_ACID_BETA_OXIDATION",
                       "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
                       "REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE",
                       "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
                       "GO_FATTY_ACID_CATABOLIC_PROCESS",
                       "KEGG_CITRATE_CYCLE_TCA_CYCLE",
                      
                       "GO_OXIDATION_REDUCTION_PROCESS",
                       "KEGG_GLUTATHIONE_METABOLISM"
                      
  ),,drop=F]
  pathways$`Classification`<-c(rep("Iron metabolism",6),rep("Fatty acid metabolism",7),rep("Glutathione",2))
  pathways$Classification<-factor(pathways$Classification,levels = c("Iron metabolism","Fatty acid metabolism","Glutathione"))
  
  
  df.plot<-gsvascore[,pathways$pathways]
  df.plot<-t(df.plot)
  cluster<-"tissue"
  df <- data.frame(group= seu.obj@meta.data[,cluster],
                   cell = colnames(seu.obj))
  data<-df.plot[pathways$pathways,df$cell]
  
  ha.t <- HeatmapAnnotation(Group = seu.obj@meta.data[,cluster], 
                            Patient=seu.obj$group_patient,
                            col = list(Group=pal_CTC,
                                       Patient=mypal.patient))
  
  df.plot <- t(scale(t(data)))
  cols = colorRamp2(seq(from = -3, to = 3, length.out = 11), 
                    c("#053061", "#2166AC", "#4393C3",
                      "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                      "#F4A582", "#D6604D","#B2182B", "#67001F"))
  ## width
  x.length <- unlist(lapply(strsplit(rownames(df.plot), split = ""), length))
  width = (max(x.length) * 0.12) + 5
  height = (0.13 * length(rownames(df.plot)) + 3)
  
  rownames(df.plot)<-gsub("GO_|KEGG_|REACTOME_","",rownames(df.plot))
  rownames(df.plot)<-gsub("_"," ",rownames(df.plot))
  rownames(df.plot)<- str_to_sentence(tolower(rownames(df.plot)))
  
  pdf(file = paste0("./13_20260612/05_Iron_Glutathione_pathways_heatmap_raster10_",t,".pdf"), height = height, width = width)
  print(
    Heatmap(as.matrix(df.plot), 
            cluster_rows = F, 
            cluster_columns = F, 
            show_column_dend = F, 
            show_column_names = F, 
            top_annotation = ha.t, 
            column_split = df$group, 
            cluster_column_slices = F, 
            row_split = pathways$Classification, 
            cluster_row_slices = pathways$Classification,
            col = cols, 
            column_title_gp = gpar(fontsize = 10), 
            row_names_gp = gpar(fontsize = 10), 
            name = "Z score", 
            use_raster = TRUE,
            
            # 2. 使用 CairoPNG 配合高 DPI（需安装 Cairo 包，无图形界面服务器极佳）
            raster_device = "CairoPNG",
            raster_device_param = list(res = 300), 
            
            # 3. 提高清晰度倍数
            raster_quality = 10,
            row_names_max_width = max_text_width(rownames(df.plot), gp = gpar(fontsize = 10)))
  )
  dev.off()
}

## ===== Fig 2F、G、H =====
if(F){
  seu.obj<-P11s_CTC_tumor_cells_seu.obj
  degene.use <- degenes[degenes$p_val_adj < 0.05, ]
  genes <- degene.use$gene
  
  geneset <- degene.use[rownames(degene.use) %in% gene2$Symbol,]
  geneset$avg_log2FC<-ifelse(geneset$cluster=="CTC",geneset$avg_log2FC,geneset$avg_log2FC*-1)
  # geneset <- geneset[!geneset$gene %in% c("TSC22D3","HBA1"),]
  geneset<-geneset[order(geneset$avg_log2FC,decreasing = T),]
  geneset$type<-ifelse(geneset$avg_log2FC>0,"Up","Down")
  geneset$type<-factor(geneset$type,levels = c("Up","Down"))
  
  genes <- geneset %>%
    arrange(desc(avg_log2FC)) %>%        # 按 logFC 从大到小排序
    slice(c(1:10, (n()-9):n()))    
  ## annotation
  
  # colors
  cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                    c("#053061", "#2166AC", "#4393C3",
                      "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                      "#F4A582", "#D6604D","#B2182B", "#67001F"))
  # annotation
  ha.t <- HeatmapAnnotation(Group = seu.obj$tissue, 
                            Patient=seu.obj$group_patient, 
                            col = list(Group = pal_CTC, 
                                       Patient = mypal.patient))
  
  df.plot <- as.matrix(GetAssayData(seu.obj, slot = "data")[genes$gene, ])
  df.plot[1:5, 1:5]
  df.plot <- t(scale(t(df.plot)))
  min(df.plot)
  
  width = 12
  height <- (0.2 * nrow(df.plot) + 6)
  
  ## 1 plot
  pdf(paste0("./12_Patients11/07_scRNA_CTC_paperFERROPTOSIS_DEgenes_heatmap_raster10_",t,".pdf"), height = height, width = width)
  Heatmap(df.plot, 
          cluster_rows = F, 
          show_row_dend = F, 
          cluster_columns = F, 
          cluster_column_slices = F,
          column_split = seu.obj$tissue, 
          row_split = genes$type,
          top_annotation = ha.t, 
          col = cols, 
          name = "Z score", 
          show_column_names = F,
          width = unit(8, "in"),
          height = unit(0.3 * nrow(df.plot), "in"),
          column_title_gp = gpar(fontsize = 10), 
          row_names_gp = gpar(fontsize = 10), 
          use_raster = TRUE,
          
          # 2. 使用 CairoPNG 配合高 DPI（需安装 Cairo 包，无图形界面服务器极佳）
          raster_device = "CairoPNG",
          raster_device_param = list(res = 300), 
          
          # 3. 提高清晰度倍数
          raster_quality = 10,
          row_names_max_width = max_text_width(rownames(df.plot), 
                                               gp = gpar(fontsize = 10)))
  dev.off()
  
  
  ## 绘制这些基因的featureplot
  seu.obj<-P11s_CTC_tumor_cells_seu.obj
  FeaturePlot(seu.obj, features = genes$gene, reduction = "tsne", raster = F, ncol = 5)
  ggsave(filename = paste0("./12_Patients11/07_scRNA_CTC_paperFERROPTOSIS_DEgenes_featureplot_",t,".pdf"), width = 20, height = 10.5/3*4)
  
  
  ## 绘制基因barplot
  groupBy="tissue"
  plot_data = FetchData(object = seu.obj, 
                        vars = c(genes$gene, groupBy), 
                        assy="RNA",
                        slot = 'data') %>% 
    dplyr::rename(group = as.name(groupBy)) %>% 
    tidyr::pivot_longer(cols = -group, names_to = 'Feat', values_to = 'Expr')
  
  fig<-list()
  for(gene in genes$gene){
    df.plot<-plot_data[plot_data$Feat==gene,]
    lable.y<-df.plot %>% 
      dplyr::group_by(group) %>%
      dplyr::summarise(mean = mean(Expr,na.rm = T),)
    
    
    
    fig[[gene]]<-ggbarplot(df.plot, x = "group", y="Expr", add="mean_se",fill="group",  
                           palette = pal_CTC,
                           position=position_dodge(0.8)
                           
    ) +
      stat_compare_means(
        ref.group = "CTC",
        label = "p.signif",
        label.y = 5.2
      )+
      coord_cartesian(ylim = c(0, 5.5)) +
      # facet_wrap(~Group)+
      xlab("") +
      ylab(gene) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title  = element_text(size = 16)) +
      guides(fill = F)
  }
  fig<-fig[!sapply(fig, is.null)]
  # w<-ceiling(sqrt(length(fig)-1))
  w<-5
  p<-do.call(grid.arrange,c(fig,ncol=w))
  ggsave(p, filename = paste0("./12_Patients11/07_scRNA_CTC_paperFERROPTOSIS_DEgenes_barplot_",t,".pdf"),height = length(fig)/w*4.5,width=w*3.5)
}

## ===== Fig 4A、C、D =====
if(F){
  genes <- degene.use$gene
  df.plot <- data.frame(LFC = degene.use$avg_log2FC, Pvalue = degene.use$p_val_adj, corscore = corsecor[genes, "GPX4"])
  rownames(df.plot) <- rownames(degene.use)
  
  # plot
  ## 3.1 all DEgenes
  df.plot <- df.plot[order(df.plot$corscore, decreasing = T), ]
  df.plot <- df.plot[-1, ]
  
  # reneme names
  df.plot$gsym <- row.names(df.plot)
  head(df.plot)
  # df.plot <- df.plot[df.plot$corscore > 0.3, ]
  
  # create ggplot layers and tweaks via ggplot style options
  df.plot$log10P<- -1*log10(df.plot$Pvalue+min(df.plot[df.plot$Pvalue != 0, "Pvalue"], na.rm = TRUE))
  df<-df.plot
  df$label<-df$gsym
  df$label[7:nrow(df)]<-""
  
  library(plotly)
  p<-plot_ly(df,x=~log10P,y=~LFC,z=~corscore,
             type      = "scatter3d",
             mode      = "markers",
             text      = df$label,
             hoverinfo = "text",
             color=~corscore, colors = c("blue", "white", "#FF4829","red"),
             size = ~log10P, sizes = c(10, 300)
  ) %>%
    add_markers() %>%
    layout(scene = list(
      xaxis = list(title = '-Log10(adj.p value)'),
      yaxis = list(title = 'Log2 (Fold change)'),
      zaxis = list(title = 'Correlation coefficient'))
    )
  htmlwidgets::saveWidget(p, "GPX4_correlation.html")
  
  ## 查看top gene的表达情况
  ## 绘制这些基因的featureplot
  head(df.plot)
  genes<-df.plot$gsym[1:6]
  seu.obj<-P11s_CTC_tumor_cells_seu.obj
  FeaturePlot(seu.obj, features = genes, reduction = "tsne", raster = F, ncol = 3)
  ggsave(filename = paste0("./12_Patients11/08_GPX4_corGene_top6_featureplot_",t,".pdf"), width = 12, height = 10.5/3*2)
  
  
  ## 绘制基因barplot
  groupBy="tissue"
  plot_data = FetchData(object = seu.obj, 
                        vars = c(genes, groupBy), 
                        assy="RNA",
                        slot = 'data') %>% 
    dplyr::rename(group = as.name(groupBy)) %>% 
    tidyr::pivot_longer(cols = -group, names_to = 'Feat', values_to = 'Expr')
  
  fig<-list()
  for(gene in genes){
    df.plot<-plot_data[plot_data$Feat==gene,]
    lable.y<-df.plot %>% 
      dplyr::group_by(group) %>%
      dplyr::summarise(mean = mean(Expr,na.rm = T),)
    
    
    fig[[gene]]<-ggbarplot(df.plot, x = "group", y="Expr", add="mean_se",fill="group",  
                           palette = pal_CTC,
                           position=position_dodge(0.8)
                           
    ) +
      stat_compare_means(
        ref.group = "CTC",
        label = "p.signif",
        label.y = 6.1
      )+
      # facet_wrap(~Group)+
      xlab("") +
      ylab(gene) +
      coord_cartesian(ylim = c(0, 6.5)) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title  = element_text(size = 16)) +
      guides(fill = F)
  }
  fig<-fig[!sapply(fig, is.null)]
  # w<-ceiling(sqrt(length(fig)-1))
  w<-3
  p<-do.call(grid.arrange,c(fig,ncol=w))
  ggsave(p, filename = paste0("./12_Patients11/08_GPX4_corGene_top6_barplot_",t,".pdf"),height = length(fig)/w*4.5,width=w*3.5)
}

## ===== Fig 6B =====
if(F){
  score2 <- fread(file ="./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/C2_CP_KEGG_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score3 <- fread(file ="./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/C2_CP_REACTOME_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score4 <- fread(file = "./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/C5_BP_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  score1 <- fread(file = "./SingleCellSignatureExplorer/SingleCellSignatureScorer/results/H_02_Pt11s_Ep_RNA_data.tsv", header = T, data.table = F)
  
  gsvascore <- cbind(score1,score2, score3, score4)
  gsvascore<-as.data.frame(gsvascore)
  rownames(gsvascore)<-gsvascore$id
  gsvascore<-gsvascore[,-1]
  gsvascore<-gsvascore[colnames(seu.obj),]
  gsvascore<-as.data.frame(t(gsvascore))
  tmp<-rownames(gsvascore)
  tmp2<-colnames(gsvascore)
  gsvascore <- data.frame(apply(gsvascore, 2, as.numeric))
  rownames(gsvascore)<-tmp
  colnames(gsvascore)<-tmp2
  
  
  pathways<-rownames(gsvascore)
  pathways<-pathways[grepl("ENDOCYTOSIS|Pinocytosis|Phagocytosis|Macropinocytosis",pathways,ignore.case = T)]
  pathways<-pathways[!grepl("REGULATION_OF|PHOSPHOLIPIDS|TRANSFERRIN|PRESYNAPTIC|CHOLESTEROL",pathways)]
  pathways
  
  ## 计算CTC与其它两类细胞之间的差异
  seu.obj<-P11s_CTC_tumor_cells_seu.obj
  
  cluster<-"tissue"
  df <- data.frame(group= seu.obj@meta.data[,cluster],
                   cell = colnames(seu.obj))
  
  group<-df
  group$cell.group<-ifelse(group$group=="CTC",1,0)
  cell.group<-group$cell.group
  cell.group <- as.numeric(cell.group)
  design <- cbind(sampleGroup1=1, sampleGroup2vs1=cell.group)
  
  data<-gsvascore[pathways,df$cell]
  fit <- lmFit(data, design)    #limma
  fit <- eBayes(fit)
  sigPathways <- topTable(fit, coef="sampleGroup2vs1", number=Inf, p.value=1, adjust="BH")
  sigPathways<-sigPathways[order(sigPathways$logFC,decreasing = T),]
  sigPathways$type<-ifelse(sigPathways$logFC>0,"Up","Down")
  sigPathways$type<-factor(sigPathways$type,levels = c("Up","Down"))
  
  df.plot <- data[rownames(sigPathways),,drop=F ]
  df.plot <- t(scale(t(df.plot)))
  cols = colorRamp2(seq(from = -2, to = 2, length.out = 11), 
                    c("#053061", "#2166AC", "#4393C3",
                      "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                      "#F4A582", "#D6604D","#B2182B", "#67001F"))
  ## width
  x.length <- unlist(lapply(strsplit(rownames(df.plot), split = ""), length))
  width = (max(x.length) * 0.1) + 10
  height = (0.13 * length(rownames(df.plot)) + 3)
  
  pdf(file = paste0("./12_Patients11/09_ENDOCYTOSIS_pathways_heatmap_",t,".pdf"), height = height, width = width)
  print(
    Heatmap(as.matrix(df.plot), 
            cluster_rows = F, 
            cluster_columns = F, 
            show_column_dend = F, 
            show_column_names = F, 
            top_annotation = ha.t, 
            column_split = group$group, 
            cluster_column_slices = F, 
            row_split = sigPathways$type, 
            col = cols, 
            column_title_gp = gpar(fontsize = 10), 
            row_names_gp = gpar(fontsize = 10), 
            name = "Z score", 
            row_names_max_width = max_text_width(rownames(df.plot), gp = gpar(fontsize = 10)))
  )
  dev.off()
}

## ===== Fig 7A、B =====
if(F){
  res<-GPX4up_lfcres.list[["2cellline"]]
  resOrdered <- res[order(res$log2FoldChange,decreasing = T),]  #order根据padj从小到大排序结果
  resOrdered <- as.data.frame(resOrdered)
  resOrdered <- na.omit(resOrdered)
  resOrdered<-as.data.frame(resOrdered)
  sigPathways<-resOrdered[resOrdered$padj<0.05,]
  
  df<-rbind(msigdb_hsa.pathway.list[["Hall"]],
            msigdb_hsa.pathway.list[["WiKi"]],
            msigdb_hsa.pathway.list[["Reactome"]],
            msigdb_hsa.pathway.list[["GObp"]],
            msigdb_hsa.pathway.list[["BIOCARTA"]])
 
  ## 进行差异表达基因的富集分析
  kk<-enricher(rownames(sigPathways),
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               TERM2GENE = rbind(msigdb_hsa.pathway.list[["KEGG"]],
                                 df[df$gs_name %in% c("GOBP_ACTIVATION_OF_PROTEIN_KINASE_C_ACTIVITY",
                                                      "GOBP_PROTEIN_KINASE_C_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",
                                                      "GOBP_PROTEIN_KINASE_C_SIGNALING",
                                                      "WP_AMPACTIVATED_PROTEIN_KINASE_AMPK_SIGNALING"),]),
               maxGSSize = 500
  )
  save(kk,file = "./02_data/RAP1B_DEG_pathway_enrichment.Rdata")
  
  spli<-function(x){
    strsplit(x[4],"/")[[1]][1]}
  
  df<-as.data.frame(kk)
  df<-df[grepl("signaling",df$ID,ignore.case = T),]
  df<-df[!grepl("NEUROTROPHIN|B_CELL|INSULIN|T_CELL|ADIPOCYTOKINE|HELICOBACTER|FC_EPSILON_RI",df$ID,ignore.case = T),]
  df<-df[order(df$pvalue),]
  df<-df[1:15,]
 
  df$term_count<- apply(df, 1,spli )
  df$rich_factor<-as.numeric(df$Count)/as.numeric(df$term_count) 
  df<-df[order(df$pvalue,decreasing = F),]
  hm.palette<-colorRampPalette(colors = c("#210B61", "#D6604D","#B2182B", "#67001F"))(100)
  df<-df[order(df$Count,decreasing = T),]
  df$Description<-factor(df$Description,levels=rev(df$Description))
  head(df)
  p<-ggplot(df,aes( -1*log10(pvalue),Description))+
    geom_point(aes(size=Count,color=-1*log10(pvalue)))+
    scale_size_area(max_size=8)+
    scale_colour_gradient(low = "#4300FF", high = "red")+
    # scale_colour_gradientn(colours = hm.palette)+
    labs(color="-log10(pvalue)",size="Gene count",x="-log10(pvalue)",
         y="Pathway name")+
    theme_pubclean() +
    theme(axis.title = element_text(size = 16), 
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10,color = "black"),
          axis.text.y = element_text(size = 10,color = "black"),
          legend.position = "right")+
    force_panelsizes(rows = unit(5, "in"),
                     cols = unit(3, "in"))
  
  pdf(file = paste0("./11_RAP1Boe_sequence_GPX4up/02_DEG_enrichKEGG_signaling_top15_",t,".pdf"), 
      width = 15, 
      height = 10)
  print(p)
  dev.off()
  
  
  ## 绘制横向的柱状图
  df$label<-gsub("WP_|KEGG_","",df$ID)
  df$label<-gsub("_"," ",df$label)
  df$label<-str_to_sentence(tolower(df$label))
  df$label<-factor(df$label,levels = rev(df$label))
  df$logP<- -1*log10(df$pvalue)
  
  max_logP  <- max(df$logP)
  max_Count <- max(df$Count)
  scale_factor <- max_logP / max_Count  # 计算系数
  
  # 创建用于绘图的缩放版 Count 列
  df$Count_scaled <- df$Count * scale_factor
  
  ggplot(df, aes(y = label)) +
    
    # ① 绘制横向柱状图：X轴映射为主轴 logP，颜色也映射给 logP
    geom_col(aes(x = logP, fill = logP), width = 0.65, alpha = 0.85) +
    
    # ② 绘制 Gene count 折线：使用缩放后的数值
    # 【核心防坑点 2】必须用 geom_path 且明确 group = 1，线条才能沿着纵向通路顺次连接
    geom_path(aes(x = Count_scaled, group = 1), color = "#2c3e50", size = 1) +
    
    # ③ 为折线图的每个节点添加实体圆点，方便对齐阅读
    geom_point(aes(x = Count_scaled), color = "#2c3e50", size = 2.5) +
    
    # ④ 颜色渐变美化：模仿原图从蓝紫色（低显著）到红粉色（高显著）的渐变
    scale_fill_gradient(low = "#4300FF", high = "red", name = "-log10(pvalue)") +
    
    # ⑤ 【关键步骤】设置双 X 轴：下方为主轴(logP)，上方为次轴(Gene count)
    scale_x_continuous(
      name = "-log10(pvalue)",
      # sec_axis 内传入逆向转换公式 (点除以系数)，将其还原为真实的 Count 标签刻度值
      sec.axis = sec_axis(~ . / scale_factor, name = "Gene count",
                          breaks = seq(0, max_Count, by = 3)) 
    ) +
    
    # ⑥ 沿用高水平期刊要求的纯净主题设置
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
      
      # 彻底去除背景网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 在图形外侧加上黑色实线边框
      panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
      
      # 让所有坐标轴（主轴、次轴、Y轴）的刻度线重新显现
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.12, "cm"),
      
      # 标签文字优化
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 12),
      
      axis.title.y = element_blank(),
      axis.title.x = element_text(color = "black", size = 14),
      # 图例调整到图表右侧，保持简洁
      legend.position = "right"
    )
  ggsave(filename = paste0("./13_20260612/02_DEG_enrichKEGG_signaling_top15_barplot_",t,".pdf"),width = 9, height = 7)
  
  
  ## 绘制GSVA score的气泡图
  sample_info<-GPX4up_sample_info
  sample_info$treatment<-factor(sample_info$treatment,levels = c("CPPT","RAP1B"))
  
  df<-df.plot[,sample_info$ID]
  design <- model.matrix(~ cellline + treatment, data = sample_info)
  print(colnames(design))
  fit <- lmFit(df, design)
  fit <- eBayes(fit)
  res_table <- topTable(fit, coef = "treatmentRAP1B", number = Inf, adjust.method = "BH")
  
  # 【A. 宽表格转长表格，并融入样本信息】
  gsva_long <- as.data.frame(df) %>%
    mutate(Pathway = rownames(.)) %>%
    pivot_longer(cols = -Pathway, names_to = "ID", values_to = "GSVA_score") %>%
    inner_join(sample_info, by = "ID")
  
  # 【B. 执行两步法统计转换】
  plot_data <- gsva_long %>%
    # 步骤①：必须先在每个 Cellline 内部，对所有原始样本（4个样本）统一计算 Z-score
    group_by(Pathway, cellline) %>%
    mutate(Z_score_raw = as.vector(scale(GSVA_score))) %>%
    ungroup() %>%
    
    # 步骤②：再按照处理组进行合并，求取两个重复样本 Z-score 的平均值
    group_by(Pathway, cellline, treatment) %>%
    summarise(
      Mean_Z_score = mean(Z_score_raw), 
      .groups = 'drop'
    )
  
  limma_logfc <- data.frame(
    Pathway = rownames(res_table),
    logFC = res_table$logFC
  )
  
  final_data <- plot_data %>%
    inner_join(limma_logfc, by = "Pathway")
  
  final_data$Pathway<-gsub("KEGG_|WP_","",final_data$Pathway)
  final_data$Pathway_clean <- str_to_sentence(tolower(gsub("_", " ", final_data$Pathway)))
  
  pathway_order <- str_to_sentence(tolower(gsub("_", " ", gsub("KEGG_|WP_","",rownames(res_table)))))
  final_data$Pathway_clean <- factor(final_data$Pathway_clean, levels = rev(pathway_order))
  
  final_data$treatment <- factor(final_data$treatment, levels = c("CPPT", "RAP1B"))
  
  
  ggplot(final_data, aes(x = treatment, y = Pathway_clean)) +
    
    # 【解决方案】：用 abs(logFC) 映射大小 (Size)，代表变化的剧烈程度
    # 用 Cellline 内部 scale 后的 Z_score 映射颜色 (Color)，代表表达的高低
    geom_point(aes(color = Mean_Z_score, size = abs(logFC)), alpha = 0.9) +
    
    # 按 Cellline 进行分面绘制 (横向并排展示不同的细胞系)
    facet_grid(~ cellline, scales = "free_x", space = "free") +
    
    # 颜色映射：经典的学术蓝白红/黄红渐变（符合原图 ae09d4af28ecbf40b655c183788d5ed9.png 风格）
    scale_color_gradientn(
      colors = c("black", "blue", "yellow", "red"),
      name = " mean Z-score"
    ) +
    
    # 气泡大小范围微调，并规范图例名称
    scale_size_continuous(
      range = c(2, 8), 
      name = "|logFC|"
    ) +
    
    # 深度视觉美化
    theme_minimal(base_size = 13) +
    labs(x = NULL, y = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      
      # 彻底去除网格线，使气泡背景干净
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 每一个分面四周加上标准的黑色实线外框
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      
      # 显现刻度线
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.12, "cm"),
      
      # 字体样式
      axis.text.x = element_text(color = "black",  size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      
      # 分面顶部条形标签 (MIApaca2 和 Su86.86 的通栏抬头) 美化
      strip.background = element_blank(),
      strip.text = element_text(color = "black", size = 14),
      legend.position = "right"
    )
  ggsave(filename = paste0("./13_20260612/02_DEG_enrichKEGG_signaling_gsva_point_",t,".pdf"),width = 8,height = 6)
}

## ===== Fig 7D、E =====
if(F){
  x<-list("MAPK downstream TF"=MAPK_downstreamTF$MAPK_downstream_TF,
          "GPX4 TF"=GPX4_TF$TF,
          "RAP1B_OE DEG"=rownames(DEGs))
  
  library(ggvenn)
  ggvenn(
    x, 
    fill_color = pal_nejm()(4),
    stroke_size = 0.5, 
    set_name_size = 4
  )
  ggsave(filename = paste0("./11_RAP1Boe_sequence_GPX4up/04_GPX4TF_ven_",t,".pdf"),width = 6,height = 6)
  
  
  ## 查看筛选后的转录因子的表达水平
  df<-GPX4_TF[GPX4_TF$TF %in% MAPK_downstreamTF$MAPK_downstream_TF,]
  df<-df[df$TF %in% rownames(DEGs),]
  df$TF
  
  fig<-list()
  for(cellline in unique(GPX4up_sample_info$cellline)){
    group<-GPX4up_sample_info[GPX4up_sample_info$cellline==cellline,]
    df<-GPX4up_FC_TPM[c("CEBPA","FOS"),group$ID]
    df.plot <- t(scale(t(df)))
    
    ha.t<-HeatmapAnnotation(Cellline = group$cellline, 
                            Treatment=group$treatment, 
                            col = list(Cellline = pal_cellline, 
                                       Treatment = pal_RAP1B),
                            show_annotation_name = F)
    
    fig[[cellline]]<-Heatmap(as.matrix(df.plot), 
                             cluster_rows = F, 
                             cluster_columns = F, 
                             show_column_dend = F, 
                             show_column_names = F, 
                             top_annotation = ha.t,
                             column_split = group$treatment,
                             cluster_column_slices = F, 
                             # row_split = factor(gene.use$celltype, levels = unique(gene.use$celltype)), 
                             col = cols, 
                             width = ncol(df.plot)*unit(6, "mm"), 
                             height = nrow(df.plot)*unit(4, "mm"),
                             column_title_gp = gpar(fontsize = 10), 
                             row_names_gp = gpar(fontsize = 10), 
                             name = "Z score", 
                             row_names_max_width = max_text_width(rownames(df.plot), gp = gpar(fontsize = 10)))
    
    
    
  }
  
  pdf(file = paste0("./12_Patients11/06_RAP1b_GPX4_TF_exp_heatmap_",t, ".pdf"), height = 4, width = 5)
  fig[["MIApaca2"]]+fig[["Su86.86"]]
  dev.off()
}