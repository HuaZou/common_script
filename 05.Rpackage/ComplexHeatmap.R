create_complexheatmap <- function(indata){
  set.seed(123)
  
  sample.annotation <- data.frame(SampleType=sample(c("groupA","groupB"),
                                                    size = ncol(indata),
                                                    replace = TRUE))
  
  topha <- HeatmapAnnotation(df = sample.annotation,
                             col = list(SampleType = c("groupA" =  "coral3",
                                                      "groupB" = "aquamarine4")),
                             height = unit(0.333, "cm"))
  
  maintitle <- "Example Heatmap"
  
  gene.annotation <- data.frame(GeneAnnot1=sample(c("Gene_groupA","Gene_groupB"),
                                                  size = nrow(indata),
                                                  replace = TRUE),
                                GeneAnnot2=sample(c("Gene_groupC","Gene_groupD"),
                                                 size = nrow(indata),
                                                 replace = TRUE))
  
  ha_row = HeatmapAnnotation(df = gene.annotation,
                                   col = list(GeneAnnot1 = c("Gene_groupA" =  "#4DAF4A",
                                                              "Gene_groupB" = "#984EA3"),
                                              GeneAnnot2 = c("Gene_groupC" =  "#FFFF33",
                                                              "Gene_groupD" = "skyblue")),
                                   which = "row",
                                   width = unit(0.666, "cm"),
                                   show_legend = TRUE)
  
                                   
  h1 <- Heatmap(indata,
                cluster_rows = T,
                cluster_columns = T,
                show_row_names = F,
                show_column_names = T,
                row_title_gp = gpar(fontsize =10),
                top_annotation = topha,
                name="Expression",
                column_title = maintitle,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(color_bar = "continuous"))
  draw(h1+ha_row,row_dend_side = "left", annotation_legend_side = "bottom")
}

library(ComplexHeatmap)

randomdata <- data.frame(matrix(rnorm(625), nrow=25))
rownames(randomdata) <- paste(rep("gene_"), 1:25, sep='')
colnames(randomdata) <- paste(rep("sample_"), 1:25, sep='')

create_complexheatmap(randomdata)