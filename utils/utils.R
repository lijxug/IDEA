#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("methods")) # Rscript for CMD doesn't load this automatically

# ---- loading dependencies ----
# utils : for tmp functions that can be used multiple time but not generalized enough
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("reshape2"))
# suppressPackageStartupMessages(library("ggfortify"))
# suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("assertthat"))

# ---- predefines ----
# bioType groups
# GENE.INFO = read_rds(WhereAmI("01.DATA/01.ASinLung/01.raw/detailed_GENE.INFO.RDS"))

bioType = list(
  protein_coding = 
    unlist(strsplit(
      "IG_C_gene, IG_D_gene, IG_J_gene, IG_LV_gene, IG_M_gene, IG_V_gene, IG_Z_gene, nonsense_mediated_decay, nontranslating_CDS, non_stop_decay, polymorphic_pseudogene, protein_coding, TR_C_gene, TR_D_gene, TR_gene, TR_J_gene, TR_V_gene",
      split = ",\\s*")),
  pseudogenes = 
    unlist(strsplit(
      "disrupted_domain, IG_C_pseudogene, IG_J_pseudogene, IG_pseudogene, IG_V_pseudogene, processed_pseudogene, pseudogene, transcribed_processed_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, translated_unprocessed_pseudogene, TR_J_pseudogene, TR_V_pseudogene, unitary_pseudogene, unprocessed_pseudogene, transcribed_unitary_pseudogene", 
      split = ",\\s*")),
  long_nc = 
    unlist(strsplit(
      "3prime_overlapping_ncrna, ambiguous_orf, antisense_RNA, antisense, lincRNA, ncrna_host, non_coding, processed_transcript, retained_intron, sense_intronic, sense_overlapping, bidirectional_promoter_lncRNA",
      split = ",\\s*")),
  short_nc = 
    unlist(strsplit(
      "miRNA, miRNA_pseudogene, misc_RNA, misc_RNA_pseudogene, Mt_rRNA, Mt_tRNA, Mt_tRNA_pseudogene, ncRNA, pre_miRNA, RNase_MRP_RNA, RNase_P_RNA, rRNA, rRNA_pseudogene, scRNA_pseudogene, snlRNA, snoRNA, snoRNA_pseudogene, snRNA, snRNA_pseudogene, SRP_RNA, tmRNA, tRNA, tRNA_pseudogene"
      , split = ",\\s*"))
  # others = 
  #   c("TEC")
)
bioType_tbl = melt(bioType) %>% as_tibble()
colnames(bioType_tbl) = c("transType", "bioType")
bioType_tbl = bioType_tbl %>% dplyr::mutate_if(is.factor, as.character)

type2bio = function(x){
  y = bioType_tbl$bioType[ match(x, bioType_tbl$transType) ]
  y[is.na(y)] = "others"
  y
}

trans2type = function(x, gene_info = GENE.INFO){
  # transform transcript names into transcript types
  t_idx = match(x, gene_info$`Transcript name`)
  y = gene_info$`Transcript type`[t_idx]
  y
}

trans2length = function(x, gene_info = GENE.INFO){
  # transform transcript names into transcript types
  t_idx = match(x, gene_info$`Transcript name`)
  y = gene_info$`Transcript length (including UTRs and CDS)`[t_idx]
  y
}

gene2type = function(x, gene_info = GENE.INFO){
  # transform transcript names into transcript types
  t_idx = match(x, gene_info$`Gene name`)
  y = gene_info$`Transcript type`[t_idx]
  y
}
trans2bio = function(x, ...) {
  # a wrapper for trans2type, this function is a wrapper which is about to deprecated
  y = trans2type(x, ...)
  y
}

transName2transid = function(tran_names, gene_info = GENE.INFO){
  t_idx = match(tran_names, gene_info$`Transcript name`)
  y = gene_info$`Transcript stable ID`[t_idx]
  y
}

trans2gene = function(trans){
  genes = stringr::str_split(trans, "-\\d{3}|_MIX|_NO_EXPR|-MIX|-NO_EXPR$", simplify = T)[, 1]
  # tmp_idx = match(trans, GENE.INFO$`Transcript name`)
  # invisible(assert_that(all(GENE.INFO$`Transcript name`[tmp_idx] == trans)))
  # genes = GENE.INFO$`Gene name`[tmp_idx]
}

cellid2cluster = function(x, cell.info = pCells, varname.cell.id = "title", varname.cluster = NULL){
  cluster.info = cell.info[,grep("[C|c]luster", names(cell.info))] %>% unlist()
  names(cluster.info) = eval(parse(text = sprintf("cell.info$`%s`", varname.cell.id)))
  y = cluster.info[x]
  return(y)
}

cellid2tissue = function(x, cell.info = pCells, varname.cell.id = "title", varname.cluster = NULL){
  cluster.info = cell.info[,grep("[T|t]issue", names(cell.info))] %>% unlist()
  names(cluster.info) = eval(parse(text = sprintf("cell.info$`%s`", varname.cell.id)))
  y = cluster.info[x]
  return(y)
}

tissues = c("N","P","T")
CDs = c("CD4","CD8")

# ---- FUNCTION ----
# ---- Extraction ----
extract_GENE.INFO = function(trans_names, GENE.INFO){
  require(dplyr)
  trans_names = trans_names %>% unique()
  ch_idx = match(trans_names, GENE.INFO$`Transcript name`)
  assert_that(all(GENE.INFO$`Transcript name`[ch_idx] == trans_names, na.rm = T))
  transcript_types = GENE.INFO$`Transcript type`[ch_idx]
  transcript_bio = type2bio(transcript_types)
  gene_names = GENE.INFO$`Gene name`[ch_idx]
  gene_types = GENE.INFO$`Gene type`[ch_idx]
  gene_bios = type2bio(gene_types)
  
  info_tbl = tibble(transcript_name = trans_names,
                    trans_type = transcript_types,
                    trans_bio = transcript_bio,
                    gene_name = gene_names,
                    gene_type = gene_types,
                    gene_bio = gene_bios)
  invisible(info_tbl)
}
# originally used in expl02_specificAS.R
extractTissue = function(Mat, spare = 6, metaData, tissue){
  require(dplyr)
  chosen_cells = (metaData %>% filter(`characteristics: tissueType` == tissue))$title
  chosen_cells = chosen_cells[chosen_cells %in% colnames(Mat)]
  tMat = bind_cols(Mat[, 1:spare], Mat[,chosen_cells])
  invisible(tMat)
}

# originally used in expl02_specificAS.R
extractCD = function(Mat, spare = 6, metaData, CD){
  require(dplyr)
  CDType = gsub(".*(CD[4|8]).+", "\\1", metaData$`characteristics: majorCluster`, perl = T)
  metaData = metaData %>% mutate(`characteristics: CDType` = CDType)
  chosen_cells = (metaData %>% filter(`characteristics: CDType` == CD))$title
  chosen_cells = chosen_cells[chosen_cells %in% colnames(Mat)]
  tMat = bind_cols(Mat[, 1:spare], Mat[,chosen_cells])
  invisible(tMat)
}

# Plot ----


# easy plot PCA
# PCbiplot <- function(PC, data, colour, ...) {
#   require(ggplot2)
#   require(ggfortify)
#   # PC being a prcomp object
#   args = list(...)
#   args$loadings = ifelse(is.null(args$loadings), T, args$loadings)
#   args$loadings.label = ifelse(is.null(args$loadings.label), T, args$loadings.label)
#   args$loadings.label.size = ifelse(is.null(args$loadings.label.size), 3, args$loadings.label.size)
#   args$loadings.colour = ifelse(is.null(args$loadings.colour), "purple", args$loadings.colour)
#   args$draw.now = ifelse(is.null(args$draw.now), T, args$draw.now)
#   
#   prop_var <- PC$sdev^2 / sum(PC$sdev^2)
#   
#   plot <-
#     autoplot(
#       PC, data = data, colour = colour,
#       loadings = args$loadings, loadings.label = args$loadings.label,
#       loadings.label.size = args$loadings.label.size,
#       loadings.colour = args$loadings.colour
#     ) + 
#     labs(
#       x = sprintf("PC1 : %.2f %% var explained", prop_var[1] * 100),
#       y = sprintf("PC2 : %.2f %% var explained", prop_var[2] * 100)) +
#     theme_gray(base_size = 20)
#   
#   if(args$draw.now) print(plot)
#   return(invisible(plot))
# }

#' ... could contain: pat CD tissue 
plotEnrichment = function(gene_name, counts_matrix, ...){
  require("pheatmap")
  
  chosen_gene = gene_name
  args = list(...)
  pat = ifelse(is.null(args$patientID), "ALL", args$patientID)
  CD = ifelse(is.null(args$CD), "CD4", args$CD)
  tissue = ifelse(is.null(args$tissue), "all", args$tissue)
  returnCountData = ifelse(is.null(args$returnCountData), F, args$returnCountData)
  printCounts = ifelse(is.null(args$printCounts), F, args$printCounts)
  TOTAL_counts = counts_matrix
  
  # tissue = "T"
  # group = "tissueType"
  group = ifelse(is.null(args$group), "MCtype", args$group)
  drawnow = ifelse(is.null(args$drawnow), T, args$drawnow)
  
  t_countsMat = eval(parse(text = paste0("TOTAL_counts$", chosen_gene)))
  calc_dt = transCountsMat(
    t_countsMat,
    CD = CD,
    pat = pat,
    group = group,
    tissue = tissue
  )
  # calc_dt = transCountsMat(t_countsMat, CD = "CD4")
  if(printCounts)
    print(calc_dt)
  chi_res = chisq.test(calc_dt + 1e-10)
  plot_dt = chi_res$residuals
  # plot_dt = (chi_res$observed - chi_res$expected) / chi_res$expected
  # plot_dt = chi_res$observed / chi_res$expected
  colors <- c(min(plot_dt), seq(-4, 4, by = 0.1), max(plot_dt))
  my_palette = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(colors))
  p = pheatmap(
    plot_dt,
    main = paste("Residuals:", group, "vs.", chosen_gene, "transcripts,",
                 "Patient:", pat, CD, tissue),
    display_numbers = T,
    cluster_rows = F,
    cluster_cols = F,
    fontsize = 12,
    fontsize_number = 14,
    breaks = colors,
    color = my_palette,
    silent = T
  )
  if(drawnow)
    plot(p$gtable)
  
  if(returnCountData){
    invisible(calc_dt)  
  }else{
    invisible(p)
  }
}


plotUsage = function(gene_name, ...){
  suppressPackageStartupMessages("dplyr")
  suppressPackageStartupMessages("ggplot2")
  usageMat = usageMat_all_trim
  
  args = list(...)
  pat = ifelse(is.null(args$patientID), "ALL", args$patientID)
  CD = ifelse(is.null(args$CD), "CD4", args$CD)
  tissue = ifelse(is.null(args$tissue), "all", args$tissue)
  # tissue = "T"
  # group = "tissueType"
  group = ifelse(is.null(args$pat), "MCtype", args$group)
  drawnow = ifelse(is.null(args$drawnow), T, args$drawnow)
  
  t_usageMat = dplyr::slice(usageMat, grep(paste0("\\b", chosen_gene, "\\b"), usageMat$`Gene name`))
  t_usageMat2 = as.data.frame(t(t_usageMat[, -1:-6]))
  colnames(t_usageMat2) = t_usageMat$`Transcript name`
  assert_that(all(metaData_trim$title == rownames(t_usageMat2)))
  CDType = metaData_trim$`characteristics: majorCluster`
  CDType[grep("CD4", CDType)] = "CD4"
  CDType[grep("CD8", CDType)] = "CD8"
  t_usageMat3 = bind_cols(
    cellID = rownames(t_usageMat2),
    CDType = CDType,
    tissueType = metaData_trim$`characteristics: tissueType`,
    patientID = metaData_trim$`characteristics: patient`,
    MCtype = metaData_trim$`characteristics: majorCluster`,
    t_usageMat2
  )
  t_usageMat4 = t_usageMat3
  if (pat != "ALL")
    t_usageMat4 = t_usageMat4 %>% filter(patientID == pat)
  t_usageMat4 = t_usageMat4 %>% filter(CDType == CD)
  if (tissue != "all")
    t_usageMat4 = t_usageMat4 %>% filter(tissueType == tissue)
  
  t_usageMat5 = t_usageMat4
  is_mix = apply(t_usageMat5[, -1:-5], 1, function(x)
    ! sum(x > 0.7))
  t_usageMat6 = t_usageMat5[is_mix, ]
  plot_tbl = melt(t_usageMat6)
  p = ggplot(plot_tbl, aes_string(color = "variable", y = "value", x = group)) +
    geom_boxplot() +
    # geom_violin() +
    labs(y = "Usage",
         title = paste("Usage:", pat, CD, tissue, "by", group)) +
    facet_grid(variable ~ .) +
    geom_hline(yintercept = 0.7, linetype = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(drawnow)
    plot(p)
  invisible(p)
  
}

plotTPM = function(gene_name , ...){
  
  suppressPackageStartupMessages("dplyr")
  suppressPackageStartupMessages("ggplot2")
  
  args = list(...)
  pat = ifelse(is.null(args$patientID), "ALL", args$patientID)
  CD = ifelse(is.null(args$CD), "CD4", args$CD)
  tissue = ifelse(is.null(args$tissue), "all", args$tissue)
  # tissue = "T"
  # group = "tissueType"
  group = ifelse(is.null(args$pat), "MCtype", args$group)
  drawnow = ifelse(is.null(args$drawnow), T, args$drawnow)
  
  t_exprMat = dplyr::slice(exprMat, grep(paste0("\\b",chosen_gene,"\\b"), exprMat$`Gene name`))
  t_exprMat2 = as.data.frame(t(t_exprMat[,-1:-6]))
  colnames(t_exprMat2) = t_exprMat$`Transcript name`
  assert_that(all(metaData_trim$title == rownames(t_exprMat2)))
  CDType = metaData_trim$`characteristics: majorCluster`
  CDType[grep("CD4", CDType)] = "CD4"
  CDType[grep("CD8", CDType)] = "CD8"
  t_exprMat3 = bind_cols(cellID = rownames(t_exprMat2), 
                         CDType = CDType,
                         tissueType = metaData_trim$`characteristics: tissueType`, 
                         patientID = metaData_trim$`characteristics: patient`,
                         MCtype = metaData_trim$`characteristics: majorCluster`,
                         t_exprMat2)
  t_exprMat4 = t_exprMat3
  if(pat != "ALL")
    t_exprMat4 = t_exprMat4 %>% filter(patientID == pat)
  t_exprMat4 = t_exprMat4 %>% filter(CDType == CD)
  if(tissue != "all")
    t_exprMat4 = t_exprMat4 %>% filter(tissueType == tissue)
  is_expr = as.logical(apply(t_exprMat4[,-1:-5],1, function(x) sum(x) > 0))
  t_exprMat5 = t_exprMat4[is_expr,]
  plot_tbl = melt(t_exprMat5)
  plot_tbl$value = log2(plot_tbl$value + 1)
  # plot_tbl$value[plot_tbl$value == 0] = NA
  p = ggplot(plot_tbl, aes_string(color = "variable", y = "value", x = group)) + 
    geom_violin() +
    # geom_violin()+
    labs(y = "log2(tpm+1)", title = paste("TPM:",pat,CD,tissue,"by",group)) +
    geom_hline(yintercept = 1, linetype = 3) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    facet_grid(variable~.)
  p
  
  if(drawnow)
    plot(p)
  invisible(p)
}




plot_bioTypes = function(transcript_names, GENE.INFO, dttype = "bio", ttl = "", mode = "transcript"){
  require(dplyr)
  require(ggplot2)
  require(gridExtra)
  # library(waffle)
  
  transcript_names = transcript_names %>% unique()
  ch_idx = match(transcript_names, GENE.INFO[[grep("transcript[_ \\.]*name", colnames(GENE.INFO), ignore.case = T)[1]]])
  assert_that(all(GENE.INFO[[grep("transcript[_ \\.]*name", colnames(GENE.INFO), ignore.case = T)[1]]][ch_idx] == transcript_names, na.rm = T))
  transcript_types = GENE.INFO[[grep("transcript[_ \\.]*type", colnames(GENE.INFO), ignore.case = T)[1]]][ch_idx]
  transcript_bio = type2bio(transcript_types)
  gene_names = GENE.INFO[[grep("gene[_ \\.]*name", colnames(GENE.INFO), ignore.case = T)[1]]][ch_idx] %>% unique()
  ch_gene_idx = match(gene_names, GENE.INFO[[grep("gene[_ \\.]*name", colnames(GENE.INFO), ignore.case = T)[1]]])
  gene_types = GENE.INFO[[grep("gene[_ \\.]*type", colnames(GENE.INFO), ignore.case = T)[1]]][ch_gene_idx]
  gene_bios = type2bio(gene_types)
  
  info_tbl = tibble(transcript_name = transcript_names,
                    trans_type = transcript_types,
                    bio_type = transcript_bio)
  gene_info_tbl = tibble(gene_name = gene_names,
                         gene_type = gene_types,
                         gene_bio = gene_bios)
  
  bio_type_table = switch(dttype,
                          "bio" = table(info_tbl$bio_type),
                          "trans" = table(info_tbl$trans_type)
  )
  
  blank_theme <- theme_minimal(base_size = 20)+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(face="bold")
    )
  bio_type_tbl = as_tibble(bio_type_table) %>% dplyr::mutate(ratio = round(n / sum(n) * 100, digits = 2))
  gene_bio_tbl = table(gene_info_tbl$gene_bio) %>% as_tibble() %>% dplyr::mutate(ratio = round(n/sum(n) *100, digits = 2))
  p.pie = 
    ggplot(bio_type_tbl, aes(x = factor(1), fill = Var1, y = n )) + 
    geom_col(width = 1) + 
    blank_theme + theme(axis.text.x=element_blank(), axis.text.y = element_blank()) + 
    labs(fill = "Transcript biotypes", title = ttl) + 
    scale_fill_discrete(labels = sprintf("%s: %s%%", bio_type_tbl$Var1, bio_type_tbl$ratio)) +
    coord_polar("y")
  
  p.gene.pie = 
    ggplot(gene_bio_tbl, aes(x = factor(1), fill = Var1, y = n )) + 
    geom_col(width = 1) + 
    blank_theme + theme(axis.text.x=element_blank(), axis.text.y = element_blank()) + 
    labs(fill = "Genes biotypes", title = "") + 
    scale_fill_discrete(labels = sprintf("%s: %s%%", gene_bio_tbl$Var1, gene_bio_tbl$ratio)) +
    coord_polar("y")
  
  p.col = 
    ggplot(bio_type_tbl, aes(x = reorder(Var1, -n), y = n, fill = Var1)) + geom_col(show.legend = F) +
    geom_text(aes(y = n + max(n)*0.04), label = bio_type_tbl$n)+
    theme_classic(base_size = 20) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  p.gene.col = 
    ggplot(gene_bio_tbl, aes(x = reorder(Var1, -n), y = n, fill = Var1)) + geom_col(show.legend = F) +
    geom_text(aes(y = n + max(n)*0.04), label = gene_bio_tbl$n) +
    theme_classic(base_size = 20) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  gObj = gridExtra::grid.arrange(
    p.pie, p.gene.pie, p.col, p.gene.col,
    nrow = 2, ncol = 2
  )
  invisible(list(trans_info = info_tbl, gene_info = gene_info_tbl, graphObj = gObj))
}



# split violin 
# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# plot protein domains ----
draw_mainChain <- function(data = data,
                           size = 0.5,
                           label_size = 4,
                           margin = 0.1,
                           ...) {
  labels = unique(data$entryName)[order(unique(data$order))]
  p <- ggplot2::ggplot(data = data) + 
    ggplot2::geom_rect(mapping=ggplot2::aes(xmin=1,
                                            xmax=total_length, 
                                            ymin=order-margin, 
                                            ymax=order+margin), 
                       ...) +
    ggplot2::scale_y_continuous(breaks = 1:length(labels), 
                                labels = labels)
  
  p <- p + ggplot2::labs(x = "Amino acid number") # label x-axis
  p <- p + ggplot2::labs(y = "") # label y-axis
  
  p <- p +
    ggplot2::theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank())
  
  return(p)
}
draw_substructure <- function(p,
                              data = data,
                              label_size = 4, 
                              margin = .1,
                              ...){
  begin=end=description=NULL
  p <- p + ggplot2::geom_rect(data= data[data$type == "substructure",],
                              mapping=ggplot2::aes(xmin=begin,
                                                   xmax=end,
                                                   ymin=order-margin,
                                                   ymax=order+margin,
                                                   fill=item), ...)
  return(p)
}
draw_domains <- function(p,
                         data = data,
                         label_domains = TRUE,
                         label_size = 4, 
                         margin = .25,
                         ...){
  begin=end=description=NULL
  p <- p + ggplot2::geom_rect(data= data[data$type == "DOMAIN",],
                              mapping=ggplot2::aes(xmin=begin,
                                                   xmax=end,
                                                   ymin=order-margin,
                                                   ymax=order+margin,
                                                   fill=description), ...)
  
  if(label_domains == TRUE){
    p <- p + ggplot2::geom_label(data = data[data$type == "DOMAIN", ],
                                 ggplot2::aes(x = begin + (end-begin)/2,
                                              y = order,
                                              label = description),
                                 size = label_size)
  }
  return(p)
}

# a wrapper for simple plots of txs 
local_drawProteinStructure = function(t_domain_info_tbl){
  # hard coded plotting
  if(!nrow(t_domain_info_tbl)){
    return(ggplot())
  }
  # data preparation
  t_domain_info_tbl_pfam = t_domain_info_tbl %>% dplyr::filter(database %in% c("Pfam", "TMHMM", "SignalP_EUK")) %>%
    dplyr::mutate(type = ifelse(database == "Pfam", "DOMAIN", "substructure"))
  
  # make plot
  p = draw_mainChain(t_domain_info_tbl_pfam, fill = "grey")
  p = draw_substructure(p, t_domain_info_tbl_pfam)
  p = draw_domains(p, t_domain_info_tbl_pfam, label_domains = F, color = "white", alpha = .8, margin = .2) + 
    theme(legend.direction = "vertical", legend.position = "bottom")
  p = p + labs(fill = "")
  return(p)
}

# 
local_splitExon = function(x, UTR_start, UTR_end){
  # x being a row of tbl
  breakpoints = table(sort(c(x$chunk_start, x$chunk_end, UTR_start, UTR_end)))
  middle = (names(breakpoints) %>% as.numeric() %>% sort())[2]
  common_end = which(breakpoints>1) %>% names() %>% as.numeric()
  x = x[c(1,1), ]
  x$chunk_start[2] = UTR_start
  x$chunk_end[2] = UTR_end
  x$chunk_coding_type[2] = F
  if(x$chunk_start[1] == common_end){ # row 1 has to be the longer one
    x$chunk_start[1] = middle
  } else {
    x$chunk_end[1] = middle
  }
  return(x)
}
# plot transcript structure ----
local_drawTranscriptStrucutre = function(exon_info_tbl,
                                         margin_MC = .05,
                                         margin_EXON = .1) {
  if(!nrow(exon_info_tbl)){
    return(ggplot())
  }
  
  # data preparation
  tx_coding_tbl = exon_info_tbl %>% dplyr::select(`Transcript name`, bioType) %>% distinct()
  tx_coding_type = tx_coding_tbl$bioType == "protein_coding" 
  names(tx_coding_type) = tx_coding_tbl$`Transcript name`
  exon_info_tbl = exon_info_tbl %>% dplyr::mutate(
    order = as.numeric(as.factor(`Transcript name`)), 
    tx_coding_type = as.logical(tx_coding_type[match(`Transcript name`, names(tx_coding_type))]),
    chunk_start = `Exon region start (bp)`, 
    chunk_end = `Exon region end (bp)`,
    chunk_coding_type = tx_coding_type
  )
  # add existed UTR
  chunk_change_codingtype = with(exon_info_tbl, `5' UTR start` == chunk_start & `5' UTR end` == chunk_end) | with(exon_info_tbl, `3' UTR start` == chunk_start & `3' UTR end` == chunk_end)
  exon_info_tbl$chunk_coding_type[chunk_change_codingtype] = F
  # modification for chunk plots
  t_UTR_ridx = (!is.na(exon_info_tbl$`5' UTR start`)) | (!is.na(exon_info_tbl$`3' UTR start`))
  UTR_info_tbl = exon_info_tbl[t_UTR_ridx, ] # extract every row that UTRs are not NA
  
  # split exon ----
  t_UTR_rec_tbl = tibble()
  for(i in 1:nrow(UTR_info_tbl)){
    t_UTR_info_row = UTR_info_tbl[i, ]
    
    
    # if 5' exists
    if(!is.na(t_UTR_info_row$`5' UTR end`)){
      UTR_start = t_UTR_info_row$`5' UTR start`
      UTR_end = t_UTR_info_row$`5' UTR end`
      if(!((t_UTR_info_row$chunk_end - t_UTR_info_row$chunk_start) == (UTR_end - UTR_start))){
        t_UTR_info_row = local_splitExon(t_UTR_info_row, UTR_start, UTR_end)
      }
    } 
    
    # if 3' exists
    if(!is.na(t_UTR_info_row[1, ]$`3' UTR end`)){
      UTR_start = t_UTR_info_row[1, ]$`3' UTR start`
      UTR_end = t_UTR_info_row[1, ]$`3' UTR end`
      if(!((t_UTR_info_row[1, ]$chunk_end - t_UTR_info_row[1, ]$chunk_start) == (UTR_end - UTR_start))){
        t_UTR_info_row = bind_rows(
          t_UTR_info_row[-1, ],
          local_splitExon(t_UTR_info_row[1, ], UTR_start, UTR_end)
        )
      }
    } 
    
    t_UTR_rec_tbl = bind_rows(t_UTR_rec_tbl, t_UTR_info_row)
  }
  
  exon_info_tbl = bind_rows(
    exon_info_tbl[!t_UTR_ridx, ], 
    t_UTR_rec_tbl)
  
  
  # draw main chain ----
  # setting direction
  labels = unique(exon_info_tbl$`Transcript name`)[order(unique(exon_info_tbl$order))]
  
  p <- ggplot(data = exon_info_tbl) + 
    geom_rect(aes(
      xmin = `Transcript start (bp)`,
      xmax = `Transcript end (bp)`,
      ymin = order - margin_MC,
      ymax = order + margin_MC), 
      alpha = .5, fill = "grey"
    ) + 
    scale_y_continuous(breaks = 1:length(labels), labels = labels) +
    # theme_classic(base_size = 16) +
    coord_cartesian(ylim = c(.5 , length(labels) + .5)) +
    theme(axis.line.y = element_blank(), 
          axis.ticks.y = element_blank())
  
  if(unique(exon_info_tbl$Strand) < 0) {
    # labels = paste0("< ", )
    p = p + scale_x_reverse()
  } 
  
  # draw chunks ----
  coding_colors = gg_color_hue(2) %>% sort()
  p = 
    p + 
    geom_rect(aes(
      xmin = chunk_start, 
      xmax = chunk_end, 
      ymin = order - margin_EXON, 
      ymax = order + margin_EXON,
      color = tx_coding_type, 
      fill = chunk_coding_type)) +
    scale_color_manual(labels = c("non-coding", "coding"), values = coding_colors, limits = c(F, T)) + 
    scale_fill_manual(labels = c("non-coding", "coding"), values = c("white", coding_colors[2]), limits = c(F, T)) 
  p = p + theme(legend.position = "bottom", legend.direction = "vertical")
  return(p)
}

# wrapper for single gene ----
f_plotGeneInspect = function(gene_name,
                             chosen_txs = NULL,
                             by = "cluster", 
                             chosen_cats = NULL,
                             annotation_path = tsv_path, 
                             return_data = F) {
  
  ## Plot expression, enrichment, raw count and transcript/protein expression
  ## 
  ## Args:
  #'  @gene_name: Target gene
  #'  @chosen_txs: Which txs should be plot
  #'  @by: "cluster" or "tissue" to plot by
  #'  @chosen_cats: Which cluster should be plot
  #'  
  #'  Returns: 
  #'  A ggplot item with combined plot for the chosen gene. 
  
  if(by == "cluster") {
    expressed_cat = rep(T, length(unique(pCells$characteristics..majorCluster)))
    names(expressed_cat) = sort(unique(pCells$characteristics..majorCluster))
    obj = all_obj_lst[[gene_name]]
  } else {
    expressed_cat = rep(T, length(unique(pCells$characteristics..tissueType)))
    names(expressed_cat) = sort(unique(pCells$characteristics..tissueType))
    obj = all_obj_lst_byTissue[[gene_name]]
  }
  
  if(is.null(chosen_cats)){
    chosen_cats = names(expressed_cat[expressed_cat])
  }
  
  
  # choose only meaningful txs to plot structure
  chn_fishers = obj$fisher.p.adj
  canonical_tx = (obj$count[grep(
    rownames(obj$count),
    pattern = "NO_EXPR|MIX",
    invert = T,
    value = T
  ),] %>% apply(1, sum) %>% sort(decreasing = T))[1] %>% names() 
  txsForStructures = apply(obj$fisher.p.adj, 1, function(x){sum(x>1.3) > 0})
  txsForStructures = names(txsForStructures)[txsForStructures]
  txsForStructures = c(txsForStructures, canonical_tx)
  txsForStructures = txsForStructures %>% grep(pattern = "NO_EXPR|MIX", value = T, invert = T)
  # choose txs, by default use the arguments
  if(is.null(chosen_txs)){
    chn_txs = 
      c(rownames(chn_fishers)[chn_fishers %>% apply(1, function(x) {
        sum(x > 1.3)}) > 0] %>% 
          grep(pattern = "NO_EXPR", invert = T, value = T), canonical_tx) %>% sort()
  } else {
    chn_txs = chosen_txs
  }
  chn_txs = chn_txs %>% unique() %>% sort()
  
  # plt1 for gene level expression ----
  p1 = exprMat_sum_tbl %>% dplyr::filter(Gene_name == gene_name) %>% dplyr::filter_(.dots = paste0(by, "%in% chosen_cats")) %>% 
    ggplot(aes_string(x = by, y = "log2(tpm+1)", color = by, fill = by)) +
    geom_bar(stat = "summary", fun.y = mean, size = 1.2, width = 0.8) +
    geom_errorbar(stat = "summary", fun.y = mean, fun.ymin = function(x){mean(x) - sd(x)/sqrt(length(x))}, fun.ymax = function(x){mean(x) + sd(x)/sqrt(length(x))}, width = .4, size = 1) +
    # geom_boxplot(width = 0.1, aes_string(color = by), fill = "white") +
    labs(title = gene_name) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      # axis.text.x = element_blank(),
      axis.title.x = element_blank())
  
  # count View with pval fill color ----
  p3 = obj$count[chn_txs, chosen_cats] %>%
    melt(varnames = c("id", "category"), value.name = "cell_count") %>% 
    left_join(
      chn_fishers[chn_txs, chosen_cats] %>% melt(varnames = c("id", "category"), value.name = "-log10(fisher.p.adj)"),
      by = c("id", "category")
    ) %>% as_tibble() %>% 
    dplyr::group_by(category) %>%
    ggplot(aes(x = category, y = id)) +
    geom_tile(aes(fill = `-log10(fisher.p.adj)`), color = "white") +
    geom_text(aes(label = `cell_count`)) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 2), na.value = "red") +
    labs(fill = "-log10(q value)") +
    theme(
      axis.text.y = element_text(face = colorado(src = chn_txs, canonical_tx)),
      axis.text.x = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  # # count View ----
  # p3 = obj$count[chn_txs, chosen_cats] %>%
  #   melt(varnames = c("id", "category"), value.name = "cell_count") %>% as_tibble() %>% 
  #   dplyr::group_by(category) %>% 
  #   ggplot(aes(x = category, y = id)) +
  #   geom_tile(aes(fill = `cell_count`), color = "white") +
  #   geom_text(aes(label = `cell_count`)) +
  #   scale_fill_gradient(low = "white") +
  #   labs(fill = "Cell count") +
  #   theme(
  #     axis.text.y = element_text(face = colorado(src = chn_txs, canonical_tx)),
  #     axis.text.x = element_blank(), 
  #     axis.line = element_blank(),
  #     axis.ticks = element_blank(),
  #     axis.title = element_blank()
  #   )
  # 
  # # enrichment Fisher View ----
  # p4 = chn_fishers[chn_txs, chosen_cats] %>% melt(varnames = c("id", "category"), value.name = "-log10(fisher.p.adj)") %>% 
  #   ggplot(aes(x = category, y = id)) +
  #   geom_tile(aes(fill = `-log10(fisher.p.adj)`), color = "white") +
  #   geom_text(aes(label = `-log10(fisher.p.adj)`)) +
  #   scale_fill_gradient(low = "white", high = "red", limits = c(0, 2), na.value = "red") +
  #   labs(fill = "-log10(q value)") +
  #   theme(
  #     axis.text.y = element_text(face = colorado(src = chn_txs, canonical_tx)),
  #     axis.text.x = element_blank(), 
  #     axis.line = element_blank(),
  #     axis.ticks = element_blank(),
  #     axis.title = element_blank()
  #   )
  
  # transcript view ----
  test_struct_tbl = STRUCT.INFO %>% dplyr::filter(`Transcript name` %in% txsForStructures) %>% dplyr::mutate(bioType = type2bio(GENE.INFO$`Transcript type`[match(`Transcript name`, GENE.INFO$`Transcript name`)]))
  p6 = local_drawTranscriptStrucutre(test_struct_tbl)
  p6 = p6 + theme(
    legend.position = "none", 
    axis.text.y = element_text(face = colorado(src = chn_txs, canonical_tx)),
    axis.text.x = element_text(angle = 45, hjust = 1))
  
  # protein view ----
  annotation_tsvs = paste0(annotation_path, txsForStructures, ".fasta.tsv.txt")
  t_domain_info_tbl = local_loadDomains(annotation_tsvs)
  if(nrow(t_domain_info_tbl)){
    t_domain_info_tbl = t_domain_info_tbl %>% dplyr::mutate(entryName = `Transcript name`)
    p5 = local_drawProteinStructure(t_domain_info_tbl) 
    p5 = p5 + theme(
      axis.text.y = element_text(face = colorado(src = chn_txs, canonical_tx)),
      legend.position = "bottom", legend.direction = "vertical", axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")
  } else{
    p5 = ggplot()
  }
  
  # plot ----
  t_nProtein = t_domain_info_tbl$`Transcript name` %>% unique() %>% length()
  require(patchwork)
  # total = p1 + p3 + p4 + p6 + p5 + patchwork::plot_layout(
  total = p1 + p3 + p6 + p5 + patchwork::plot_layout(
    ncol = 1, 
    heights = c(1,1,1,length(chn_txs) * 0.3,t_nProtein * 0.3))
  if(return_data){
    return(list(plot = total, count = obj$count[chn_txs, chosen_cats]))
  }
  return(total)
}
# total

# Calculation ----
DoubleMAD <- function(x, zero.mad.action="warn"){
  # http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/
  # (1) the median absolute deviation from the median of all points less than or equal to the median and 
  # (2) the median absolute deviation from the median of all points greater than or equal to the median. 
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  x         <- x[!is.na(x)]
  m         <- median(x)
  abs.dev   <- abs(x - m)
  left.mad  <- median(abs.dev[x<=m])
  right.mad <- median(abs.dev[x>=m])
  if (left.mad == 0 || right.mad == 0){
    if (zero.mad.action == "stop") stop("MAD is 0")
    if (zero.mad.action %in% c("warn", "warn and na")) warning("MAD is 0")
    if (zero.mad.action %in% c(  "na", "warn and na")){
      if (left.mad  == 0) left.mad  <- NA
      if (right.mad == 0) right.mad <- NA
    }
  }
  return(c(left.mad, right.mad))
}

DoubleMADsFromMedian <- function(x, zero.mad.action="warn"){
  # http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  two.sided.mad <- DoubleMAD(x, zero.mad.action)
  m <- median(x, na.rm=TRUE)
  x.mad <- rep(two.sided.mad[1], length(x))
  x.mad[x > m] <- two.sided.mad[2]
  mad.distance <- abs(x - m) / x.mad
  mad.distance[x==m] <- 0
  return(mad.distance)
}

# testset_ans4 = c(3,5,4,8,10)
# testset_ans3 = c(25,8,5,3,3)
# test_Hindex = function(x){
#   x = x[!is.nan(x)]
#   x = x[!is.na(x)]
#   if(length(x) == 0)
#     return(NA)
#   n = length(x)
#   x1 = x / max(x)
#   y = (x1 * n) %>% sort(decreasing = T) %>% bind_cols(sn = seq(n), value = .)
#   z = (y[(which(y$sn > y$value) %>% head(1))-1,])$sn / n * max(x)
#   return(z)
# }

# return duplicated items
duplicated.item = function(ori, return.idx = F){
  tmp_idx1 = duplicated(ori)
  tmp_idx2 = duplicated(ori, fromLast = T)
  chosen_ones = ori[tmp_idx1 | tmp_idx2]
  if(return.idx)
    return((1:length(ori))[tmp_idx1 | tmp_idx2])
  else
    return(chosen_ones)
}

Counts2TPM <- function(countsMat, feature_lengths, FragmentLength = 100){
  # each length correspond to a transcript name
  # identify transcripts & samples specific fragmentlength are too much works... 
  # Just use a empirical mean value ...
  # Transcripts shorter than FragmentLength would be remove
  # So, WARN!! nrow would change!
  stopifnot(all(names(feature_lengths) == rownames(countsMat)))
  stopifnot(!anyNA(FragmentLength))
  if(length(FragmentLength) == 1){
    FragmentLength = rep(FragmentLength, ncol(countsMat))
  }
  n_cells = ncol(countsMat)
  n_tx = nrow(countsMat)
  
  # calculate effective length
  stopifnot(length(FragmentLength) == 1 | length(FragmentLength) == n_cells)
  feature_length_mt = matrix(rep(feature_lengths, n_cells), ncol = n_cells)
  FragmentLength_mt = t(matrix(rep(FragmentLength, n_tx), ncol = n_tx))
  eff_lengths_mt = feature_length_mt - FragmentLength_mt + 1
  
  # remove short transcripts ( < 0) 
  r_idx = !as.logical(rowSums(eff_lengths_mt <= 0))
  
  countsMat = countsMat[r_idx, ]
  eff_lengths_mt = eff_lengths_mt[r_idx, ]
  # duplicate feature length to matrix for easier calculation
  rownames(eff_lengths_mt) = rownames(countsMat)
  colnames(eff_lengths_mt) = colnames(countsMat)
  
  # calculate TPM
  if(all(dim(countsMat) == dim(eff_lengths_mt)) ){
    scaled.countsMat <- divMatrix(countsMat, eff_lengths_mt)
    tpm <- 10^6 * t(t(scaled.countsMat) / colSums(scaled.countsMat))
    rownames(tpm) = rownames(countsMat)
    colnames(tpm) = colnames(countsMat)
    tpm = round(tpm, digits = 2)
    return(tpm)
  } else{
    warning("The dimensions of countsMat and featurelength matrix are different")
    return(NULL)
  }
}

divMatrix <- function(m1,m2){
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( all(dim_m1 == dim_m2) ){
    # div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow=dim_m1[1] )
    # row.names(div.result) <- row.names(m1)
    # colnames(div.result) <- colnames(m1)
    div.result = m1 / m2
    # for(i in 1:dim_m1[1]){
    #   for(j in 1:dim_m1[2]){
    #     div.result[i,j] <- m1[i,j] / m2[i,j]
    #   }
    # }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

dea.ttest = function(x, ...){
  # x being the list with two named expression matrixes,
  # ... whose rownames correspond to tx_names and colnames to sample_names
  # other augments: 
  mt1 = x[[1]]
  mt2 = x[[2]]
  
  args = list(...)
  progress = ifelse(is.null(args$progress), T, args$progress)
  
  stopifnot(nrow(mt1) == nrow(mt2))
  
  n_tx = nrow(mt1)
  p_vals = rep(NA, n_tx)
  means_mt1 = rep(NA, n_tx)
  means_mt2 = rep(NA, n_tx)
  iter01 = ".I1_fun"
  initiatePB(iter01)
  for(i in 1:n_tx){
    t_res = t.test(mt1[i, ], mt2[i, ])
    p_vals[i] = t_res$p.value
    means_mt1[i] = mean(mt1[i,])
    means_mt2[i] = mean(mt2[i,])
    if(progress) processBar(iter01, i, n_tx, tail = "ETA")
  }
  FC = means_mt2 / (means_mt1 + .0001)
  LFC = log2(FC)
  p_adjs = p.adjust(p_vals)
  
  res_tbl = tibble(id = rownames(mt1), 
                   means_mt1 = means_mt1,
                   means_mt2 = means_mt2,
                   FoldChange = FC,
                   LogFC = LFC,
                   P.Value = p_vals,
                   adj.P.Val = p_adjs)
  
  return(res_tbl)
}

dea.wilcox = function(x, ...){
  # x being the list with two named expression matrixes,
  # ... whose rownames correspond to tx_names and colnames to sample_names
  # other augments: 
  mt1 = x[[1]]
  mt2 = x[[2]]
  
  args = list(...)
  progress = ifelse(is.null(args$progress), T, args$progress)
  
  stopifnot(nrow(mt1) == nrow(mt2))
  
  n_tx = nrow(mt1)
  p_vals = rep(NA, n_tx)
  means_mt1 = rep(NA, n_tx)
  means_mt2 = rep(NA, n_tx)
  iter01 = ".I1_fun"
  initiatePB(iter01)
  for(i in 1:n_tx){
    t_res = wilcox.test(mt1[i, ], mt2[i, ])
    p_vals[i] = t_res$p.value
    means_mt1[i] = mean(mt1[i,])
    means_mt2[i] = mean(mt2[i,])
    if(progress) processBar(iter01, i, n_tx, tail = "ETA")
  }
  FC = means_mt2 / (means_mt1 + .0001)
  LFC = log2(FC)
  p_adjs = p.adjust(p_vals)
  
  res_tbl = tibble(id = rownames(mt1), 
                   means_mt1 = means_mt1,
                   means_mt2 = means_mt2,
                   FoldChange = FC,
                   LogFC = LFC,
                   P.Value = p_vals,
                   adj.P.Val = p_adjs)
  
  return(res_tbl)
}

dea.limma <- function(data_set = NULL, ...){
  # only support 2 groups in the list
  # data_set = list(a = dataframe1, b = dataframe2), ... is for limma::decideTests and limma::topTable
  stopifnot(length(data_set) ==2)
  
  args = list(...)
  
  suppressMessages(require(limma))
  suppressMessages(require(dplyr))
  
  nrows <- as.numeric(lapply(
    data_set,
    FUN = function(x)
      nrow(x)
  ))
  stopifnot(all((nrows - median(nrows)) == 0))
  
  exprMat <- do.call(cbind, data_set)
  group <- c()
  # for(i in 1:length(data_set)) group <- c(group, rep(names(data_set)[i], length(data_set[[i]])))
  group = c(rep(names(data_set)[1], ncol(data_set[[1]])), rep(names(data_set)[2], ncol(data_set[[2]])))
  design <- model.matrix(~0+factor(group))
  colnames(design) <- c("STC","CTL")
  fit <- lmFit(exprMat, design)
  contrast.matrix <- makeContrasts(STC-CTL, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  # results <- decideTests(fit2, method = "global", adjust.method = "BH", p.value=0.01, lfc = 2)
  results <- decideTests(fit2)
  x <- topTable(fit2, adjust.method = "BH", sort.by = "B", resort.by = "M", ...)
  
  x = dplyr::bind_cols(id = rownames(x),x )
  invisible(x)
}

dea.kruskal = function(x, progress = T, ...){
  # to be modified or deprecated, for the tests should be done by multi-groups in moust cases.
  # x being the list with two named expression matrixes,
  # ... for lower level arguments
  mt1 = x[[1]]
  mt2 = x[[2]]
  
  # args = list(...)
  # progress = ifelse(is.null(args$progress), T, args$progress)
  
  stopifnot(nrow(mt1) == nrow(mt2))
  
  n_tx = nrow(mt1)
  p_vals = rep(NA, n_tx)
  means_mt1 = rep(NA, n_tx)
  means_mt2 = rep(NA, n_tx)
  iter01 = ".I1_fun"
  initiatePB(iter01)
  for(i in 1:n_tx){
    t_res = kruskal.test(list(mt1[i, ], mt2[i, ]), ...)
    p_vals[i] = t_res$p.value
    means_mt1[i] = mean(mt1[i,])
    means_mt2[i] = mean(mt2[i,])
    if(progress) processBar(iter01, i, n_tx, tail = "ETA")
  }
  FC = means_mt2 / means_mt1 + 0.001
  LFC = log2(FC)
  p_adjs = p.adjust(p_vals)
  
  res_tbl = tibble(id = rownames(mt1), 
                   means_mt1 = means_mt1,
                   means_mt2 = means_mt2,
                   FoldChange = FC,
                   LogFC = LFC,
                   P.Value = p_vals,
                   adj.P.Val = p_adjs)
  
  return(res_tbl)
}

dea.limma.voom = function(data_set = NULL, ...){
  # only support 2 groups in the list
  # data_set = list(a = dataframe1, b = dataframe2), ... is for limma::decideTests and limma::topTable
  
  stopifnot(length(data_set) == 2) 
  args = list(...)
  
  suppressMessages(require(limma))
  suppressMessages(require(dplyr))
  
  nrows <- as.numeric(lapply(
    data_set,
    FUN = function(x)
      nrow(x)
  ))
  stopifnot(all((nrows - median(nrows)) == 0))
  
  exprMat <- do.call(cbind, data_set)
  group <- c()
  for(i in 1:length(data_set)) group <- c(group, rep(names(data_set)[i], ncol(data_set[[i]])))
  design <- model.matrix(~0+factor(group))
  colnames(design) <- c("STC","CTL")
  v <- voom(exprMat, design, plot=TRUE)
  fit <- lmFit(v, design)
  contrast.matrix <- makeContrasts(STC-CTL, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  results <- decideTests(fit2, ...)
  x <- topTable(fit2, coef = 1, number = 30000, adjust.method = "BH", sort.by = "B", resort.by = "M")
  
  x = dplyr::bind_cols(id = rownames(x),x )
  invisible(x)
}

Remove_Batch <- function(raw_counts, batch){
  #' Remove batch method adapted from Pagoda2
  # raw_counts being a matrix with row corresponding to transcripts and cols to cells
  # batch being a factor with length of cols
  require(bit64)
  
  batch = as.factor(batch)
  row_names = rownames(raw_counts)
  raw_counts = raw_counts %>% as.data.frame() %>% lapply(as.integer64) %>% do.call(cbind,.)
  rownames(raw_counts) = row_names
  #Depth for each cell and each batch
  # cell.depth <- colSums(raw_counts)
  # batch.depth <- tapply(cell.depth,batch,sum,simplify = T)
  # batch.depth <- setNames(as.integer(batch.depth), names(batch.depth))
  cell.depth = lapply(raw_counts %>% as.data.frame, sum.integer64) %>% do.call(c, .)
  batch.depth.lst = list()
  for(i in levels(batch)){
    batch.depth.lst[[i]] = sum.integer64(cell.depth[batch==i])
  }
  batch.depth = do.call(c, batch.depth.lst)
  
  #Average counts for each gene
  # gene.average <-  (rowSums(raw_counts) + length(levels(batch)))/(sum(cell.depth) + length(levels(batch)))
  row_sum_raw_counts = raw_counts %>% t() %>% as.data.frame() %>% lapply(sum.integer64) %>% do.call(c,.) 
  gene.average = (row_sum_raw_counts + length(levels(batch)))/(sum.integer64(cell.depth) + length(levels(batch)))
  
  #Total counts percentage for each gene in each batch
  gene_counts.batch.lst = list()
  for(i in levels(batch)){
    gene_counts.batch.lst[[i]] = lapply(raw_counts[, batch == i] %>% t() %>% as.data.frame() , 
                                        function(x) {y = sum.integer64(x); log(y+1)}) %>% unlist()
  }
  gene_counts.batch_log = do.call(rbind, gene_counts.batch.lst) 
  gene_percentage.batch = t(gene_counts.batch_log - log(batch.depth+1))
  
  # gene_counts.batch <- log(rowsum(t(raw_counts), batch) +1 )
  # gene_percentage.batch <- t(gene_counts.batch - log(batch.depth+1))
  
  #Batch_factor for each gene in each batch
  batch_factor <- exp(gene_percentage.batch - log(gene.average)) 
  
  #Split all cells by batch
  cells.batch = split(x = colnames(raw_counts), f = batch)
  
  #Calculated the scaled counts
  scaled_counts <- lapply(seq_along(cells.batch), function(idx) {
    batch.name = names(cells.batch)[idx]
    cells_in_batch = cells.batch[[idx]]
    return(raw_counts[, cells_in_batch, drop = F]/batch_factor[,batch.name])
  })
  scaled_counts = do.call("cbind", scaled_counts)
  scaled_counts = scaled_counts[, colnames(raw_counts)]
  
  return(scaled_counts)
}

# Others ----
# originally used in expl01_specificAS.R
transCountsMat = function(a, CD = NA, pat = NA, tissue = NA, group = "tissueType"){
  require(dplyr)
  
  CDType = a$MCtype
  CDType[grep("CD4", a$MCtype)] = "CD4"
  CDType[grep("CD8", a$MCtype)] = "CD8"
  b = bind_cols(tissueType = substr(a$sampleType, 1, 1), CDType = CDType, a)
  if(!is.na(CD))
    b = b %>% dplyr::filter(CDType == CD)
  if(!is.na(pat) && pat != "ALL")
    b = b %>% dplyr::filter(patientID == pat)
  if(!is.na(tissue) && tissue != "all")
    b = b %>% dplyr::filter(tissueType == tissue)
  
  c = b %>% dplyr::group_by_(group) %>% dplyr::summarise_if(is.numeric, sum)
  
  d = t(as.data.frame(c[, -1]))
  colnames(d) = eval(parse(text = sprintf("c$%s", group)))
  invisible(d)
}

colorado = function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}