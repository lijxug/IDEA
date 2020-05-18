# ---- loading dependencies ----
suppressPackageStartupMessages(library("assertthat"))

# ---- ToolKit ----
loginfo <- function(..., printnow = T) {
  msg = paste0(list(...), collapse = "")
  msg <- paste0("[",format(Sys.time()), "] ", msg,"\n")
  if(printnow)
    cat(msg)
  invisible(msg)
}

msgASAP <- function(ttl, msg = ""){
  # warning: authorized access only
  SCKEY="SCU11537T7fea873af4eba0e95c9ef25143da731159b38f97bd60d"
  website="https://sc.ftqq.com/"
  url = paste0(website, SCKEY, ".send?")
  
  assert_that(is.character(ttl))
  ttl = gsub("\\s", "%20", ttl)
  
  url = paste0(url, "text=", ttl)
  if(nchar(msg) != 0){
    msg = gsub("\\s", "%20", msg)
    url = paste0(url, "&desp=", msg)
  }
  x = readLines(url,encoding="UTF-8", warn = F)
  invisible(x)
}

# process bar ----
# a terminal process Bar, use initiatePB to initialize at the beginning of the forloop,
# and insert processBar with objName, i, cycles inside the for loop
processBar = function(objName,
                      i,
                      cycles,
                      title = "Process",
                      scale = 40,
                      sign = "#", 
                      tail = "",
                      terminal = "R", # terminal could be R/log, others default to shell
                      final = "Work done!") {
  suppressPackageStartupMessages(require("iterators"))
  if (!exists(objName)) {
    if (terminal != "R")
      words_list = unlist(lapply(1:cycles, function(x) {
        sprintf(
          paste0("\033[?25l\r%s %5.1f%% | %-", scale, "s | "),
          title,
          x * 100 / cycles ,
          paste0(rep(sign, ceiling(x * scale / cycles)), collapse = "")
        )
      }))#\033[?25l hide the cursor - linux control code
    else
      words_list = unlist(lapply(1:cycles, function(x) {
        sprintf(
          paste0("\r%s %5.1f%% | %-", scale, "s | "),
          title,
          x * 100 / cycles ,
          paste0(rep(sign, ceiling(x * scale / cycles)), collapse = "")
        )
      }))
    eval(parse(text = sprintf("%s <<- iter(words_list)", objName)))
    eval(parse(text = sprintf("%s <<- Sys.time()", paste0(".TIC_",objName))))
    # if i didn't start at 1
    times = i
    while (times > 1) {
      msg = eval(parse(text = sprintf("nextElem(%s)", objName)))
      times = times - 1
    }
  }
  
  msg = eval(parse(text = sprintf("nextElem(%s)", objName)))
  if(tail == "ETA"){
    .tic = eval(parse(text = sprintf("%s", paste0(".TIC_",objName))))
    if(terminal != "R"){
      tail = paste0("ETA: ", format(round((Sys.time() - .tic) / i * (cycles - i), digits = 2)), "\033[K")
    }
    else {
      tail = paste0("ETA: ", format(round((Sys.time() - .tic) / i * (cycles - i), digits = 2)), "   ")
    }
  }
  if(terminal == "log")
    tail = paste0(tail, "\n")
  cat(paste0(msg, tail))
  if(i == cycles){
    if(nchar(final)) final = loginfo(final, printnow = F)
    if(terminal != "R") cat("\033[?25h")
    cat(paste0("\n", final))
    rm(list = objName, inherits = T)
  }
}

initiatePB = function(iterOBJ){
  .tic = sprintf("%s", paste0(".TIC_", iterOBJ))
  rm_list = c(iterOBJ, .tic)
  if(any(exists(rm_list, inherits = T)))
    rm(list = rm_list, inherits = T)
}

# entrez symbol translattion ----
entrezToXXX<-function(x,type="SYMBOL",species="human")
{
  ret=c()
  x=as.character(x)
  if(species=="human")
  {
      suppressPackageStartupMessages(require("org.Hs.eg.db"))
  }else if(species=="mouse")
  {
      suppressPackageStartupMessages(require("org.Mm.eg.db"))
  }else
  {
    stop("species must be human or mouse\n")
  }
  if(species=="human")
  {
    if(type=="SYMBOL")
    {
      ret=unlist(as.list(org.Hs.egSYMBOL))[x]
    }else if(type=="ENSG")
    {
      ret=unlist(as.list(org.Hs.egENSEMBL))[x]
    }else if(type=="GENENAME")
    {
      ret=unlist(as.list(org.Hs.egGENENAME))[x]
    }
  }else if(species=="mouse")
  {
    if(type=="SYMBOL")
    {
      ret=unlist(as.list(org.Mm.egSYMBOL))[x]
    }else if(type=="ENSG")
    {
      ret=unlist(as.list(org.Mm.egENSEMBL))[x]
    }
  }
  return(as.character(ret))
}

XXXToEntrez<-function(x,type="SYMBOL",species="human")
{
  ret=c()
  x=as.character(x)
  if(species=="human")
  {
    suppressPackageStartupMessages(require("org.Hs.eg.db"))
  }else if(species=="mouse")
  {
    suppressPackageStartupMessages(require("org.Mm.eg.db"))
  }else
  {
    stop("species must be human or mouse\n")
  }
  if(species=="human")
  {
    if(type=="SYMBOL")
    {
      ret=unlist(as.list(org.Hs.egSYMBOL2EG))[x]  
    }else if(type=="ENSG")
    {
      ret=unlist(as.list(org.Hs.egENSEMBL2EG))[x]
    }
  }else if(species=="mouse")
  {
    if(type=="SYMBOL")
    {
      ret=unlist(as.list(org.Mm.egSYMBOL2EG))[x]  
    }else if(type=="ENSG")
    {
      ret=unlist(as.list(org.Mm.egENSEMBL2EG))[x]
    }
  }
  return(as.character(ret))
}




# gene enrichment
geneEnricher_bycP = function(geneList, db = "gsea", pvalCutoff = 0.05){
  # db = "gsea" | "kegg" | "go" 
  require("clusterProfiler")
  require(org.Hs.eg.db)
  require("dplyr")
  
  h_all_MSigDB = read.gmt(paste0(getwd(), "/00.ref_tools/h.all.v6.1.symbols.gmt"))
  geneList_id = XXXToEntrez(geneList)
  enRes = switch(
    db,
    "gsea" = enricher( gene = geneList, TERM2GENE = h_all_MSigDB, pvalueCutoff = pvalCutoff ),
    "kegg" = enrichKEGG( gene = geneList_id, organism = 'hsa', pvalueCutoff = pvalCutoff ),
    "go" =  enrichGO(gene = geneList_id, OrgDb="org.Hs.eg.db", ont = "BP", pvalueCutoff = pvalCutoff), 
    default = stop("Invalid db choice.")
  )
  enRes = as_tibble(enRes)
  if(nrow(enRes) > 0)
    enRes_tbl = enRes %>% mutate(geneSymbol = unlist(lapply(strsplit(enRes$geneID, "/"), function(x)
      paste(entrezToXXX(x), collapse = "/"))))
  else
    enRes_tbl = enRes
  return(enRes_tbl)
}
# plot ----
extractLegend<-function(a.gplot){ 
  # a function to extract legends from ggplot object
  # http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

getDensity = function(x, y, n = 100) {
  require(MASS)
  tryCatch({
    dens <- MASS::kde2d(x = x, y = y, n = n)
  }, error = function(e) {
    print(e)
    warning("Swith bandwidth to h = 1")
    dens <<- MASS::kde2d(x = x,
                         y = y,
                         n = n,
                         h = 1)
  })
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## ggvenn from GuangchuangYu/yyplot
##' venn plot using ggplot2
##'
##'
##' @title ggvenn
##' @param x data
##' @param alpha transparency of color
##' @return ggplot object
##' @importFrom rvcheck get_fun_from_pkg
## @importFrom venneuler venneuler
##' @importFrom ggforce geom_circle
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 coord_fixed
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_
##' @export
##' @author guangchuang yu; altered by Jason Li
##' @examples
##' \dontrun{
##' set.seed(2017-11-08)
##' x <- matrix(sample(0:4, 40, TRUE, c(.5, .1, .1, .1, .1)), ncol=4)
##' colnames(x) <- LETTERS[1:4]
##' ggvenn(x)
##' }
ggvenn <- function(x, alpha = 0.5, return_data = F) {
  ## the venneuler depends on rJava
  ## maybe wrapping VennDiagram (output gList) as geom layer
  # add myself ---
  require(rvcheck)
  require(venneuler)
  require(ggforce)
  # ---
  venneuler <- get_fun_from_pkg("venneuler", "venneuler") ## for easy installation, as many users have issue of rJava installation
  y <- venneuler(x)
  d <- data.frame(y$centers,
                  diameters = y$diameters,
                  labels = y$labels,
                  stringsAsFactors = FALSE)
  if(return_data) return(d)
  ggplot(d) +
    geom_circle(aes_(x0 = ~x, y0 = ~y,
                     r = ~diameters/2, fill = ~labels),
                alpha = alpha) +
    geom_text(aes_(x = ~x, y = ~y, label = ~labels)) +
    coord_fixed()
}

# Venn.plot = function(vennList, ...){
#   suppressPackageStartupMessages(require(Vennerable))
#   suppressPackageStartupMessages(require(grid))
#   args = list(...)
#   doWeights = ifelse(is.null(args$doWeights), F, args$doWeights)
#   drawNow = ifelse(is.null(args$drawNow), T, args$drawNow)
#   type = ifelse(is.null(args$type), "circles", args$type)
#   
#   Vstem = Venn(vennList)
#   p = compute.Venn(Vstem, doWeights = doWeights, type = type)
#   if(drawNow){
#     grid.newpage()
#     plot(p)
#   }
#   invisible(p)
# }
# 
# VennDiagram.plot = function(vennList, ...){
#   suppressPackageStartupMessages(require(VennDiagram))
#   args = list(...)
#   
#   invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
#   # doWeights = ifelse(is.null(args$doWeights), T, args$doWeights)
#   # drawNow = ifelse(is.null(args$drawNow), T, args$drawNow)
#   # type = ifelse(is.null(args$type), "circles", args$type)
#   
#   # Vstem = Venn(vennList)
#   # p = compute.Venn(Vstem, doWeights = doWeights, type = type)
#   p = venn.diagram(vennList, filename = NULL)
#   if(drawNow){
#     grid.newpage()
#     grid.draw(p)
#     # plot(p)
#   }
#   invisible(p)
# }

# saving ----
save.png = function(graphObj, path, ...){
  args = list(...)
  args$width = ifelse(is.null(args$width), 800, args$width)
  args$height = ifelse(is.null(args$height), 600, args$height)
  png(filename = path,
      width = args$width, height = args$height) 
  plot(graphObj)
  invisible(dev.off())
}

EXWB.initiate = function(...) {
  suppressPackageStartupMessages(require(openxlsx))
  args = list(...)
  
  wb <- createWorkbook()
  options("openxlsx.borderColour" = "#4F80BD")
  options("openxlsx.borderStyle" = "thin")
  modifyBaseFont(wb, fontSize = 16, fontName = "Arial")
  invisible(wb)
}

EXWB.writeSheet = function(wb, sheetName, sheetData, overwrite = F,...){
  suppressPackageStartupMessages(require(openxlsx))
  suppressPackageStartupMessages(require(dplyr))
  args = list(...)
  
  stopifnot(nchar(sheetName) <= 31)
  # sheetName cannot have strange characters like !@#$%^&*()
  sheetIndex = ifelse(length(which(wb$sheet_names == sheetName)) == 0, -1, which(wb$sheet_names == sheetName))
  if(sheetIndex>0){ # exists
    if(overwrite){
      openxlsx::removeWorksheet(wb, sheetIndex)
    }else{
      return(invisible(wb))
    }
  }
  addWorksheet(wb, sheetName = sheetName, gridLines = T)
  sheetData = sheetData %>% mutate_if(is.double, round, digits = 4)
  writeDataTable(wb, sheet = sheetName, x = sheetData, colNames = TRUE, rowNames = TRUE,tableStyle = "TableStyleLight9")
  invisible(wb)
}


#' x      numeric vector for each slice
#' group  vector identifying the group for each slice
#' labels vector of labels for individual slices
#' col    colors for each group
#' radius radius for inner and outer pie (usually in [0,1])
donuts <- function(x, group = 1, labels = NA, col = NULL, radius = c(.7, 1)) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]
  
  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    Vectorize(adjustcolor)(x, alpha.f = al)
  })
  plot.new()
  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],
      col = unlist(col.sub), labels = labels)
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],
      col = unlist(col.main), labels = NA)
}

# for sashimi plot
write_sashimi_job = function(sample_ids, groups, event_name, path){
  require(dplyr)
  require(readr)
  output_tbl = tibble(sample_ids = sample_ids, group = groups, event_name = event_name)
  write_delim(output_tbl, path = path, delim = "\t", col_names = F)
}

# for transcripts typing
write_interpro_job = function(transcript_ids, output_path){
  require(dplyr)
  require(readr)
  transcript_ids = tibble(id = transcript_ids)
  write_delim(transcript_ids, path = output_path, col_names = F)
}

# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}