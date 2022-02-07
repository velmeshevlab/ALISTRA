#' @export
annotate_gene_peaks_sub <- function(gene, d, cells, age){
exp_age = cbind(d[gene,cells], age[cells,])
colnames(exp_age)[1] <- "exp"
res = exp_age[exp_age$exp >= quantile(exp_age$exp, 0.99),]
res = c(Mode(res$age_range), mean(res$age_num))
names(res) <- c("age_range", "age_num")
res
}

#' @export
annotate_gene_peaks <- function(data, cds, genes, lineage, age){
d = GetAssayData(object = data, assay = "RNA", slot = "data")
cds_name = deparse(substitute(cds))
input = paste0("cells = ",cds_name,"@lineages$", lineage)
eval(parse(text=input))
res = pbsapply(genes, annotate_gene_peaks_sub, d = d, cells = cells, age = age)
res
}

#' @export
AUC_window_sub <- function(gene, cds, lineage, comp_lineages, factor, window_ratio){
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  fit.sel = fit[,gene]
  N = length(fit.sel)
  window = N*window_ratio
  out = c()
  top = 0
  mat = matrix(nrow = N, ncol = 0, )
  for(lin in comp_lineages){
    cds_name = deparse(substitute(cds))
    input = paste0("fit = ",cds_name,"@expectation$", lin)
    eval(parse(text=input))
    if(gene %in% colnames(fit) & !(is.na(fit[,gene]))){
    fit = fit[,gene]
    mat = cbind(mat, fit)
    top = max(fit, top)
    }
  }
      res = c()
      score = c()
      for(i in 1:N){
          if(fit.sel[i] > top*factor){
            res = append(res, 1)
          }
          else{
            res = append(res, 0)
          }
          score = append(score, fit.sel[i]-max(mat[i,]))
        }
      target = paste(res,collapse="")
      res = rle(res)$lengths[rle(res)$values==1]
      if(any(res > window))
      {
        search = paste(rep(1, max(rle(res)$values)),collapse="")
        index = gregexpr(pattern = search, target)[[1]][1]
        score = sum(score[index:(index+max(rle(res)$values))-1])/top
        score
      }
      else{
        FALSE
      }
    }

#' @export
AUC_window_sub_complex <- function(gene, cds, lineage, comp_lineages, factor, window_ratio){
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  fit.sel = fit[,gene]
  N = length(fit.sel)
  window = N*window_ratio
  out = c()
  res = c()
      for(i in 1:N){
        top = 0
        for(lin in comp_lineages){
          cds_name = deparse(substitute(cds))
          input = paste0("fit = ",cds_name,"@expectation$", lin)
          eval(parse(text=input))
          if(gene %in% colnames(fit) & !(is.na(fit[,gene]))){
            fit = fit[,gene]
            top = max(fit[i], top)
          }
        }
          if(fit.sel[i] > top*factor){
            res = append(res, TRUE)
          }
          else{
            res = append(res, FALSE)
          }
        }
      res = rle(res)$lengths[rle(res)$values==TRUE]
      if(any(res > window))
      {
        TRUE
      }
      else{
        FALSE
      }
    }

#' @export
get_dist <- function(genes, lineage, lineages, dist)
out = sapply(genes, get_dist_sub, lineage = lineage, lineages = lineages, dist = dist)

get_dist_sub <- function(gene, lineage, lineages, dist){
gene.name = paste0(lineage, "__", gene)
dist.sel = dist[gene.name,]
dist.sel = dist.sel[names(dist.sel) %in% paste0(lineages, "__", gene)]
names(dist.sel) <- gsub(paste0("__", gene), "", names(dist.sel))
if(length(dist.sel) != length(lineages)){
dist.sel = rep(0, length(lineages))
names(dist.sel) <- lineages
}
dist.sel
}

#' @export
get_Is <-  function(res, lineages, action){
out = matrix(nrow = length(res), ncol = 0,)
for(lin in lineages){
load(file = paste("Moran_", lin,".R", sep = ""))
Is = cds_pr_test_res$morans_I
names(Is) <- rownames(cds_pr_test_res)
I = sapply(res, lookup_I, I = Is, action = 0)
out = cbind(out, I)
}
colnames(out) <- lineages
rownames(out) <- res
out
}

#' @export
get_lineage_genes <- function(lineage, lineages, exlude_lineages = F, factor = 4, high_I = 0.1){
load(file = paste("Moran_", lineage,".R", sep = ""))
lin.genes = cds_pr_test_res[cds_pr_test_res$morans_I >= high_I & cds_pr_test_res$q_value < 0.05, ]
lin.genes = lin.genes[,c("morans_I", "q_value")]
test_lineages = lineages[!(lineages == lineage)]
if(exlude_lineages != F){
test_lineages = test_lineages[!(test_lineages %in% exlude_lineages)]
}
Is = get_Is(rownames(lin.genes), test_lineages, action = 0)
out = c()
for(i in 1:nrow(lin.genes)){
test = lin.genes$morans_I[i]
d = Is[i,]
res = (d* factor) < test
res = sum(res == TRUE)
if(res == length(test_lineages)){
out = append(out, rownames(Is)[i])
}
}
out
}

#' @export
get_lineage_genes_multiple <- function(test_lineages, lineages, factor = 4, high_I = 0.1){
res = get_lineage_genes(test_lineages[1], lineages, exlude_lineages = test_lineages, high_I = high_I, factor = factor)
for(test_lineage in test_lineages[2:length(test_lineages)]){
local_res = get_lineage_genes(test_lineage, lineages, exlude_lineages = test_lineages, high_I = high_I, factor = factor)
res = intersect(res, local_res)
}
out = get_Is(res, lineages, action = 0)
out
}

prep_mat <- function(gene, mat){
name = paste0(lineage, "__", gene)
cols = colnames(mat)[grepl(paste0("__", gene), colnames(mat))]
res = mat[name,cols]
unlist(res)
}

#' @export
get_distance <- function(lineage, genes){
d = read.table("adult_dist_mat.txt", sep="\t",header=TRUE,row.names=1,check.names=FALSE)
res = sapply(genes, prep_mat, mat = d)
res = t(res)
colnames(res) <- gsub(paste0("__", genes[1]), "", colnames(res))
res
}

lookup_I <- function(gene, I, action = "NA"){
if(gene %in% names(I)){
res = I[gene]
}
else{
res = action
}
res
}

#' @export
get_Moran_I <- function(lineages, genes){
out = matrix(nrow = length(genes), ncol = 0)
for(lineage in lineages){
load(file = paste("Moran_", lineage,".R", sep = ""))
I = cds_pr_test_res$morans_I
names(I) <- rownames(cds_pr_test_res)
res = sapply(genes, lookup_I, I = I)
out = cbind(out, res)
}
colnames(out) <- lineages
out
}

#' @export
calc_dist <- function(lineage1, lineage2, pt_genes = FALSE, dist_method = "CORT", q = 0.05, I = 0.15, adjust = F, cores = 1){
if(pt_genes == FALSE){
pt_genes = c()
for(lineage in c(lineage1, lineage2)){
d = read.table(paste0("pt_DGE_", lineage,".txt"), sep = "\t", header = T)
d = d[d$q_value < q & d$morans_I >= I,]
if(length(pt_genes) == 0){
pt_genes = toupper(d$gene_short_name)
}
else{
pt_genes = union(pt_genes, toupper(d$gene_short_name))
}
}
}
input = paste0("fit = cds@expectation$", lineage1)
eval(parse(text=input))
colnames(fit) <- toupper(colnames(fit))
fit1 = fit[,colnames(fit)%in%pt_genes]
input = paste0("fit = cds@expectation$", lineage2)
eval(parse(text=input))
colnames(fit) <- toupper(colnames(fit))
fit2 = fit[,colnames(fit)%in%pt_genes]
pt_genes = intersect(pt_genes, intersect(colnames(fit1), colnames(fit2)))
fit1 = fit1[,pt_genes]
fit2 = fit2[,pt_genes]
fit1 = t(fit1)
fit2 = t(fit2)
fit1 = log10(fit1+1)
fit2 = log10(fit2+1)
fit1.list <- split(fit1, seq(nrow(fit1)))
names(fit1.list) <- rownames(fit1)
fit2.list <- split(fit2, seq(nrow(fit2)))
names(fit2.list) <- rownames(fit2)
if(dist_method == "CORT"){
out = pbmapply(CortDistance, fit1.list, fit2.list, deltamethod="DTW")
}
if(dist_method == "ACF"){
out = pbmapply(ACFDistance, fit1.list, fit2.list)
}
if(dist_method == "DTW"){
out = pbmapply(ACFDistance, fit1.list, fit2.list)
}
return(out)
}

#' @export
calc_dist_lineages <- function(lineages, pt_genes = FALSE, dist_method = "CORT", suffix = "", q = 0.05, I = 0.15){
if(pt_genes == FALSE){
pt_genes = c()
for(lineage in lineages){
d = read.table(paste0("pt_DGE_", lineage,".txt"), sep = "\t", header = T)
d = d[d$q_value < q & d$morans_I >= I,]
if(length(pt_genes) == 0){
pt_genes = toupper(d$gene_short_name)
}
else{
pt_genes = union(pt_genes, toupper(d$gene_short_name))
}
}
}
input = paste0("fit = cds@expectation$", lineages[1])
eval(parse(text=input))
fit.comb = matrix(ncol = 0, nrow = nrow(fit), 0)
for(lineage in lineages){
input = paste0("fit = cds@expectation$", lineage)
eval(parse(text=input))
colnames(fit) <- toupper(colnames(fit))
fit = fit[,colnames(fit)%in%pt_genes]
colnames(fit) <- paste(lineage, colnames(fit), sep = "_")
fit.comb = cbind(fit.comb, fit)
}
fit.comb = apply(fit.comb, 2, as.numeric)
fit.comb = as.data.frame(fit.comb)
d = log10(t(fit.comb)+1)
if(dist_method == "CORT"){
dis = TSclust::diss(d, dist_method, deltamethod="DTW")
}
else{
dis = TSclust::diss(d, dist_method)
}
return(list("dis" = dis, "d" = d))
}

#' @export
get_pt_exp <-function(cds, lineages, pt_genes = FALSE, q = 0.05, I = 0.15){
if(pt_genes == FALSE){
pt_genes = c()
for(lineage in lineages){
d = read.table(paste0("pt_DGE_", lineage,".txt"), sep = "\t", header = T)
d = d[d$q_value < q & d$morans_I >= I,]
if(length(pt_genes) == 0){
pt_genes = toupper(d$gene_short_name)
}
else{
pt_genes = union(pt_genes, toupper(d$gene_short_name))
}
}
}
cds_name = deparse(substitute(cds))
input = paste0("fit = ",cds_name,"@expectation$", lineages[1])
eval(parse(text=input))
fit.comb = matrix(ncol = 0, nrow = nrow(fit), 0)
for(lineage in lineages){
input = paste0("fit = ",cds_name,"@expectation$", lineage)
eval(parse(text=input))
colnames(fit) <- toupper(colnames(fit))
fit = fit[,colnames(fit)%in%pt_genes]
colnames(fit) <- paste(lineage, colnames(fit), sep = "__")
fit.comb = cbind(fit.comb, fit)
}
fit.comb = apply(fit.comb, 2, as.numeric)
fit.comb = as.data.frame(fit.comb)
d = log10(t(fit.comb)+1)
}

#' @export
clust_lineages <- function(d, dis, k, lineage_names, colors, suffix = "", myCol = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",  "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",  "orange", "tomato", "coral2", "khaki1", "violetred", "red2",  "springgreen2", "yellowgreen", "palegreen4",  "wheat2", "tan", "tan2", "tan3", "brown",  "grey70", "grey50", "grey30", "aquamarine", "bisque3", "cornflowerblue", "darkseagreen1", "darkred", "lightgreen","hotpink"), method = "ward.D2"){
tree = hclust(dis, method = method)
png(filename = paste0("tree_",suffix,".png"), width = 3500, height = 2080, bg = "white",  res = 300);
plot(tree)
rect.hclust(tree, k = k, border = 2:10)
dev.off();
clust = cutree(tree, k = k)
out = matrix(nrow = length(clust)/length(lineage_names),ncol = 0, 0)
for(lineage_name in lineage_names){
clust.lin = clust[grepl(paste0(lineage_name,"_"), names(clust))]
out = cbind(out, as.character(clust.lin))
}
rownames(out) = gsub(paste0(lineage_name, "_"), "", names(clust.lin))
colnames(out) <- lineage_names
look = cbind(unique(clust[tree$order]), myCol[1:k])
new <- out
new[] <- look[,2][match(unlist(out), look[,1])]
write.table(new, paste0("dynamic_clust_",suffix,".txt"), sep = "\t", quote = F)
hcd = as.dendrogram(tree)
cols = rownames(d)
N = 1
for(lineage_name in lineage_names){
cols = gsub(paste0(lineage_name, "_.+"), colors[N], cols)
N = N+1
}
png(filename = paste0("dendo_L_",suffix,".png"), width = 3500, height = 2080, bg = "white",  res = 300);
hcd = as.dendrogram(tree)
labels(hcd) <- c(1:nrow(d))
hcd %>% set("branches_lwd", 5) %>% plot(horiz = F, leaflab = "none")
colored_bars(cols, hcd, horiz = F, y_shift = -40, y_scale = 200)
hcd %>% rect.dendrogram(cluster = clust, k = k, horiz = F, lwd = 5, lower_rect = -280, text = unique(clust[tree$order]), border = myCol[1:k])
dev.off()
png(filename = paste0("dendo_",suffix,".png"), width = 3500, height = 2080, bg = "white",  res = 300);
hcd = as.dendrogram(tree)
labels(hcd) <- c(1:nrow(d))
hcd %>% set("branches_lwd", 5) %>% plot(horiz = F, leaflab = "none")
colored_bars(cols, hcd, horiz = F, y_shift = -40, y_scale = 200)
hcd %>% rect.dendrogram(cluster = clust, k = k, horiz = F, lwd = 5, lower_rect = -280, border = myCol[1:k])
dev.off()
}

#' @export
average_curves <- function(d, clust, colors = c("red", "green", "blue", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan", "blanchedalmond", "coral")){
clusters = unique(clust)
N = ncol(d)
M = 1
for(cluster in clusters){
color = colors[M]
print(cluster)
M = M+1
sel = d[names(clust)[clust == cluster],]
pt = seq(from=0, to=25, by = (25/N)+0.0001)
fit = fit.m3(colMeans(sel), pt, max(pt))
dd = cbind(pt, fit)
dd = as.data.frame(dd)
colnames(dd) <- c("pseudotime", "average")
q <- ggplot(data = dd)
q <- q + geom_line(aes(x = pseudotime, y = average, size = I(5)), color = color)
q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime")
q <- q + ylim(y = c(0,max(fit))) + theme(plot.title = element_text(size = 28, face="bold", hjust = 0.5), axis.text=element_text(size=36), axis.title=element_text(size=40,face="bold"), legend.position = "none")
q <- q + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
q
ggsave(file=paste(cluster, ".png",sep=""),width = 8, height = 6, units = "in",  dpi = 600);
}
}

