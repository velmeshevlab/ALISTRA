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
clust_lineages_DTW <- function(cds, lineages, k, pt_genes = FALSE, q = 0.05, I = 0.15){
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
colnames(fit) <- paste(lineage, colnames(fit), sep = "_")
fit.comb = cbind(fit.comb, fit)
}
fit.comb = apply(fit.comb, 2, as.numeric)
fit.comb = as.data.frame(fit.comb)
d = log10(t(fit.comb)+1)
print(paste0("Clustering ", length(pt_genes), " genes across ", length(lineages), " lineages"))
res = tsclust(d, k = k)
res
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
average_curves <- function(d, dis, k, colors = c("red", "green", "blue", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), N = 500){
tree = hclust(dis, method = "ward.D2")
clust = cutree(tree, k = k)
clusters = unique(clust[tree$order])
M = 1
for(cluster in clusters){
color = colors[M]
print(color)
M = M+1
sel = d[names(clust[clust == cluster]),]
sel = (10^sel)-1
dd = cbind(seq(from=0, to=25, by = (25/N)+0.0001), log10(colMeans(sel)+1))
dd = as.data.frame(dd)
colnames(dd) <- c("pseudotime", "average")
q <- ggplot(data = dd)
q <- q + geom_line(aes(x = pseudotime, y = average, size = I(5)), color = color)
q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime")
q <- q + ylim(-0.02, max(dd$average)) + theme(plot.title = element_text(size = 28, face="bold", hjust = 0.5), axis.text=element_text(size=36), axis.title=element_text(size=40,face="bold"), legend.position = "none")
q
ggsave(file=paste(color, ".png",sep=""),width = 8, height = 6, units = "in",  dpi = 600);
}
}
