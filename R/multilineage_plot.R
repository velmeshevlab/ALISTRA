monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(plot.title = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.y = element_line(size=1, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' @export
compress <- function(df, window = 3, step){
df.comp = SlidingWindow("mean", df, window, step)
}

#' @export
fit.m3 <- function(exp.sel, pt, max.pt, model = "expression ~ splines::ns(pseudotime, df=3)", N = 500){
require(speedglm)
family = stats::quasipoisson()
exp_data.sel = cbind(pt, exp.sel)
colnames(exp_data.sel) <- c("pseudotime","expression")
exp_data.sel = as.data.frame(exp_data.sel)
exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
colnames(d) <- c("pseudotime")
fit = stats::predict(fit, newdata=d, type="response")
return(fit)
}, error=function(cond) {return(rep("NA", N))})
}

#' @export
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#' @export
compress_lineage <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F){
exp = compress_expression(cds, lineage, start=start, gene = gene, N = N, cores = cores)
input = paste0("cds@expression$", lineage, " <- exp$expression")
eval(parse(text=input))
input = paste0("cds@expectation$", lineage, " <- exp$expectation")
eval(parse(text=input))
return(cds)
}

#' @export
compress_expression <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F){
if(cores != F){
cl <- makeCluster(cores)
clusterEvalQ(cl, c(library(evobiR)))
}
input = paste0("get_lineage_object(cds, '", lineage, "',", start, ")")
cds_subset = eval(parse(text=input))
family = stats::quasipoisson()
model = "expression ~ splines::ns(pseudotime, df=3)"
names(cds_subset) <- rowData(cds_subset)$gene_short_name
exp = as.data.frame(as_matrix(exprs(cds_subset)))
exp = t(t(exp) /  pData(cds_subset)[, 'Size_Factor'])
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt <- as.data.frame(pt)
colnames(pt) <- c("pseudotime")
#exp = cbind(pt, round(t(exp)))
exp = cbind(pt, t(exp))
exp = exp[order(exp$pseudotime),]
step = ((nrow(exp)-3)/N)
pt = exp[,"pseudotime"]
#use sliding window to compress expression values and pseudotime
pt.comp = SlidingWindow("mean", pt, 3, step)
max.pt = max(pt.comp)
if(gene != F){
exp = exp[,c("pseudotime", gene)]
exp.comp = compress(exp[,gene], 3, step = step)
}
else{
print(paste0("Compressing lineage ", lineage, " and fitting curves"))
mat <- as.data.frame(exp[,2:ncol(exp)])
if(cores != F){
exp.comp = pbsapply(mat, compress, step = step, cl = cl)
}
else{
exp.comp = pbsapply(mat, compress, step = step)
}
}
if(gene != F){
exp_data.sel = cbind(pt.comp, exp.comp)
exp_data.sel = as.data.frame(exp_data.sel)
colnames(exp_data.sel) <- c("pseudotime", "expression")
exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
colnames(d) <- c("pseudotime")
tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
fit = predict(fit, newdata = d, type='response')
}, error=function(cond) {fit = as.data.frame(rep(0, N))})
exp = as.data.frame(cbind(exp.comp, fit, d))
colnames(exp) <- c("expression", "expectation", "pseudotime")
exp$expression[exp$expression < 0] <- 0
exp$expectation[exp$expectation < 0] <- 0
}
else{
mat <- as.data.frame(exp.comp)
if(cores != F){
fit = pbsapply(mat, fit.m3, pt = pt.comp, max.pt = max.pt, N = N, cl = cl)
}
else{
fit = pbsapply(mat, fit.m3, pt = pt.comp, max.pt = max.pt, N = N)
}
fit = apply(fit, 2, as.numeric)
return(list("expression" = exp.comp, "expectation" = fit))
}
exp$expression[exp$expression < 0] <- 0
exp$expectation[exp$expectation < 0] <- 0
if(cores != F){
stopCluster(cl)
}
return(exp)
}

#' @export
ridge_plot <- function(cds, cds1, gene, lineages_1, lineages_2, N = 500, colors1 = c("red", "blue", "green"), colors2 = c("coral", "cornflowerblue", "darkseagreen2"), label = F){
exps_1 = list()
means_1 = c()
for(lineage in lineages_1){
exp = compress_expression(cds, lineage, gene = toupper(gene), N)
exp = exp$expectation
mean = mean(quantile(exp, probs = seq(0, 1, by= 0.01))[6:96])
means_1 = append(means_1, mean)
exps_1[[lineage]] = exp
}
exps_2 = list()
means_2 = c()
for(lineage in lineages_2){
exp = compress_expression(cds1, lineage, gene = gene, N)
exp = exp$expectation
mean = mean(quantile(exp, probs = seq(0, 1, by= 0.01))[6:96])
means_2 = append(means_2, mean)
exps_2[[lineage]] = exp
}
max1 = max(means_1)
max2 = max(means_2)
factor = max1/max2
ymax1 = max(unlist(exps_1))
ymax2 = max(unlist(lapply(exps_2, "*" , factor)))
ymax = log10(max(ymax1, ymax2)+1)
dd = cbind(seq(from=0, to=25, by = (25/N)+0.0001), log10(as.data.frame(exps_1[1])+1))
dd = as.data.frame(dd)
colnames(dd) <- c("pseudotime", "exp")
q1 <- ggplot(data = dd) + geom_ridgeline(aes(x = pseudotime, y = rep(0, N), height = exp), fill = colors1[1], color = "black", size = I(1.2)) + ylab(names(exps_1)[1]) + ylim(0, ymax) + theme_opts() + theme(axis.line.x = element_line(size=1, color="black"))
qs = list()
if(length(exps_1) > 1){
M = 1
rest_exp1 = exps_1[2:length(exps_1)]
for(exp in rest_exp1){
dd = cbind(seq(from=0, to=25, by = (25/N)+0.0001), log10(as.data.frame(exp)+1))
dd = as.data.frame(dd)
colnames(dd) <- c("pseudotime", "exp") 
q <- ggplot(data = dd) + geom_ridgeline(aes(x = pseudotime, y = rep(0, N), height = exp), fill = colors1[M+1], color = "black", size = I(1.2)) + ylab(names(rest_exp1)[M]) + ylim(0, ymax) + theme_opts()
qs[[names(rest_exp1)[M]]] <- q
M = M+1
}
}
M = 1
for(exp in exps_2){
dd = cbind(seq(from=0, to=25, by = (25/N)+0.0001), log10(as.data.frame(exp*factor)+1))
dd = as.data.frame(dd)
colnames(dd) <- c("pseudotime", "exp") 
q <- ggplot(data = dd) + geom_ridgeline(aes(x = pseudotime, y = rep(0, N), height = exp), fill = colors2[M], color = "black", size = I(1.2)) + ylab(names(exps_2)[M]) + ylim(0, ymax) + theme_opts()
qs[[names(exps_2)[M]]] <- q
M = M+1
}
input = paste0("ggarrange(")
for(M in 0:(length(qs)-1)){
input = paste0(input, 'qs[["', names(qs)[length(qs)-M],'"]], ')
}
if(label == TRUE){
input = paste0(input, 'q1, labels = c(rev(names(qs)), names(exps_1)[1]), ncol = 1)')
}
else{
input = paste0(input, 'q1, ncol = 1)')
}
input = noquote(input)
figure <- eval(parse(text=input))
figure
}

#' @export
get_pt <- function(lineage){
load(file = paste0(lineage, ".R"))
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt = pt[order(pt)]
return(pt)
}

#' @export
plot_multiple <- function(cds, gene, lineages, start, colors = c("red", "blue", "green", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), colors.mids = c("yellow", "deepskyblue", "darkolivegreen1"), N = 500, legend_position = "right"){
cds_name = deparse(substitute(cds))
input = paste0("get_lineage_object(",cds_name,", '", lineage, "',", start, ")")
cds_subset = eval(parse(text=input))
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt <- as.data.frame(pt)
pt = pt[,1]
step = ((length(pt)-3)/N)
pt.comp = SlidingWindow("mean", pt, 3, step)
max.pt = max(pt.comp)
dd = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
cols = c("pseudotime")
fits = c()
for(lineage in lineages){
input = paste0("exp = ",cds_name,"@expression$", lineage,"[,'",gene,"']")
eval(parse(text=input))
input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
eval(parse(text=input))
dd = cbind(dd, exp, fit)
cols = append(cols, paste0("exp_", lineage))
cols = append(cols, paste0("fit_", lineage))
fits = c(fits, fit)
}
colnames(dd) <- cols
ymax = max(fits)
q <- ggplot(data = dd)
for(N in 1:length(lineages)){
loop_input1 = paste0("geom_point(aes_string(x='pseudotime',y = '", paste0('exp_', lineages[N]), "',color='pseudotime'),size=I(1))")
loop_input2 = paste0("scale_color_gradient2(lineages[N],low='grey',mid='",colors.mids[N],"',", "high='",colors[N],"')")
loop_input3 = "new_scale_color()"
loop_input4 = paste0("geom_line(aes_string(x='pseudotime', y = '", paste0('fit_', lineages[N]), "',size = I(1.2)), color = '", colors[N],"')")
q <- q + eval(parse(text=loop_input1)) + eval(parse(text=loop_input2)) + eval(parse(text=loop_input3)) + eval(parse(text=loop_input4))
}
q <- q + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
q <- q + ylim(y = c(0,ymax))
q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) + theme(plot.title = element_text(size = 36, face="bold", hjust = 0.5), axis.text=element_text(size=28), axis.title=element_text(size=28,face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=24, face = "bold"), legend.position = legend_position)
q
}
