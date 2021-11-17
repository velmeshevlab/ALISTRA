setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list")) -> cell_data_set_ext

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

combine_lineages <- function(cds, start){
lineage = names(cds@lineages)[1]
input = paste0("principal_graph(cds@lineages$", lineage,")[['UMAP']]")
for(lineage in names(cds@lineages)[2:length(names(cds@lineages))]){
input = paste0(input,",principal_graph(cds@lineages$", lineage,")[['UMAP']]")
}
input = paste0("union(", input, ")")
g = eval(parse(text=input))
nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
principal_graph(cds)[["UMAP"]] <- g
cds@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(g))]
cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(g))]))
closest_vertex = as.data.frame(closest_vertex)
cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
cds <- project2MST(cds, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(g))])
cds <- order_cells(cds, root_pr_nodes = as.character(start))
return(cds)
}

node_plot <- function(cds, filter = F, N = 50){
Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
d = as.data.frame(t(Y))
if(filter == T){
g = principal_graph(cds)[["UMAP"]]
dd = degree(g)
names1 = names(dd[dd > 2 | dd == 1])
names2 = names(dd[dd == 2])
names2 = sample(names2, length(names2)/N, replace = F)
d.f = d[c(names1, names2),]
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d.f, aes(x=UMAP_1, y=UMAP_2), label=rownames(d.f), size=0.3, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
}
else{
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d, aes(x=UMAP_1, y=UMAP_2), label=rownames(d), size=0.3, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
}
}

path.distance <- function(path){
dists=c()
for(i in 2:nrow(path)){
x1 = path[i-1,1]
y1 = path[i-1,2]
x2 = path[i,1]
y2 = path[i,2]
d.x = x2 - x1
d.y = y2 - y1
dist = sqrt(d.x*d.x + d.y*d.y)
dists = append(dist,dists)
}
return(mean(dists))
}

cell.selector_sub2 <- function(cell, coords, r){
x2 = cell[1]
y2 = cell[2]
d.x = x2 - coords[1]
d.y = y2 - coords[2]
dist = sqrt(d.x*d.x + d.y*d.y)
if(dist <= r){
return(TRUE)
}
else{
return(FALSE)
}
}

selector_sub <- function(node, cells, r){
x1 = node[1]
y1 = node[2]
res = apply(cells, 1, cell.selector_sub2, coords = c(x1, y1), r = r, simplify = T)
res = names(res[res == TRUE])
return(res)
}

cell.selector <- function(path, cells, r, cl){
sel.cells = c()
sel.cells = pbapply(path, 1, selector_sub, cells = cells, r = r, cl = cl, simplify = T)
return(unique(unlist(sel.cells)))
}

make_graph <- function(sub.graph){
edges = names(sub.graph)
start.edges = c()
end.edges = c()
for(i in 1:(length(edges)-1)){
start.edges = append(start.edges, edges[i])
end.edges = append(end.edges, edges[i+1])
}
d = cbind(start.edges, end.edges)
g = graph_from_data_frame(d, directed = F)
return(g)
}

included <- function(graph, include_nodes){
all(include_nodes %in% names(graph))
}

isolate_graph <- function(cds, start, end, lineage, include_nodes = F){
#get lineage graph
reduction_method = "UMAP"
graph = cds@principal_graph[[reduction_method]]
#select cells that are 1) progenitor cells from the region of interest (MGE, CGE) or 2) lineage-committed cells
sub.graph = all_simple_paths(graph, start, end)
if(include_nodes != F){
sub.graph = sub.graph[sapply(sub.graph, included, include_nodes = include_nodes)]
}
lengths = lengths(sub.graph)
#get the shortest path
n = which(lengths==min(lengths))[1]
sub.graph = sub.graph[[n]]
input = paste0("cds@graphs$", lineage, " <- make_graph(sub.graph)")
eval(parse(text=input))
return(cds)
}

#this is an "improved" lineage isolation function that recalculates pseudotime based on selected branch
isolate_lineage <- function(cds, lineage, sel_clusters = F, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
input = paste0("sub.graph = cds@graphs$", lineage)
eval(parse(text=input))
nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
if(subset == F){
nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
}
else{
g = principal_graph(cds)[["UMAP"]]
dd = degree(g)
names1 = names(dd[dd > 2 | dd == 1])
names2 = names(dd[dd == 2])
names2 = sample(names2, length(names2)/subset, replace = F)
names = c(names1, names2)
names = intersect(names(V(sub.graph)), names)
nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
}
#select cells along the graph
mean.dist = path.distance(nodes_UMAP.sub)
r = mean.dist*N
cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
#only keep cells in the progenitor and lineage-specific clusters
sel.cells1 = c()
sel.cells2 = sel.cells
if(starting_clusters != F){
sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
}
if(start_regions != F){
sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
}
if(sel_clusters != F){
sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
}
cells = unique(c(sel.cells1, sel.cells2))
sel.cells = sel.cells[sel.cells %in% cells]
#subset the moncole object
cds_subset = cds[,sel.cells]
#set the graph, node and cell UMAP coordinates
principal_graph(cds_subset)[["UMAP"]] <- sub.graph
cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
#recalculate closest vertex for the selected cells
cells_UMAP = as.data.frame(reducedDims(cds_subset)["UMAP"])
closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(sub.graph))]))
closest_vertex = as.data.frame(closest_vertex)
cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(sub.graph))])
cds_subset <- order_cells(cds_subset, root_pr_nodes = c(paste0("Y_", as.character(start))))
input = paste0("cds@lineages$", lineage, " <- cds_subset")
eval(parse(text=input))
return(cds)
}

calculate_closest_vertex <- function(cells, nodes){
new.pos = as.numeric(cells)
nearest.idx <- which.min(colSums((nodes - new.pos)^2))
out = as.integer(gsub("Y_", "", names(nearest.idx)))
}

connect_nodes <- function(cds, node1, node2){
graph.old = cds@principal_graph[["UMAP"]]
graph.new <- add_edges(graph.old, c(node1, node2))
cds@principal_graph[["UMAP"]] <- graph.new
return(cds)
}

import_monocle <-function(cds){
cds <- as(cds,"cell_data_set_ext")
return(cds)
}
