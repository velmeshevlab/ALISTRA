# ANMOLIB
Analysis of Monocle Lineage Branches
Analysis of Monocle Lineage Branches

This package is an extension of monocle3 (https://github.com/cole-trapnell-lab/monocle3) that allows isolating trajecory branches corresonding to individual lineages, identifying lineage-specific dynamically expressed genes, classifying genes based on their expression pattern and plotting gene expression along multiple lineages.

# Tutorial
## Monocle object import and isolating lineages.
1) After building trajectories, import Monocle object:
> cds = import_monocle(cds)

2) Displate nodes of the trajectory graph:
> node_plot(cds)

If the nodes are too dense to see, filter = T option can be used to skip non-terminal and non-branching nodes, and the node density can be set with the N. Larger N skips more nodes.

3) If necessary, connect nodes that you believe should be part of one branch based on biological assumptions, but were not joined by Monocle:

> cds = connect_nodes(cds, "Y_790", "Y_550")

4) 
