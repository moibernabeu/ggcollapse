library(treeio)
library(ggtree)

get_sp <- function(header) {
  return(header)
}

set.seed(2024-11-04)
tree <- rtree(20)

tree_data <- data.frame(tip = tree$tip.label,
                        group = c(rep('A', 10), rep('B', 10)))
atree <- annotate_tree(tree = tree,
                       tree_data = tree_data,
                       get_sp = get_sp,
                       data_sp_column = 1)

# Collapsing monophyletics of the group
mphy <- get_monophyletics(atree, 'group')
cols <- c('A' = 'steelblue', 'B' = 'darkorange')

ggcollapse(tree, mphy, collapse_mode = 'mixed', node_colours = cols)

# Collapsing sisters
sist <- get_sisters(atree, 36, 'group')

ggcollapse(tree, sist, collapse_mode = 'mixed') +
  hexpand(0.15)