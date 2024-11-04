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

mphy <- get_monophyletics(atree, 'group')

mphy

ggtree(atree) +
  geom_nodelab(aes(subset = node %in% mphy, label = node), geom = 'label') +
  geom_tiplab(aes(label = group))