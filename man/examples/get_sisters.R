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

ggtree(atree) +
  geom_nodelab(aes(label = node), geom = 'label') +
  geom_tiplab()

# Getting sisters to node 36
get_sisters(tree, 36)

# Getting sisters to t7
get_sisters(atree, 't7')

# Getting sisters to node 36 annotated with the percentage of tips with a
# certain value of the given naming_column
get_sisters(atree, 36, naming_column = 'group')