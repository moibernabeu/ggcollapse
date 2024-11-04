#' Plot a collapsed tree
#' 
#' This function takes the nodes to collapse, prunes all the leaves except the
#' longest and the shortest, collapse them and assign some colours and names
#' according to the arguments given.
#' 
#' @param tree Phylogenetic tree.
#' @param nodes Named vector conaining the nodes to collapse.
#' @param collapse_mode String to select the aesthetics of the collapsing triangle, 'mixed', 'min' or 'max'.
#' @param node_colours Vector named with the same names as the nodes with the colour codes.
#' @param get_sp Function to extract the species from the tip.label.
#' @param tree_data Dataframe with the data linked to each species or sequence.
#' @param data_sp_column Column of the tree_data containing the species code.
#' 
#' @return Collapsed tree ggplot plot object.
#' @export
#' @author Moisès Bernabeu
#' @example man/examples/ggcollapse.R
ggcollapse <- function(tree, nodes, collapse_mode='mixed', node_colours=NULL,
                       get_sp=NULL, tree_data=NULL, data_sp_column=NULL) {
  
  if (!methods::is(tree, 'phylo')) {
    stop('The tree must be a phylo object.')
  }
  
  # Get distances from the tree
  dists <- ape::dist.nodes(tree)
  
  # Annotate tree
  to_drop <- c()
  kept <- c()
  ntips <- c()
  for (node in 1:length(nodes)) {
    # Getting descendants
    desc <- phytools::getDescendants(tree, nodes[node])
    
    # Getting just descendant tips
    desc <- desc[which(desc <= ape::Ntip(tree))]
    
    # Retrieving the number of tips in the collapsed node
    ntips[node] <- length(desc)
    
    # Getting the distances of the node to be collapsed
    ndists <- dists[nodes[node], desc]
    
    # Getting nodes to drop and to keep
    to_drop_nnode <- as.numeric(names(ndists[-c(which.min(ndists), which.max(ndists))]))
    to_keep_nnode <- as.numeric(names(ndists[c(which.min(ndists), which.max(ndists))]))
    
    # Getting tip labels of the nodes to drop
    to_drop <- c(to_drop, tree$tip.label[to_drop_nnode])
    
    # Getting the tip labels 
    kept[[node]] <- c(tree$tip.label[to_keep_nnode])
  }
  
  # Prunning the tree to get the vertically symmetric tree
  symm_tree <- treeio::drop.tip(tree, to_drop)
  
  # Adjusting the branch lengths to be able to scale the tree horizontally
  # and fit it in the final plot.
  symm_tree$edge.length <- symm_tree$edge.length / max(diag(ape::vcv(symm_tree)))
  
  # Annotating if all the proper data to annotate the tree is given
  if (!is.null(tree_data) & !is.null(get_sp) & !is.null(data_sp_column)) {
    symm_tree <- annotate_tree(symm_tree, tree_data, get_sp, data_sp_column)
    symm_tree_phylo <- symm_tree@phylo
  } else {
    symm_tree_phylo <- symm_tree
  }
  
  new_nodes <- c()
  offset <- c()
  vjust <- c()
  hjust <- c()
  dists <- ape::dist.nodes(symm_tree_phylo)
  for (node in 1:length(nodes)) {
    new_nodes[node] <- ape::getMRCA(symm_tree_phylo, kept[[node]])
    desc <- phytools::getDescendants(symm_tree_phylo, new_nodes[node])
    desc <- desc[which(desc <= ape::Ntip(symm_tree_phylo))]
    ndists <- dists[new_nodes[node], desc]
    if (collapse_mode == 'max') {
      offset[node] <- max(ndists)
      vjust[node] <- 0.5
      hjust[node] <- -0.015
    } else if (collapse_mode == 'min') {
      offset[node] <- min(ndists)
      vjust[node] <- 0.5
      hjust[node] <- -0.015
    } else if (length(unique(round(ndists, digits = 1))) == 1) {
      offset[node] <- mean(ndists)
      vjust[node] <- 1
      hjust[node] <- -0.015
    } else if (collapse_mode == 'mixed') {
      offset[node] <- 0.85 * max(ndists)
      vjust[node] <- 1
      hjust[node] <- 0
    }
  }
  names(new_nodes) <- names(nodes)
  names(offset) <- names(nodes)
  names(vjust) <- names(nodes)
  names(hjust) <- names(nodes)
  
  p <- ggtree::ggtree(symm_tree, ladderize = TRUE, right = TRUE)
  
  if (is.null(node_colours)) {
    node_colours <- rep('grey80', length(nodes))
    names(node_colours) <- names(nodes)
    nodelab_colours <- rep('black', length(nodes))
    names(nodelab_colours) <- names(nodes)
  } else {
    nodelab_colours <- node_colours
  }
  
  for (node in 1:length(new_nodes)) {
    p <- ggtree::collapse(p, new_nodes[node], mode = collapse_mode,
                          fill = node_colours[names(nodes)[node]],
                          colour = 'black')
  }
  
  for (node in 1:length(new_nodes)) {
    lab = paste0(names(nodes)[node], ' (', ntips[node], ')')
    p <- p + ggtree::geom_cladelab(new_nodes[node],
                                   label = lab,
                                   offset.text = offset[node],
                                   hjust = hjust[node],
                                   vjust = vjust[node],
                                   textcolour = nodelab_colours[names(nodes)[node]])
  }
  
  return(p)
}
