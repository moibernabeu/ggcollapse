#' Annotate the tree
#'
#' Function to annotate the tree with the data given in the tree information
#' dataframe 'data'. This dataframe must contain one column with the species
#' code matching the species code in the headers which is extracted by the
#' 'get_sp' function
#'
#' @param tree Phylogenetic tree imported with the treeio::read.tree().
#' @param tree_data Dataframe containing information for the species in the tree.
#' @param get_sp function to extract the species code from the tip.label.
#' @param data_sp_column Column containing the species code.
#'
#' @return A phylo4 tree object with a data frame with the annotations.
#' @export
#' @author Moisès Bernabeu
#' @example man/examples/annotate_tree.R
annotate_tree <- function(tree, tree_data, get_sp, data_sp_column) {
  if (is_tibble(tree_data)) {
    tree_data <- data.frame(tree_data)
  }
  row.names(tree_data) <- tree_data[, data_sp_column]

  ttree <- treeio::as_tibble(tree)
  seqids <- ttree$label[1:ape::Ntip(ttree)]
  taxids <- get_sp(seqids)

  tree_dat <- treeio::tibble(label = seqids, tree_data[taxids, ])
  otree <- treeio::full_join(ttree, tree_dat, by = 'label')
  otree <- treeio::as.treedata(otree)

  return(otree)
}

#' Check if a node is son
#'
#' A function to check wether a tip belongs to a specific node.
#'
#'
#' @param tree Phylogenetic tree.
#' @param parents A vector with a set of parent nodes where to check whether the node belongs to.
#' @param son Node to be tested.
#'
#' @return A boolean value about whether the son is in parents.
#' @export
#' @author Moisès Bernabeu
is_son <- function(tree, parents, son) {
  for (parent in parents) {
    if (son %in% phytools::getDescendants(treeio::as.phylo(tree), parent)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Get sisters
#'
#' Get the number of the ancestor of the sister nodes to a specified node.
#'
#' Args:
#' @param tree Phylogenetic tree.
#' @param node The number of name (in case of a tip) to start extracting sister nodes.
#' @param naming_column Column to name the sister groups based on the proportion of the column values among the children tips.
#'
#' @return integer: vector with the sister nodes.
#' @export
#' @author Moisès Bernabeu
#' @example man/examples/get_sisters.R
get_sisters <- function(tree, node, naming_column=NULL) {
  if (typeof(tree) == 'S4') {
    if (!is.null(naming_column)) {
      tree_data_condition <- TRUE
    } else {
      tree_data_condition <- FALSE
    }
    tree_data <- data.frame(treeio::as_tibble(tree)[, -c(1:3)])
    tree <- tree@phylo
  } else {
    tree_data_condition <- FALSE
    tree <- tree
  }
  if (is.character(node)) {
    node <- which(tree$tip.label == node)
  }
  curr_node <- 0
  sisters <- c()
  prev_sist <- phytools::getSisters(tree, node)
  curr_node <- ape::getMRCA(tree, c(prev_sist, node))
  while (length(prev_sist) == 1 & curr_node != treeio::rootnode(tree)) {
    curr_node <- ape::getMRCA(tree, c(prev_sist, node))
    if (!treeio::isTip(tree, prev_sist)) {
      sisters <- c(sisters, prev_sist)
      descendants <- phytools::getDescendants(tree, prev_sist)
      descendants <- tree$tip.label[descendants[which(descendants <= ape::Ntip(tree))]]
      if (tree_data_condition) {
        desc_summary <- table(tree_data[which(tree_data[, 1] %in% descendants), naming_column]) / length(descendants)
        clade_name <- names(which.max(desc_summary))
        clade_prop <- desc_summary[which.max(desc_summary)]
        names(sisters)[length(sisters)] <- paste(clade_name, ' ', round(clade_prop * 100, 0), '%', sep = '')
      }
    }
    prev_sist <- phytools::getSisters(tree, curr_node)
  }
  return(sisters)
}

#' Get monophyletic groups
#'
#' Get the monophyletic groups assigned to a specific column.
#'
#' @param tree Phylogenetic tree.
#' @param group_column Column to extract the monophyletic groups.
#' @param tree_data Dataframe containing information for the species in the tree.
#' @param get_sp Function to extract the species code from the tip.label.
#' @param data_sp_column Column containing the species code.
#'
#' @return Vector with the sister nodes.
#' @export
#' @author Moisès Bernabeu
#' @example man/examples/get_monophyletics.R
get_monophyletics <- function(tree, group_column, tree_data=NULL,
                              get_sp=NULL, data_sp_column=NULL) {
  if (typeof(tree) != 'S4') {
    if (!is.null(tree_data)) {
      tree <- ggcollapse::annotate_tree(tree, tree_data, get_sp, data_sp_column)
    } else {
      quit('No correct treedata specified', status = 1)
    }
  } else {
    ntips <- ape::Ntip(treeio::as.phylo(tree))
    tree_data <- as.data.frame(treeio::as_tibble(tree))
    tree_data_rows <- tree_data[, 'node'] <= ntips
    tree_data_names <- tree_data[tree_data_rows, 4]
    tree_data <- tree_data[tree_data_rows, -c(1:4)]
    row.names(tree_data) <- tree_data_names
  }

  ntips <- ape::Ntip(treeio::as.phylo(tree))
  int_nodes <- (ntips + 1):(2 * ntips - 1)

  monophyletic_groups <- c()
  j <- 1
  for (i in int_nodes) {
    desc <- phytools::getDescendants(treeio::as.phylo(tree), node = i)
    desc_tips <- desc[which(desc <= ntips)]

    desc_groups <- unique(tree@data[desc_tips, group_column][[group_column]])

    if (length(desc_groups) == 1 & !is_son(tree, monophyletic_groups, i)) {
      monophyletic_groups[j] <- i
      names(monophyletic_groups)[j] <- desc_groups
      j <- j + 1
    }
  }

  return(monophyletic_groups)
}
