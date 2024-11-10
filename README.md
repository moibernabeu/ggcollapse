# ggcollapse: an R package to collapse nodes using ggtree
This package contains a set of functions that expand the `ggtree` package. The
main function of this set, `ggcollapse` allows the user to collapse nodes with
the same height for the collapsing triangles. The package also contains some
functions that help in the task of finding monophyletic groups of a certain
taxonomic rank and the full set of sister clades to a given node.

# Installation

To install `ggcollapse`, you just needs to execute the following
command:

```r
if (!require('devtools', quietly = TRUE)) {
    install.packages('devtools')
}

devtools::install_github('moibernabeu/ggcollapse')
```

# Tutorial

## Data for `ggcollapse`

Although the functions are versatile, the main idea of them is to use
three data types, a `phylo` object for the tree, a `data.frame` for the
data, and named vectors for the nodes (`numeric`) to collapse and the
colours of the collapsed clades (`character`).

### Tree and their tip labels
In this tutorial we will simulate our data, first, we simulate a tree:


``` r
library(treeio)

set.seed(2024-11-04)
tree <- rtree(20)
```

To see the whole potential of the package, we will assign species to
each tip, the tip label format chosen for this tutorial will be
`tip_SPECIES`, therefore,


``` r
tree$tip.label[3:20] <- paste(tree$tip.label[3:20], LETTERS[2:19], sep = '_')
# Adding a species duplication in the tree
tree$tip.label[1:2] <- paste(tree$tip.label[1:2], 'A', sep = '_')
```

The final tree, with the number of the internal nodes, would be:


``` r
library(ggtree)

ggtree(tree) +
  geom_tiplab() +
  geom_nodelab(aes(label = node), geom = 'label') +
  xlim(0, 5.5)
```

<img src="tutorial_figures/unnamed-chunk-4-1.svg" width="60%" style="display: block; margin: auto;" />

### Function to get the species from the tip label
To extract the species from each tip, we will design a function, which
will be added to the `ggcollapse` functions later on:


``` r
library(stringr)
get_sp <- function(tip.label) {
  return(str_split(tip.label, '_', n = 2, simplify = TRUE)[, 2])
}

tree$tip.label[1]
```

```
## [1] "t14_A"
```

``` r
get_sp(tree$tip.label[1])
```

```
## [1] "A"
```

Other options, and advantages, of using a species-getting function are the
search of the species in a table relating the tip value with the species, a
vector, etc.

### Data linked to the tree
Imagine that the species of your tree are distributed into two groups. To
incorporate this information to the tree we can take two strategies:

1. generate a `data.frame` where each tip has the information, in this case,
the species-getting function must return the header, or

2. generate a `data.frame` relating the species with its information, in this
case, the species-generating function must return the species.

In this tutorial, we will use the second strategy to generate de `data.frame`,
we also will generate a palette for the groups:


``` r
groups <- c(rep('G1', 7), rep('G2', 3), rep('G3', 9))
tree_data <- data.frame(species = LETTERS[1:19],
                        group = groups)
head(tree_data, n = 4)
```

```
##   species group
## 1       A    G1
## 2       B    G1
## 3       C    G1
## 4       D    G1
```

``` r
group_colours <- c('G1' = 'coral3',
                   'G2' = 'steelblue3',
                   'G3' = 'darkolivegreen4')
```

## `ggcollapse` functions

### Annotating the tree with the data using `annotate_tree`
To incorporate the tree data (`tree_data`) into the tree, we can use the 
`annotate_tree` function. This function associates the tip labels with the
information of the tree data table by using a linking function, in this case,
the species-getting function. We added the argument `data_sp_column` to force
the user to incorporate the value (as character or numeric) of the column where
the species code is. The resulting tree is a `treedata S4` object.


``` r
library(ggcollapse)

annot_tree <- annotate_tree(tree = tree,
                            tree_data = tree_data,
                            get_sp = get_sp,
                            data_sp_column = 'species')
```

Now, we can plot data across the tree, we can plot the group they belong using
points on the tips:


``` r
ggtree(annot_tree) +
  geom_tiplab(offset = 0.035) +
  geom_tippoint(aes(colour = group)) +
  scale_colour_manual(values = group_colours) +
  xlim(0, 5.5)
```

<img src="tutorial_figures/unnamed-chunk-8-1.svg" width="60%" style="display: block; margin: auto;" />

### Obtaining monophyletic clades of the groups with `get_monophyletics` 
Using the function `get_monophyletics` we can get the value of the nodes whose
all their descendants belong to the group specified in the `group_column`
argument. Here we show the example of a tree that has been already annotated,
however, the function can get the `tree_data`, `get_sp` and `data_sp_column` in
the case the tree is not annotated.


``` r
gr_mphy <- get_monophyletics(tree = annot_tree,
                             group_column = 'group')
gr_mphy
```

```
## G1 G1 G2 G3 
## 22 26 30 32
```

### Getting sister groups with `get_sisters`
An interesting function of this package, is the function `get_sisters`. This
function allows the user to obtain the vector of sister clades to a specific
tip or internal node. The function will return a vector with the number of the
sister clade nodes sorted by increasing topological distance to the node. To
get the sisters of a tip, we can specify the tree and the tip label, as well
as the tip number:


``` r
tip_sisters <- get_sisters(tree = tree,
                           node = 't8_O')
tip_sisters
```

```
## [1] 38 34 39 30 26 22
```

``` r
tip_num <- which(tree$tip.label == 't8_O')
tip_sisters_num <- get_sisters(tree = tree,
                               node = tip_num)
tip_sisters_num == tip_sisters
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

In case we want to name the vector of sisters with the most abundant value of
a column in the descendants a specific clade, we can use the annotated tree,
the node, and the naming column, in this case, the group. Using these arguments,
the function will return a named vector with the number of the sister clade
nodes named with the most abundant value for the column and its percentage
(e.g., if in a clade there are 3 tips from the `group` `G1` and 1 tip from the
`group` `G2` the clade would be annotated as `G1 75%`):


``` r
tip_sisters_nmd <- get_sisters(tree = annot_tree,
                               node = 't8_O',
                               naming_column = 'group')
tip_sisters_nmd
```

```
## G3 100% G3 100% G3 100% G2 100% G1 100% G1 100% 
##      38      34      39      30      26      22
```

We can follow the same procedure for an internal node:


``` r
node_sisters <- get_sisters(tree = tree,
                            node = 26)
node_sisters
```

```
## [1] 29 22
```

``` r
node_sisters_nmd <- get_sisters(tree = annot_tree,
                                node = 26,
                                naming_column = 'group')
node_sisters_nmd
```

```
##  G3 75% G1 100% 
##      29      22
```

### Plotting a collapsed tree with `ggcollapse`
The main point of this package is to allow the user to collapse some clades
easily. The `ggcollapse` function just requires a `phylo` tree object and
a named vector with the nodes to be collapsed. Imagine now that we want to
collapse the sister groups to `t8_O` that we computed previously:


``` r
# Naming the sisters
names(tip_sisters) <- paste('S', 1:length(tip_sisters), sep = '')

# Collapsing the sister nodes
ggcollapse(tree = tree,
           nodes = tip_sisters) +
  geom_tiplab() +
  xlim(0, 1.1)
```

<img src="tutorial_figures/unnamed-chunk-13-1.svg" width="60%" style="display: block; margin: auto;" />

We can also annotate the tree with the `ggcollapse` function, and it would
return a `ggtree` object with the data. In the following example, we will plot
the collapsed sisters with the uncollapsed nodes showing the species, rather
than the tip, label.


``` r
ggcollapse(tree = tree,
           nodes = tip_sisters,
           get_sp = get_sp,
           tree_data = tree_data,
           data_sp_column = 'species') +
  geom_tiplab(aes(label = species)) +
  xlim(0, 1.1)
```

<img src="tutorial_figures/unnamed-chunk-14-1.svg" width="60%" style="display: block; margin: auto;" />

Using the same logic, we can apply this to the monophyletic groups, for which
we have the node colours.


``` r
# Collapsing the monophyletic groups
ggcollapse(tree = tree,
           nodes = gr_mphy,
           node_colours = group_colours)
```

<img src="tutorial_figures/unnamed-chunk-15-1.svg" width="60%" style="display: block; margin: auto;" />

#### Collapsing modes
We can also change the visualisation of the collapsed triangle using the same
argument values as in `ggtree` (`min`, `max`, `mixed`). Our default is `mixed`.


``` r
ggcollapse(tree = tree,
           nodes = gr_mphy,
           collapse_mode = 'mixed',
           node_colours = group_colours) +
  xlim(0, 1.15) +

ggcollapse(tree = tree,
           nodes = gr_mphy,
           collapse_mode = 'max',
           node_colours = group_colours) +
  xlim(0, 1.15) +

ggcollapse(tree = tree,
           nodes = gr_mphy,
           collapse_mode = 'min',
           node_colours = group_colours) +
  xlim(0, 0.8)
```

<img src="tutorial_figures/unnamed-chunk-16-1.svg" width="100%" style="display: block; margin: auto;" />

#### Adding bootstrap values
As the `ggcollapse` function returns a `ggtree` object, the user can use all the
functionalities of this package, such as adding collapsing functions.


``` r
set.seed(09-08-1997)
supports <- round(runif(length(tree$tip.label) - 2, 25, 100))
tree$node.label <- c(NA, supports)

ggcollapse(tree = tree,
           nodes = gr_mphy,
           node_colours = group_colours) +
  geom_nodelab(aes(x = branch), nudge_y = 0.15)
```

<img src="tutorial_figures/unnamed-chunk-17-1.svg" width="60%" style="display: block; margin: auto;" />
