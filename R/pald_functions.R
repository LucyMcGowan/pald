#' @import igraph graphics
NULL
library(igraph, graphics)

#' Cohesion Matrix
#'
#'
#' Creates a matrix of (pairwise) cohesion values from a matrix of pairwise distances or a 'dist' object.
#'
#'
#' Computes the matrix of (pairwise) cohesion values, C_xw, from a matrix of pairwise distances or 'dist' object. Cohesion is an interpretable probability that reflects the strength of alignment of a point, w, to another point, x.
#' The rows of the cohesion matrix can be seen as providing neighborhood weights.  These values may be used for defining associated weighted graphs (for the purposes of community analysis) as in BMM22.
#'
#' Given an n x n distance matrix, the sum of the entries in the resulting cohesion matrix is always equal to n/2.
#' Cohesion is partitioned local depth (see `local_depths`) and thus the row sums of the cohesion matrix provide a measure of local depth centrality.
#'
#' @param D A distance matrix or `dist' object.
#' @return The matrix, C, of cohesion values.
#' @examples
#' plot(exdata1)
#' text(exdata1 + .08, lab = 1:8)
#' D <- dist(exdata1)
#' C<-cohesion_matrix(D)
#' C
#' #neighbor weights (provided by cohesion) for the 8th point in exdata1
#' C[8, ]
#' localdepths<-rowSums(C)
#' @references K. S. Berenhaut, K. E. Moore, R. L. Melvin, A social perspective on perceived distances reveals deep community structure. Proc. Natl. Acad. Sci., 119(3), 2022.
#'
#' @export
cohesion_matrix <- function(D) {
  D <- round(as.matrix(D), 15) ## Rounds distances to fifteen decimal places
  if (dim(D)[1] != dim(D)[2]) {
    return(print("D is not square, please provide a distance matrix or 'dist' object"))
  } else {
    n <- dim(D)[1]
    if (is.null(rownames(D)[1])) {
      rownames(D) <- 1:n
    }
    C <- matrix(0, n, n)
    for (x in 1:(n - 1)) {
      for (y in (x + 1):n) {
        dx <- D[, x]
        dy <- D[, y]
        Uxy <- (which((dx <= D[y, x]) | (dy <= D[x, y])))
        wx <- 1 * (dx[Uxy] < dy[Uxy]) + .5 * ((dx[Uxy] == dy[Uxy]))
        C[x, Uxy] <- C[x, Uxy] + 1 / (length(Uxy)) * wx
        C[y, Uxy] <- C[y, Uxy] + 1 / (length(Uxy)) * (1 - wx)
      }
    }
    rownames(C) <- rownames(D)
    colnames(C) <- rownames(D)
    return(C/(n-1))
  }
}



#' Local (Community) Depths
#'
#' Creates a vector of local depths from a matrix of distances (or 'dist' object).
#'
#' Local depth is an interpretable probability which reflects aspects of relative position and centrality via distance comparasions (i.e., d(z, x) < d(z, y)).
#' The average of the local depth values is always 1/2.  Cohesion is partitioned local depth (see `cohesion_matrix`); the row-sums of the cohesion matrix are the values of local depth.
#' One can optionally provide a pre-computed cohesion matrix using the optional input 'is.cohesion = TRUE'.
#'
#' @param D A matrix of pairwise distances or 'dist' object.
#' @param is.cohesion Set to TRUE when D is a pre-computed cohesion matrix; the default is FALSE.
#' @return A vector of local depths.
#' @examples
#' D<-dist(exdata1)
#' local_depths(D)
#' C<-cohesion_matrix(D)
#' local_depths(C, is.cohesion =  TRUE)
#' #local depths are the row sums of the cohesion matrix
#' rowSums(C)
#' #cognate distance data
#' ld_lang<-sort(local_depths(cognate_dist))
#' @export
local_depths <- function(D, is.cohesion = FALSE) {
  if(is.cohesion==FALSE){
    C<-cohesion_matrix(D)
  }else{C<-as.matrix(D)}
  return(apply(C, 1, sum))
}


#' Cohesion Threshold for Strong Ties
#'
#' Given a cohesion matrix, provides the value of the threshold above which values of cohesion are considered "particularly strong".
#'
#' The threshold considered in BMM22 which may be used for distinguishing between strong and weak ties.
#' The threshold is equal to half the average of the diagonal of the cohesion matrix, see BMM22.
#'
#' @param C A matrix of cohesion values (see `cohesion_matrix`).
#' @return The value of the threshold.
#' @examples
#' C <- cohesion_matrix(dist(exdata1))
#' strong_threshold(C)
#' mean(diag(C))/2
#' #points whose cohesion are greater than the threshold may be considered (strong) neighbors
#' which(C[3, ] > strong_threshold(C))
#' #note that the number of (strongly-cohesive) neighbors varies across the space
#' which(C[4, ] > strong_threshold(C))
#' C[4, c(2, 3, 4, 6)] #cohesion values can provide neighbor weights
#'
#' @references K. S. Berenhaut, K. E. Moore, R. L. Melvin, A social perspective on perceived distances reveals deep community structure. Proc. Natl. Acad. Sci., 119(3), 2022.
#'
#' @export
strong_threshold <- function(C) {
  return(mean(diag(C)) / 2)
}



#' Cohesion Matrix: Strong Ties
#'
#' Provides the symmetrized and thresholded matrix of cohesion values.
#'
#' The threshold is that provided by strong_threshold (and is equal to half of the average of the diagonal of C).
#' Values of the cohesion matrix which are less than the threshold are set to zero.
#' The symmetrization, if desired, is computed using the entry-wise minimum of C and the transpose of C (i.e., min(C_ij, C_ji)).
#' The matrix provided by cohesion_strong (with default symmetric = TRUE) is the adjacency matrix for the graph of strong ties (the cluster graph), see `community_graphs` and `pald`.
#'
#' @param C A matrix of cohesion values (see `cohesion_matrix`).
#' @param symmetric Whether the returned matrix should be made symmetric (using the minimum); the default is TRUE.
#' @return The symmetrized cohesion matrix in which all entries corresponding to weak ties are set to zero.
#' @examples
#' C <- cohesion_matrix(dist(exdata2))
#' strong_threshold(C)
#' cohesion_strong(C)
#' #To illustrate the calculation performed
#' C_strong<-C
#' C_strong[C < strong_threshold(C)]<-0 #C_strong is equal to cohesion_strong(C, symmetric = FALSE)
#' C_strong_sym<-pmin(C_strong, t(C_strong)) #C_strong_sym is equal to cohesion_strong(C)
#' #The (cluster) graph whose adjacency matrix, CS, is the matrix of strong ties
#' CS<-cohesion_strong(C)
#' library(igraph)
#' G_strong <- simplify(graph.adjacency(CS, weighted = TRUE, mode = "undirected"))
#' plot(G_strong)
#' @export
cohesion_strong <- function(C, symmetric = TRUE) {
  threshold <- strong_threshold(C)
  C[C < threshold] <- 0
  if (symmetric == TRUE) {
    C <- pmin(C, t(C)) # uses minimum for mutual relationship, min(Cxw, Cwx)
  }
  return(C)
}



#' Community Graphs
#'
#' Provides the graphs whose edge weights are (mutual) cohesion together with a graph layout.
#'
#' Constructs the graph objects whose edge weights are (mutual) cohesion (see `cohesion_matrix`), self-loops are removed.
#' The graph G has adjacency matrix equal to the symmetrized cohesion matrix (using the entry-wise minimum of C and its transpose).
#' The graph G_strong has adjacency matrix equal to the thresholded and symmetrized cohesion matrix (see `cohesion_strong`).  The threshold is equal to half of the average of the diagonal of the
#' cohesion matrix (see `strong_threshold`).
#'
#' A layout is also computed using the Fructerman Reingold (FR) force-directed graph drawing algorithm.  As a result, it will provide a somewhat different layout each time it is run.
#'
#' @param C A matrix of cohesion values (see `cohesion_matrix`).
#' @return A list consisting of:
#'  \itemize{
#'        \item G: the weighted (community) graph whose edge weights are mutual cohesion
#'        \item  G_strong: the weighted (community) graph consisting of edges for which mutual cohesion is greater than the threshold for strong ties (see strong_threshold)
#'        \item  layout: the layout, using the Fruchterman Reingold (FR) force-directed graph drawing for the graph G
#' }
#' @examples
#' C <- cohesion_matrix(dist(exdata2))
#' plot(community_graphs(C)$G_strong)
#' plot(community_graphs(C)$G_strong, layout = community_graphs(C)$layout)
#' @export
community_graphs <- function(C) {
  C_symmetric <- pmin(C, t(C)) # relationship strength between x and w is mutual, i.e., min(Cxw, Cwx)
  G <- simplify(graph.adjacency(C_symmetric, weighted = TRUE, mode = "undirected"))
  G_strong <- simplify(graph.adjacency(cohesion_strong(C, symmetric = TRUE), weighted = TRUE, mode = "undirected"))
  lay <- layout_with_fr(G)

  # Points, x, for which mutual cohesion with all others are zero are omitted from the graph
  # This prints an alert that such points are not included in G
  Cdiagz <- C_symmetric
  diag(Cdiagz) <- 0
  isolated <- rownames(C)[apply(Cdiagz == 0, 1, all)]
  if (length(isolated) > 0) {
    print(paste("These points have no (strong nor weak) ties:", isolated))
  }
  return(list(G = G, G_strong = G_strong, layout = lay))
}





#' Plot Community Graphs
#'
#' Provides a plot of the community graphs, with connected components of the graph of strong ties colored by connected component.
#'
#' Plots the community graph, G, with the sub-graph of strong ties emphasized and colored by connected component.  If no layout is provided, the Fructerman-Reingold (FR) graph drawing algorithm is used.
#' Note that the FR graph drawing algorithm provides a somewhat different layout each time it is run.  You can also access and save a given graph layout using `community_graphs(C)$layout`.
#' The example below shows how to display only a subset of vertex labels.
#'
#' Note that the parameter `emph_strong` is for visualization purposes only and does not influence the network layout.
#'
#' @param C A matrix of cohesion values (see `cohesion_matrix`).
#' @param layout A layout for the graph.  If none is specified, FR-graph drawing algorithm is used.
#' @param show.labels Set to FALSE to omit vertex labels (to display a subset of labels, use optional input "vertex.lab" and modify the label list).
#' @param only_strong Set to TRUE if only strong ties, G_strong, should be displayed; the default FALSE will show both strong (colored by connected component) and weak ties (in gray).
#' @param emph_strong The numeric factor by which the edge widths of strong ties are emphasized in the display; the default is 2.
#' @param edge.width.factor Modify to change displayed edge widths.
#' @param vertex.lab A vector containing label names.
#' @param colors A vector of display colors, if none is given a default is provided.
#' @param vertex.size A numeric value for vertex size.
#' @param vertex.color.vec A vector of colors names for coloring vertices.
#' @param vertex.label.cex A numeric value for modifying the vertex label size.
#' @return NULL
#' @examples
#' C <- cohesion_matrix(dist(exdata1))
#' plot_community_graphs(C, emph_strong = 1, layout = exdata1)
#' plot_community_graphs(C, only_strong = TRUE)
#'
#' C2<-cohesion_matrix(cognate_dist)
#' subset_lang_names<-rownames(C2)
#' subset_lang_names[sample(1:87, 60)]<-''
#' plot_community_graphs(C2, vertex.lab = subset_lang_names, vertex.size = 3)
#' @export
plot_community_graphs <- function(C, layout = NULL, show.labels = TRUE, only_strong = FALSE,
                                  emph_strong = 2, edge.width.factor = 50, vertex.lab = rownames(C), vertex.label.cex = 1, colors = NULL, vertex.color.vec = NULL, vertex.size = 1) {
  # Hide vertex labels if show.labels=FALSE.  Otherwise, displays vertex labels.
  if (show.labels == FALSE) {
    vertex.lab <- ""
  }
  # Call function community_graphs to obtain weighted graphs G and G_strong, as well as the FR layout of G (for display purposes)
  c_graphs <- community_graphs(C)
  # Color the edges of G_strong according to the connected components of the graph
  cluster_labels <- clusters(c_graphs$G_strong)$membership
  # if no layout is given, the FR layout for G is used
  if (length(layout) == 0) {
    layout <- c_graphs$layout
  }
  # if no vector of color names is given, a default list is given
  if (length(colors) == 0) {
    colors <- c("#5F4690", "#73AF48", "#1D6996", "#CC503E", "#38A6A5", "#EDAD08", "#994E95","#0F8554",   "#CC6677", "#E17C05", "#94346E",  "#666666", "#88CCEE", "#AA4499","#117733", "#332288", "#44AA99","#6F4070", "#999933",  "#DDCC77","#882255", "#661100", "#6699CC", "#888888")
  }
  # Note: the default layout is still determined by G (consisting of both strong and weak ties)
  if (only_strong == FALSE) {
    plot(c_graphs$G, vertex.label = rep('', dim(C)[[1]]), vertex.size = vertex.size,  edge.color = "gray", edge.width = E(c_graphs$G)$weight * edge.width.factor, layout = layout, asp = 0, xlim = c(-1, 1), ylim = c(-1, 1))
    par(new = TRUE)
  }
  if(length(vertex.color.vec)==0){
    vertex.color.vec <- colors[cluster_labels]
  }
  # plot the graph, G_strong, color edges by connected component
  plot(c_graphs$G_strong,
       vertex.label = vertex.lab, vertex.size = vertex.size, vertex.color = vertex.color.vec, vertex.label.cex = vertex.label.cex, vertex.label.color="black", edge.color = colors[cluster_labels[get.edgelist(c_graphs$G_strong)[, 1]]],
       edge.width = emph_strong * E(c_graphs$G_strong)$weight * edge.width.factor, layout = layout, asp = 0, xlim = c(-1, 1), ylim = c(-1, 1)
  )
}


#' Partitioned Local Depths (PaLD)
#'
#' A wrapper function which computes the cohesion matrix, local depths, community graphs and provides a plot of the community graphs with connected components of the graph of strong ties colored by connected component.
#'
#' This function re-computes the cohesion matrix each time it is run.
#' To avoid unnecessary computation when creating visualizations, use the function `cohesion_matrix` to
#' compute the cohesion matrix which may then be taken as input for `local_depths`,
#' `strong_threshold`, `cohesion_strong`, `community_graphs`,
#' and `plot_community_graphs`.  For further detail about each component, see the documentation for each of the above functions.
#'
#'
#' @param D A dist object or matrix of pairwise distances.
#' @param layout A layout for the graph.  If none is specified, FR-graph drawing algorithm is used.
#' @param show.labels Set to FALSE to omit vertex labels (to display a subset of labels, use optional input "vertex.lab" and modify the label list).
#' @param only_strong Set to TRUE if only strong ties, G_strong, should be displayed; the default FALSE will show both strong (colored by connected component) and weak ties (in gray).
#' @param emph_strong The numeric factor by which the edge widths of strong ties are emphasized in the display; the default is 2.
#' @param edge.width.factor Modify to change displayed edge widths.
#' @param vertex.lab A vector containing label names.
#' @param colors A vector of display colors, if none is given a default is provided.
#' @param vertex.size A numeric value for vertex size.
#' @param vertex.color.vec A vector of colors names for coloring vertices.
#' @param vertex.label.cex A numeric value for modifying the vertex label size.
#' @param show.plot Set to TRUE to display plot; the default is TRUE.
#' @return A list consisting of:
#' \itemize{
#'         \item C: the matrix of cohesion values
#'        \item  local_depths: a vector of local depths
#'        \item clusters: a vectors of (community) cluster labels
#'       \item  threshold: the threshold above which cohesion is consider
#'       particularly strong
#'        \item  C_strong: the thresholded matrix of cohesion values
#'      \item   G: the graph whose edges weights are mutual cohesion
#'       \item  G_strong: the weighted graph whose edges are those for
#'       which cohesion is particularly strong
#'       \item  layout: a FR force-directedlayout associated with G
#' }
#' @examples
#' D <- dist(exdata2)
#' pald_results<-pald(D)
#' pald_results$local_depths
#' pald(D, layout=exdata2, show.labels = FALSE)
#'
#' C<-cohesion_matrix(D)
#' local_depths(C)
#' plot_community_graphs(C, layout = exdata2, show.labels = FALSE)
#'
#' pald_languages<-pald(cognate_dist)
#' head(pald_languages$local_depths)
#'
#' @references K. S. Berenhaut, K. E. Moore, R. L. Melvin, A social perspective on perceived distances reveals deep community structure. Proc. Natl. Acad. Sci., 119(3), 2022.
#'
#' @export
pald <- function(D, show.plot = TRUE, layout = NULL, show.labels = TRUE, only_strong = FALSE,
                 emph_strong = 2, edge.width.factor = 50, vertex.lab = rownames(C), colors = NULL, vertex.color.vec = NULL, vertex.label.cex = 1, vertex.size = 1) {
  C <- cohesion_matrix(D)
  C_graphs <- community_graphs(C)
  if (show.plot == TRUE) {
    plot_community_graphs(C,
                          layout = layout, show.labels = show.labels, only_strong = only_strong,
                          emph_strong = emph_strong, edge.width.factor = edge.width.factor, vertex.lab = vertex.lab, colors = colors, vertex.color.vec = vertex.color.vec, vertex.size = vertex.size
    )
  }

  return(list(C = C, local_depths = local_depths(C), clusters= clusters(C_graphs$G_strong)$membership, threshold = strong_threshold(C), C_strong = cohesion_strong(C), G = C_graphs$G, G_strong = C_graphs$G_strong, layout = C_graphs$layout))
}


#' Distance Cohesion Plot
#'
#' Provides a plot of cohesion against distance with the threshold indicated by a horizontal line.
#'
#' The plot of cohesion against distance provides a visualization for the manner in which distance is transformed.
#' The threshold distinguishing strong and weak ties is indicated by a horizontal line.
#' When there are separated regions with different density, one observes vertical bands of color, see example below and BMM22.  For each distance pair in D,
#' the corresponding value of cohesion is computed.  If the pair is within a single cluster, the point is colored
#' (with the same color provided by the pald and plot_community_graph functions).  Weak ties appear below the threshold.
#'
#' Note that cohesion is not symmetric, and so all n^2 points are plotted.  A gray point above the threshold corresponds to a pair in which the value of cohesion is greater than
#' the threshold in only one direction.  If one only wants to observe mutual cohesion (i.e., cohesion made symmetric using the minimum),
#' set "mutual=TRUE".
#'
#' @param D A matrix of distance values.
#' @param colors A vector of color names, if none is given a default is provided.
#' @param mutual Set to TRUE to consider mutual cohesion (i.e., symmetrized using the minimum); the default is FALSE.
#' @param xlim_max If desired, set the maximum value of distance which is displayed on the x-axis.
#' @param cex Factor by which points should be scaled relative to the default.
#' @param weak_gray Set to TRUE to display the plot with all weak ties plotted in gray; the default is FALSE.
#' @return A plot of cohesion against distance with threshold indicated by a horiztonal line.
#' @examples
#' D <- dist(exdata2)
#' dist_cohesion_plot(D)
#' dist_cohesion_plot(D, mutual = TRUE)
#' C<-cohesion_matrix(D)
#' threshold<-strong_threshold(C) #the horizontal line
#' dist_cohesion_plot(D, mutual = TRUE, weak_gray = TRUE)
#' @export
dist_cohesion_plot<-function(D, xlim_max = NULL, colors = NULL, mutual = FALSE, cex = 1, weak_gray = FALSE){
  par(xpd=FALSE)
  if(length(xlim_max)==0){xlim_max = max(D)}
  D<-as.matrix(D)
  if (dim(D)[1] != dim(D)[2]) {
    return(print("D is not square, please provide a distance matrix or 'dist' object"))
  }else{
  C<-cohesion_matrix(D)
  if(mutual ==TRUE){C <- pmin(C, t(C))}
  if (length(colors) == 0) {
    colors <- c("#5F4690", "#73AF48", "#1D6996", "#CC503E", "#38A6A5", "#EDAD08", "#994E95","#0F8554",   "#CC6677", "#E17C05", "#94346E",  "#666666", "#88CCEE", "#AA4499","#117733", "#332288", "#44AA99","#6F4070", "#999933",  "#DDCC77","#882255", "#661100", "#6699CC", "#888888")
  }
  n<-dim(C)[[1]]
  C_graphs<-community_graphs(C)
  pald_clusters<-clusters(C_graphs$G_strong)$membership
  C_colors<-matrix(colors[pald_clusters], n, n)
  C_colors[matrix(pald_clusters, n, n) != t(matrix(pald_clusters, n, n))]<-"gray"
  if(weak_gray ==TRUE){C_colors[C < strong_threshold(C)]<-"gray" }
  diag(C_colors)<-"black"
  C_pch<-matrix(16, n, n)
  C_pch[C < strong_threshold(C)]<-1
  if(mutual == TRUE){ylable = "mutual cohesion"}else{ylable = "cohesion"}
  cex1<-cex
  plot(D, C, col=C_colors, pch = C_pch, ylab = ylable, xlab = "distance", xlim = c(0, xlim_max), cex = cex1)
  abline(h = strong_threshold(C))}
}


