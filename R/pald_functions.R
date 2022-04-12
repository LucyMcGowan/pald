#' Cohesion Matrix
#'
#'
#' Creates a matrix of (pairwise) cohesion values from a matrix of pairwise
#' distances or a [`dist`] object.
#'
#'
#' Computes the matrix of (pairwise) cohesion values, C_xw, from a matrix of
#' pairwise distances or a [`dist`] object. Cohesion is an interpretable probability
#' that reflects the strength of alignment of a point, `w`, to another point, `x`.
#' The rows of the cohesion matrix can be seen as providing neighborhood
#' weights.  These values may be used for defining associated weighted graphs
#' (for the purpose of community analysis) as in  Berenhaut, Moore, and
#' Melvin (2022).
#'
#' Given an n x n distance matrix, the sum of the entries in the resulting
#' cohesion matrix is always equal to n/2.
#' Cohesion is partitioned local depth (see [`local_depths`]) and thus the row
#' sums of the cohesion matrix provide a measure of local depth centrality.
#'
#' If you have a matrix that is already a cohesion matrix and you would like to
#' add the class, see [`as_cohesion_matrix()`].
#'
#' @param d A matrix of pairwise distances or a [`dist`] object.
#' @return The matrix of cohesion values. An object of class `cohesion_matrix`.
#' @examples
#'
#' plot(exdata1)
#' text(exdata1 + .08, lab = 1:8)
#'
#' D <- dist(exdata1)
#' C <- cohesion_matrix(D)
#' C
#'
#' ## neighbor weights (provided by cohesion) for the 8th point in exdata1
#' C[8, ]
#' localdepths <- rowSums(C)
#'
#' @references K. S. Berenhaut, K. E. Moore, R. L. Melvin, A social perspective
#' on perceived distances reveals deep community structure. Proc. Natl.
#' Acad. Sci., 119(4), 2022.
#'
#' @export
cohesion_matrix <- function(d) {
  d <- check_dist(d)
  n <- dim(d)[1]
  if (is.null(rownames(d)[1])) {
    rownames(d) <- 1:n
  }
  c <- matrix(0, n, n)
  for (x in 1:(n - 1)) {
    for (y in (x + 1):n) {
      dx <- d[, x]
      dy <- d[, y]
      uxy <- (which((dx <= d[y, x]) | (dy <= d[x, y])))
      wx <- 1 * (dx[uxy] < dy[uxy]) + .5 * ((dx[uxy] == dy[uxy]))
      c[x, uxy] <- c[x, uxy] + 1 / (length(uxy)) * wx
      c[y, uxy] <- c[y, uxy] + 1 / (length(uxy)) * (1 - wx)
    }
  }
  rownames(c) <- rownames(d)
  colnames(c) <- rownames(d)

  c <- as_cohesion_matrix(c / (n - 1))

  return(c)
}


#' Coerce a matrix to a cohesion matrix object
#'
#' `as_cohesion_matrix()` converts an existing matrix into an object of class
#' `cohesion_matrix`.
#'
#' @param c A matrix of cohesion values (see [`cohesion_matrix`]).
#'
#' @return Object of class `cohesion_matrix`
#' @export
#'
#' @examples
#' C <- matrix(
#'   c(0.25, 0.125, 0.125, 0,
#'    0.125, 0.25, 0, 0.125,
#'    0.125, 0, 0.25, 0.125,
#'    0, 0.125, 0.125, 0.25
#' ), nrow = 4, byrow = TRUE)
#'
#' class(C)
#'
#' C <- as_cohesion_matrix(C)
#' class(C)
as_cohesion_matrix <- function(c) {
  if (!is.matrix(c)) {
    stop_glue("The cohesion matrix input must be a matrix.\n * You input an",
              "object of class `{class(c)}`")
  }

  if (dim(c)[1] != dim(c)[2]) {
    stop_glue("The cohesion matrix must be a square matrix.")
  }
  if (is.null(rownames(c)[1])) {
    rownames(c) <- 1:nrow(c)
    colnames(c) <- 1:nrow(c)
  }
  cl <- class(c)
  structure(c, class = c("cohesion_matrix", cl))
}

#' Local (Community) Depths
#'
#' Creates a vector of local depths from a matrix of distances (or `dist`
#' object).
#'
#' Local depth is an interpretable probability which reflects aspects of
#' relative position and centrality via distance comparisons
#' (i.e., d(z, x) < d(z, y)).
#'
#' The average of the local depth values is always 1/2.  Cohesion is
#' partitioned local depth (see [`cohesion_matrix`]); the row-sums of the
#' cohesion matrix are the values of local depth.
#'
#'
#' @param d A matrix of pairwise distances, a [`dist`] object, or a
#'   [`cohesion_matrix`] object.
#'
#' @return A vector of local depths.
#' @examples
#'
#' D <- dist(exdata1)
#' local_depths(D)
#' C <- cohesion_matrix(D)
#' local_depths(C)
#'
#' ## local depths are the row sums of the cohesion matrix
#' rowSums(C)
#'
#' ## cognate distance data
#'
#' ld_lang <- sort(local_depths(cognate_dist))
#' @export
local_depths <- function(d) {
  if (!is_cohesion_matrix(d)) {
    d <- cohesion_matrix(d)
  }
  return(apply(d, 1, sum))
}


#' Cohesion Threshold for Strong Ties
#'
#' Given a cohesion matrix, provides the value of the threshold above which
#' values of cohesion are considered "particularly strong".
#'
#' The threshold considered in Berenhaut, Moore, and Melvin (2022) which may be
#' used for distinguishing between strong and weak ties.
#' The threshold is equal to half the average of the diagonal of the cohesion
#' matrix, see Berenhaut, Moore, and Melvin (2022).
#'
#' @param c A `cohesion_matrix` object, a matrix of cohesion values
#'   (see [`cohesion_matrix`]).
#' @return The value of the threshold.
#'
#' @examples
#' C <- cohesion_matrix(dist(exdata1))
#' strong_threshold(C)
#' mean(diag(C)) / 2
#'
#' ## points whose cohesion are greater than the threshold may be considered
#' ## (strong) neighbors
#' which(C[3, ] > strong_threshold(C))
#'
#' ## note that the number of (strongly-cohesive) neighbors varies across the
#' ## space
#' which(C[4, ] > strong_threshold(C))
#' C[4, c(2, 3, 4, 6)] # cohesion values can provide neighbor weights
#'
#' @references K. S. Berenhaut, K. E. Moore, R. L. Melvin, A social perspective
#' on perceived distances reveals deep community structure. Proc. Natl. Acad.
#' Sci., 119(4), 2022.
#'
#' @export
strong_threshold <- function(c) {
  check_cohesion_matrix(c)
  return(mean(diag(c)) / 2)
}



#' Cohesion Matrix: Strong Ties
#'
#' Provides the symmetrized and thresholded matrix of cohesion values.
#'
#' The threshold is that provided by strong_threshold (and is equal to half of
#' the average of the diagonal of `c`).
#' Values of the cohesion matrix which are less than the threshold are set to
#' zero.
#' The symmetrization, if desired, is computed using the entry-wise (parallel)
#' minimum of C and ts transpose (i.e., `min(C_ij, C_ji)`).
#' The matrix provided by cohesion_strong (with default `symmetric = TRUE`) is
#' the adjacency matrix for the graph of strong ties (the cluster graph), see
#' [`community_graphs`] and [`pald`].
#'
#' @param c A `cohesion_matrix` object, a matrix of cohesion values
#'   (see [`cohesion_matrix`]).
#' @param symmetric Logical. Whether the returned matrix should be made
#' symmetric (using the minimum); the default is `TRUE`.
#' @return The symmetrized cohesion matrix in which all entries corresponding
#' to weak ties are set to zero.
#'
#' @examples
#' C <- cohesion_matrix(dist(exdata2))
#' strong_threshold(C)
#' cohesion_strong(C)
#'
#' ## To illustrate the calculation performed
#' C_strong <- C
#'
#' ## C_strong is equal to cohesion_strong(C, symmetric = FALSE)
#' C_strong[C < strong_threshold(C)] <- 0
#'
#' ## C_strong_sym is equal to cohesion_strong(C)
#' C_strong_sym <- pmin(C_strong, t(C_strong))
#'
#' ## The (cluster) graph whose adjacency matrix, CS,
#' ## is the matrix of strong ties
#' CS <- cohesion_strong(C)
#'
#' if (requireNamespace("igraph", quietly = TRUE)) {
#' G_strong <- igraph::simplify(
#'   igraph::graph.adjacency(CS, weighted = TRUE, mode = "undirected")
#'   )
#' plot(G_strong)
#' }
#' @export
cohesion_strong <- function(c, symmetric = TRUE) {
  check_cohesion_matrix(c)
  threshold <- strong_threshold(c)
  c[c < threshold] <- 0
  if (symmetric == TRUE) {
    c <- pmin(c, t(c)) # uses minimum for mutual relationship, min(Cxw, Cwx)
  }
  return(c)
}



#' Community Graphs
#'
#' Provides the graphs whose edge weights are (mutual) cohesion, together with
#' a graph layout.
#'
#' Constructs the graphs whose edge weights are (mutual) cohesion
#' (see [`cohesion_matrix`]), self-loops are removed.
#' The graph G has adjacency matrix equal to the symmetrized cohesion matrix
#' (using the entry-wise parallel minimum of C and its transpose).
#' The graph G_strong has adjacency matrix equal to the thresholded and
#' symmetrized cohesion matrix (see [`cohesion_strong`]).  The threshold is
#' equal to half of the average of the diagonal of the
#' cohesion matrix (see [`strong_threshold`]).
#'
#' A layout is also computed using the Fruchterman-Reingold (FR)
#' force-directed graph drawing algorithm.  As a result, it may provide a
#' somewhat different layout each time it is run.
#'
#' @param c A `cohesion_matrix` object, a matrix of cohesion values
#'   (see [`cohesion_matrix`]).
#' @return A list consisting of:
#'  \itemize{
#'        \item `G`: the weighted (community) graph whose edge weights are mutual
#'        cohesion
#'        \item `G_strong`: the weighted (community) graph consisting of edges
#'        for which mutual cohesion is greater than the threshold for strong
#'        ties (see [`strong_threshold`])
#'        \item `layout`: the layout, using the Fruchterman Reingold (FR)
#'        force-directed graph drawing for the graph `G`
#' }
#' @examples
#' C <- cohesion_matrix(dist(exdata2))
#' plot(community_graphs(C)$G_strong)
#' plot(community_graphs(C)$G_strong, layout = community_graphs(C)$layout)
#' @export
community_graphs <- function(c) {
  check_cohesion_matrix(c)
  # relationship strength between x and w is mutual, i.e., min(Cxw, Cwx)
  c_symmetric <- pmin(c, t(c))
  g <- igraph::simplify(
    igraph::graph.adjacency(c_symmetric, weighted = TRUE, mode = "undirected")
  )
  g_strong <- igraph::simplify(
    igraph::graph.adjacency(
      cohesion_strong(c, symmetric = TRUE),
      weighted = TRUE, mode = "undirected")
  )
  lay <- igraph::layout_with_fr(g)

  # Points, x, for which mutual cohesion with all others are zero
  # This prints an alert
  any_isolated(c)
  return(list(G = g, G_strong = g_strong, layout = lay))
}


#' Any isolated
#'
#' Checks for isolated points.
#'
#' @param c A `cohesion_matrix` object, a matrix of cohesion values
#'   (see [`cohesion_matrix`]).
#'
#' @return Logical, indicating whether any points are isolated.
#' @examples
#' d <- data.frame(
#'   x1 = c(1, 2, 3, 6),
#'   x2 = c(2, 1, 3, 10)
#'   )
#' D <- dist(d)
#' C <- cohesion_matrix(D)
#' any_isolated(C)
#' @export
any_isolated <- function(c) {
  check_cohesion_matrix(c)
  c_symmetric <- pmin(c, t(c))
  cdiagz <- c_symmetric
  diag(cdiagz) <- 0
  isolated <- rownames(c)[apply(cdiagz == 0, 1, all)]
  is_isolated <- length(isolated) > 0
  if (is_isolated) {
    message_glue(
      "These points have no (strong nor weak) ties: \n * {isolated}"
    )
  }
  return(invisible(is_isolated))
}


#' Plot Community Graphs
#'
#' Provides a plot of the community graphs, with connected components of the
#' graph of strong ties colored by connected component.
#'
#' Plots the community graph, G, with the sub-graph of strong ties emphasized
#' and colored by connected component.  If no layout is provided, the
#' Fruchterman-Reingold (FR) graph drawing algorithm is used.
#' Note that the FR graph drawing algorithm may provide a somewhat different
#' layout each time it is run.  You can also access and save a given graph
#' layout using `community_graphs(C)$layout`.
#' The example below shows how to display only a subset of vertex labels.
#'
#' Note that the parameter `emph_strong` is for visualization purposes
#' only and does not influence the network layout.
#'
#' @param c A `cohesion_matrix` object, a matrix of cohesion values
#'   (see [`cohesion_matrix`]).
#' @param show_labels Set to `FALSE` to omit vertex labels (to display a subset
#'   of labels, use optional parameter `vertex.label` to modify the label list).
#'   Default: `TRUE`.
#' @param only_strong Set to `TRUE` if only strong ties, G_strong, should be
#'   displayed; the default `FALSE` will show both strong (colored by connected
#'   component) and weak ties (in gray).
#' @param emph_strong Numeric. The numeric factor by which the edge widths of
#'    strong ties are emphasized in the display; the default is `2`.
#' @param edge_width_factor Numeric. Modify to change displayed edge widths.
#'   Default: `50`.
#' @param colors A vector of display colors, if none is given a default list
#'   (of length 24) is provided.
#' @param ... Optional parameters to pass to the [`igraph::plot.igraph`].
#'  function. Some commonly passed arguments include:
#'    * `layout`  A layout for the graph.  If none is specified, FR-graph
#'   drawing algorithm is used.
#'    * `vertex.label` A vector containing label names. If none is given,
#'    the rownames of `c` are used
#'    * `vertex.size` A numeric value for vertex size (default = `1`)
#'    * `vertex.color.vec` A vector of color names for coloring the vertices
#'    * `vertex.label.cex` A numeric value for modifying the vertex label size.
#'    (default = `1`)
#'
#' @return A plot of the community graphs.
#' @examples
#' C <- cohesion_matrix(dist(exdata1))
#' plot_community_graphs(C, emph_strong = 1, layout = as.matrix(exdata1))
#' plot_community_graphs(C, only_strong = TRUE)
#'
#' C2 <- cohesion_matrix(cognate_dist)
#' subset_lang_names <- rownames(C2)
#' subset_lang_names[sample(1:87, 60)] <- ""
#' plot_community_graphs(C2, vertex.label = subset_lang_names, vertex.size = 3)
#' @export
plot_community_graphs <- function(c,
                                  show_labels = TRUE,
                                  only_strong = FALSE,
                                  emph_strong = 2,
                                  edge_width_factor = 50,
                                  colors = NULL,
                                  ...) {

  check_cohesion_matrix(c)
  dots <- list(...)

  dots[["vertex.size"]] <- dots[["vertex.size"]] %||% 1
  dots[["vertex.label.cex"]] <- dots[["vertex.label.cex"]] %||% 1
  dots[["asp"]] <- dots[["asp"]] %||% 0
  dots[["xlim"]] <- dots[["xlim"]] %||% c(-1, 1)
  dots[["ylim"]] <- dots[["ylim"]] %||% c(-1, 1)

  # Hide vertex labels if show_labels = FALSE. Otherwise, displays vertex labels.
  if (!show_labels) {
    dots[["vertex.label"]] <- NA
  } else {
    dots[["vertex.label"]] <- dots[["vertex.label"]] %||% rownames(c)
  }
  ## Call function community_graphs to obtain weighted graphs G and G_strong,
  ## as well as the FR layout of G (for display purposes)
  c_graphs <- community_graphs(c)

  # Color the edges of G_strong according to the connected components of the
  # graph
  cluster_labels <- igraph::clusters(c_graphs$G_strong)$membership

  # if no layout is given, the FR layout for G is used
  dots[["layout"]] <- dots[["layout"]] %||% c_graphs$layout

  # if no vector of color names is given, a default list (of length 24) is given
  if (is.null(colors)) {
    colors <- pald::pald_colors
  }

  old <- list(egde.color = dots[["edge.color"]],
              edge.width = dots[["edge.width"]],
              vertex.label = dots[["vertex.label"]])

  ## Note: the default layout is still determined by G (consisting of both
  ## strong and weak ties)
  if (!only_strong) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    dots[["edge.color"]] <- "grey"
    dots[["edge.width"]] <- igraph::E(c_graphs$G)$weight * edge_width_factor
    dots[["x"]] <- c_graphs$G
    dots[["vertex.label"]] <- dots[["vertex.label"]] %||% rep("", dim(c)[[1]])
    do.call(igraph::plot.igraph, dots)
    graphics::par(new = TRUE)
  }

  dots[["vertex.label"]] <- old[["vertex.label"]]
  dots[["vertex.color"]] <- dots[["vertex.color"]] %||% colors[cluster_labels]
  dots[["edge.color"]] <- old[["edge.color"]] %||% colors[
    cluster_labels[igraph::get.edgelist(c_graphs$G_strong)[, 1]]
  ]
  dots[["vertex.label.color"]] <- dots[["vertex.label.color"]] %||% "black"
  dots[["edge.width"]] <- old[["edge.width"]] %||%
    emph_strong * igraph::E(c_graphs$G_strong)$weight * edge_width_factor
  dots[["x"]] <- c_graphs$G_strong

  # plot the graph, G_strong, color edges by connected component
  do.call(igraph::plot.igraph, dots)
}

#' Community clusters
#'
#' @param c A `cohesion_matrix` object, a matrix of cohesion values
#'   (see [`cohesion_matrix`]).
#'
#' @return A data frame with two columns:
#'  * `point`: The points from cohesion matrix `c`
#'  * `cluster`: The (community) cluster labels
#'
#' @examples
#' D <- dist(exdata2)
#' C <- cohesion_matrix(D)
#' community_clusters(C)
#' @export
community_clusters <- function(c) {
  check_cohesion_matrix(c)
  c_graphs <- community_graphs(c)
  cl <- igraph::clusters(c_graphs$G_strong)$membership
  data.frame(
    point = names(cl),
    cluster = cl
  )
}
#' Partitioned Local Depths (PaLD)
#'
#' A wrapper function which computes the cohesion matrix, local depths,
#' community graphs and provides a plot of the community graphs with connected
#' components of the graph of strong ties colored by connected component.
#'
#' This function re-computes the cohesion matrix each time it is run.
#' To avoid unnecessary computation when creating visualizations, use the
#' function [`cohesion_matrix`] to compute the cohesion matrix which may then
#' be taken as input for [`local_depths`], [`strong_threshold`],
#' [`cohesion_strong`], [`community_graphs`], and [`plot_community_graphs`].
#' For further details regarding each component, see the documentation for
#' each of the above functions.
#'
#'
#' @param d A matrix of pairwise distances or a [`dist`] object.
#' @param show_plot Set to `TRUE` to display plot; the default is `TRUE`.
#' @param show_labels Set to `FALSE` to omit vertex labels (to display a subset
#'   of labels, use optional parameter `vertex.label` to modify the label list).
#'   Default: `TRUE`.
#' @param only_strong Set to `TRUE` if only strong ties, G_strong, should be
#'   displayed; the default `FALSE` will show both strong (colored by connected
#'   component) and weak ties (in gray).
#' @param emph_strong Numeric. The numeric factor by which the edge widths of
#'    strong ties are emphasized in the display; the default is `2`.
#' @param edge_width_factor Numeric. Modify to change displayed edge widths.
#'   Default: `50`.
#' @param colors A vector of display colors, if none is given a default list
#'   (of length 24) is provided.
#' @param ... Optional parameters to pass to the [`igraph::plot.igraph`].
#'  function. Some commonly passed arguments include:
#'    * `layout`  A layout for the graph.  If none is specified, FR-graph
#'   drawing algorithm is used.
#'    * `vertex.label` A vector containing label names. If none is given,
#'    the rownames of `c` are used
#'    * `vertex.size` A numeric value for vertex size (default = `1`)
#'    * `vertex.color.vec` A vector of color names for coloring the vertices
#'    * `vertex.label.cex` A numeric value for modifying the vertex label size.
#'    (default = `1`)
#' @return A list consisting of:
#' \itemize{
#'    \item `C`: the matrix of cohesion values
#'    \item  `local_depths`: a vector of local depths
#'    \item `clusters`: a vector of (community) cluster labels
#'    \item  `threshold`: the threshold above which cohesion is considered
#'       particularly strong
#'    \item  `C_strong`: the thresholded matrix of cohesion values
#'    \item   `G`: the graph whose edges weights are mutual cohesion
#'    \item  `G_strong`: the weighted graph whose edges are those for
#'       which cohesion is particularly strong
#'    \item  `layout`: a FR force-directed layout associated with G
#' }
#' @examples
#' D <- dist(exdata2)
#' pald_results <- pald(D)
#' pald_results$local_depths
#' pald(D, layout = as.matrix(exdata2), show_labels = FALSE)
#'
#' C <- cohesion_matrix(D)
#' local_depths(C)
#' plot_community_graphs(C, layout = as.matrix(exdata2), show_labels = FALSE)
#'
#' pald_languages <- pald(cognate_dist)
#' head(pald_languages$local_depths)
#'
#' @references K. S. Berenhaut, K. E. Moore, R. L. Melvin, A social perspective
#' on perceived distances reveals deep community structure. Proc. Natl. Acad.
#' Sci., 119(4), 2022.
#'
#' @export
pald <- function(d,
                 show_plot = TRUE,
                 show_labels = TRUE,
                 only_strong = FALSE,
                 emph_strong = 2,
                 edge_width_factor = 50,
                 colors = NULL,
                 ...) {

  d <- check_dist(d)
  c <- cohesion_matrix(d)
  c_graphs <- community_graphs(c)

  if (show_plot) {
    plot_community_graphs(c,
                          show_labels = show_labels,
                          only_strong = only_strong,
                          emph_strong = emph_strong,
                          edge_width_factor = edge_width_factor,
                          colors = colors,
                          ...
    )
  }

  return(
    invisible(
      list(C = c,
           local_depths = local_depths(c),
           clusters = igraph::clusters(c_graphs$G_strong)$membership,
           threshold = strong_threshold(c),
           C_strong = cohesion_strong(c),
           G = c_graphs$G,
           G_strong = c_graphs$G_strong,
           layout = c_graphs$layout
      )
    )
  )
}


#' Distance Cohesion Plot
#'
#' Provides a plot of cohesion against distance, with the threshold indicated
#' by a horizontal line.
#'
#' The plot of cohesion against distance provides a visualization for the
#' manner in which distance is transformed.
#' The threshold distinguishing strong and weak ties is indicated by a
#' horizontal line.
#' When there are separated regions with different density, one can often
#' observe vertical bands of color, see example below and Berenhaut, Moore, and
#' Melvin (2022).  For each distance pair in `d`, the corresponding value of
#' cohesion is computed.  If the pair is within a single cluster, the point is
#' colored (with the same color provided by the [`pald`] and
#' [`plot_community_graphs`] functions).  Weak ties appear below the threshold.
#'
#' Note that cohesion is not symmetric, and so all `n^2` points are plotted.
#' A gray point above the threshold corresponds to a pair in which the value
#' of cohesion is greater than the threshold in only one direction.  If one
#' only wants to observe mutual cohesion (i.e., cohesion made symmetric via
#' the minimum), set `mutual = TRUE`.
#'
#' @param d A matrix of pairwise distances or a [`dist`] object.
#' @param mutual Set to `TRUE` to consider mutual cohesion (i.e., symmetrized
#'   using the minimum); the default is `FALSE`.
#' @param xlim_max If desired, set the maximum value of distance which
#'  is displayed on the x-axis.
#' @param cex Factor by which points should be scaled relative to the default.
#' @param colors A vector of color names, if none is given a default is
#'   provided.
#' @param weak_gray Set to `TRUE` to display the plot with all weak ties plotted
#'   in gray; the default is `FALSE`.
#' @return A plot of cohesion against distance with threshold indicated by a
#'   horizontal line.
#' @examples
#' D <- dist(exdata2)
#' dist_cohesion_plot(D)
#' dist_cohesion_plot(D, mutual = TRUE)
#' C <- cohesion_matrix(D)
#' threshold <- strong_threshold(C) #the horizontal line
#' dist_cohesion_plot(D, mutual = TRUE, weak_gray = TRUE)
#' @export
dist_cohesion_plot <- function(d,
                               mutual = FALSE,
                               xlim_max = NULL,
                               cex = 1,
                               colors = NULL,
                               weak_gray = FALSE) {

  d <- check_dist(d)

  if (is.null(xlim_max)) {
    xlim_max <- max(d)
  }

  c <- cohesion_matrix(d)

  if (mutual) {
    c <- pmin(c, t(c))
  }
  if (is.null(colors)) {
    colors <- c(
      "#5F4690", "#73AF48", "#1D6996", "#CC503E", "#38A6A5", "#EDAD08",
      "#994E95", "#0F8554", "#CC6677", "#E17C05", "#94346E",  "#666666",
      "#88CCEE", "#AA4499", "#117733", "#332288", "#44AA99", "#6F4070",
      "#999933",  "#DDCC77", "#882255", "#661100", "#6699CC", "#888888")
  }
  n <- dim(c)[[1]]

  c_graphs <- community_graphs(c)

  pald_clusters <- igraph::clusters(c_graphs$G_strong)$membership

  c_colors <- matrix(colors[pald_clusters], n, n)

  c_colors[matrix(pald_clusters, n, n) !=
             t(matrix(pald_clusters, n, n))] <- "gray"

  if (weak_gray) {
    c_colors[c < strong_threshold(c)] <- "gray"
  }
  diag(c_colors) <- "black"
  c_pch <- matrix(16, n, n)
  c_pch[c < strong_threshold(c)] <- 1
  if (mutual) {
    ylab <- "mutual cohesion"
  } else {
    ylab <- "cohesion"
  }
  graphics::plot(d,
                 c,
                 col = c_colors,
                 pch = c_pch,
                 ylab = ylab,
                 xlab = "distance",
                 xlim = c(0, xlim_max),
                 cex = cex)
  graphics::abline(h = strong_threshold(c))

}
