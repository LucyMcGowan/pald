stop_glue <- function(..., .sep = "", .envir = parent.frame(),
                      call. = FALSE, .domain = NULL) {
  stop(
    glue::glue(..., .sep = .sep, .envir = .envir),
    call. = call., domain = .domain
  )
}

warning_glue <- function(..., .sep = "", .envir = parent.frame(),
                         call. = FALSE, .domain = NULL) {
  warning(
    glue::glue(..., .sep = .sep, .envir = .envir),
    call. = call., domain = .domain
  )
}

message_glue <- function(..., .sep = "", .envir = parent.frame(),
                         .domain = NULL) {
  message(
    glue::glue(..., .sep = .sep, .envir = .envir),
    domain = .domain
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else x
}

is_cohesion_matrix <- function(x) inherits(x, "cohesion_matrix")

check_cohesion_matrix <- function(x) {
  if (!is_cohesion_matrix(x)) {
    stop_glue("You input an object of type:\n * {class(x)}\nThis funtion is ",
    "expecting an input of type `cohesion_matrix`. See the `cohesion_matrix` ",
    "function or use `as_cohesion_matrix()` to convert an object of type `matrix`",
    "to a `cohesion_matrix`.")
  }
}

check_dist <- function(x) {
  if (inherits(x, "dist")) {
    return(round(as.matrix(x), 15))
  }

  if (!is.matrix(x)) {
    stop_glue("`d` is not a distance matrix or `dist` object. You input an ",
    "object of type:\n * {class(x)}\nPlease provide a distance matrix or ",
    "`dist` object")
  }

  if (!isSymmetric.matrix(x)) {
    stop_glue("`d` is not a symmetric square matrix. Please provide a distance ",
    "matrix or `dist` object")
  }

  return(round(x, 15))
}
