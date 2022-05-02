test_that("check_dist works", {
  expect_error(check_dist("c"),  "`d` is not a distance matrix")
  expect_error(check_dist(1), "`d` is not a distance matrix")
  non_sym_mat <- matrix(c(1, 2, 3))
  d <- dist(non_sym_mat)
  expect_true(is.matrix(check_dist(d)))
})

test_that("coehsion_matrix gives expected matrix", {
  d <- matrix(rep(c(0, 1.5, 0, 1.5,
                    1.5, 0, 1.5, 0), 2), nrow = 4, byrow = TRUE)
  c <- matrix(rep(c(0.25, 0, 0.25, 0,
                    0, 0.25, 0, 0.25), 2), nrow = 4, byrow = TRUE)
  c <- as_cohesion_matrix(c)
  expect_equal(c, cohesion_matrix(d))
})

test_that("local_depths works", {
  D <- dist(exdata1)
  d_d <- local_depths(D)
  C <- cohesion_matrix(D)
  d_c <- local_depths(C)

  expect_equal(d_d, rowSums(C))
  expect_equal(d_c, rowSums(C))
})

test_that("strong_threshold works", {
  C <- cohesion_matrix(dist(exdata1))
  expect_equal(strong_threshold(C), mean(diag(C)) / 2)
})

test_that("cohesion_strong works", {
  C <- cohesion_matrix(dist(exdata2))
  C_strong <- C
  C_strong[C < strong_threshold(C)] <- 0
  expect_equal(C_strong,
               cohesion_strong(C, symmetric = FALSE))

  C_strong_sym <- pmin(C_strong, t(C_strong))
  expect_equal(C_strong_sym,
               cohesion_strong(C))
})

test_that("any_isolated works", {
  d <- data.frame(
    x1 = c(1, 2, 3, 6),
    x2 = c(2, 1, 3, 10)
    )
  D <- dist(d)
  C <- cohesion_matrix(D)
  expect_message(any_isolated(C), "These points")
  expect_true(any_isolated(C))
})

test_that("community_clusters works", {
  D <- dist(exdata1)
  C <- cohesion_matrix(D)
  cc <- community_clusters(C)
  expect_equal(cc$cluster, c(1, 1, 1, 1, 2, 2, 2, 3))
})
