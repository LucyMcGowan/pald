## code to prepare `DATASET` dataset goes here



exdata1 <- cbind(c(-2, -1, -1, 0, 1.4, 1.4, 1.4, 0.5), c(0, 1.28, -1.28, 0, -.9, 0, .9, 3))
exdata2 <- cbind(c(48.00, 22.83, 85.77, 62.84, 39.97, 446.23, 149.83, 20.00, 314.12, 549.29, 349.56, 249.00, 399.49, 466.71, 640.56, 787.00), c(18.00, 29.41, 3.27, 37.59, 0.00, 73.54, 69.56, 117.16, 114.39, 633.39, 670.00, 450.33, 533.55, 350.54, 411.00, 613.42))
exdata3<-as.matrix(read.csv("~/Dropbox (Amherst College)/git/pald/pald/data-raw/fig4f_data_eight.csv", header=TRUE))

cognate_data <- read.csv("~/Dropbox (Amherst College)/git/pald/pald/data-raw/cognate_data.csv", header=TRUE)
cognate<-cognate_data[, 2:2666]
rownames(cognate)<-cognate_data[, 1]
cognate_dist<-dist(cognate)

usethis::use_data(exdata1, overwrite = TRUE)
usethis::use_data(exdata2, overwrite = TRUE)
usethis::use_data(exdata3, overwrite = TRUE)
usethis::use_data(cognate_dist, overwrite = TRUE)
