# Test combine Functions
library(microbenchmark)
library(rstan)
library(doParallel)
library(irtpar)


# Increasing number of iterations
n <- 1e3
sub_post_1 <- list(matrix(rnorm(n), nc = 10),
                  matrix(rnorm(n), nc = 10),
                  matrix(rnorm(n), nc = 10),
                  matrix(rnorm(n), nc = 10)
                 )

n <- 5e3
sub_post_2 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 1e4
sub_post_3 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 5e4
sub_post_4 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 1e5
sub_post_5 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 5e5
sub_post_6 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 1e6
sub_post_7 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 5e6
sub_post_8 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

n <- 1e7
sub_post_9 <- list(matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10),
                   matrix(rnorm(n), nc = 10)
)

out <- NULL
combine_p(sub_post_1, out)

microbenchmark(
  combine_p(sub_post_1),
  combine_sp(sub_post_1),
  combine_np(sub_post_1),
  times = 2
  )

microbenchmark(
  combine_p(sub_post_2),
  combine_sp(sub_post_2),
  combine_np(sub_post_2),
  times = 10
)

microbenchmark(
  combine_p(sub_post_3),
  combine_sp(sub_post_3),
  combine_np(sub_post_3),
  times = 2
)

microbenchmark(
  combine_p(sub_post_4),
  combine_sp(sub_post_4),
  combine_np(sub_post_4),
  times = 2
)

# Increasing number of sub posteriors
sub_post_2 <- vector(mode = "list", length = 2)
sub_post_2 <- lapply(sub_post_2, function(x) matrix(rnorm(n), nc = 10))

sub_post_5 <- vector(mode = "list", length = 5)
sub_post_5 <- lapply(sub_post_5, function(x) matrix(rnorm(n), nc = 10))

sub_post_10 <- vector(mode = "list", length = 10)
sub_post_10 <- lapply(sub_post_10, function(x) matrix(rnorm(n), nc = 10))

sub_post_20 <- vector(mode = "list", length = 20)
sub_post_20 <- lapply(sub_post_20, function(x) matrix(rnorm(n), nc = 10))

sub_post_50 <- vector(mode = "list", length = 50)
sub_post_50 <- lapply(sub_post_50, function(x) matrix(rnorm(n), nc = 10))

sub_post_100 <- vector(mode = "list", length = 100)
sub_post_100 <- lapply(sub_post_100, function(x) matrix(rnorm(n), nc = 10))

microbenchmark(
  combine_p(sub_post_2),
  combine_p(sub_post_5),
  combine_p(sub_post_10),
  combine_p(sub_post_20),
  combine_p(sub_post_50),
  combine_p(sub_post_100),
  times = 2
  )

