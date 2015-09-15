library(testthat)
library(aeep)

test_that("combine_p produces output of correct dimensionality", {
    set.seed(123)
    post_list <- list(matrix(rnorm(1000), nr = 100, nc = 10),
                      matrix(rnorm(1000), nr = 100, nc = 10))
    combine_p_obj <- combine_p(post_list)
    expect_equal(dim(combine_p_obj$post_means), c(10, 1))
    expect_equal(dim(combine_p_obj$post_variance), c(10, 10))
})