# load the source code of the functions to be tested
source("ars.R")

# context("Tests make_h_x")
context("Tests make_h_x")

test_that("make_hx creates proper log(function)", {
  h <- make_hx(dnorm)
  
  expect_is(h,"function")
  expect_equal(h(1), log(dnorm(1)))
})

# context("Tests calc_z_j")
context("Tests calc_z_j")

test_that("calc_z_j returns expected values for normal distribution", {
  t_k = c(-5, 0.0, 5)
  j <- 2
  h_x <- make_hx(dnorm)
  expect_type(calc_z_j(t_k, j, h_x), 'double')
  expect_equal(calc_z_j(t_k, j, h_x), 2.5)

  t_k = c(-5, 0.1, 5)
  expect_gt(calc_z_j(t_k, j, h_x), 2.5)
})

# context("Tests get_xj_uk")
context("Tests get_xj_uk")

test_that("get_xj_uk returns expected index for normal distribution", {
  t_k = c(-5, 0.1, 5)
  x_star <- 1
  h_x <- make_hx(dnorm)
  
  expect_type(get_xj_uk(t_k, x_star, h_x), 'double')
  expect_equal(get_xj_uk(t_k, x_star, h_x), 2)
})

# context("Tests u_k")
context("Tests u_k")

test_that("u_k returns expected mx+b line for normal distribution", {
  t_k = c(-5, 0, 5)
  x_star <- -2.9
  h_x <- make_hx(dnorm)
  line <- u_k(t_k, x_star, h_x)
  
  expect_is(line, 'function')
  expect_equal(line(1)-line(0),5)
  
  x_star <- 1
  line2 <- u_k(t_k, x_star, h_x)
  expect_equal(line2(0),log(dnorm(0)))
})

# context("Tests get_xj_lk")
context("Tests get_xj_lk")

test_that("get_xj_lk returns expected index for normal distribution", {
  t_k = c(-5, 0.1, 5)
  x_star <- 1
  
  expect_type(get_xj_lk(t_k, x_star), 'double')
  expect_equal(get_xj_lk(t_k, x_star), 2)
})

# context("Tests s_k)
context("Tests s_k")

test_that("s_k returns proper function and integrals that sum to 1",{
  t_k = c(-5, 0.1, 5)
  h_x <- make_hx(dnorm)
  
  s <- s_k(t_k, h_x, -10, 10)
  expect_is(s[[1]],"function")
  expect_lt(s[[1]](-4),s[[1]](-1))
  expect_equal(sum(s[[2]]),1)
})

# context("Tests is_concave)
context("Tests is_concave")

test_that("is_concave properly tests for concavity",{
  h_x <- make_hx(dnorm)
  bool <- is_concave(-10, h_x, 10)
  
  expect_is(bool,"logical")
  expect_equal(bool, TRUE)
  
  f <- function(x) {
    f <- x^2
  }
  h_x2 <- make_hx(f)
  bool2 <- is_concave(-10, h_x2, 10)
  expect_equal(bool2, FALSE)
})

# context("Tests find_sp")
context("Tests find_sp")

test_that("find_sp returns valid starting values",{
  h_x <- make_hx(dnorm)
  sps <- find_sp(h_x, -10, 10) 
  
  expect_is(sps, "list")
  expect_lte(sps[[1]],sps[[2]])
})


# context("Tests check_boundary")
context("Tests check_boundary")

test_that("check_boundary changes boundaries in appropriate cases", {
  h_x <- make_hx(dnorm)
  cb <- check_boundary(-10, 10, h_x)
  
  expect_is(cb, "list")
  expect_equal(cb[[4]], FALSE)
  
  cb2 <- check_boundary(-200, 200, h_x)
  
  expect_equal(cb2[[4]], TRUE)
})

# context("Tests is_unif_exp")
context("Tests is_unif_exp")

test_that("is_unif_exp catches uniform or exponential distributions", {
  h_x <- make_hx(dnorm)
  check <- is_unif_exp(-2, h_x, 8) 
  
  expect_equal(check[[2]],FALSE)
  expect_equal(check[[4]],FALSE)
  
  h_x2 <- make_hx(dexp)
  check2 <- is_unif_exp(2, h_x2, 8) 
  expect_equal(check2[[2]],TRUE)
  expect_equal(check2[[4]],FALSE)
  
  h_x3 <- make_hx(dunif)
  check3 <- is_unif_exp(0.1, h_x3, 0.8) 
  expect_equal(check3[[2]],FALSE)
  expect_equal(check3[[4]],TRUE)
})


# context("Tests ARS")
context("Tests ARS")

test_that("ARS returns accurate normal distributions", {
  x <- ARS(dnorm, 100, sp = c(-5, 5), l_bound = -10, u_bound = 10)
  x_true <- rnorm(100)
  ks.norm <- ks.test(x, x_true)$p.value
  
  expect_type(x, 'double')
  expect_gt(ks.norm, 0.01)
})

test_that("ARS returns accurate chisquared distributions", {
  f <- function(x) {
    f <- dchisq(x,df=5)
  }
  x <- ARS(f, 300, sp = c(2, 8), l_bound = 0, u_bound = 30)
  x_true <- rchisq(300, df = 5)
  ks.chisq <- ks.test(x, x_true)$p.value
  
  expect_gt(ks.chisq, 0.01)
})

test_that("ARS returns accurate gamma distributions", {
  g <- function(x) {
    g <- dgamma(x, shape = 2, rate = 2)
  }
  x <- ARS(g, 300, sp = c(2, 8), l_bound = 0, u_bound = 30)
  x_true <- rgamma(300, shape = 2, rate = 2)
  ks.gamma <- ks.test(x, x_true)$p.value

  expect_gt(ks.gamma, 0.01)
})



