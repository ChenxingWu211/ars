library(pracma)
library(RGeode)

# This function takes a function and wraps it to be the log of the original function
# Input: f_x: A function to be calculated
# A wrapper function to access the logged function
make_hx = function(f_x){
  h_x = function(x){
    log(f_x(x))
  }
  
  return(h_x)
}

calc_z_j = function(t_k, j, h_x){
  x_j = t_k[j]
  x_j_plus_1 = t_k[j + 1]
  
  up = h_x(x_j_plus_1) - h_x(x_j) - 
    x_j_plus_1*(grad(h_x, x_j_plus_1)) +
    x_j*(grad(h_x, x_j))
  down = grad(h_x, x_j) - grad(h_x, x_j_plus_1)
  
  z_j = up/down
  return(z_j)
}

get_xj_uk = function(t_k, x_star, h_x){
  curr_idx = 1
  
  # Check which interval in z_j x_star belongs to 
  while (curr_idx < length(t_k)) {
    if (x_star < calc_z_j(t_k, curr_idx, h_x)) {
      return(curr_idx)
    }
    else{
      curr_idx = curr_idx + 1
    }
  }
  return(length(t_k))
}

u_k = function(t_k, x_star, h_x){
  x_j = t_k[get_xj_uk(t_k, x_star, h_x)]
  
  m = grad(h_x, x_j)
  b = h_x(x_j) - x_j*grad(h_x, x_j)
  
  u_k_line = function(x){
    x_line = m*x + b
    return(x_line)
  }
  return(u_k_line)
}

get_xj_lk = function(t_k, x_star){
  curr_idx = 2
  
  # Check with interval in t_k x_star belongs to 
  while (curr_idx <= length(t_k)) {
    if (x_star < t_k[curr_idx]) {
      return(curr_idx - 1)
    }
    else{
      curr_idx = curr_idx + 1
    }
  }
  return(curr_idx - 1)
}

l_k = function(t_k, x_star, h_x){
  if (x_star < t_k[1]) {
    result <- list(is_finite = FALSE)
    return(result)
  }
  
  if (x_star > t_k[length(t_k)]) {
    result = list(is_finite = FALSE)
    return(result)
  }
  
  j = get_xj_lk(t_k, x_star)
  x_j = t_k[j]
  x_j_plus_1 = t_k[j + 1]
  
  m2 = (h_x(x_j_plus_1) - h_x(x_j))/(x_j_plus_1 - x_j)
  b2 = (h_x(x_j)*x_j_plus_1 - h_x(x_j_plus_1)*x_j)/(x_j_plus_1 - x_j)
  
  l_k_line = function(x){
    x2_line = m2*x + b2
    return(x2_line)
  }
  result = list(is_finite = TRUE, line = l_k_line)
  return(result)
}

s_k = function(t_k, h_x, lower_bound, upper_bound){
  
  k = length(t_k)
  
  # calculate integration constant for denominator
  # integral of each s_k
  integrals <- matrix(0,k)
  # entire sum of integrals
  integ_const <- 0
  
  # loop through j many sections of u_k, and integrate u_j from z_(j-1) to z_j
  # sum them to get entire integral of exp(u_k)
  for (j in 1:k) {
    
    # find proper integration bounds of section u_j
    if (j == 1) {
      lower_limit = lower_bound
    } else {
      lower_limit = calc_z_j(t_k, j - 1, h_x)
    }
    if (j == k) {
      upper_limit = upper_bound
    } else {
      upper_limit = calc_z_j(t_k, j, h_x)
    }
    
    # get sections line
    func_line = u_k(t_k, t_k[j], h_x)
    
    # pull out slope m
    m = func_line(1) - func_line(0)
    # pull out y-int b
    b = func_line(0)
    
    # exp(u_j)
    f = function(x) {
      return(exp(m*x + b))
    }
    
    # integrate jth section, save to vector
    integrals[j] = integrate(f, lower_limit, upper_limit)$value
    # sum to total integ constant
    integ_const = integ_const + integrals[j]
    
  }
  
  # returns this s_k function that can take in x value
  s_k_func = function(x){
    # u_k at x
    u_k_line = u_k(t_k, x, h_x)
    
    # pull out slope m
    m2 = u_k_line(1) - u_k_line(0)
    # pull out y-int b
    b2 = u_k_line(0)
    
    # exponentiate
    f2 = function(x2) {
      return(exp(m2*x2 + b2))
    }
    
    # return calculated value at x
    s_k_of_x = f2(x) / integ_const
    return(s_k_of_x)
  }
  
  # normalized integral of each j section
  integrals = integrals/integ_const
  
  # return s_k function and normalize integrals
  return(list(s_k_func,integrals))
}

sampler = function(t_k, h_x, lower_bound, upper_bound){
  k = length(t_k)
  s = s_k(t_k, h_x, lower_bound, upper_bound)
  
  # determine which jth section
  u1 = runif(1, 0, 1 - 2*.Machine$double.eps)
  partition = sum(cumsum(s[[2]]) < u1) + 1
  
  # u_k at x
  u_k_line = u_k(t_k, t_k[partition], h_x)
  # pull out slope m
  m = u_k_line(1) - u_k_line(0)
  # pull out z_(j-1)
  z1 = NA
  if (partition == 1) {
    z1 = lower_bound
  } else{
    z1 = calc_z_j(t_k, partition - 1, h_x)
  }
  # if tangent is straight line, sample uniformly
  x = NA
  if (m == 0) {
    if (partition == k) {
      z2 = upper_bound
    } else{
      z2 = calc_z_j(t_k, partition, h_x)
    }
    x = runif(1,z1,z2)
  } else{
    # pull out y-int b
    b = u_k_line(0)
    
    integ_const_j = s[[2]][partition]
    
    u2 = runif(1)
    x = (log(integ_const_j*m*u2 + exp(m*z1 + b)) - b)/m
  }  
  return(x)
}


find_sp = function(h_x, l_bound, u_bound) {
  # Start finding from the mid point
  start_val = (u_bound - l_bound)/2 + l_bound
  step = (u_bound - l_bound)/10
  x1 = start_val - step
  xk = start_val + step
  
  # Find a place where h`(x1) > 0 and h`(xk) < 0
  while ((x1 >= l_bound) & (grad(h_x, x1) <= 0)) {
    x1 = x1 - step
  }
  
  while ((xk <= u_bound) & (grad(h_x, xk) >= 0)) {
    xk = xk + step
  }
  
  retval = list(x1 = x1, xk = xk)
  return(retval)
}

check_boundary = function(in_l, in_u, h_x) {
  step = (in_u - in_l)/1000
  new_l = in_l
  changed_bound = FALSE
  
  # Check if upper and lower bound return finite values for h_x
  while (!is.finite(h_x(new_l)) & !is.finite(grad(h_x, new_l))) {
    changed_bound = TRUE
    if (new_l >= in_u) {
      retval = list(reach_bound = TRUE, new_u = NA, new_l = NA, changed_bound = changed_bound)
      return(retval)
    }
    new_l = new_l + step
  }
  
  new_u = in_u
  while (!is.finite(h_x(new_u)) & !is.finite(grad(h_x, new_u))) {
    changed_bound = TRUE
    if (new_u <= new_l) {
      retval = list(reach_bound = TRUE, new_u = NA, new_l = NA, changed_bound = changed_bound)
      return(retval)
    }
    new_u = new_u - step
  }
  
  # Return new bounds 
  retval = list(reach_bound = FALSE, new_u = new_u, new_l = new_l, changed_bound = changed_bound)
  return(retval)
}

is_unif_exp = function(start_point, h_x, end_point) {
  step = (end_point - start_point) / 1000
  x1 = start_point
  x2 = x1 + step
  x3 = x2 + step
  
  first_slope = (h_x(x2) - h_x(x1)) / step
  second_slope = (h_x(x3) - h_x(x2)) / step
  
  check_unif = FALSE
  check_exp = FALSE
  #first check is to see if the first two slopes are equal
  #to the degree of machine accuracy
  if (isTRUE(all.equal(first_slope, second_slope, 1e-08))) {
    #if basically equal to 0 in terms of machine accuracy,
    #then test if uniform
    if (isTRUE(all.equal(first_slope, 0, 1e-08)))
    {
      check_unif = TRUE
      while (check_unif && x2 < end_point) {
        this_slope = (h_x(x2) - h_x(x1)) / step
        check_unif = check_unif && isTRUE(all.equal(this_slope, 0))
        x1 = x1 + step
        x2 = x2 + step
      }
      return_list =
        list(
          is_exp = check_exp,
          lambda = FALSE,
          is_unif = check_unif
        )
      return(return_list)
    }
    else{
      #if first 2 slopes are equal to each other and not 0
      #check to see if exponential to see if all the other
      # slopes are equal
      check_exp = TRUE
      lambda = first_slope
      end_point = min(end_point, 700/(-lambda))
      while (check_exp && x2 < end_point) {
        this_slope = (h_x(x2) - h_x(x1)) / step
        check_exp = check_exp &&
          isTRUE(all.equal(this_slope, lambda, 1e-08))
        x1 = x1 + step
        x2 = x2 + step
      }
      return_list =
        list(
          is_exp = check_exp,
          lambda = -lambda,
          is_unif = check_unif
        )
      return(return_list)
    }
  }
  else{
    return_list =
      list(
        is_exp = FALSE,
        lambda = NA,
        is_unif = FALSE
      )
    return(return_list)
  }
}

is_concave = function(start_point, h_x, end_point) {
  step = (end_point - start_point)/100
  x1 = start_point
  x2 = x1 + step
  x3 = x2 + step
  
  while (x3 < end_point) {
    slope1 = (h_x(x2) - h_x(x1))/step
    slope2 = (h_x(x3) - h_x(x2))/step
    
    if (slope2 > slope1) {
      return(FALSE)
    }
    
    x1 = x2
    x2 = x3
    x3 = x3 + step
  }
  
  return(TRUE)
}

#' The Adaptive Rejection Sampling Function
#'
#' This function allows you to sample from any univariate log-concave probability density function f(x).
#' @param fx The log-concave function to be sampled from. Needs to be in the format of a function.
#' @param n Number of samples. Needs to be an integer.
#' @param sp Vector of Starting points. Needs to be at least length 2. Default to NA.
#' @param l_bound Value for lower bound. Needs to be a number. Defaults to -Inf.
#' @param u_bound Value for upper bound. Needs to be a number. Defaults to Inf.
#' @return A vector of samples from f(x).
#' @source Gilks, W. R. and Wild, P. (1992) Adaptive Rejection Sampling for Gibbs Sampling. Appl.Statist., 41 337 - 348.
#' @export
#' @examples
#' ars(dnorm, 100, l_bound = -10, u_bound = 10)
#' 
ars = function(fx, n, sp = NA, l_bound = -Inf, u_bound = Inf) {
  # Check if inputs are valid
  
  if (!is.function(fx)) {
    print("Please input a valid function.")
    return()
  }
  
  if (!is.numeric(n)) {
    print("Please input a number for n.")
    return()
  }
  
  if (n %% 1 != 0) {
    print("n has to be an integer.")
    return()
  }
  
  if (!is.numeric(l_bound)) {
    print("Please input a number for l_bound.")
    return()
  }
  
  if (!is.numeric(u_bound)) {
    print("Please input a number for u_bound.")
    return()
  }
  
  
  if ((!is.na(sp)) & (length(sp) < 2)) {
    print("Please input vector with length larger than or equal to 2 as starting points. Or don't input sp if you don't want to specify the starting points.")
    return()
  }
  
  max_iter = n * 100
  h_x = make_hx(fx)
  
  #Check if the bound is reasonable, shrink it if not
  if (!is.finite(l_bound)) {
    l_bound = -2000
  }
  
  if (!is.finite(u_bound)) {
    u_bound = 1500
  }
  
  check_out = check_boundary(l_bound, u_bound, h_x)
  if (check_out$reach_bound) {
    print("The boundary is not valid. Try again :)")
    return()
  }
  
  if (check_out$changed_bound) {
    l_bound = check_out$new_l
    u_bound = check_out$new_u
    print(paste("The boundary was too wide. It has been shrunk to (", l_bound, ", ", u_bound, ")", sep = ""))
  }
  
  # Catch uniform and exponential distributions
  which_type = is_unif_exp(l_bound, h_x, u_bound)
  
  if (which_type$is_unif) {
    return(runif(n, min = l_bound, max = u_bound))
  }
  
  if (which_type$is_exp) {
    return(rexptr(n = n, lambda = which_type$lambda, range = c(l_bound, u_bound)))
  }
  
  # Check if fx is log concave
  if (!is_concave(l_bound, h_x, u_bound)) {
    print("The input function should be log concave. Try again :)")
    return()
  }
  
  # Give starting points if the user doesn't, check if they are valid if the user gives starting points
  if (is.na(sp)) {
    x_found = find_sp(h_x, l_bound, u_bound)
    sp = c(x_found$x1, x_found$xk)
    
  }
  else{
    # check if h_prime(sp[1]) is > 0
    if ((!is.finite(l_bound)) & (grad(h_x, sp[1]) <= 0)) {
      print("h'(x_1) should be larger than 0. Try again :)")
      return()
    }
    
    # check if h_prime(sp[length(sp)]) is < 0
    if ((!is.finite(u_bound)) & (grad(h_x, sp[length(sp)]) >= 0)) {
      print("h'(x_k) should be smaller than 0. Try again :)")
      return()
    }
  }
  
  t_k = sort(sp)
  sampled_x = c()
  
  num_iter = 0
  
  # Do the actual sampling
  while ((length(sampled_x) < n) & (num_iter < max_iter)) {
    need_sample = min((n - length(sampled_x)), n/20)
    sampled_temp = c()
    
    for (i in seq(1, need_sample)) {
      sampled_temp[i] = sampler(t_k, h_x, l_bound, u_bound)
    }
    
    # Check if x_star can use or we need to update functions
    for (x_star in sampled_temp) {
      sampled_w = runif(1)
      
      u_k_func = u_k(t_k, x_star, h_x)
      l_k_func = l_k(t_k, x_star, h_x)
      
      # l_k is infinite if x_star is smaller than x1 or larger than xk
      l_k_func_val = 0
      if (!(l_k_func$is_finite)) {
        l_k_func_val = -Inf
      }
      else{
        l_k_func_val = l_k_func$line(x_star)
      }
      
      # Also check for log concave along theway
      if ((h_x(x_star) < l_k_func_val) | (h_x(x_star) > u_k_func(x_star))) {
        print("The input function should be log concave. Try again :)")
        return()
      }
      
      # Accept, reject, or update functions
      if (sampled_w <= exp(l_k_func_val - u_k_func(x_star))) {
        sampled_x[length(sampled_x) + 1] = x_star
      }
      else if (sampled_w <= exp((h_x(x_star) - u_k_func(x_star)))) {
        sampled_x[length(sampled_x) + 1] = x_star
        
        t_k[length(t_k) + 1] = x_star
        t_k = sort(t_k)
        break
      }
      else{
        
      }
    }
    
    # In case it takes too many interations
    num_iter = num_iter + 1
  }
  
  return(sampled_x)
}
