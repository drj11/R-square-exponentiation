#!/usr/bin/Rscript

# Response to
# https://eric.netlify.com/2017/08/08/taking-powers-of-a-matrix-in-r/

library(tidyverse)
library(microbenchmark)

same = function(e, w) {
  eps = max(abs(e - w))
  eps == 0
}

test_properties = function(f, A) {
  list(
    "A^4 == AAAA"       = same(f(A, 4), A %*% A %*% A %*% A),
    "A^1 == A"          = same(f(A, 1), A),
    "A^0 == I"          = same(f(A, 0), diag(1, nrow(A)))
  )
}

set.seed(1)
n = 3
A = matrix(sample(10, n^2, TRUE), n, n)

pow_iter = function (A, k) {
  if(k==0) {
    diag(1, nrow(A))
  } else {
    A %*% pow_iter(A, (k-1))
  }
}

# Calculate exponentatiation by repeated squaring.
# See SICP 2nd edition, 1.2.4
# https://mitpress.mit.edu/sicp/chapter1/node15.html

pow_sq_recursive = function (A, k) {
  if(k == 0) {
    # identity
    diag(1, nrow(A))
  } else if(k %% 2 == 1) {
    A %*% pow_sq_recursive(A, k-1)
  } else {
    pow_sq_recursive(A %*% A, k %/% 2)
  }
}

# Iterative version of pow_sq_recursive
pow_sq_ = function (A, k, f) {
  if(k == 0) {
    f
  } else if(k %% 2 == 1) {
    pow_sq_(A, k-1, f %*% A)
  } else {
    pow_sq_(A %*% A, k %/% 2, f)
  }
}

pow_sq_iter = function (A, k) {
  I = diag(1, nrow(A))
  pow_sq_(A, k, I)
}

`%^%` = pow_sq_recursive

# Eric's Eigenvalue method.
# https://eric.netlify.com/2017/08/08/taking-powers-of-a-matrix-in-r/
pow_eigen = function(A, k) {
  eig = eigen(A)
  stopifnot(length(unique(eig$values)) == nrow(A))

  P = eig$vectors
  D = diag(eig$values ^ k)
  Ak = P %*% D %*% solve(P)
  Re(Ak)
}

`%^^%` = pow_eigen

print("pow_iter")
test_properties(pow_iter, A)

print("pow_sq")
test_properties(pow_sq_recursive, A)

print("pow_eigen")
test_properties(pow_eigen, A)

times = 1:200 %>%
  map_df(~
    microbenchmark(
      "repeated squaring" = A %^% .,
      "diagonalization" = A %^^% .,
      times = 50
    ) %>%
    as_data_frame() %>%
    mutate(k = .x)
  ) %>%
  filter(time < quantile(time, 0.99))

ggplot(times, aes(x = k, y = time / 1e3, group = expr)) +
  geom_jitter(aes(color = expr), alpha = 0.1) +
  geom_smooth(method = "lm", color = "black", size = 0.1) +
# geom_point(data = df_intersect) +
# geom_text(aes(label = sprintf("(%0.0f, %0.0f)", k, time)),
#           data = df_intersect, 
#           vjust = 2, hjust = 0, size = 3) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_light() +
  theme(legend.position = c(1, 0), 
        legend.justification = c(1, 0),
        legend.background = element_rect(color = "grey70")) +
  labs(y = "Time (microseconds)", color = NULL,
       title = "Powers of a matrix",
       subtitle = paste0(
         "Time taken to raise a 3x3 matrix to the", 
         " kth power using different approaches")
       )
