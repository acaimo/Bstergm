random_signed <- function(n, p.active = 0.4, p.pos = 0.5) {
  A <- matrix(stats::rbinom(n * n, 1, p.active), n)
  A[lower.tri(A, diag = TRUE)] <- 0
  A <- A + t(A)
  S <- matrix(sample(c(-1L, 1L), n * n, TRUE, prob = c(1 - p.pos, p.pos)), n)
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  m <- A * S
  storage.mode(m) <- "integer"
  m
}

gw <- function(s, a) exp(a) * (1 - (1 - exp(-a))^s)

brute_signed_stats <- function(m, a = 0.6) {
  P <- (m > 0) * 1; N <- (m < 0) * 1
  ut <- upper.tri(m)
  spPP <- P %*% P; spNN <- N %*% N; spPN <- P %*% N + N %*% P
  c(gwesf_pos = sum(gw(spPP[ut & P == 1], a)),
    gwesf_neg = sum(gw(spPP[ut & N == 1], a)),
    gwese_pos = sum(gw(spNN[ut & P == 1], a)),
    gwese_neg = sum(gw(spNN[ut & N == 1], a)),
    mixed_pos = sum(gw(spPN[ut & P == 1], a)),
    mixed_neg = sum(gw(spPN[ut & N == 1], a)),
    gwdegree_neg = sum(gw(rowSums(N), a)),
    esp2 = sum(spNN[ut & P == 1] == 2),
    degree_neg1 = sum(rowSums(N) == 1))
}

native_formula <- ~ gwesf_pos(0.6) + gwesf_neg(0.6) + gwese_pos(0.6) + gwese_neg(0.6) +
  gwesp_signed(0.6, "+", "+-") + gwesp_signed(0.6, "-", "+-") + gwdegree_neg(0.6) +
  esp_signed(2, "+", "--") + degree_neg(1)

small_senate <- function(n = 30, times = 1:3) {
  data(senators70_am, package = "Bstergm", envir = environment())
  keep <- seq_len(n)
  va <- attr(senators70_am$S99, "senators")[keep, ]
  signed_series(lapply(senators70_am[times], function(m) m[keep, keep]), vertex.attr = va)
}
