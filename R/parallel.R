.rng_streams <- function(seed, njobs) {
  old.kind <- RNGkind()
  old.seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old.kind))
    if (is.null(old.seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old.seed, envir = globalenv())
  })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- get(".Random.seed", envir = globalenv())
  out <- vector("list", njobs)
  for (j in seq_len(njobs)) {
    out[[j]] <- s
    s <- parallel::nextRNGStream(s)
  }
  out
}

.run_jobs <- function(jobs, FUN, cores = 1L) {
  run1 <- function(job) {
    old.kind <- RNGkind()
    old.seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
      do.call(RNGkind, as.list(old.kind))
      if (is.null(old.seed)) {
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
          rm(".Random.seed", envir = globalenv())
      } else assign(".Random.seed", old.seed, envir = globalenv())
    })
    assign(".Random.seed", job$rng, envir = globalenv())
    FUN(job)
  }
  cores <- max(1L, min(as.integer(cores), length(jobs)))
  if (cores == 1L) return(lapply(jobs, run1))

  if (.Platform$OS.type == "unix" && !isTRUE(getOption("Bstergm.psock", FALSE))) {
    res <- parallel::mclapply(jobs, run1, mc.cores = cores, mc.preschedule = FALSE,
                              mc.set.seed = FALSE, mc.silent = FALSE)
  } else {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(p) invisible(.libPaths(p)), .libPaths())
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(Bstergm)))
    res <- parallel::parLapply(cl, jobs, function(job) tryCatch(run1(job), error = function(e) {
      structure(conditionMessage(e), class = "try-error", condition = e)
    }))
  }
  failed <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), TRUE)
  if (any(failed)) {
    msg <- unique(vapply(res[failed], function(r)
      if (is.null(r)) "a worker process died (out of memory?)" else conditionMessage(attr(r, "condition")), ""))
    stop("Sampling failed in ", sum(failed), " job(s): ", paste(msg, collapse = "; "), call. = FALSE)
  }
  res
}

.chunks <- function(n, k) {
  k <- max(1L, min(as.integer(k), n))
  split(seq_len(n), ceiling(seq_len(n) * k / n))
}

.sim_jobs <- function(groups, nsim, cores, seed) {
  streams <- .rng_streams(seed, length(groups) * nsim)
  jobs <- list()
  for (g in seq_along(groups)) {
    for (b in .chunks(nsim, cores)) {
      rngs <- streams[(g - 1L) * nsim + b]
      jobs[[length(jobs) + 1L]] <- list(component = groups[[g]], rows = b, rngs = rngs, rng = rngs[[1L]])
    }
  }
  jobs
}

.use_stream <- function(s) assign(".Random.seed", s, envir = globalenv())
