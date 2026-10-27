.layer_kind <- function(spec, env) {
  val <- tryCatch(eval(spec, env), error = function(e) NULL)
  if (is.list(val) && length(val) == 1) val <- val[[1]]
  if (!inherits(val, "formula")) return(NA_character_)
  rhs <- gsub("\\s", "", paste(deparse(val[[length(val)]]), collapse = ""))
  repeat {
    if (startsWith(rhs, "(") && endsWith(rhs, ")")) rhs <- substr(rhs, 2, nchar(rhs) - 1) else break
  }
  if (rhs == "positive") return("positive")
  if (rhs == "active") return("active")
  if (rhs %in% c("active&!positive", "!positive&active", "active&!(positive)", "!(positive)&active"))
    return("negative")
  NA_character_
}

.arg <- function(cl, name, pos = NULL) {
  a <- as.list(cl)[-1]
  nm <- names(a)
  if (is.null(nm)) nm <- rep("", length(a))
  if (name %in% nm) return(a[[which(nm == name)[1]]])
  if (!is.null(pos)) {
    unnamed <- which(nm == "")
    if (length(unnamed) >= pos) return(a[[unnamed[pos]]])
  }
  NULL
}

.translate_formula <- function(formula) {
  env <- environment(formula)
  terms <- tryCatch(list_rhs.formula(formula), error = function(e) NULL)
  if (is.null(terms) || any(attr(terms, "sign") < 0)) return(NULL)
  out <- list()
  dropped <- 0L
  for (tm in terms) {
    if (!is.call(tm)) return(NULL)
    fn <- as.character(tm[[1]])
    if (fn == "L") {
      inner <- .arg(tm, "formula", 1)
      kind <- .layer_kind(.arg(tm, "Ls", 2), env)
      inner <- tryCatch(eval(inner, env), error = function(e) NULL)
      if (!inherits(inner, "formula") || is.na(kind)) return(NULL)
      if (kind == "active") { dropped <- dropped + 1L; next }
      if (kind != "positive") return(NULL)
      it <- list_rhs.formula(inner)
      if (any(attr(it, "sign") < 0)) return(NULL)
      if (any(unlist(lapply(it, all.names)) %in% .multilayer_term_names())) return(NULL)
      out <- c(out, unclass(it))
    } else if (fn %in% c("gwespL", "dgwespL")) {
      decay <- tryCatch(eval(.arg(tm, "decay", 1), env), error = function(e) NULL)
      fixed <- tryCatch(eval(.arg(tm, "fixed", 2), env), error = function(e) NULL)
      inord <- tryCatch(eval(.arg(tm, "L.in_order"), env), error = function(e) NULL)
      if (!is.numeric(decay) || length(decay) != 1 || !isTRUE(fixed) || isTRUE(inord)) return(NULL)
      known <- c("decay", "fixed", "L.base", "Ls.path", "L.in_order", "")
      if (!all(names(as.list(tm)[-1]) %in% known)) return(NULL)
      base <- .layer_kind(.arg(tm, "L.base"), env)
      path <- tryCatch(eval(.arg(tm, "Ls.path"), env), error = function(e) NULL)
      if (!is.list(path) || length(path) != 2 || !(base %in% c("positive", "negative"))) return(NULL)
      pk <- vapply(path, function(p) .layer_kind(p, environment()), "")
      if (any(is.na(pk)) || !all(pk %in% c("positive", "negative"))) return(NULL)
      b <- if (base == "positive") "+" else "-"
      p <- if (all(pk == "positive")) "++" else if (all(pk == "negative")) "--" else "+-"
      fname <- c("++" = "gwesf", "--" = "gwese", "+-" = NA)[[p]]
      out[[length(out) + 1]] <- if (is.na(fname)) call("gwesp_signed", decay, base = b, path = p)
                                else call(paste0(fname, if (b == "+") "_pos" else "_neg"), decay)
    } else return(NULL)
  }
  if (!length(out)) return(NULL)
  rhs <- Reduce(function(a, b) call("+", a, b), out)
  f <- eval(call("~", rhs))
  environment(f) <- env
  attr(f, "dropped") <- dropped
  f
}

.translate_formulas <- function(formulas, data) {
  tr <- lapply(formulas, .translate_formula)
  if (any(vapply(tr, is.null, TRUE))) return(NULL)
  tn <- transition_networks(data)
  ok <- tryCatch({
    all(vapply(seq_along(formulas), function(t) {
      all(vapply(c("F", "P"), function(U) {
        sl <- .target_spec(tn, t, "sign", U, "layer", data$vertex.attr)
        sn <- .target_spec(tn, t, "sign", U, "native", data$vertex.attr)
        a <- suppressMessages(summary(formulas[[t]], basis = sl$basis))
        b <- suppressMessages(summary(tr[[t]], basis = sn$basis))
        keep <- !grepl("^L\\(active\\)", names(a))
        length(b) == sum(keep) && isTRUE(all.equal(unname(a[keep]), unname(b), tolerance = 1e-8))
      }, TRUE))
    }, TRUE))
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  lapply(tr, function(f) { attr(f, "dropped") <- NULL; f })
}
