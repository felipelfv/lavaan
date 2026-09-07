# check if a fitted model is admissible
lav_object_post_check <- function(object) {
  stopifnot(inherits(object, "lavaan"))
  lavpartable <- object@ParTable
  lavmodel <- object@Model
  lavdata <- object@Data

  var_ov_ok <- var_lv_ok <- result_ok <- TRUE
  var_na <- FALSE

  # 1a. check for negative variances ov
  var_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lav_object_vnames(object, "ov") &
    lavpartable$lhs == lavpartable$rhs)
  if (any(is.na(lavpartable$est[var_idx]))) {
    # perhaps estimator = "IV" + stage 1 only
    var_na <- TRUE
  } else if (length(var_idx) > 0L && any(lavpartable$est[var_idx] < 0.0)) {
    result_ok <- var_ov_ok <- FALSE
    lav_msg_warn(gettext("some estimated ov variances are negative"))
  }

  # 1b. check for negative variances lv
  var_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lav_object_vnames(object, "lv") &
    lavpartable$lhs == lavpartable$rhs)
  if (any(is.na(lavpartable$est[var_idx]))) {
    # perhaps estimator = "IV" + stage 1 only
    var_na <- TRUE
  } else if (length(var_idx) > 0L && any(lavpartable$est[var_idx] < 0.0)) {
    result_ok <- var_lv_ok <- FALSE
    lav_msg_warn(gettext("some estimated lv variances are negative"))
  }

  # 2. is cov.lv (PSI) positive definite? (only if we did not already warn
  # for negative variances)
  if (!var_na && var_ov_ok && var_lv_ok) {
    # keep the dummy lv's: residual covariances among observed endogenous
    # variables live in psi, not in theta, and are dropped from cov.lv
    eta <- lav_model_veta(lavmodel, remove_dummy_lv = FALSE)
    for (b in seq_along(eta)) {
      if (nrow(eta[[b]]) == 0L) next
      txt_block <- if (length(eta) > 1L) {
        gettextf("in block %s", lavdata@block.label[b])
      } else ""
      eigvals <- eigen(eta[[b]], symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
        lav_msg_warn(gettextf(
          "covariance matrix of latent variables is not positive definite %s;
          use lavInspect(fit, \"est\")$psi to investigate.", txt_block
        ))
        result_ok <- FALSE
      }
    }
  }

  # 3. is THETA positive definite (but only for numeric variables)
  # and if we have not already warned for negative ov variances
  if (!var_na && var_ov_ok) {
    mm_theta <- lavTech(object, "theta")
    for (b in seq_along(mm_theta)) {
      num_idx <- lavmodel@num.idx[[b]]
      if (length(num_idx) > 0L) {
        txt_block <- if (length(mm_theta) > 1L) {
          gettextf("in block %s", lavdata@block.label[b])
        } else ""
        eigvals <- eigen(mm_theta[[b]][num_idx, num_idx, drop = FALSE],
          symmetric = TRUE,
          only.values = TRUE
        )$values
        if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
          lav_msg_warn(gettextf(
            "the covariance matrix of the residuals of the observed variables
            (theta) is not positive definite %s; use lavInspect(fit, \"theta\")
            to investigate.", txt_block))
          result_ok <- FALSE
        }
      }
    }
  }

  result_ok
}
