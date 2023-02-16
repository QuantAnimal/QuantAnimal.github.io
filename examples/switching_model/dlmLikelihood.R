#Modified from dlmLL function from dlm package to calculate likelihood of each observation
dlmLikelihood <- function(y, mod) 
{
    lls =rep(0.0, length(y))

    eps <- .Machine$double.eps^0.3
    y <- as.matrix(y)
    n <- nrow(y)
    ll <- 0
    if (is.null(mod$JFF)) 
      tvFF <- FALSE
    else {
      tvFF <- TRUE
      nz <- mod$JFF != 0
      mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], 
                       mod$JFF[nz])
    }
    if (is.null(mod$JV)) 
      tvV <- FALSE
    else {
      tvV <- TRUE
      nz <- mod$JV != 0
      mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], 
                      mod$JV[nz])
    }
    if (is.null(mod$JGG)) 
      tvGG <- FALSE
    else {
      tvGG <- TRUE
      nz <- mod$JGG != 0
      mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], 
                       mod$JGG[nz])
    }
    if (is.null(mod$JW)) 
      tvW <- FALSE
    else {
      tvW <- TRUE
      nz <- mod$JW != 0
      mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], 
                      mod$JW[nz])
    }
    tvFV <- tvFF || tvV
    if (!tvV) {
      tmp <- La.svd(mod$V, nu = 0)
      Dv <- sqrt(tmp$d)
      if (any(Dv < eps)) {
        Dv <- pmax(Dv, eps)
        warning("a numerically singular 'V' has been slightly perturbed to make it nonsingular")
      }
      Dv.inv <- 1/Dv
      sqrtVinv <- Dv.inv * tmp$vt
      sqrtV <- Dv * tmp$vt
      if (!tvFF) 
        tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
    }
    if (!tvW) {
      svdW <- La.svd(mod$W, nu = 0)
      sqrtW <- sqrt(svdW$d) * svdW$vt
    }
    tmp <- La.svd(mod$C0, nu = 0)
    Ux <- t(tmp$vt)
    Dx <- sqrt(tmp$d)
    for (i in seq(length = n)) {
      if (tvFF) 
        mod$FF[mod$JFF[, -3, drop = FALSE]] <- mod$X[i, 
                                                     mod$JFF[, 3]]
      if (tvV) {
        mod$V[mod$JV[, -3, drop = FALSE]] <- mod$X[i, 
                                                   mod$JV[, 3]]
        tmp <- La.svd(mod$V, nu = 0)
        Dv <- sqrt(tmp$d)
        Dv.inv <- 1/Dv
        Dv.inv[abs(Dv.inv) == Inf] <- 0
        sqrtVinv <- Dv.inv * tmp$vt
        sqrtV <- sqrt(tmp$d) * tmp$vt
      }
      if (tvGG) 
        mod$GG[mod$JGG[, -3, drop = FALSE]] <- mod$X[i, 
                                                     mod$JGG[, 3]]
      if (tvW) {
        mod$W[mod$JW[, -3, drop = FALSE]] <- mod$X[i, 
                                                   mod$JW[, 3]]
        svdW <- La.svd(mod$W, nu = 0)
        sqrtW <- sqrt(svdW$d) * svdW$vt
      }
      if (tvFV) 
        tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
      if (!any(whereNA <- is.na(y[i, ]))) {
        a <- mod$GG %*% mod$m0
        tmp <- La.svd(rbind(Dx * t(mod$GG %*% Ux), sqrtW), 
                      nu = 0)
        Ux.prior <- t(tmp$vt)
        Dx.prior <- tmp$d
        f <- mod$FF %*% a
        tmp <- La.svd(rbind(Dx.prior * t(mod$FF %*% Ux.prior), 
                            sqrtV), nu = 0)
        Uy <- t(tmp$vt)
        Dy <- tmp$d
        D.inv <- 1/Dx.prior
        D.inv[abs(D.inv) == Inf] <- 0
        tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% Ux.prior, 
                            diag(x = D.inv, nrow = length(D.inv))), nu = 0)
        Ux <- Ux.prior %*% t(tmp$vt)
        Dx <- 1/tmp$d
        Dx[abs(Dx) == Inf] <- 0
        e <- as.matrix(y[i, ] - f)
        mod$m0 <- a + crossprod(Dx * t(Ux)) %*% tF.Vinv %*% 
          e
        ll <- ll + 2 * sum(log(Dy)) + crossprod(crossprod(Uy, 
                                                           e)/Dy)
        lls[i] <- 2 * sum(log(Dy)) + crossprod(crossprod(Uy, 
                                                         e)/Dy)
      }
      else {
        if (all(whereNA)) {
          mod$m0 <- mod$GG %*% mod$m0
          tmp <- La.svd(rbind(Dx * t(mod$GG %*% Ux), 
                              sqrtW), nu = 0)
          Ux <- t(tmp$vt)
          Dx <- tmp$d
        }
        else {
          good <- !whereNA
          tmp <- La.svd(mod$V[good, good], nu = 0)
          Dv <- sqrt(tmp$d)
          Dv.inv <- 1/Dv
          Dv.inv[abs(Dv.inv) == Inf] <- 0
          sqrtVinvTMP <- Dv.inv * tmp$vt
          tF.VinvTMP <- t(mod$FF[good, , drop = FALSE]) %*% 
            crossprod(sqrtVinvTMP)
          sqrtVTMP <- Dv * tmp$vt
          a <- mod$GG %*% mod$m0
          tmp <- La.svd(rbind(Dx * t(mod$GG %*% Ux), 
                              sqrtW), nu = 0)
          Ux.prior <- t(tmp$vt)
          Dx.prior <- tmp$d
          f <- mod$FF[good, , drop = FALSE] %*% a
          tmp <- La.svd(rbind(Dx.prior * t(mod$FF[good, 
                                                  , drop = FALSE] %*% Ux.prior), sqrtVTMP), 
                        nu = 0)
          Uy <- t(tmp$vt)
          Dy <- tmp$d
          D.inv <- 1/Dx.prior
          D.inv[abs(D.inv) == Inf] <- 0
          tmp <- La.svd(rbind(sqrtVinvTMP %*% mod$FF[good, 
                                                     , drop = FALSE] %*% Ux.prior, diag(x = D.inv, 
                                                                                        nrow = length(D.inv))), nu = 0)
          Ux <- Ux.prior %*% t(tmp$vt)
          Dx <- 1/tmp$d
          Dx[abs(Dx) == Inf] <- 0
          e <- as.matrix(y[i, good] - f)
          mod$m0 <- a + crossprod(Dx * t(Ux)) %*% tF.VinvTMP %*% 
            e
          ll <- ll + 2 * sum(log(Dy)) + crossprod(crossprod(Uy, 
                                                            e)/Dy)
          lls[i] <- 2 * sum(log(Dy)) + crossprod(crossprod(Uy, 
                                                           e)/Dy)
        }
      }
    }
  # Return likelihood instead of negative loglikelihood return by dlmLL
  return(lls = exp(-lls*0.5)) 
}