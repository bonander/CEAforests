#' @title Internal function.
#' @description \code{dr.cate.regress} Fits a linear regression to the estimated heterogeneous effects to assess the determinants of heterogeneity. Uses augmented inverse probability weighting (AIPW) to de-bias the scores, with sandwich standard errors.
#'
#' @param dr.scores Doubly robust scores for tau.
#' @param X The variables to include in the model. If NULL, and intercept only model is run.
#' @param alpha The desired significance level, defaults to 0.05.
#' @param robust.se Whether or not to compute robust (sandwich) standard errors. Defaults to TRUE.
#' @param subset A specified subset.
#'
#' @return Returns the results from linear regression on the estimated heterogeneous effects with supplied covariates with asymptotically valid standard errors. To be used with causal_forest object directly. For cea_forests or nmb_forests, use the explain_forest() function.
#' @keywords internal
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export

dr.cate.regress = function(dr.scores, X=NULL, alpha=0.05, robust.se=TRUE, subset=NULL) {# Fix output with only one variable!

  #Get regression coefficients for OLS on dr.scores

  if (is.null(X)) { #Intercept only if no X is provided
    form = as.formula(paste("dr.scores~1", sep=""))
    mo = lm(form)
  } else {X = data.frame(X)
    if (is.null(subset)) {
      subset <- 1:nrow(X)
    }
    if (class(subset) == "logical" & length(subset) ==
        nrow(X)) {
      subset <- which(subset)
    }
    if (!all(subset %in% 1:nrow(X))) {
      stop(paste("If specified, subset must be a vector contained in 1:n,",
                 "or a boolean vector of length n."))
    }
    M=data.frame(X[subset,])
    colnames(M)=colnames(X)
    regdf = data.frame(cbind(dr.scores, M))
    b=paste(colnames(M), collapse="+")
    form=as.formula(paste("dr.scores~", b, sep=""))
    mo=lm(form, data=regdf)
  }

  #Get coefficients
  coefs = coef(mo)

  #Get standard errors
  if (isTRUE(robust.se)) {
  se <- sqrt(diag(sandwich::vcovHC(mo)))} else {se = sqrt(diag(vcov(mo)))}

  #Confidence intervals
  lowers = coefs-se*qt(1-(alpha/2), df=mo$df.residual)
  uppers = coefs+se*qt(1-(alpha/2), df=mo$df.residual)

  #Get p-values
  p_value <- qnorm(1-(0.05/2))*pt(abs(coefs/se), df=mo$df.residual,
                                  lower.tail= FALSE)

  #Return
  return(cbind(Estimates = coefs, StdError = se, LowerCI = lowers, UpperCI = uppers, P.val = round(p_value, 4)))
}

#' @title Internal function.
#' @description \code{cate.prepare} Provides de-biased tau estimates based on AIPW (for use in custom explainers, e.g., regression models or GAMs).
#'
#' @param forest A trained causal forest or instrumental forest object.
#'
#' @return Returns a vector of transformed outcomes that can be used in parametric or semi-parametric regressions to explain the results from the forest.
#' @keywords internal
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export

cate.prepare = function(forest, subset=NULL, compliance.score=NULL) {

  if (is.null(subset)) {
    subset <- 1:length(forest$Y.hat)
  }
  if (class(subset) == "logical" & length(subset) ==
      length(forest$Y.hat)) {
    subset <- which(subset)
  }
  if (!all(subset %in% 1:length(forest$Y.hat))) {
    stop(paste("If specified, subset must be a vector contained in 1:n,",
               "or a boolean vector of length n."))
  }

  observation.weight <- 1:length(forest$Y.hat)
  subset.W.orig <- forest$W.orig[subset]
  subset.W.hat <- forest$W.hat[subset]

  subset.Y.orig <- forest$Y.orig[subset]
  subset.Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]
  subset.weights <- observation.weight[subset]

  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)

  Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise

  if (("causal_forest" %in% class(forest))) {

    dr.scores = tau.hat.pointwise+(((subset.Y.orig-Y.hat.1)*subset.W.orig)/subset.W.hat)-((subset.Y.orig-Y.hat.0)*(1-subset.W.orig))/(1-subset.W.hat)


  } else if (("instrumental_forest" %in% class(forest))) {
    if (is.null(compliance.score)) {
    compliance.forest <- grf::causal_forest(forest$X.orig, Y = forest$W.orig,
                                       W = forest$Z.orig, Y.hat = forest$W.hat, W.hat = forest$Z.hat)
    compliance.score <- predict(compliance.forest)$predictions
    }
    subset.compliance.score <- compliance.score[subset]
    subset.Z.orig <- forest$Z.orig[subset]
    subset.Z.hat <- forest$Z.hat[subset]
    subset.g.hat = (subset.Z.orig - subset.Z.hat)/(subset.Z.hat *
                                                     (1 - subset.Z.hat))/subset.compliance.score
    dr.scores = tau.hat.pointwise + subset.g.hat * (subset.Y.orig -
                                                      subset.Y.hat - (subset.W.orig - subset.W.hat) * tau.hat.pointwise)
  } else {
    stop("Unrecognized class. Please supply a causal_forest or instrumental_forest object.")
  }
  return(dr.scores)
}


#' @title Explain the results from CEA forests using best linear projection.
#' @description \code{explain_forest} Fits a linear regression to the estimated heterogeneous effects to assess the determinants of heterogeneity. Uses augmented inverse probability weighting (AIPW) to de-bias the scores.
#'
#' @param forest A trained CEA forest.
#' @param X The variables to include in the model. If NULL, and intercept only model is run, which gives the ATE.
#' @param alpha The desired significance level, defaults to 0.05.
#' @param WTP The WTP for the net monetary benefit forest. Uses the WTP supplied to the CEA forest if NULL.
#' @param robust.se If robust (sandwich) standard errors are desired. Defaults to TRUE.
#' @param subset A specified subset of the data to compute the regression on.
#' @references Chernozhukov, Victor, and Vira Semenova. "Simultaneous inference for Best Linear Predictor of the Conditional Average Treatment Effect and other structural functions." arXiv preprint arXiv:1702.06240 (2017).
#' @return Returns the results from linear regression(s) on the estimated heterogeneous effects with supplied covariates with asymptotically valid standard errors.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export

explain_forest = function(forest, X=NULL, alpha=0.05, WTP=NULL, robust.se=TRUE, subset=NULL) {
  if (is.null(WTP)) {WTP=forest[["WTP"]]}
  if (isTRUE(any(class(forest) %in% c("CEAforests")))) {
   res = list()
   dr.y = cate.prepare(forest[["outcome.forest"]], subset=subset)
   dr.cost = cate.prepare(forest[["cost.forest"]], subset=subset)
   dr.nmb = dr.y*WTP-dr.cost
   res[["Outcome forest"]] = dr.cate.regress(dr.y, X, alpha, robust.se=robust.se, subset=subset)
   res[["Cost forest"]] = dr.cate.regress(dr.cost, X, alpha, robust.se=robust.se, subset=subset)
   res[[paste("Net monetary benefit forest, WTP = ", WTP, sep="")]] = dr.cate.regress(dr.nmb, X, alpha, robust.se=robust.se, subset=subset)

  } else {stop("Unrecognized or unsupported forest object. Please supply a CEAforests object.")}
  return(res)
}

#' @title Get de-biased heterogeneous effect estimates (for custom explainers)
#' @description \code{debias_effects} Provides de-biased tau estimates based on AIPW (for use in custom explainers, e.g., regression models or GAMs).
#'
#' @param forest A trained CEAforests object.
#' @param WTP The willingness to pay for a one-unit increase in Y. Defaults to the WTP supplied to the CEAforests object unless manually specified.
#'
#' @return Returns a matrix of transformed outcomes that can be used in parametric or semi-parametric regressions to explain the results from the forest.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export

debias_effects = function(forest, WTP=NULL) {
  if (is.null(WTP)) {WTP = forest[["WTP"]]}
  if (isTRUE(any(class(forest) %in% c("CEAforests")))) {
    debiased.delta_outcome = cate.prepare(forest[["outcome.forest"]])
    debiased.delta_cost = cate.prepare(forest[["cost.forest"]])
    debiased.nmb = debiased.delta_outcome*WTP-debiased.delta_cost
    res = cbind(debiased.delta_outcome, debiased.delta_cost, debiased.nmb)
    return(res)
  } else {stop("Unrecognized or unsupported forest object. Please supply CEAforests object.")}
}



