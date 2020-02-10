#' @title Estimate average effects after CEA forests.
#' @description \code{avg_effects} Estimate de-biased average causal effects for the entire sample or for a specified subpopulation.
#'
#' @param forest A trained CEA forest.
#' @param WTP Willingness to pay per one-unit increase in Y.
#' @param subset A specified subpopulation to obtain average effects for (optional).
#' @param robust.se Whether or not robust (sandwich) standard errors are desired. Defaults to FALSE. Bootstrapped CIs or ICER CIs are not affected by this setting.
#' @param ci.level The desired confidence level.
#' @param compliance.scores An optional two-column matrix containing pre-fitted compliance scores for instrumental forests (col 1 = outcomes, col 2 = costs).
#' @param icer.ci Logical, if confidence intervals for ICERs are desired. Uses Fieller's method if boot.ci=FALSE. Defaults to TRUE.
#' @param boot.ci Logical, if bootstrapped confidence intervals (BCa) are desired. Defaults to FALSE.
#' @param R Number of bootstrap replicates for bootstrapped CIs. Defaults to 999.
#'
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., & Newey, W. (2017). Double/debiased/neyman machine learning of treatment effects. American Economic Review, 107(5), 261-65.
#'
#' @return Returns de-biased, pooled effect estimates with asymptotic variance estimates or accelerated percentile bootstrap confidence intervals if boot.ci=TRUE. Confidence intervals for ICERs are estimated using Fieller's method unless boot.ci=TRUE.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export

avg_effects = function(forest, WTP=NULL, subset=NULL, robust.se=FALSE, ci.level=0.95, compliance.scores=NULL, icer.ci=TRUE, boot.ci=FALSE, R=999) {

  if (!("CEAforests" %in% class(forest))) {
    stop("This function only works with CEAforests objects.")
  }

  if (is.null(subset)) {
    subset <- 1:length(forest[["outcome.forest"]]$Y.hat)
  }
  if (class(subset) == "logical" & length(subset) ==
      length(forest[[1]]$Y.hat)) {
    subset <- which(subset)
  }
  if (!all(subset %in% 1:length(forest[["outcome.forest"]]$Y.hat))) {
    stop(paste("If specified, subset must be a vector contained in 1:n,",
               "or a boolean vector of length n."))
  }
  if (is.null(WTP)) {
    #If another WTP is not supplied, take the WTP from the forest object.
    WTP = forest[["WTP"]]
  }

  alpha=1-ci.level

  observation.weight <- 1:length(forest[[1]]$Y.hat)
  subset.W.orig <- forest[["outcome.forest"]]$W.orig[subset]
  subset.W.hat <- forest[["outcome.forest"]]$W.hat[subset]

  subset.Y.orig <- forest[["outcome.forest"]]$Y.orig[subset]
  subset.Y.hat <- forest[["outcome.forest"]]$Y.hat[subset]
  tau.hat.pointwise <- predict(forest[[1]])$predictions[subset]
  subset.weights <- observation.weight[subset]

  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)

  Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise

  if (("causal_forest" %in% class(forest[["outcome.forest"]]))) {
  raw.y = tau.hat.pointwise
  cate.y = raw.y+(((subset.Y.orig-Y.hat.1)*subset.W.orig)/subset.W.hat)-((subset.Y.orig-Y.hat.0)*(1-subset.W.orig))/(1-subset.W.hat)
  #gamma.y1 = Y.hat.1+((subset.Y.orig-Y.hat.1)*subset.W.orig)/subset.W.hat
  #gamma.y0 = Y.hat.0+((subset.Y.orig-Y.hat.0)*(1-subset.W.orig))/(1-subset.W.hat)
  cate.nmb.temp = cate.y*WTP
  } else if (("instrumental_forest" %in% class(forest[["outcome.forest"]]))) {
    if (is.null(compliance.scores)){
    compliance.forest <- grf::causal_forest(forest[["outcome.forest"]]$X.orig, Y = forest[["outcome.forest"]]$W.orig,
                                       W = forest[["outcome.forest"]]$Z.orig, Y.hat = forest[["outcome.forest"]]$W.hat, W.hat = forest[["outcome.forest"]]$Z.hat)
    compliance.score <- predict(compliance.forest)$predictions
    } else {compliance.score <- compliance.scores[,1]}
    subset.compliance.score <- compliance.score[subset]
    subset.Z.orig <- forest[["outcome.forest"]]$Z.orig[subset]
    subset.Z.hat <- forest[["outcome.forest"]]$Z.hat[subset]
    subset.g.hat = (subset.Z.orig - subset.Z.hat)/(subset.Z.hat *
                                                     (1 - subset.Z.hat))/subset.compliance.score
    cate.y = tau.hat.pointwise + subset.g.hat * (subset.Y.orig -
                                                             subset.Y.hat - (subset.W.orig - subset.W.hat) * tau.hat.pointwise)
    cate.nmb.temp = cate.y*WTP

  }

  if (is.null(subset)) {
    subset <- 1:length(forest[["cost.forest"]]$Y.hat)
  }
  if (class(subset) == "logical" & length(subset) ==
      length(forest[["cost.forest"]]$Y.hat)) {
    subset <- which(subset)
  }
  if (!all(subset %in% 1:length(forest[["cost.forest"]]$Y.hat))) {
    stop(paste("If specified, subset must be a vector contained in 1:n,",
               "or a boolean vector of length n."))
  }
  observation.weight <- 1:length(forest[["cost.forest"]]$Y.hat)
  subset.W.orig <- forest[["cost.forest"]]$W.orig[subset]
  subset.W.hat <- forest[["cost.forest"]]$W.hat[subset]
  subset.Y.orig <- forest[["cost.forest"]]$Y.orig[subset]
  subset.Y.hat <- forest[["cost.forest"]]$Y.hat[subset]
  tau.hat.pointwise <- predict(forest[["cost.forest"]])$predictions[subset]
  subset.weights <- observation.weight[subset]

  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)

  Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise

  if (("causal_forest" %in% class(forest[["cost.forest"]]))) {

    raw.cost = tau.hat.pointwise
    cate.cost = raw.cost+(((subset.Y.orig-Y.hat.1)*subset.W.orig)/subset.W.hat)-((subset.Y.orig-Y.hat.0)*(1-subset.W.orig))/(1-subset.W.hat)
    #gamma.c1 = Y.hat.1+((subset.Y.orig-Y.hat.1)*subset.W.orig)/subset.W.hat
    #gamma.c0 = Y.hat.0+((subset.Y.orig-Y.hat.0)*(1-subset.W.orig))/(1-subset.W.hat)


  } else if (("instrumental_forest" %in% class(forest[["cost.forest"]]))) {
    if (is.null(compliance.scores)){
    compliance.forest <- grf::causal_forest(forest[["cost.forest"]]$X.orig, Y = forest[["cost.forest"]]$W.orig,
                                       W = forest[["cost.forest"]]$Z.orig, Y.hat = forest[["cost.forest"]]$W.hat, W.hat = forest[["cost.forest"]]$Z.hat)
    compliance.score <- predict(compliance.forest)$predictions
    } else {compliance.score <- compliance.scores[,2]}
    subset.compliance.score <- compliance.score[subset]
    subset.Z.orig <- forest[["cost.forest"]]$Z.orig[subset]
    subset.Z.hat <- forest[["cost.forest"]]$Z.hat[subset]
    subset.g.hat = (subset.Z.orig - subset.Z.hat)/(subset.Z.hat *
                                                     (1 - subset.Z.hat))/subset.compliance.score
    cate.cost = tau.hat.pointwise + subset.g.hat * (subset.Y.orig -
                                                   subset.Y.hat - (subset.W.orig - subset.W.hat) * tau.hat.pointwise)
  }

  #Combine
  cate.nmb = cate.nmb.temp-cate.cost

  #Estimate means and standard errors

  form.y = as.formula(paste("cate.y~1", sep=""))
  form.c = as.formula(paste("cate.cost~1", sep=""))
  form.nb = as.formula(paste("cate.nmb~1", sep=""))
  mo.y = lm(form.y)
  mo.c = lm(form.c)
  mo.nb = lm(form.nb)

  #Get coefficients
  average.delta_y = as.vector(coef(mo.y))
  average.delta_cost = as.vector(coef(mo.c))
  average.nmb = as.vector(coef(mo.nb))

  if (!isTRUE(boot.ci)) {

  #Get standard errors

  if (isTRUE(robust.se)) {
  se.delta_y <- as.vector(sqrt(diag(sandwich::vcovHC(mo.y))))
  se.delta_cost = as.vector(sqrt(diag(sandwich::vcovHC(mo.c))))
  se.nmb = as.vector(sqrt(diag(sandwich::vcovHC(mo.nb))))
  } else {
    se.delta_y <- as.vector(sqrt(diag(vcov(mo.y))))
    se.delta_cost = as.vector(sqrt(diag(vcov(mo.c))))
    se.nmb = as.vector(sqrt(diag(vcov(mo.nb))))
  }

  #Confidence intervals
  lower.y = average.delta_y-se.delta_y*qt(1-(alpha/2), df=length(cate.nmb)-1)
  upper.y = average.delta_y+se.delta_y*qt(1-(alpha/2), df=length(cate.nmb)-1)
  lower.c = average.delta_cost-se.delta_cost*qt(1-(alpha/2), df=length(cate.nmb)-1)
  upper.c = average.delta_cost+se.delta_cost*qt(1-(alpha/2), df=length(cate.nmb)-1)
  lower.n = average.nmb-se.nmb*qt(1-(alpha/2), df=length(cate.nmb)-1)
  upper.n = average.nmb+se.nmb*qt(1-(alpha/2), df=length(cate.nmb)-1)

  # Estimate ICER

  if (isTRUE(icer.ci)) {
    icers = icer.ci(cate.cost, cate.y, alpha=1-ci.level)
  } else { icers = c(NA,NA,NA)}

  } else if (isTRUE(boot.ci)) {

   boot.inference = boot.dr_scores(dy=cate.y, dc=cate.cost,
                                   WTP=WTP, alpha=1-ci.level, R=R)[[1]]

   se.delta_y = boot.inference[1,1]
   se.delta_cost = boot.inference[2,1]
   se.nmb = boot.inference[3,1]
   lower.y = boot.inference[1,2]
   lower.c = boot.inference[2,2]
   lower.n = boot.inference[3,2]
   upper.y = boot.inference[1,3]
   upper.c = boot.inference[2,3]
   upper.n = boot.inference[3,3]

   icers = c(mean(cate.cost)/mean(cate.y), boot.inference[4,2], boot.inference[4,3])

  }

  accept.normal = 1-stats::pnorm(0, average.nmb, se.nmb)

  #Store results
  res.y = c(estimate=average.delta_y, std.err = se.delta_y, lower=lower.y, upper=upper.y, accept.prob.normal=NA)
  res.c = c(estimate=average.delta_cost, std.err = se.delta_cost, lower=lower.c, upper=upper.c, accept.prob.normal=NA)
  res.nb = c(estimate=average.nmb, std.err = se.nmb, lower=lower.n, upper=upper.n, accept.prob.normal=accept.normal)
  res.icer = c(estimate=icers[[1]], std.err = NA, lower=icers[[2]], upper=icers[[3]], accept.prob.normal=NA)

  res = rbind(res.y, res.c, res.nb, res.icer)
  rownames(res) = c("Average effect on outcome", "Average effect on costs", paste("Average net monetary benefit (WTP: ", WTP,")", sep=""), "Incremental cost-effectiveness ratio")
  return(round(res, 4))
}

#' @title ICERs with Fieller's method CIs
#' @description \code{fieller} Fieller's method confidence intervals for ICERs.
#'
#' @param dy delta_y
#' @param dc delta_cost
#' @param alpha Desired confidence level.
#' @keywords internal
#' @return Returns ICER with confidence intervals based on Fieller's method assuming a bivariate normal distribution for the averages of delta costs and delta y.
#' @export

icer.ci <- function(dc,dy,alpha){
    a = mean(dc); b = mean(dy)
    theta <- a/b
    v11 <- var(dc)/length(dc)
    v12 <- cov(dc,dy)/length(dc)
    v22 <- var(dy)/length(dy)
    z <- qt(1-alpha/2, df=length(dc)-1)
    g <- z*v22/b^2
    c = v11 - 2*theta*v12 + theta^2 * v22 - g*(v11-v12^2/v22)
    if (!isTRUE(c<=0)) {
    C <- sqrt(v11 - 2*theta*v12 + theta^2 * v22 - g*(v11-v12^2/v22))
    lower <- (1/(1-g))*(theta- g*v12/v22 - z/b * C)
    upper <- (1/(1-g))*(theta- g*v12/v22 + z/b * C)
    } else {
      message("Fieller's method CIs failed to compute (NaN). You may want to try setting boot.ci=TRUE to obtain bootstrapped (BCa) confidence intervals.")
      lower <- NA
      upper <- NA
    }
    if (isTRUE(upper<lower)){
      message("Fieller's method CIs appear to be unbounded and may seem strange, but this behavior is to be expected if the incremental outcome does not differ significantly from zero. You may want to try setting boot.ci=TRUE to obtain bootstrapped (BCa) confidence intervals, but some authors recommend against this (see documentaion for details).")
    }
  return(c(ICER=theta,lower,upper))
}

#' @title Bootstrap average effects
#' @description \code{boot.dr_scores} Bootstraps doubly robust scores and obtains accelerated bootstrap confidence intervals (BCa).
#'
#' @param dy delta_y
#' @param dc delta_cost
#' @param alpha Desired confidence level.
#' @keywords internal
#' @return Returns a matrix with estimated standard errors and BCa confidence intervals.
#' @export
#'
boot.dr_scores <- function(dc, dy, R, WTP, alpha) {
  x = dc
  y = dy
  df = as.data.frame(cbind(x,y))
  n = nrow(df)
  bfun = function(data, indices, x,y, WTP) {
    d=data[indices,]
    DY = mean(d[,y])
    DCOST = mean(d[,x])
    NMB=DY*WTP-DCOST
    ICER=DCOST/DY
    res=c(DY,DCOST,NMB,ICER)
  return(res)}
  b=boot::boot(data=df, bfun, R=R, x="x", y="y", WTP=WTP)
  dy_se=sd(b$t[,1])
  dc_se=sd(b$t[,2])
  nmb_se=sd(b$t[,3])
  icer_se=NA
  res = list()
  if (R<=n) {
    warning("Number of bootstrap replicates R is smaller than the number of rows in the data. BCa confidence intervals cannot not be computed. This warning can be ignored for CE planes. For other applications, please increase R.")
    res[[1]] = cbind(c(NA,NA,NA,NA), c(NA,NA,NA,NA), c(NA,NA,NA,NA))}
  if (R>n) {
  bci_y=boot::boot.ci(b, index=1, conf=1-alpha, type="bca")
  bci_c=boot::boot.ci(b, index=2, conf=1-alpha, type="bca")
  bci_nmb=boot::boot.ci(b, index=3, conf=1-alpha, type="bca")
  bci_icer=boot::boot.ci(b, index=4, conf=1-alpha, type="bca")
  lower_y = bci_y$bca[,4]
  upper_y = bci_y$bca[,5]
  lower_c = bci_c$bca[,4]
  upper_c = bci_c$bca[,5]
  lower_nmb = bci_nmb$bca[,4]
  upper_nmb = bci_nmb$bca[,5]
  lower_icer = bci_icer$bca[,4]
  upper_icer = bci_icer$bca[,5]
  ses = c(dy_se,dc_se,nmb_se,icer_se)
  lowers = c(lower_y,lower_c,lower_nmb,lower_icer)
  uppers = c(upper_y, upper_c, upper_nmb, upper_icer)
  res[[1]] = cbind(ses, lowers, uppers)
  }
  res[[2]] = b$t
  return(res)
}


