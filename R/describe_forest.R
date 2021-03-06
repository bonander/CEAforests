#' @title Create a descriptive table of the results.
#' @description \code{describe_forest} Generates a table of descriptives after stratifying the data into cost-effectiveness certainty groups after CEA forests.
#'
#' @param forest A trained CEA forest.
#' @param X A matrix of variables to include in the table. If NULL, the X matrix from the forest object is used.
#' @param WTP The willingness to pay threshold. Defaults to WTP supplied to CEA forest object. Ignored when certainty_groups is TRUE.
#' @param alpha The desired significance level. Defaults to 0.05. Only used when certainty_groups is TRUE.
#' @param certainty_groups Divide the sample into groups based on sampling uncertainty? Defaults to FALSE.
#' @param ... Other options to be passed to tableone::CreateTableOne (e.g., which variables are to be treated as factors, see examples).
#'
#' @return Returns a table object from the tableone package. See the print function from the tableone package for additional options.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export
describe_forest = function(forest, X=NULL, WTP=NULL, alpha=0.05, certainty_groups=FALSE, ...) {

  if (isTRUE(any(class(forest) %in% c("CEAforests")))==FALSE) {stop("Unrecognized class. Please supply a CEAforests object.")}


  if (isTRUE(certainty_groups)==TRUE){

  preds = predict(forest)
  dy = preds$predicted.delta_y
  lowerx = dy-sqrt(preds$variance.delta_y)*qnorm(1-alpha/2)
  upperx = dy+sqrt(preds$variance.delta_y)*qnorm(1-alpha/2)
  dc = preds$predicted.delta_cost
  lowery = dc-sqrt(preds$variance.delta_cost)*qnorm(1-alpha/2)
  uppery = dc+sqrt(preds$variance.delta_cost)*qnorm(1-alpha/2)

  if (is.null(X)) {
    X=forest[["outcome.forest"]]$X.orig
  }

    nmb.point = preds$predicted.nmb

    nmb.lower = nmb.point-sqrt(preds$variance.nmb)*qnorm(1-alpha/2)
    nmb.upper = nmb.point+sqrt(preds$variance.nmb)*qnorm(1-alpha/2)

    certain=as.numeric(rep(0, length(nmb.point)))
    certain[nmb.point>0 & nmb.lower>0 & nmb.upper>0] = 1
    certain[nmb.point<0 & nmb.lower<0 & nmb.upper<0] = -1

    if (isTRUE(length(unique(certain))==3)){
      certain = factor(certain, labels=c("Not cost-effective",
                                         "Undetermined",
                                         "Cost-effective"))

    } else if (isTRUE(identical(base::sort(unique(certain)), c(-1,0)))) {
      certain = factor(certain, labels=c("Not cost-effective", "Undetermined"))
    } else if (isTRUE(identical(base::sort(unique(certain)), c(0,1)))) {
      certain = factor(certain, labels=c("Undetermined", "Cost-effective"))
    } else if (isTRUE(identical(unique(certain), c(0)))) {
      certain = factor(certain, labels="Undetermined")
    } else if (isTRUE(identical(unique(certain), c(-1)))) {
      certain = factor(certain, labels="Not cost effective")
    } else if (isTRUE(identical(unique(certain), c(1)))) {
      certain = factor(certain, labels="Cost-effective")
    }

  tdat = as.data.frame(X)
  tdat$Certainty_Group = certain
  if (isTRUE(length(unique(certain))>1)) {
  tab = tableone::CreateTableOne(vars=colnames(tdat)[!colnames(tdat) %in% c("Certainty_Group")], data=tdat, strata="Certainty_Group", ...)
  } else {tab = tableone::CreateTableOne(vars=colnames(tdat)[!colnames(tdat) %in% c("Certainty_Group")], data=tdat, ...)
  warning(paste("All observations belong to a single certainty group: ", levels(certain)),". Table reflects the entire sample.", sep="")}
  }

  if (isTRUE(certainty_groups)==FALSE) {
    if (is.null(WTP)) {WTP = forest[["WTP"]]} else {WTP=WTP}
    preds = predict(forest)
    dy = preds$predicted.delta_y
    dc = preds$predicted.delta_cost
    nmb.point = preds$predicted.nmb
    below_threshold = as.numeric((dy*WTP-dc)>0)

    if (is.null(X)) {
      X=forest[["outcome.forest"]]$X.orig
    }

    tdat = as.data.frame(X)
    Cost_Effective = ifelse(below_threshold==1, "Yes","No")
    tdat$Cost_Effective = Cost_Effective

    if (isTRUE(length(unique(Cost_Effective))>1)) {

      tab = tableone::CreateTableOne(vars=colnames(tdat)[!colnames(tdat) %in% c("Cost_Effective")], data=tdat, strata="Cost_Effective", ...)

    } else {

      tab = tableone::CreateTableOne(vars=colnames(tdat)[!colnames(tdat) %in% c("Cost_Effective")], data=tdat, ...)
      warning(paste("All observations are either cost-effective or not cost-effective. Table reflects the entire sample.", sep=""))

    }

  }

  return(tab)
}
