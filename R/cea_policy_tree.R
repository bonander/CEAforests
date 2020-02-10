#' @title Train a policy tree after a CEA forest.
#' @description \code{cea_policy_tree} Trains an efficient policy decision tree given a CEA forest (a wrapper for policytree::policy_tree).
#'
#' @param forest A trained CEA forest.
#' @param X A covariate matrix containing variables that are to be used in the policy tree.
#' @param WTP Willingness to pay for a one unit increase in the outcome. If NULL, the WTP supplied to the CEA forest is used.
#' @param depth The desired depth for the decision tree.
#' @param ci.level Desired significance level (for confidence intervals).
#' @param robust.se Whether or not robust (sandwich) standard errors are desired. Defaults to FALSE.
#'
#'
#' @references Athey, S., & Wager, S. (2017). Efficient policy learning. arXiv preprint arXiv:1702.02896.
#'
#' @return Returns a trained policy tree.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @importFrom utils installed.packages
#' @import stats
#' @export

cea_policy_tree = function(forest, X, WTP=NULL, depth=2, ci.level=0.95, robust.se=FALSE) {
  if (isTRUE("policytree" %in% rownames(installed.packages())==FALSE)) {
    stop("Package \"policytree\" must be installed to estimate policy trees.")
  }
  if (isTRUE(any(class(forest) %in% c("CEAforests")))) {

    if (is.null(WTP)==TRUE) {
      WTP = forest[["WTP"]]
      warning(paste("WTP not specified, assuming WTP = ", WTP, ", as supplied to the CEAforests object.", sep=""))
    }

    gamma1 = cate.prepare(forest[["outcome.forest"]])*WTP
    gamma2 = cate.prepare(forest[["cost.forest"]])
    gamma = gamma1-gamma2
    Gamma = cbind(control=-gamma, treated=gamma)

    treefit = policytree::policy_tree(X, Gamma, depth=depth) }  else {
      stop("Unrecognized or unsupported forest object. Please supply a CEAforests object.")}

  results = list()
  treefit$n.sample = nrow(Gamma)
  results[["tree"]] = treefit
  results[["X"]] = X
  class(results) = c("cea_policy_tree", "CEAforests")

  return(results)
}

#' @title Conduct inference for a personalized treatment policy.
#' @description Conduct inference for a personalized treatment policy, either using a manually specified policy or a learned policy.
#'
#' @param forest A trained CEA forest.
#' @param treat.policy A logical vector or cea policy tree defining the subset covered by the policy.
#' @param WTP Willingness to pay for a one unit increase in the outcome. If NULL, the WTP supplied to the CEA forest is used.
#' @param ci.level Desired significance level (for confidence intervals).
#' @param robust.se Whether or not robust (sandwich) standard errors are desired. Defaults to FALSE. Ignored when boot.ci=TRUE.
#' @param boot.ci Whether or not bootstrapped confidence intervals are desired. Defaults to FALSE.
#' @param R The number of bootstrap replications. Defaults to 999. Ignored when boot.ci=FALSE.
#'
#' @return Returns a matrix containing estimates for the average welfare gain per population member under various treatment policies (treat everyone vs. treat no one; treat suggested subset vs. treat no one; treat suggested subset vs. treat everyone). Also outputs the share of the popuation covered by the policy.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @import stats
#' @import boot
#' @export
infer_policy = function(forest, treat.policy, WTP=NULL, ci.level=0.95, robust.se=FALSE, boot.ci=FALSE, R=999) {

  subset = treat.policy

  if (isTRUE(any(class(forest) %in% c("CEAforests")))) {#Check if forest class OK

    if (is.null(WTP)==TRUE) {
      WTP = forest[["WTP"]]
      warning(paste("WTP not specified, assuming WTP = ", WTP, ", as supplied to the CEAforests object.", sep=""))
    }

    #Estimate doubly robust scores

    gamma1 = cate.prepare(forest[["outcome.forest"]])*WTP
    gamma2 = cate.prepare(forest[["cost.forest"]])
    gamma = gamma1-gamma2

  if (isTRUE("cea_policy_tree" %in% class(subset))) {
  #Predict the suggested policy using the supplied X
  X = subset[["X"]]
  predicted.action = predict(subset[["tree"]], newdata=X)
  P = as.numeric(predicted.action==2)
  } else {

    if (class(subset) == "logical" & length(subset) ==
        length(forest[["outcome.forest"]]$Y.hat)) {
      subset <- which(subset)
    }
    if (!all(subset %in% 1:length(forest[["outcome.forest"]]$Y.hat))) {
      stop(paste("treat.policy must be a vector contained in 1:n,",
                 "a boolean vector of length n or a trained CEA policy tree."))
    }
    P = rep(0, length(gamma)); P[subset] = 1
  }

  tau_tr = gamma #Scores for treat everyone
  policy_suggested = tau_tr*P #Scores for suggested policy
  policy_diff = policy_suggested-tau_tr #Scores for difference between suggested and treat everyone

  #Fit intercept only models to get mean and variance
  all.m = lm(tau_tr~1)
  sugg.m = lm(policy_suggested~1)
  vs.m = lm(policy_diff~1)

  est.tr.all = as.vector(coef(all.m))
  est.suggested = as.vector(coef(sugg.m))
  est.diff = as.vector(coef(vs.m))

  ests = c(est.tr.all, est.suggested, est.diff)

  if (!isTRUE(boot.ci)) {#Asymptotic variance estimates

  #Extract estimates and standard errors
  if (isTRUE(robust.se)) {
    se.tr.all = as.vector(sqrt(diag(sandwich::vcovHC(all.m))))
    se.suggested = as.vector(sqrt(diag(sandwich::vcovHC(sugg.m))))
    se.diff = as.vector(sqrt(diag(sandwich::vcovHC(vs.m))))
  } else {
    se.tr.all = as.vector(sqrt(diag(vcov(all.m))))
    se.suggested = as.vector(sqrt(diag(vcov(sugg.m))))
    se.diff = as.vector(sqrt(diag(vcov(vs.m)))) }

  ses = c(se.tr.all, se.suggested, se.diff)

  #Get confidence intervals

  lowers = ests-ses*qt(1-(1-ci.level)/2, df=length(gamma)-1)
  uppers = ests+ses*qt(1-(1-ci.level)/2, df=length(gamma)-1)

  } else {#Else bootstrap scores
    bootres=boot.policy_scores(tau_tr, policy_suggested, R, 1-ci.level)[[1]]
    ses = bootres[,1]
    lowers = bootres[,2]
    uppers = bootres[,3]
  }

  #Share of population who are treated (in suggested policy)
  tr.sugg.share = mean(P)

  #Tidy up and print results
  res = as.data.frame(cbind(ests,ses,lowers,uppers))
  res = rbind(res, c(tr.sugg.share,NA,NA,NA))
  colnames(res) = c("Estimate", "Std.Err", "Lower.CI", "Upper.CI")
  rownames(res) = c("Average welfare gain per population member, treat everyone vs treat no one",
                    "Average welfare gain per population member, suggested policy vs treat no one",
                    "Difference in welfare gain (suggested vs. treat everyone)",
                    "Share of the population treated, suggested policy")

  return(res)

  } else (stop("Unrecognized forest object."))
}



#' Writes each node information
#' If it is a leaf node: show it in different color, show number of samples, show leaf id
#' If it is a non-leaf node: show its splitting variable and splitting value
#' @param tree the tree to convert
#' @param index the index of the current node
#' @keywords internal
cea_create_dot_body <- function(tree, index = 1) {

  #n = tree$n.sample

  node <- tree$nodes[[index]]

  # Leaf case: print label only
  if (node$is_leaf) {
    action <- node$action
    action <- ifelse(action==1, "Do not treat", "Treat")
    line_label <- paste(index - 1, ' [shape=box,style=filled,color="White", height=0.2, label="', action, "\n", '"];', sep="")
    return(line_label)
  }

  # Non-leaf case: print label, child edges
  if (!is.null(node$left_child)) {
    edge <- paste(index - 1, "->", node$left_child - 1)
    if (index == 1) {
      edge_info_left <- paste(edge, '[labeldistance=2.5, labelangle=45, headlabel="Yes"];')
    }
    else {
      edge_info_left <- paste(edge, " ;")
    }
  }
  else {
    edge_info_right <- NULL
  }

  if (!is.null(node$right_child)) {
    edge <- paste(index - 1, "->", node$right_child - 1)
    if (index == 1) {
      edge_info_right <- paste(edge, '[labeldistance=2.5, labelangle=-45, headlabel="No"]')
    } else {
      edge_info_right <- paste(edge, " ;")
    }
  } else {
    edge_info_right <- NULL
  }

  variable_name <- tree$columns[node$split_variable]
  node_info <- paste(index - 1, '[label="', variable_name, "<=", round(node$split_value, 2), '"] ;')

  this_lines <- paste(node_info,
                      edge_info_left,
                      edge_info_right,
                      sep = "\n"
  )

  left_child_lines <- ifelse(!is.null(node$left_child),
                             cea_create_dot_body(tree, index = node$left_child),
                             NULL
  )

  right_child_lines <- ifelse(!is.null(node$right_child),
                              cea_create_dot_body(tree, index = node$right_child),
                              NULL
  )

  lines <- paste(this_lines, left_child_lines, right_child_lines, sep = "\n")

  return(lines)
}

#' Export a tree in DOT format.
#' This function generates a GraphViz representation of the tree,
#' which is then written into `dot_string`.
#' @param tree the tree to convert
#' @keywords internal
cea_export_graphviz <- function(tree) {
  header <- "digraph nodes { \n node [shape=box] ;"
  footer <- "}"
  body <- cea_create_dot_body(tree)

  dot_string <- paste(header, body, footer, sep = "\n")

  return(dot_string)
}

#' Plot a cea_policy_tree tree object.
#' @param x The tree to plot
#' @param ... Additional options (currently ignored).
#'
#' @method plot cea_policy_tree
#' @export
plot.cea_policy_tree <- function(x, ...) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package \"DiagrammeR\" must be installed to plot trees.")
  }

  dot_file <- cea_export_graphviz(x[["tree"]])
  DiagrammeR::grViz(dot_file)
}

#' @title Bootstrap average effects
#' @description \code{boot.dr_scores} Bootstraps doubly robust scores and obtains accelerated bootstrap confidence intervals (BCa).
#'
#' @param Gamma_all Scores for treating everyone vs treating no-one.
#' @param Gamma_policy Scores for suggested policy vs treating no-one.
#' @param R Number of bootstrap replicates.
#' @param alpha Desired confidence level.
#' @keywords internal
#' @return Returns a matrix with estimated standard errors and BCa confidence intervals.
#' @export
#'
boot.policy_scores <- function(Gamma_all, Gamma_policy, R, alpha) {
  df = as.data.frame(cbind(Gamma_all, Gamma_policy))
  n = nrow(df)
  bfun = function(data, indices, Gamma_all, Gamma_policy) {
    d=data[indices,]
    tr_all=mean(d[,Gamma_all])
    tr_policy=mean(d[,Gamma_policy])
    tr_diff=tr_policy-tr_all
    res=c(tr_all,tr_policy,tr_diff)
    return(res)}
  b=boot::boot(data=df, bfun, R=R, Gamma_all="Gamma_all", Gamma_policy="Gamma_policy")
  all_se=sd(b$t[,1])
  policy_se=sd(b$t[,2])
  diff_se=sd(b$t[,3])
  res = list()
  if (R<=n) {
    warning("Number of bootstrap replicates R is smaller than the number of rows in the data. BCa confidence intervals cannot not be computed. Please increase R.")
    res[[1]] = cbind(c(NA,NA,NA), c(NA,NA,NA), c(NA,NA,NA))}
  if (R>n) {
  bci_all=boot::boot.ci(b, index=1, conf=1-alpha, type="bca")
  bci_policy=boot::boot.ci(b, index=2, conf=1-alpha, type="bca")
  bci_diff=boot::boot.ci(b, index=3, conf=1-alpha, type="bca")
  lower_all = bci_all$bca[,4]
  upper_all = bci_all$bca[,5]
  lower_policy = bci_policy$bca[,4]
  upper_policy = bci_policy$bca[,5]
  lower_diff = bci_diff$bca[,4]
  upper_diff = bci_diff$bca[,5]
  ses = c(all_se,policy_se,diff_se)
  lowers = c(lower_all,lower_policy,lower_diff)
  uppers = c(upper_all, upper_policy, upper_diff)
  res = list()
  res[[1]] = cbind(ses, lowers, uppers)
  }
  res[[2]] = b$t
  return(res)
}
