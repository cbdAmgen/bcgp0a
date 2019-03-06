##' Class \code{bcgpfit}: fitted bcgp model
##'
##' @slot sim A list containing simulation results including the
##' posterior draws as well as various pieces of metadata used by many of the
##' methods for \code{bcgpfit} objects.
##' @slot priors A list containing the values of the prior parameters.
##' @slot inits The initial values (either user-specified or generated
##' randomly) for all chains. This is a list with one component per chain. Each
##' component is a named list containing the initial values for each parameter
##' for the corresponding chain.
##' @slot algorithm Either "Stan" or "M-H and Gibbs". Tells the method of sampling from the
##' posterior
##' @section TODO: Decide whether to put in \code{prototype} and \code{validity}
##' arguments.
setClass("bcgpfit",
         slots = c(model_pars = "character",
                   par_dims = "list",
                   sim = "list",
                   priors = "list",
                   inits = "list",
                   args = "list",
                   algorithm = "character",
                   scale = "matrix"))


# setMethod("show", "bcgpfit",
#           function(object) {
#             print.bcgpfit(x = object, pars = object@sim$pars_oi)
#           })
#
# print.bcgpfit <- function(x, pars = x@sim$pars_oi,
#                           probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
#                           digits_summary = 2, include = TRUE, ...) {
#
#   if(!include) pars <- setdiff(x@sim$pars_oi, pars)
#   s <- summary(x, pars, probs, ...)
#   if (is.null(s)) return(invisible(NULL))
#   n_kept <- x@sim$n_save - x@sim$warmup2
#   cat("Inference for Stan model: ", x@model_name, '.\n', sep = '')
#   cat(x@sim$chains, " chains, each with iter=", x@sim$iter,
#       "; warmup=", x@sim$warmup, "; thin=", x@sim$thin, "; \n",
#       "post-warmup draws per chain=", n_kept[1], ", ",
#       "total post-warmup draws=", sum(n_kept), ".\n\n", sep = '')
#
#   if (!is.null(x@stan_args[[1]]$method) &&
#       x@stan_args[[1]]$method == "variational") {
#     print(round(s$summary, digits_summary), ...)
#     cat("\nApproximate samples were drawn using VB(", x@stan_args[[1]]$algorithm, ") at ", x@date,
#         ".\n", sep = '')
#     message("We recommend genuine 'sampling' from the posterior distribution for final inferences!")
#     return(invisible(NULL))
#   }
#
#   # round n_eff to integers
#   s$summary[, 'n_eff'] <- round(s$summary[, 'n_eff'], 0)
#
#   print(round(s$summary, digits_summary), ...)
#
#   sampler <- attr(x@sim$samples[[1]], "args")$sampler_t
#
#   cat("\nSamples were drawn using ", sampler, " at ", x@date, ".\n",
#       "For each parameter, n_eff is a crude measure of effective sample size,\n",
#       "and Rhat is the potential scale reduction factor on split chains (at \n",
#       "convergence, Rhat=1).\n", sep = '')
#   return(invisible(NULL))
# }
