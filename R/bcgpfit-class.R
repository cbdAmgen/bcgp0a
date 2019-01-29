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
