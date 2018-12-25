.onAttach <- function(...) {
  bcgpLib <- dirname(system.file(package = "bcgp0a"))
  pkgdesc <- packageDescription("bcgp0a", lib.loc = bcgpLib)
  packageStartupMessage(paste("bcgp (Version ", pkgdesc$Version,")", sep = ""))
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM call\n",
                        "options(mc.cores = parallel::detectCores()).")
}
