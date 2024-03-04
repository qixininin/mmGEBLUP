#' qtxnetwork.perform function
#' 1D scanning for additive and additive-by-environment qtls
#' @param qtxnetwork_path path to the executable QTXNetwork file
#' @param gen_file (input file) path and name for .gen file
#' @param phe_file (input file) path and name for .phe file
#' @param pre_file (output file) path and name for .pre file
#'
#' @return NULL
#' @export
#'
#' @examples \dontrun{qtxnetwork.perform(
#'                    qtxnetwork_path = "/public3/zqx/QTXNetwork_4.0/build/QTXNetwork"
#'                    gen_file = "./inst/Simulation.gen"
#'                    phe_file = "./inst/Simulation_SimTrait.phe"
#'                    pre_file = "./inst/Simulation_SimTrait.pre")}

qtxnetwork.perform <- function(qtxnetwork_path, gen_file, phe_file, pre_file)
{

  command = paste0(qtxnetwork_path,
                   " --map ", gen_file,
                   " --txt ", phe_file,
                   " --out ", pre_file,
                   " --QTX 1 --only-1D 1")
  cat(paste0("*** QTXNetwork command:\n", command))

  system(command)

}
