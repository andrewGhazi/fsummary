.onAttach = function(libname, pkgname) {
  ver = utils::packageVersion("fsummary")

  cli::cli_inform("{cli::symbol$pointer} Use {.code mirai::daemons(4)} to start background daemons for parallelization.",
                  class = "packageStartupMessage")
}
