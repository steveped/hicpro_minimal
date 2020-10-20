local({
  conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
  if (!is.na(conda) && startsWith(R.home(), conda)) {
    .libPaths(R.home("library"))
  }
})

