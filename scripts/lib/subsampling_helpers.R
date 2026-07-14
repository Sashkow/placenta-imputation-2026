#!/usr/bin/env Rscript

library(yaml)

parse_config_arg <- function(default_config = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- default_config
  for (arg in args) {
    if (grepl("^--config=", arg)) {
      config_path <- sub("^--config=", "", arg)
    }
  }
  if (is.null(config_path))
    stop("No config specified. Use --config=PATH")
  yaml::read_yaml(config_path)
}

parse_int_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^--", name, "=")
  for (arg in args) {
    if (grepl(pat, arg)) return(as.integer(sub(pat, "", arg)))
  }
  default
}
