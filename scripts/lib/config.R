# Configuration and command-line helpers shared by every script.
# All scripts are run from the repository root.

`%||%` <- function(a, b) if (is.null(a)) b else a

# --name=value and --flag arguments
.cli_args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), .cli_args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
has_flag  <- function(name) any(.cli_args == paste0("--", name))
split_arg <- function(name, default) {
  v <- get_arg(name)
  if (is.null(v)) default else strsplit(v, ",")[[1]]
}

# Overlay `override` on `base`: named lists merge key by key, anything else
# (scalars, vectors, unnamed lists) is replaced; a NULL value removes the key.
merge_config <- function(base, override) {
  for (key in names(override)) {
    both_maps <- is.list(base[[key]]) && !is.null(names(base[[key]])) &&
                 is.list(override[[key]]) && !is.null(names(override[[key]]))
    if (both_maps) base[[key]] <- merge_config(base[[key]], override[[key]])
    else base[key] <- list(override[[key]])
    if (is.null(override[[key]])) base[[key]] <- NULL
  }
  base
}

# Read a YAML config. A file may start from another one with
#   base: other.yaml        (path relative to the file itself)
# and then list only the settings that differ.
load_config <- function(path) {
  if (!file.exists(path)) stop("Config file not found: ", path)
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$base)) return(cfg)
  base <- load_config(file.path(dirname(path), cfg$base))
  cfg$base <- NULL
  merge_config(base, cfg)
}

# Config named by --config=PATH, else the script's default
config_from_cli <- function(default) load_config(get_arg("config", default))
int_arg <- function(name, default) as.integer(get_arg(name, default))
