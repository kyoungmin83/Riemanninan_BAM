local_lib <- normalizePath(file.path(getwd(), ".r-lib"), winslash = "/", mustWork = FALSE)
if (!dir.exists(local_lib)) {
  dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(unique(c(local_lib, .libPaths())))

options(repos = c(CRAN = "https://cloud.r-project.org"))

required <- c("Seurat", "remotes", "reticulate")
missing_required <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  install.packages(missing_required, lib = local_lib)
}

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  remotes::install_github("mojaveazure/seurat-disk", lib = local_lib, upgrade = "never")
}

message("Library paths:")
message(paste0("  - ", .libPaths(), collapse = "\n"))
message("Seurat setup completed.")
