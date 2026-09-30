#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(Matrix)
  library(SingleCellExperiment)
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
  library(ggplot2)
  library(patchwork)
  library(zellkonverter)
})

# ============================================================
# Logging / small helpers
# ============================================================
log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] ", ts), sprintf(...), "\n", sep = "")
  flush.console()
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

as_char_vec <- function(x) {
  if (is.null(x)) return(NULL)
  as.character(unlist(x, use.names = FALSE))
}

as_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  as.numeric(unlist(x, use.names = FALSE))
}

timed_step <- function(label, expr) {
  start_time <- Sys.time()
  log_msg("START: %s", label)

  result <- tryCatch(
    force(expr),
    error = function(e) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      log_msg(
        "FAILED: %s | elapsed %.2f min | error: %s",
        label, elapsed, conditionMessage(e)
      )
      stop(e)
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  log_msg("DONE: %s | elapsed %.2f min", label, elapsed)
  result
}

safe_timed_step <- function(label, expr, default = NULL) {
  start_time <- Sys.time()
  log_msg("START: %s", label)

  result <- tryCatch(
    force(expr),
    error = function(e) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      log_msg(
        "WARNING: %s failed | elapsed %.2f min | error: %s",
        label, elapsed, conditionMessage(e)
      )
      return(default)
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  log_msg("DONE: %s | elapsed %.2f min", label, elapsed)
  result
}

is_drawable_plot <- function(x) {
  inherits(x, c("ggplot", "grob", "gtable", "patchwork"))
}

# ============================================================
# Config loader / defaults
# ============================================================
default_config <- function() {
  list(
    input = list(
      input_dir = NULL,
      h5ad_paths = NULL,
      file_pattern = "_hvg\\d+\\.h5ad$"
    ),
    output = list(
      out_dir = "hdwgcna_outputs",
      save_seurat_rds = TRUE,
      save_plots = TRUE,
      save_tom = TRUE,
      save_me_csv = TRUE,
      save_power_table_csv = TRUE,
      save_summary_json = TRUE
    ),
    data = list(
      assay_name = "RNA",
      # Do not include logcounts as a counts source by default.
      # NormalizeData expects raw-like counts.
      counts_layer_candidates = c("counts", "X"),
      subclass_key = "Subclass",
      donor_key = "donor_id",
      harmony_key = "donor_id",
      idents_key = "Subclass",
      min_cells = 200,
      min_genes = 50,
      use_all_genes_in_file = TRUE,
      custom_gene_list = NULL,
      transpose_if_needed = FALSE
    ),
    preprocess = list(
      normalize_method = "LogNormalize",
      scale_factor = 10000,
      nfeatures_variable = 2000,
      run_find_variable_features = FALSE,
      run_scale_data = TRUE,
      vars_to_regress = NULL,
      run_pca = TRUE,
      npcs = 30,
      # Conservative default for donor-disease-confounded data.
      run_harmony = FALSE,
      harmony_theta = NULL,
      run_umap = FALSE,
      umap_dims = 1:30,
      seed = 42
    ),
    metacell = list(
      group_by = c("Subclass", "donor_id"),
      ident_group = "Subclass",
      reduction = "pca",
      dims = 1:25,
      k = 25,
      max_shared = 10,
      min_cells = 100,
      target_metacells = 700,
      max_iter = 3000,
      mode = "average"
    ),
    network = list(
      network_type = "signed",
      cor_fnc = "bicor",
      soft_powers = c(1:10, seq(12, 30, by = 2)),
      min_power = 6,
      deepSplit = 4,
      pamStage = FALSE,
      detectCutHeight = 0.995,
      minModuleSize = 50,
      mergeCutHeight = 0.2,
      sampleForCalibration = TRUE,
      sampleForCalibrationFactor = 1000,
      useDiskCache = TRUE,
      chunkSize = NULL
    ),
    eigengenes = list(
      compute = TRUE,
      harmony_group_by_vars = NULL,
      vars_to_regress = NULL,
      scale_model_use = "linear",
      pc_dim = 1,
      compute_connectivity = TRUE,
      connectivity_group_by = "Subclass"
    ),
    runtime = list(
      threads = 8,
      stop_on_error = FALSE,
      verbose = TRUE
    )
  )
}

merge_lists <- function(base, override) {
  if (is.null(override)) return(base)
  out <- base
  for (nm in names(override)) {
    if (is.list(out[[nm]]) && is.list(override[[nm]])) {
      out[[nm]] <- merge_lists(out[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

load_config <- function(config_path) {
  user_cfg <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  merge_lists(default_config(), user_cfg)
}

# ============================================================
# Input discovery
# ============================================================
discover_h5ad_files <- function(cfg) {
  input_dir <- cfg$input$input_dir
  h5ad_paths <- cfg$input$h5ad_paths
  pattern <- cfg$input$file_pattern %||% "\\.h5ad$"

  files <- character(0)

  if (!is.null(input_dir)) {
    files <- c(files, list.files(
      path = input_dir,
      pattern = pattern,
      full.names = TRUE,
      recursive = TRUE
    ))
  }

  if (!is.null(h5ad_paths)) {
    files <- c(files, unlist(h5ad_paths, use.names = FALSE))
  }

  files <- unique(normalizePath(files, winslash = "/", mustWork = FALSE))
  files <- files[file.exists(files)]

  if (length(files) == 0) {
    stop(paste0("No input .h5ad files were found in: ", input_dir))
  }

  log_msg("Found %d h5ad files.", length(files))
  files
}

# ============================================================
# H5AD -> Seurat
# ============================================================
choose_assay_matrix <- function(sce, candidates) {
  anames <- assayNames(sce)
  log_msg("Available assays: %s", paste(anames, collapse = ", "))

  for (cand in candidates) {
    if (cand %in% anames) {
      return(list(name = cand, mat = assay(sce, cand)))
    }
  }

  if (length(anames) == 0) {
    stop("No assays found in H5AD/SCE object.")
  }

  list(name = anames[[1]], mat = assay(sce, anames[[1]]))
}

make_sparse_if_possible <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  if (inherits(x, "matrix")) return(Matrix::Matrix(x, sparse = TRUE))
  if (inherits(x, "DelayedMatrix")) return(as(x, "dgCMatrix"))
  as(x, "dgCMatrix")
}

sce_to_seurat_simple <- function(sce, cfg) {
  mat_info <- choose_assay_matrix(sce, as_char_vec(cfg$data$counts_layer_candidates))
  counts <- mat_info$mat

  if (!identical(mat_info$name, "counts")) {
    log_msg(
      "[WARN] Using assay '%s' as counts source. Ideally the H5AD should contain a true counts layer.",
      mat_info$name
    )
  }

  # Correct SCE dimensions:
  # nrow(sce) = genes/features, ncol(sce) = cells/samples.
  # Do NOT use ncol(colData(sce)); that is the number of metadata columns.
  n_genes <- nrow(sce)
  n_cells <- ncol(sce)

  log_msg("Raw matrix class: %s", paste(class(counts), collapse = ", "))
  log_msg("Raw matrix dim: %d x %d", nrow(counts), ncol(counts))
  log_msg("SCE dim: %d genes x %d cells", n_genes, n_cells)
  log_msg("rowData rows: %d | colData rows: %d | colData columns: %d",
          nrow(rowData(sce)), nrow(colData(sce)), ncol(colData(sce)))

  if (nrow(counts) == n_genes && ncol(counts) == n_cells) {
    log_msg("Matrix orientation detected: genes x cells")
  } else if (nrow(counts) == n_cells && ncol(counts) == n_genes) {
    log_msg("Matrix orientation detected: cells x genes -> transposing")
    counts <- t(counts)
  } else if (isTRUE(cfg$data$transpose_if_needed)) {
    log_msg("transpose_if_needed=TRUE -> transposing matrix as fallback")
    counts <- t(counts)
  } else {
    stop(sprintf(
      paste0(
        "Matrix dimension does not match SCE metadata. ",
        "counts=%d x %d, SCE=%d genes x %d cells. ",
        "If this is intentional, set transpose_if_needed=true, but inspect dimensions first."
      ),
      nrow(counts), ncol(counts), n_genes, n_cells
    ))
  }

  counts <- make_sparse_if_possible(counts)

  if (nrow(counts) != n_genes) {
    stop(sprintf(
      "Gene dimension mismatch after orientation fix: counts rows=%d, genes=%d",
      nrow(counts), n_genes
    ))
  }

  if (ncol(counts) != n_cells) {
    stop(sprintf(
      "Cell dimension mismatch after orientation fix: counts cols=%d, cells=%d",
      ncol(counts), n_cells
    ))
  }

  rownames(counts) <- rownames(sce)
  colnames(counts) <- colnames(sce)

  meta <- as.data.frame(colData(sce))
  rownames(meta) <- colnames(sce)

  seu <- CreateSeuratObject(
    counts = counts,
    assay = cfg$data$assay_name,
    meta.data = meta,
    min.cells = 0,
    min.features = 0
  )

  if (!is.null(reducedDimNames(sce)) && length(reducedDimNames(sce)) > 0) {
    for (rd in reducedDimNames(sce)) {
      emb <- reducedDim(sce, rd)
      if (!is.null(emb) && nrow(emb) == ncol(seu)) {
        nm <- tolower(rd)
        key <- paste0(toupper(substr(nm, 1, 1)), "_")
        seu[[nm]] <- CreateDimReducObject(
          embeddings = as.matrix(emb),
          key = key,
          assay = DefaultAssay(seu)
        )
      }
    }
  }

  log_msg("Created Seurat object: %d genes x %d cells", nrow(seu), ncol(seu))
  seu
}

# ============================================================
# Preprocessing + SetupForWGCNA
# ============================================================
prepare_seurat_for_hdwgcna <- function(seu, cfg, subclass_name, stem) {
  set.seed(cfg$preprocess$seed %||% 42)
  assay_name <- cfg$data$assay_name
  DefaultAssay(seu) <- assay_name

  if (!(cfg$data$idents_key %in% colnames(seu@meta.data))) {
    stop(sprintf("idents_key '%s' not found in metadata.", cfg$data$idents_key))
  }
  Idents(seu) <- seu@meta.data[[cfg$data$idents_key]]

  if (ncol(seu) < cfg$data$min_cells) {
    stop(sprintf("Too few cells (%d < %d).", ncol(seu), cfg$data$min_cells))
  }
  if (nrow(seu) < cfg$data$min_genes) {
    stop(sprintf("Too few genes (%d < %d).", nrow(seu), cfg$data$min_genes))
  }

  seu <- timed_step(sprintf("[%s] NormalizeData", stem), {
    NormalizeData(
      seu,
      normalization.method = cfg$preprocess$normalize_method,
      scale.factor = cfg$preprocess$scale_factor,
      verbose = FALSE
    )
  })

  if (isTRUE(cfg$preprocess$run_find_variable_features)) {
    seu <- timed_step(sprintf("[%s] FindVariableFeatures", stem), {
      FindVariableFeatures(
        seu,
        selection.method = "vst",
        nfeatures = cfg$preprocess$nfeatures_variable,
        verbose = FALSE
      )
    })
  } else {
    VariableFeatures(seu) <- rownames(seu)
    log_msg("[%s] VariableFeatures set to all genes in file: %d", stem, length(VariableFeatures(seu)))
  }

  if (isTRUE(cfg$preprocess$run_scale_data)) {
    seu <- timed_step(sprintf("[%s] ScaleData", stem), {
      ScaleData(
        seu,
        features = rownames(seu),
        vars.to.regress = as_char_vec(cfg$preprocess$vars_to_regress),
        verbose = FALSE
      )
    })
  }

  if (isTRUE(cfg$preprocess$run_pca)) {
    seu <- timed_step(sprintf("[%s] RunPCA", stem), {
      RunPCA(
        seu,
        features = VariableFeatures(seu),
        npcs = cfg$preprocess$npcs,
        verbose = FALSE
      )
    })
  }

  harmony_key <- cfg$data$harmony_key
  can_run_harmony <- isTRUE(cfg$preprocess$run_harmony) &&
    !is.null(harmony_key) &&
    harmony_key %in% colnames(seu@meta.data) &&
    length(unique(as.character(seu@meta.data[[harmony_key]]))) > 1 &&
    requireNamespace("harmony", quietly = TRUE)

  if (can_run_harmony) {
    seu <- timed_step(sprintf("[%s] RunHarmony by %s", stem, harmony_key), {
      theta <- cfg$preprocess$harmony_theta
      if (is.null(theta)) {
        harmony::RunHarmony(
          object = seu,
          group.by.vars = harmony_key,
          verbose = TRUE
        )
      } else {
        harmony::RunHarmony(
          object = seu,
          group.by.vars = harmony_key,
          theta = theta,
          verbose = TRUE
        )
      }
    })
  }

  if (isTRUE(cfg$preprocess$run_umap)) {
    red_use <- if ("harmony" %in% names(seu@reductions)) {
      "harmony"
    } else if ("pca" %in% names(seu@reductions)) {
      "pca"
    } else {
      NULL
    }

    if (!is.null(red_use)) {
      max_dim <- ncol(Embeddings(seu, red_use))
      dims_use <- as.integer(as_num_vec(cfg$preprocess$umap_dims))
      dims_use <- dims_use[dims_use <= max_dim]

      if (length(dims_use) >= 2) {
        seu <- timed_step(sprintf("[%s] RunUMAP on %s", stem, red_use), {
          RunUMAP(seu, reduction = red_use, dims = dims_use, verbose = FALSE)
        })
      }
    }
  }

  gene_features <- if (isTRUE(cfg$data$use_all_genes_in_file)) {
    rownames(seu)
  } else {
    as_char_vec(cfg$data$custom_gene_list)
  }

  if (is.null(gene_features) || length(gene_features) == 0) {
    stop("No genes available for SetupForWGCNA.")
  }
  gene_features <- intersect(as.character(gene_features), rownames(seu))
  if (length(gene_features) == 0) {
    stop("Configured gene set does not overlap with Seurat genes.")
  }

  wgcna_name <- paste0("hdWGCNA_", stem)

  seu <- timed_step(sprintf("[%s] SetupForWGCNA (%d genes)", stem, length(gene_features)), {
    SetupForWGCNA(
      seurat_obj = seu,
      wgcna_name = wgcna_name,
      features = gene_features
    )
  })

  list(seu = seu, wgcna_name = wgcna_name)
}

# ============================================================
# Plot helpers
# ============================================================
# ============================================================
# Soft-power helper
# ============================================================
choose_soft_power_for_network <- function(power_df, fallback_power = 6) {
  chosen_power <- NA_real_

  if (
    !is.null(power_df) &&
    nrow(power_df) > 0 &&
    "SFT.R.sq" %in% colnames(power_df) &&
    "Power" %in% colnames(power_df)
  ) {
    valid <- power_df[
      is.finite(power_df$SFT.R.sq) &
        is.finite(power_df$Power) &
        power_df$SFT.R.sq >= 0.8 &
        power_df$Power >= 1 &
        power_df$Power <= 50,
      ,
      drop = FALSE
    ]

    if (nrow(valid) > 0) {
      chosen_power <- min(valid$Power, na.rm = TRUE)
    }
  }

  if (!is.finite(chosen_power) || is.na(chosen_power) || chosen_power < 1 || chosen_power > 50) {
    chosen_power <- fallback_power
  }

  chosen_power <- as.numeric(chosen_power)

  if (!is.finite(chosen_power) || is.na(chosen_power) || chosen_power < 1 || chosen_power > 50) {
    chosen_power <- 6
  }

  chosen_power
}

save_soft_power_plot <- function(seu, wgcna_name, out_path) {
  plot_list <- PlotSoftPowers(seu, wgcna_name = wgcna_name)
  soft_power_plot <- patchwork::wrap_plots(plot_list, ncol = 2)
  ggplot2::ggsave(
    filename = out_path,
    plot = soft_power_plot,
    width = 10,
    height = 8
  )
}

save_dendrogram_plot <- function(seu, wgcna_name, out_path) {
  # PlotDendrogram may draw directly to the active graphics device and may not
  # return a ggplot object. Use pdf() rather than ggsave().
  opened_device <- FALSE

  tryCatch({
    grDevices::pdf(out_path, width = 12, height = 6)
    opened_device <- TRUE

    p <- PlotDendrogram(seu, wgcna_name = wgcna_name)

    # If this hdWGCNA version returns a drawable object, print it.
    # If it returns numeric/NULL after drawing directly, do nothing.
    if (is_drawable_plot(p)) {
      print(p)
    }

    grDevices::dev.off()
    opened_device <- FALSE
  }, error = function(e) {
    if (opened_device && grDevices::dev.cur() > 1) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    stop(e)
  })
}

# ============================================================
# One h5ad file
# ============================================================
run_hdwgcna_one <- function(h5ad_path, cfg) {
  stem <- tools::file_path_sans_ext(basename(h5ad_path))
  subclass_folder <- basename(dirname(h5ad_path))

  subclass_outdir <- file.path(cfg$output$out_dir, subclass_folder)
  safe_dir_create(subclass_outdir)

  tom_outdir <- file.path(subclass_outdir, "TOM")
  if (isTRUE(cfg$output$save_tom)) safe_dir_create(tom_outdir)

  log_msg("============================================================")
  log_msg("START: %s", stem)
  log_msg("Input: %s", h5ad_path)
  log_msg("Output: %s", subclass_outdir)

  sce <- timed_step(sprintf("[%s] readH5AD", stem), {
    zellkonverter::readH5AD(h5ad_path, reader = "R")
  })

  seu <- timed_step(sprintf("[%s] SCE -> Seurat", stem), {
    sce_to_seurat_simple(sce, cfg)
  })

  rm(sce)
  invisible(gc())

  subclass_key <- cfg$data$subclass_key
  subclass_name <- NA_character_
  if (subclass_key %in% colnames(seu@meta.data)) {
    subclass_vals <- unique(as.character(seu@meta.data[[subclass_key]]))
    subclass_vals <- subclass_vals[!is.na(subclass_vals)]
    if (length(subclass_vals) == 1) subclass_name <- subclass_vals[[1]]
  }
  if (is.na(subclass_name)) subclass_name <- stem

  prep <- timed_step(sprintf("[%s] preprocessing + WGCNA setup", stem), {
    prepare_seurat_for_hdwgcna(seu, cfg, subclass_name, stem)
  })

  seu <- prep$seu
  wgcna_name <- prep$wgcna_name
  seu <- SetActiveWGCNA(seu, wgcna_name)

  red_for_metacells <- cfg$metacell$reduction
  if (!(red_for_metacells %in% names(seu@reductions))) {
    if ("harmony" %in% names(seu@reductions)) {
      red_for_metacells <- "harmony"
    } else if ("pca" %in% names(seu@reductions)) {
      red_for_metacells <- "pca"
    } else {
      stop("No reduction available for MetacellsByGroups.")
    }
  }

  max_dim_red <- ncol(Embeddings(seu, red_for_metacells))
  mc_dims <- as.integer(as_num_vec(cfg$metacell$dims))
  mc_dims <- mc_dims[mc_dims <= max_dim_red]
  if (length(mc_dims) < 2) {
    stop("Too few dimensions available for metacell construction.")
  }

  group_by <- as_char_vec(cfg$metacell$group_by)
  missing_group_cols <- setdiff(group_by, colnames(seu@meta.data))
  if (length(missing_group_cols) > 0) {
    stop(sprintf(
      "Missing group_by metadata columns: %s",
      paste(missing_group_cols, collapse = ", ")
    ))
  }

  log_msg(
    "[%s] MetacellsByGroups config: reduction=%s, dims=%s, k=%d, min_cells=%d, target_metacells=%d, max_iter=%d",
    stem,
    red_for_metacells,
    paste(mc_dims, collapse = ","),
    cfg$metacell$k,
    cfg$metacell$min_cells,
    cfg$metacell$target_metacells,
    cfg$metacell$max_iter
  )

  seu <- timed_step(sprintf("[%s] MetacellsByGroups", stem), {
    MetacellsByGroups(
      seurat_obj = seu,
      group.by = group_by,
      ident.group = cfg$metacell$ident_group,
      k = cfg$metacell$k,
      reduction = red_for_metacells,
      dims = mc_dims,
      assay = cfg$data$assay_name,
      layer = "counts",
      mode = cfg$metacell$mode,
      min_cells = cfg$metacell$min_cells,
      max_shared = cfg$metacell$max_shared,
      target_metacells = cfg$metacell$target_metacells,
      max_iter = cfg$metacell$max_iter,
      verbose = cfg$runtime$verbose,
      wgcna_name = wgcna_name
    )
  })

  seu <- timed_step(sprintf("[%s] NormalizeMetacells", stem), {
    NormalizeMetacells(seu, wgcna_name = wgcna_name)
  })

  group_name <- subclass_name
  group_by_expr <- cfg$data$subclass_key
  if (!(group_by_expr %in% colnames(seu@meta.data))) {
    group_by_expr <- NULL
    group_name <- as.character(unique(Idents(seu)))[1]
  }

  seu <- timed_step(sprintf("[%s] SetDatExpr", stem), {
    log_msg("[%s] SetDatExpr: group.by=%s, group_name=%s",
            stem, group_by_expr %||% "<Idents>", group_name)
    SetDatExpr(
      seurat_obj = seu,
      group_name = group_name,
      use_metacells = TRUE,
      group.by = group_by_expr,
      assay = cfg$data$assay_name,
      layer = "data",
      wgcna_name = wgcna_name
    )
  })

  seu <- timed_step(sprintf("[%s] TestSoftPowers", stem), {
    TestSoftPowers(
      seurat_obj = seu,
      powers = as_num_vec(cfg$network$soft_powers),
      networkType = cfg$network$network_type,
      corFnc = cfg$network$cor_fnc,
      wgcna_name = wgcna_name
    )
  })

  if (isTRUE(cfg$output$save_plots)) {
    safe_timed_step(sprintf("[%s] Save soft-power plot", stem), {
      save_soft_power_plot(
        seu,
        wgcna_name,
        file.path(subclass_outdir, sprintf("%s_soft_powers.pdf", stem))
      )
    })
  }

  power_df_for_network <- GetPowerTable(seu, wgcna_name = wgcna_name)
  chosen_power <- choose_soft_power_for_network(power_df_for_network, cfg$network$min_power)
  log_msg("[%s] Chosen soft_power for ConstructNetwork: %s", stem, chosen_power)

  seu <- timed_step(sprintf("[%s] ConstructNetwork", stem), {
    ConstructNetwork(
      seurat_obj = seu,
      tom_outdir = if (isTRUE(cfg$output$save_tom)) tom_outdir else tempdir(),
      tom_name = stem,
      soft_power = chosen_power,
      networkType = cfg$network$network_type,
      deepSplit = cfg$network$deepSplit,
      pamStage = cfg$network$pamStage,
      detectCutHeight = cfg$network$detectCutHeight,
      minModuleSize = cfg$network$minModuleSize,
      mergeCutHeight = cfg$network$mergeCutHeight,
      sampleForCalibration = cfg$network$sampleForCalibration,
      sampleForCalibrationFactor = cfg$network$sampleForCalibrationFactor,
      useDiskCache = cfg$network$useDiskCache,
      chunkSize = cfg$network$chunkSize,
      wgcna_name = wgcna_name
    )
  })

  if (isTRUE(cfg$output$save_plots)) {
    safe_timed_step(sprintf("[%s] Save dendrogram plot", stem), {
      save_dendrogram_plot(
        seu,
        wgcna_name,
        file.path(subclass_outdir, sprintf("%s_dendrogram.pdf", stem))
      )
    })
  }

  if (isTRUE(cfg$eigengenes$compute)) {
    harmony_vars <- as_char_vec(cfg$eigengenes$harmony_group_by_vars)
    if (!is.null(harmony_vars)) {
      harmony_vars <- harmony_vars[harmony_vars %in% colnames(seu@meta.data)]
      if (length(harmony_vars) == 0) harmony_vars <- NULL
    }

    seu <- timed_step(sprintf("[%s] ModuleEigengenes", stem), {
      ModuleEigengenes(
        seurat_obj = seu,
        group.by.vars = harmony_vars,
        vars.to.regress = as_char_vec(cfg$eigengenes$vars_to_regress),
        scale.model.use = cfg$eigengenes$scale_model_use,
        pc_dim = cfg$eigengenes$pc_dim,
        wgcna_name = wgcna_name,
        verbose = cfg$runtime$verbose
      )
    })

    if (isTRUE(cfg$eigengenes$compute_connectivity)) {
      conn_group_by <- cfg$eigengenes$connectivity_group_by
      if (is.null(conn_group_by) || !(conn_group_by %in% colnames(seu@meta.data))) {
        conn_group_by <- group_by_expr
      }

      # Connectivity can warn or fail if TOM is not discoverable.
      # Keep the run alive; kME columns are still often computed successfully.
      seu <- safe_timed_step(sprintf("[%s] ModuleConnectivity", stem), {
        ModuleConnectivity(
          seurat_obj = seu,
          group.by = conn_group_by,
          group_name = group_name,
          wgcna_name = wgcna_name
        )
      }, default = seu)
    }
  }

  modules_df <- GetModules(seu, wgcna_name = wgcna_name)
  modules_df$kME_self <- NA_real_

  for (m in unique(modules_df$module)) {
    kme_col <- paste0("kME_", m)
    if (kme_col %in% colnames(modules_df)) {
      idx <- modules_df$module == m
      modules_df$kME_self[idx] <- modules_df[[kme_col]][idx]
    }
  }

  modules_non_grey <- modules_df[modules_df$module != "grey", , drop = FALSE]
  if (nrow(modules_non_grey) > 0) {
    top_hubs <- do.call(
      rbind,
      lapply(split(modules_non_grey, modules_non_grey$module), function(df) {
        df <- df[order(-df$kME_self, na.last = TRUE), , drop = FALSE]
        head(df, 20)
      })
    )
  } else {
    top_hubs <- modules_non_grey
  }

  power_df <- GetPowerTable(seu, wgcna_name = wgcna_name)
  me_mat <- NULL
  if (isTRUE(cfg$eigengenes$compute)) {
    me_mat <- GetMEs(seu, wgcna_name = wgcna_name)
  }

  timed_step(sprintf("[%s] Save outputs", stem), {
    write.csv(
      modules_df,
      file.path(subclass_outdir, sprintf("%s_modules.csv", stem)),
      row.names = FALSE
    )

    write.csv(
      modules_df,
      file.path(subclass_outdir, sprintf("%s_modules_with_kME_self.csv", stem)),
      row.names = FALSE
    )

    write.csv(
      top_hubs,
      file.path(subclass_outdir, sprintf("%s_top20_hub_genes_by_module.csv", stem)),
      row.names = FALSE
    )

    if (isTRUE(cfg$output$save_power_table_csv)) {
      write.csv(
        power_df,
        file.path(subclass_outdir, sprintf("%s_power_table.csv", stem)),
        row.names = FALSE
      )
    }

    if (isTRUE(cfg$output$save_me_csv) && !is.null(me_mat)) {
      me_df <- as.data.frame(me_mat)
      me_df$barcode <- rownames(me_df)
      write.csv(
        me_df,
        file.path(subclass_outdir, sprintf("%s_hMEs.csv", stem)),
        row.names = FALSE
      )
    }

    if (isTRUE(cfg$output$save_seurat_rds)) {
      saveRDS(seu, file.path(subclass_outdir, sprintf("%s_hdWGCNA.rds", stem)))
    }
  })

  module_sizes <- sort(table(as.character(modules_df$module)), decreasing = TRUE)
  module_sizes <- module_sizes[names(module_sizes) != "grey"]

  soft_power_selected <- NA_real_
  if (nrow(power_df) > 0 && "SFT.R.sq" %in% colnames(power_df) && "Power" %in% colnames(power_df)) {
    valid <- power_df[power_df$SFT.R.sq >= 0.8, , drop = FALSE]
    if (nrow(valid) > 0) {
      soft_power_selected <- min(valid$Power, na.rm = TRUE)
    }
  }

  summary_list <- list(
    input_h5ad = normalizePath(h5ad_path, winslash = "/", mustWork = FALSE),
    stem = stem,
    subclass_folder = subclass_folder,
    subclass_name = subclass_name,
    output_dir = subclass_outdir,
    n_cells = ncol(seu),
    n_genes = nrow(seu),
    n_modules_non_grey = as.integer(length(module_sizes)),
    module_sizes_non_grey = as.list(as.integer(module_sizes)),
    module_names_non_grey = names(module_sizes),
    soft_power_selected_rule_ge_0_8 = soft_power_selected,
    wgcna_name = wgcna_name
  )
  names(summary_list$module_sizes_non_grey) <- names(module_sizes)

  if (isTRUE(cfg$output$save_summary_json)) {
    jsonlite::write_json(
      summary_list,
      path = file.path(subclass_outdir, sprintf("%s_summary.json", stem)),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }

  log_msg(
    "DONE: %s | cells=%d genes=%d non-grey-modules=%d",
    stem, ncol(seu), nrow(seu), length(module_sizes)
  )

  invisible(list(
    stem = stem,
    out_dir = subclass_outdir,
    summary = summary_list
  ))
}

# ============================================================
# All files
# ============================================================
run_all <- function(cfg) {
  safe_dir_create(cfg$output$out_dir)
  set.seed(cfg$preprocess$seed %||% 42)

  WGCNA::enableWGCNAThreads(nThreads = cfg$runtime$threads %||% 1)

  files <- discover_h5ad_files(cfg)
  log_msg("Discovered %d input file(s).", length(files))

  results <- vector("list", length(files))
  failures <- list()

  for (i in seq_along(files)) {
    f <- files[[i]]
    stem <- tools::file_path_sans_ext(basename(f))
    log_msg("[%d/%d] Processing %s", i, length(files), stem)

    res <- tryCatch(
      run_hdwgcna_one(f, cfg),
      error = function(e) {
        msg <- conditionMessage(e)
        log_msg("ERROR in %s: %s", stem, msg)
        list(stem = stem, error = msg)
      }
    )

    results[[i]] <- res

    if (!is.null(res$error)) {
      failures[[length(failures) + 1]] <- list(stem = stem, error = res$error)
      if (isTRUE(cfg$runtime$stop_on_error)) {
        stop(sprintf("Stopping on first error: %s", res$error))
      }
    }

    invisible(gc())
  }

  summary_rows <- lapply(results, function(x) {
    if (!is.null(x$error)) {
      data.frame(
        stem = x$stem,
        status = "error",
        out_dir = NA_character_,
        n_cells = NA_integer_,
        n_genes = NA_integer_,
        n_modules_non_grey = NA_integer_,
        error = x$error,
        stringsAsFactors = FALSE
      )
    } else {
      s <- x$summary
      data.frame(
        stem = x$stem,
        status = "ok",
        out_dir = x$out_dir,
        n_cells = s$n_cells,
        n_genes = s$n_genes,
        n_modules_non_grey = s$n_modules_non_grey,
        error = "",
        stringsAsFactors = FALSE
      )
    }
  })

  summary_df <- do.call(rbind, summary_rows)
  write.csv(
    summary_df,
    file.path(cfg$output$out_dir, "hdwgcna_batch_summary.csv"),
    row.names = FALSE
  )

  log_msg(
    "Batch finished: %d success / %d failure",
    sum(summary_df$status == "ok"),
    sum(summary_df$status == "error")
  )

  invisible(summary_df)
}

# ============================================================
# Main
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || !args[[1]] %in% c("--config", "-c")) {
  cat(
    "Usage:\n",
    "  Rscript making_hdWGCNA_from_residual_genes.r --config hdwgcna_cfg.json\n",
    sep = ""
  )
  quit(status = 1)
}

cfg_path <- args[[2]]
if (!file.exists(cfg_path)) {
  stop(sprintf("Config file not found: %s", cfg_path))
}

cfg <- load_config(cfg_path)
run_all(cfg)

