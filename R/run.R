#' @importFrom utils modifyList
NULL

PTL_PRS_defaults <- function() {
    #' PTL-PRS default configuration
    #' @keywords internal
    #' @noRd
  list(
    ref_file               = NULL,
    sum_stats_file         = NULL,
    target_sumstats_file   = NULL,
    subprop                = NULL,
    ref_file_ps            = NULL,
    LDblocks               = "EUR.hg19",
    outfile                = NULL,
    num_cores              = 1L,
    target_sumstats_train_file = NULL,
    target_sumstats_val_file   = NULL,
    ps                     = TRUE,
    pseudo_test            = FALSE,    
    random_seed            = 42L,
    lr_list                = c(1,10,100,1000),
    iter                   = 100,
    patience               = 3,
    trace                  = FALSE,
    random_seed_test       = 70L
  )
}

#' Build a PTL-PRS configuration list
#'
#' Constructs a validated configuration list for \code{PTL_PRS_train} / \code{PTL_PRS_run}.
#' You pass only the fields you want to override; unspecified fields fall back to
#' \code{PTL_PRS_defaults()}.
#' @param ... Named fields to override in the defaults.
#' @export
PTL_PRS_config <- function(...) {
  cfg <- modifyList(PTL_PRS_defaults(), list(...), keep.null=TRUE)

  required <- c("ref_file","sum_stats_file","target_sumstats_file",
                "subprop","ref_file_ps","outfile")

  missing <- required[vapply(required, function(k) is.null(cfg[[k]]), logical(1))]

  if (length(missing)) {
    stop("Missing required config fields: ", paste(missing, collapse = ", "))
  }
  
  cfg
}

.PTL_PRS_run_once <- function(cfg) {
  # Runs a single seed with the current cfg (expects cfg$outfile/seed already set)
  do.call(PTL_PRS_train, cfg)
}

#' Repeat the whole PTL-PRS pipeline across seeds (minimal args)
#'
#' @param cfg A list created by PTL_PRS_config(...).
#' @param num_repeat_ps Integer; repeats when ps=TRUE and >1, else runs once.
#' @param random_seeds_repeat Optional integer vector; if NULL, sampled randomly.
#' @export
#' @return List: results_by_seed, seeds, outfile_base, summary
PTL_PRS_run <- function(cfg,
                        num_repeat_ps     = 1L,
                        random_seeds_repeat = NULL) {

  do_repeat <- isTRUE(cfg$ps) && num_repeat_ps > 1L

  # Seeds
  if (do_repeat) {
    if (is.null(random_seeds_repeat)) {
      set.seed(2025)
      seeds <- sample.int(1e9, num_repeat_ps)
    } else {
        seeds <- as.integer(random_seeds_repeat)               
        seeds <- seeds[!is.na(seeds) & seeds > 0L]            
        if (any(duplicated(seeds))) {
        warning("Duplicate seeds provided; using unique set.") 
        seeds <- unique(seeds)                                 
        }
      if (length(seeds) < num_repeat_ps) {
        stop("Not enough seeds: got ", length(seeds),
            ", need ", num_repeat_ps,
            ". Supply `random_seeds_repeat` with at least ", num_repeat_ps, " unique seeds.") 
      } else if (length(seeds) > num_repeat_ps) {
        warning("Got larger number of seeds; using first ", num_repeat_ps, " seeds.") 
        seeds <- seeds[seq_len(num_repeat_ps)]
      }
    }
  } else {
    seeds <- as.integer(cfg$random_seed)
  }

  outfile_base <- cfg$outfile
  results_by_seed <- vector("list", length(seeds))
  names(results_by_seed) <- paste0("seed_", seeds)

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    cfg_i <- cfg
    cfg_i$random_seed <- s
    cfg_i$outfile <- if (do_repeat) paste0(outfile_base, "_seed", s) else outfile_base

    message(sprintf("\n [PTL-PRS] seed=%d  outfile='%s'", s, cfg_i$outfile))
    results_by_seed[[i]] <- .PTL_PRS_run_once(cfg_i)
  }

    # -------------- light numeric summary --------------
  summary_rows <- lapply(seq_along(results_by_seed), function(i) {
    s <- seeds[i]; res <- results_by_seed[[i]]
    best_param <- tryCatch(as.numeric(res$best.param), error = function(e) NA_real_)  
    best_R2    <- tryCatch(max(as.numeric(res$R2.list), na.rm = TRUE), error = function(e) NA_real_) 
    data.frame(seed = s, best_param = best_param, best_R2 = best_R2)
  })
  summary_df <- do.call(rbind, summary_rows)

  list(
    results_by_seed = results_by_seed,
    seeds           = seeds,
    outfile_base    = outfile_base,
    summary         = summary_df
  )
}
