#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <sstream>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]

// ----------------------------------------------------------------------
// Custom struct to hold updateCoeff results (pure C++ types)
// ----------------------------------------------------------------------
struct CoeffUpdateResult {
  std::vector<double> betatemp1; // final updated coefficients (after re-scaling)
  double R_num;                  // final numerator value
  std::vector<double> PRS;       // final polygenic risk score
};

// ----------------------------------------------------------------------
// Modified: Thread-safe internal update function using RMatrix and RVector,
//           and returning plain C++ containers (no R API objects).
// ----------------------------------------------------------------------
CoeffUpdateResult updateCoeff(const RMatrix<double>& beta_all,
                              const RMatrix<double>& GG2,
                              const RMatrix<double>& geno_ref,
                              double learningrate,
                              int maxiter,
                              const RVector<double>& betaSD,
                              const RVector<double>& sum_stats_target_val_cor,
                              int patience,
                              bool trace) {
  int n = beta_all.nrow();
  std::vector<double> u0(n), betatemp(n);
  if(beta_all.nrow() == 1) {
    // If beta_all is a one‐column matrix, assume first element is u0 and second is betatemp.
    u0[0] = beta_all(0, 0);
    betatemp[0] = beta_all(1, 0);
  } else {
    for (int i = 0; i < n; i++) {
      u0[i] = beta_all(i, 0);
      betatemp[i] = beta_all(i, 1);
    }
  }
  
  double prev_R2 = 0.0;
  int cnt = 0;
  int k = 1;
  int p = betatemp.size();
  
  int nrows = geno_ref.nrow();
  int ncols = geno_ref.ncol();  // should equal p

  std::vector<double> betatemp1(p);
  double R_num = 0.0;
  std::vector<double> PRS(nrows, 0.0);

  // Main iterative loop
  while (k <= maxiter) {
    R_num = 0.0;

    for (int j = 0; j < p; j++) {
      double beta_old = betatemp[j];
      betatemp[j] = beta_old + learningrate * u0[j];
      double diff = betatemp[j] - beta_old; // equals learningrate * u0[j]
      u0[j] = u0[j] - GG2(j, j) * diff;
      betatemp1[j] = betatemp[j] / (betaSD[j] + 1e-12); // rescaling

      // Accumulate R_num using the newly computed betatemp1[j]
      R_num += betatemp1[j] * sum_stats_target_val_cor[j];
    }

    // std::vector<double> betatemp1(p);
    // for (int i = 0; i < p; i++) {
    //   betatemp1[i] = betatemp[i] / betaSD[i];
    // }
  
    // int pp = sum_stats_target_val_cor.size();

    // if(ncols != p || p != pp) {
    //   Rcpp::stop("Dimension mismatch: geno_ref.ncol() (", ncols, 
    //             "), number of coefficients (p) (", p, 
    //             ") and sum_stats_target_val_cor.size() (", pp, 
    //             ") must be equal.");
    // }

    // double R_num = 0.0;
    // for (int i = 0; i < p; i++) {
    //   R_num += betatemp1[i] * sum_stats_target_val_cor[i];
    // }

    PRS.assign(nrows, 0.0);

    double sumPRS2 = 0.0;
    for (int i = 0; i < nrows; i++) {
        double sum = 0.0;
        for (int j = 0; j < p; j++) {  // assuming p equals ncols
            sum += geno_ref(i, j) * betatemp1[j];
        }
        PRS[i] = sum;
        sumPRS2 += sum * sum;
    }

    // Compute current R^2 while avoiding division by zero.
    double curr_R2 = (R_num * R_num) / (sumPRS2 + 1e-12);

    // std::vector<double> PRS(nrows, 0.0);
    // for (int i = 0; i < nrows; i++) {
    //   double sum = 0.0;
    //   for (int j = 0; j < ncols; j++) {
    //     sum += geno_ref(i, j) * betatemp1[j];
    //   }
    //   PRS[i] = sum;
    // }
    
    // double sumPRS2 = 0.0;
    // for (int i = 0; i < nrows; i++) {
    //   sumPRS2 += PRS[i] * PRS[i];
    // }
    // double curr_R2 = (R_num * R_num) / (sumPRS2 + 1e-12); // avoid division by zero
    
    if (curr_R2 < prev_R2) {
      cnt++;
    }
    if (cnt > patience) break;
    prev_R2 = curr_R2;
    k++;
  }
  
  // Final re-scaling of betas.
  // std::vector<double> betatemp1_final(p);
  // for (int i = 0; i < p; i++) {
  //   betatemp1_final[i] = betatemp[i] / betaSD[i];
  // }
  
  // double R_num_final = 0.0;
  // for (int i = 0; i < p; i++) {
  //   R_num_final += betatemp1_final[i] * sum_stats_target_val_cor[i];
  // }
  
  // int nrows = geno_ref.nrow();
  // int ncols = geno_ref.ncol();
  // std::vector<double> PRS_final(nrows, 0.0);
  // for (int i = 0; i < nrows; i++) {
  //   double sum = 0.0;
  //   for (int j = 0; j < ncols; j++) {
  //     sum += geno_ref(i, j) * betatemp1_final[j];
  //   }
  //   PRS_final[i] = sum;
  // }
  
  CoeffUpdateResult res;
  res.betatemp1 = betatemp1;
  res.R_num = R_num;
  res.PRS = PRS;
  
  if(trace) {
    // Note: std::cout is thread-safe for simple logging.
    std::ostringstream oss;
    oss << "stopped iteration: " << k << std::endl;
    Rcpp::Rcout << oss.str();
  }
  
  return res;
}

// ----------------------------------------------------------------------
// Custom struct to hold each block's results (pure C++ types)
// ----------------------------------------------------------------------
struct BlockResult {
  std::vector< std::vector<double> > extra_beta;  // each element is a beta vector
  std::vector<double> beta_rho;                     // corresponding R_num values
  std::vector< std::vector<double> > beta_g;        // corresponding PRS vectors
};

// ----------------------------------------------------------------------
// Modified: BlockData uses RcppParallel containers for thread safety.
// ----------------------------------------------------------------------
struct BlockData {
  RMatrix<double> beta_all;
  RMatrix<double> GG2;
  RMatrix<double> geno_ref;
  RVector<double> lr_list;
  int maxiter;
  RVector<double> sum_stats_target_val_cor;
  int patience;
  bool trace;
  RVector<double> sd;         // Declare sd before geno_info2 if used in its initialization.
  List geno_info2;            // Declare geno_info2 after sd.
  RVector<double> Beta2;

  // Constructor: Initialize in the same order as declared.
  BlockData(List block)
    : beta_all(as<NumericMatrix>(block["beta_all"])),
      GG2(as<NumericMatrix>(block["GG2"])),
      geno_ref(as<NumericMatrix>(block["geno_ref"])),
      lr_list(as<NumericVector>(block["lr_list"])),
      maxiter(as<int>(block["maxiter"])),
      sum_stats_target_val_cor(as<NumericVector>(block["sum_stats_target_val_cor"])),
      patience(as<int>(block["patience"])),
      trace(as<bool>(block["trace"])),
      sd(as<NumericVector>(as<List>(block["geno_info2"])["sd"])),
      geno_info2(as<List>(block["geno_info2"])),
      Beta2(as<NumericVector>(as<List>(block["geno_info2"])["Beta2"]))
  { }
};

BlockData convertBlockData(List block) {
  return BlockData(block);
}

// // ----------------------------------------------------------------------
// // Conversion function: convert an Rcpp List block to a BlockData struct.
// // ----------------------------------------------------------------------
// BlockData convertBlockData(List block) {
//   BlockData bd;
//   bd.beta_all = RMatrix<double>( as<NumericMatrix>(block["beta_all"]) ); // **Modified**
//   bd.GG2 = RMatrix<double>( as<NumericMatrix>(block["GG2"]) );           // **Modified**
//   bd.geno_ref = RMatrix<double>( as<NumericMatrix>(block["geno_ref"]) );   // **Modified**

//   NumericVector lr_vec = as<NumericVector>(block["lr_list"]);
//   bd.lr_list = RVector<double>(lr_vec);                                  // **Modified**

//   bd.maxiter = as<int>(block["maxiter"]);
//   bd.sum_stats_target_val_cor = RVector<double>( as<NumericVector>(block["sum_stats_target_val_cor"]) ); // **Modified**
//   bd.patience = as<int>(block["patience"]);
//   bd.trace = as<bool>(block["trace"]);

//   bd.geno_info2 = as<List>(block["geno_info2"]);                         // kept as List
//   bd.sd = RVector<double>( as<NumericVector>(bd.geno_info2["sd"]) );       // **Modified**
//   bd.Beta2 = RVector<double>( as<NumericVector>(bd.geno_info2["Beta2"]) ); // **Modified**
  
//   return bd;
// }

// ----------------------------------------------------------------------
// Modified: Parallel worker that uses only thread-safe, non-R API types.
// ----------------------------------------------------------------------
struct BlockWorker : public Worker {
  const std::vector<BlockData>& blocks;
  std::vector<BlockResult>& local_results;  // using plain C++ struct for results
  
  BlockWorker(const std::vector<BlockData>& blocks, std::vector<BlockResult>& local_results)
    : blocks(blocks), local_results(local_results) { }
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      const BlockData &bd = blocks[i];
      BlockResult blockResult;
      
      // **Modified:** Copy Beta2 (RVector) into a plain std::vector
      std::vector<double> baseBeta(bd.Beta2.size());
      for (std::size_t j = 0; j < bd.Beta2.size(); j++) {
        baseBeta[j] = bd.Beta2[j];
      }
      
      // Compute beta_rho_base = sum(baseBeta * sum_stats_target_val_cor)
      double beta_rho_base = 0.0;
      for (std::size_t j = 0; j < baseBeta.size(); j++) {
        beta_rho_base += baseBeta[j] * bd.sum_stats_target_val_cor[j];
      }
      
      // Compute beta_g_base = geno_ref %*% baseBeta
      int nrows = bd.geno_ref.nrow();

      std::vector<double> beta_g_base(nrows, 0.0);
      for (int r = 0; r < nrows; r++) {
        double s = 0.0;
        for (std::size_t j = 0; j < baseBeta.size(); j++) {
          s += bd.geno_ref(r, j) * baseBeta[j];
        }
        beta_g_base[r] = s;
      }
      
      // Store the base results.
      blockResult.extra_beta.push_back(baseBeta);
      blockResult.beta_rho.push_back(beta_rho_base);
      blockResult.beta_g.push_back(beta_g_base);
      
      // Loop over learning rates.
      for (std::size_t k = 0; k < bd.lr_list.size(); k++) {
        double cur_lr = bd.lr_list[k];
        if (cur_lr > 1.0) cur_lr = 1.0;
        // Call the internal update function (pure C++ version).
        CoeffUpdateResult res = updateCoeff(bd.beta_all, bd.GG2, bd.geno_ref,
                                              cur_lr, bd.maxiter, bd.sd,
                                              bd.sum_stats_target_val_cor, bd.patience,
                                              bd.trace);
        blockResult.extra_beta.push_back(res.betatemp1);
        blockResult.beta_rho.push_back(res.R_num);
        blockResult.beta_g.push_back(res.PRS);
      }
      
      local_results[i] = blockResult;
    }
  }
};

// ----------------------------------------------------------------------
// Exported Function: Pre-convert input blocks, run parallel worker,
// then merge results converting plain C++ types to R objects.
// ----------------------------------------------------------------------
// [[Rcpp::export]]
List block_calculation_parallel(List blocks) {
  int n = blocks.size();
  std::vector<BlockData> blockDataVec;
  blockDataVec.reserve(n);
  
  for (int i = 0; i < n; i++) {
    List blk = blocks[i];
    blockDataVec.push_back(convertBlockData(blk));
  }
  
  // Prepare a vector to hold results (one BlockResult per block)
  std::vector<BlockResult> local_results(n);
  
  // Run the parallel worker.
  BlockWorker worker(blockDataVec, local_results);
  parallelFor(0, n, worker);
  
  // Now, outside the parallel region, convert the plain C++ results into R objects.
  List out(n);
  for (int i = 0; i < n; i++) {
    // Retrieve the original geno_info2 list for output merging.
    List baseList = as<List>(blockDataVec[i].geno_info2);
    BlockResult& blockResult = local_results[i];
    
    // Add each extra_beta vector as a new column in baseList.
    for (std::size_t j = 0; j < blockResult.extra_beta.size(); j++) {
      std::vector<double>& col = blockResult.extra_beta[j];
      // Create a column name "beta_all_#" where # is j+1.
      std::ostringstream oss;
      oss << "beta_all_" << (j + 1);
      std::string colName = oss.str();
      baseList[colName] = wrap(col); // wrap std::vector<double> as R vector
    }
    
    // Create final output list for this block.
    List blockOut = List::create(
      Named("beta_byL") = DataFrame(baseList),
      Named("beta_rho") = wrap(blockResult.beta_rho),
      Named("beta_g") = wrap(blockResult.beta_g)
    );
    out[i] = blockOut;
  }
  
  return out;
}
