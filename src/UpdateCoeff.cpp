#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <sstream>
// #include <sys/resource.h>
// #include <iostream>

#include <mutex>
#include "input_queue.h"

using namespace Rcpp;
using namespace RcppParallel;

typedef InputQueue<List> QueueType;

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
// Thread-safe internal update function using RMatrix and RVector,
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
  // int ncols = geno_ref.ncol();  // should equal p

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

    std::fill(PRS.begin(), PRS.end(), 0.0);

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

    if (curr_R2 < prev_R2) {
      cnt++;
    }
    if (cnt > patience) break;
    prev_R2 = curr_R2;
    k++;
  }
  
  CoeffUpdateResult res;
  res.betatemp1 = std::move(betatemp1);
  res.R_num = R_num;
  res.PRS = std::move(PRS);
  
  if(trace) {
    // Note: std::cout is thread-safe for simple logging.
    std::ostringstream oss;
    oss << "stopped iteration: " << k << std::endl;
    Rcpp::Rcout << oss.str();
  }

  // u0.clear();
  // u0.shrink_to_fit();
  // betatemp.clear();
  // betatemp.shrink_to_fit();
  // geno_ref.clear();
  // geno_ref.shrink_to_fit();
  // beta_all.clear();
  // beta_all.shrink_to_fit();

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
// BlockData uses RcppParallel containers for thread safety.
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

// ----------------------------------------------------------------------
// Parallel worker that uses only thread-safe, non-R API types.
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
      
      // // **Modified:** Copy Beta2 (RVector) into a plain std::vector
      // std::vector<double> baseBeta(bd.Beta2.size());
      // for (std::size_t j = 0; j < bd.Beta2.size(); j++) {
      //   baseBeta[j] = bd.Beta2[j];
      // }

      std::vector<double> baseBeta(bd.Beta2.begin(), bd.Beta2.end());      

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
      // Move the computed base results into blockResult to avoid copying.
      blockResult.extra_beta.push_back(std::move(baseBeta));
      blockResult.beta_rho.push_back(beta_rho_base);
      blockResult.beta_g.push_back(std::move(beta_g_base));
      
      // Loop over learning rates.
      for (std::size_t k = 0; k < bd.lr_list.size(); k++) {
        double cur_lr = bd.lr_list[k];
        if (cur_lr > 1.0) cur_lr = 1.0;
        // Call the internal update function (pure C++ version).
        CoeffUpdateResult res = updateCoeff(bd.beta_all, bd.GG2, bd.geno_ref,
                                              cur_lr, bd.maxiter, bd.sd,
                                              bd.sum_stats_target_val_cor, bd.patience,
                                              bd.trace);
        blockResult.extra_beta.push_back(std::move(res.betatemp1));
        blockResult.beta_rho.push_back(res.R_num);
        blockResult.beta_g.push_back(std::move(res.PRS));
      }
      
      // Move blockResult into local_results to avoid an extra copy.
      local_results[i] = std::move(blockResult);
    }
  }
};

// void printMemoryUsage(const std::string& label) {
//   struct rusage usage;
//   if(getrusage(RUSAGE_SELF, &usage) == 0) {
//     std::cout << label << " memory usage: " 
//               << usage.ru_maxrss << " kilobytes" << std::endl;
//   } else {
//     std::cerr << "Failed to get memory usage" << std::endl;
//   }
// }

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
  
  // printMemoryUsage("After converting blocks");

  // Prepare a vector to hold results (one BlockResult per block)
  std::vector<BlockResult> local_results(n);
  
  // Run the parallel worker.
  BlockWorker worker(blockDataVec, local_results);
  parallelFor(0, n, worker);
  // printMemoryUsage("After parallel worker");

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
  
  // printMemoryUsage("Before returning out");

  return out;
}

// [[Rcpp::export]]
SEXP create_input_queue() {
  Rcpp::XPtr<QueueType> ptr(new QueueType(), true);
  return ptr;
}

// [[Rcpp::export]]
void push_input(SEXP queue_ptr, Rcpp::List input) {
  Rcpp::XPtr<QueueType> queue(queue_ptr);
  queue->push(input);
}

// [[Rcpp::export]]
void finish_queue(SEXP queue_ptr) {
  Rcpp::XPtr<QueueType> queue(queue_ptr);
  queue->set_finished();
}

// [[Rcpp::export]]
List block_calculation_parallel_streamed(SEXP queue_ptr) {
  Rcpp::XPtr<QueueType> queue(queue_ptr);
  std::vector<BlockData> blockDataVec;

  // ------------------------------
  // 1. Pull blocks one-by-one from R-prepared queue
  // ------------------------------
  Rcout << "Receiving block data from R (via queue)..." << std::endl;
  while (true) {
    List blk = queue->pop();
    if (Rf_isNull(blk) || blk.size() == 0) break; // end-of-queue signal
    blockDataVec.push_back(convertBlockData(blk));
  }

  int n = blockDataVec.size();
  Rcout << "Received " << n << " blocks from R" << std::endl;

  // printMemoryUsage("After converting blocks");

  // ------------------------------
  // 2. Allocate result container
  // ------------------------------
  std::vector<BlockResult> local_results(n);

  // ------------------------------
  // 3. Run parallel worker
  // ------------------------------
  BlockWorker worker(blockDataVec, local_results);
  parallelFor(0, n, worker);

  // printMemoryUsage("After parallel worker");

  // ------------------------------
  // 4. Convert results back to R format
  // ------------------------------
  List out(n);
  for (int i = 0; i < n; i++) {
    List baseList(blockDataVec[i].geno_info2);
    BlockResult& blockResult = local_results[i];

    for (std::size_t j = 0; j < blockResult.extra_beta.size(); j++) {
      std::ostringstream oss;
      oss << "beta_all_" << (j + 1);
      baseList[oss.str()] = wrap(blockResult.extra_beta[j]);
    }

    List blockOut = List::create(
      Named("beta_byL") = DataFrame(baseList),
      Named("beta_rho") = wrap(blockResult.beta_rho),
      Named("beta_g") = wrap(blockResult.beta_g)
    );
    out[i] = blockOut;
  }

  // printMemoryUsage("Before returning result");

  return out;
}
