// #include <Rcpp.h>
// #include <vector>
// #include <string>
// #include <sstream>

// using namespace Rcpp;

// // [[Rcpp::plugins(cpp11)]]

// // ----------------------------------------------------------------------
// // Modified: Thread-safe internal update function using RMatrix and RVector,
// //           and returning plain C++ containers (no R API objects).
// // ----------------------------------------------------------------------
// List updateCoeff(NumericMatrix beta_all,
//                      NumericMatrix GG2,
//                      NumericMatrix geno_ref,
//                      double learningrate,
//                      int maxiter,
//                      NumericVector betaSD,
//                      NumericVector sum_stats_target_val_cor,
//                      int patience,
//                      bool trace) {

//   int n = beta_all.nrow();
//   NumericVector u0(n), betatemp(n);
//   if(beta_all.ncol() == 1) {
//     // If beta_all is a one‐column matrix, assume first element is u0 and second is betatemp.
//     u0[0] = beta_all(0, 0);
//     betatemp[0] = beta_all(1, 0);
//   } else {
//     for (int i = 0; i < n; i++) {
//       u0[i] = beta_all(i, 0);
//       betatemp[i] = beta_all(i, 1);
//     }
//   }
  
//   double prev_R2 = 0.0;
//   int cnt = 0;
//   int k = 1;
//   int p = betatemp.size();
  
//   int nrows = geno_ref.nrow();
//   int ncols = geno_ref.ncol();  // should equal p
  
//   NumericVector betatemp1(p);
//   double R_num = 0.0;
//   NumericVector PRS(nrows);

//   // Main iterative loop
//   while (k <= maxiter) {
//     R_num = 0.0;

//     for (int j = 0; j < p; j++) {
//       double beta_old = betatemp[j];
//       betatemp[j] = beta_old + learningrate * u0[j];
//       double diff = betatemp[j] - beta_old; // equals learningrate * u0[j]
//       u0[j] = u0[j] - GG2(j, j) * diff;
//       betatemp1[j] = betatemp[j] / (betaSD[j] + 1e-12); // rescaling

//       // Accumulate R_num using the newly computed betatemp1[j]
//       R_num += betatemp1[j] * sum_stats_target_val_cor[j];
//     }

//     // PRS.assign(nrows, 0.0);

//     double sumPRS2 = 0.0;
//     for (int i = 0; i < nrows; i++) {
//         double sum = 0.0;
//         for (int j = 0; j < p; j++) {  // assuming p equals ncols
//             sum += geno_ref(i, j) * betatemp1[j];
//         }
//         PRS[i] = sum;
//         sumPRS2 += sum * sum;
//     }

//     // Compute current R^2 while avoiding division by zero.
//     double curr_R2 = (R_num * R_num) / (sumPRS2 + 1e-12);

//     if (curr_R2 < prev_R2) {
//       cnt++;
//     }
//     if (cnt > patience) break;
//     prev_R2 = curr_R2;
//     k++;
//   }
  
//   if(trace) {
//     // Note: std::cout is thread-safe for simple logging.
//     std::ostringstream oss;
//     oss << "stopped iteration: " << k << std::endl;
//     Rcpp::Rcout << oss.str();
//   }
  
//   return List::create(Named("betatemp1") = betatemp1,
//                       Named("R_num") = R_num,
//                       Named("PRS") = PRS);
// }

