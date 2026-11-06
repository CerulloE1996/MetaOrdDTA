// NicoSIMD version of Stan_fns_ordinal.stan written by NicoStan (included file md5 c1ca9159a8874883399702da0ea64f50; NicoSIMD declarations md5 e3c567fd115d2ccf49ec504688009e6b; writer implementation md5 aef1406a74d6694a1d260c4cd31c7a35)


////
//// Custom Stan functions - ordinal regression + cutpoint functions  --------------------------------------------------------------
////




real compute_AUC(  vector se_vec, 
                   vector sp_vec, 
                   int n_thr, 
                   real missing_value_marker) {
                       
      // Calculate AUC for a single test
      // Inputs:
      //   se_vec: vector of sensitivities (length n_thr)
      //   sp_vec: vector of specificities (length n_thr)
      //   n_thr: number of thresholds for this test
      //   invalid_threshold: value below which Se/Sp are considered invalid (typically -0.5)
      // Returns:
      //   AUC value
      
      // Count valid (non-missing) thresholds
      int n_valid = 0;
      for (k in 1:n_thr) {
          // Check if this threshold has actual data (not missing)
          if (se_vec[k] != missing_value_marker && sp_vec[k] != missing_value_marker) {
              n_valid += 1;
          }
      }
      ////
      // If no valid thresholds, return NaN or -1.0
      if (n_valid == 0) {
          return -1.0;
      }
      ////
      // We need to add points (0,0) and (1,1) to complete the ROC curve
      int n_points = n_valid + 2;
      vector[n_points] FPR;
      vector[n_points] TPR;
      ////
      // Add (0,0) point
      FPR[1] = 0.0;
      TPR[1] = 0.0;
      ////
      // Add only the VALID threshold points
      int idx = 2;
      for (k in 1:n_thr) {
          if (se_vec[k] != missing_value_marker && sp_vec[k] != missing_value_marker) {
              FPR[idx] = 1.0 - sp_vec[k];
              TPR[idx] = se_vec[k];
              idx += 1;
          }
      }
      ////
      // Add (1,1) point
      FPR[n_points] = 1.0;
      TPR[n_points] = 1.0;
      ////
      // Sort by FPR using bubble sort
      for (i in 1:(n_points-1)) {
          for (j in 1:(n_points-i)) {
              if (FPR[j] > FPR[j+1]) {
                  // Swap FPR
                  real temp_fpr = FPR[j];
                  FPR[j] = FPR[j+1];
                  FPR[j+1] = temp_fpr;
                  // Swap TPR
                  real temp_tpr = TPR[j];
                  TPR[j] = TPR[j+1];
                  TPR[j+1] = temp_tpr;
              }
          }
      }
      ////
      // Calculate AUC using trapezoidal rule
      real auc = 0.0;
      for (i in 1:(n_points-1)) {
          auc += 0.5 * (TPR[i] + TPR[i+1]) * (FPR[i+1] - FPR[i]);
      }
      ////
      return auc;
    
}


real SD_approx_ID_ord_prob_to_C_probit( real mean_C_cutpoint_scale, 
                                        real SD_prob_scale) {
                                 
    real SD_probit_scale = (1.0 / std_normal_pdf(mean_C_cutpoint_scale)) * SD_prob_scale;
    return SD_probit_scale;
                            
}
vector SD_approx_ID_ord_prob_to_C_probit( vector mean_C_cutpoint_scale, 
                                          vector SD_prob_scale) {
                                            
    int dim = num_elements(SD_prob_scale);
    vector[dim] SD_probit_scale = (1.0 / std_normal_pdf(mean_C_cutpoint_scale)) .* SD_prob_scale;
    return SD_probit_scale;
                            
}
////
//// ---- logit scale:
////
real SD_approx_ID_ord_prob_to_C_logit(  real mean_prob_scale, 
                                        real SD_prob_scale) {
                                 
    real SD_logit_scale = (1.0 / (mean_prob_scale * (1.0 - mean_prob_scale))) * SD_prob_scale;
    return SD_logit_scale;
                            
}
vector SD_approx_ID_ord_prob_to_C_logit( vector mean_prob_scale, 
                                         vector SD_prob_scale) {
                                           
    int dim = num_elements(SD_prob_scale);
    vector[dim] SD_logit_scale = (1.0 ./ (mean_prob_scale .* (1.0 - mean_prob_scale))) .* SD_prob_scale;
    return SD_logit_scale;
                            
}
 
  
  
////
//// Function to manually construct a cutpoint vector from raw unconstrained parameters - using exp() / log-differences. 
////
vector construct_C(  vector C_raw_vec, 
                     int softplus) {
  
  
          int n_total_cutpoints = num_elements(C_raw_vec);
          vector[n_total_cutpoints] C_vec;  
          //// ----  1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          C_vec[1] = C_raw_vec[1];
          //// ---- Rest of cutpoints are made using LOG-differences (OR softplus-differences):
          if (softplus == 1) {
                  vector[n_total_cutpoints - 1] softplus_C_vec = nico_log1p_exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + softplus_C_vec[k - 1];
                             // jacobian += log_inv_logit(C_raw_vec[k]);   //// Jacobian for trans. C_raw -> C:
                  }
          } else { 
                  vector[n_total_cutpoints - 1] exp_C_vec = nico_exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + exp_C_vec[k - 1];
                             // jacobian += C_raw_vec[k];  //// Jacobian for trans. C_raw -> C:
                  }
          }
          //// ---- Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}
row_vector construct_C(  row_vector C_raw_vec, 
                         int softplus) {
  
          int n_total_cutpoints = num_elements(C_raw_vec);
          row_vector[n_total_cutpoints] C_vec;  
          //// ----  1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          C_vec[1] = C_raw_vec[1];
          //// ---- Rest of cutpoints are made using LOG-differences (OR softplus-differences):
          if (softplus == 1) {
                  row_vector[n_total_cutpoints - 1] softplus_C_vec = nico_log1p_exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + softplus_C_vec[k - 1];
                             // jacobian += log_inv_logit(C_raw_vec[k]);   //// Jacobian for trans. C_raw -> C:
                  }
          } else { 
                  row_vector[n_total_cutpoints - 1] exp_C_vec = nico_exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + exp_C_vec[k - 1];
                             // jacobian += C_raw_vec[k];  //// Jacobian for trans. C_raw -> C:
                  }
          }
          //// ---- Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}




//// ---------------------------------------------------------------------------------------
//// Induced-Dirichlet ("ind_dir") log-density function:
//// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
//// NOTE: adapted from: Betancourt et al (see: https://betanalpha.github.io/assets/case_studies/ordinal_regression.html),
//// HOWEVER my version has a (much) more computationally efficient (lower-trianglar) Jacobian computation, which is 
//// mathematically still valid. 
////
real induced_dirichlet_given_C_lpdf(    vector p_ord,
                                        vector C,
                                        vector alpha,
                                        int use_probit_link) {
          
        int n_cat = num_elements(p_ord);
        int n_thr = n_cat - 1;
        
        real log_prob = 0.0;

        for (k in 1:n_thr) {
          if (use_probit_link == 1)  log_prob += nico_normal_lupdf(C[k] | 0.0, 1.0);
          else                       log_prob += nico_log_inv_logit(C[k]) + nico_log1m_inv_logit(C[k]);
        }
        log_prob += nico_dirichlet_lupdf(p_ord | alpha);
        
        return log_prob;
  
}
////
//// Row-vector overload for "induced_dirichlet_given_rho_lpdf":
////
real induced_dirichlet_given_C_lpdf(  row_vector p_ord,
                                      row_vector C,
                                      row_vector alpha,
                                      int use_probit_link) {
          

        int n_cat = num_elements(p_ord);
        int n_thr = n_cat - 1;
        vector[n_cat] p_ord_col_vec = to_vector(p_ord);
        vector[n_thr] C_col_vec     = to_vector(C);
        vector[n_cat] alpha_col_vec = to_vector(alpha);
        
        return induced_dirichlet_given_C_lupdf(p_ord_col_vec | C_col_vec, alpha_col_vec, use_probit_link);
  
}





//// ---------------------------------------------------------------------------------------
//// Induced-Dirichlet ("ind_dir") log-density function:
//// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
//// NOTE: adapted from: Betancourt et al (see: https://betanalpha.github.io/assets/case_studies/ordinal_regression.html)
////
// real induced_dirichlet_given_rho_lpdf(  vector p_ord,
//                                         vector rho,
//                                         vector alpha) {
          
//         int n_cat = num_elements(p_ord);
//         matrix[n_cat, n_cat] J = rep_matrix(0.0, n_cat, n_cat);
         
//         //// Jacobian computation:
//         for (k in 1:n_cat) {
//                 J[k, 1] = 1.0;
//         }
//         // // // for (k in 2:(n_cat - 1)) {
//         // // //    J[k, k] = +rho[k];
//         // // // }
//         // for (k in 2:n_cat) { // p_ord[k] = Phi(C_{k} - anchor) - Phi(C_{k - 1} - anchor)
//         //          J[k, k]     = - rho[k - 1];
//         //          J[k - 1, k] = + rho[k - 1];
//         // }
                
//         // First row (p_1 depends on c_1)
//         J[1, 2] = -rho[1];  // ∂p_1/∂c_1 = -ρ(φ-c_1)  [SIGN FLIPPED]
        
//         // Middle rows
//         for (k in 2:(n_cat-1)) {
//             J[k, k]     = + rho[k - 1];     // ∂p_k/∂c_{k-1} = ρ(φ-c_{k-1})  [SIGN FLIPPED] 
//             J[k, k - 1] = - rho[k];    // ∂p_k/∂c_k = -ρ(φ-c_k)  [SIGN FLIPPED]
//         }
        
//         // Last row (p_n_cat depends on c_{n_cat-1})
//         J[n_cat, n_cat] = rho[n_cat-1];  // ∂p_n_cat/∂c_{n_cat-1} = ρ(φ-c_{n_cat-1})  [SIGN FLIPPED]

//         return dirichlet_lpdf(p_ord | alpha) + log_determinant(J);
  
// }
// ////
// //// Row-vector overload for "induced_dirichlet_given_rho_lpdf":
// ////
// real induced_dirichlet_given_rho_lpdf(  row_vector p_ord,
//                                  row_vector rho,
//                                  row_vector alpha) {
          

//         int n_cat = num_elements(p_ord);
//         int n_thr = n_cat - 1;
//         vector[n_cat] p_ord_col_vec = to_vector(p_ord);
//         vector[n_thr] rho_col_vec   = to_vector(rho);
//         vector[n_cat] alpha_col_vec = to_vector(alpha);
        
//         return induced_dirichlet_given_rho_lpdf(p_ord_col_vec | rho_col_vec, alpha_col_vec);
  
// }





//// Induced dirichlet from Betancourt, 2019. 
//// See https://betanalpha.github.io/assets/case_studies/ordinal_regression.html#3_cut_to_the_chase
// real induced_dirichlet_lpdf( vector c, 
//                              vector alpha, 
//                              real phi) {
                                   
//       int K = num_elements(c) + 1;
//       vector[K - 1] anchoredcutoffs = c - phi;
      
//       vector[K] sigma;
//       vector[K] p;
//       matrix[K, K] J = rep_matrix(0, K, K);
      
//       sigma[1:(K-1)] = Phi(anchoredcutoffs); 
//       sigma[K] = 1;
      
//       p[1] =  sigma[1];
//       for (k in 2:(K - 1))
//         p[k] = sigma[k] - sigma[k - 1];
//       p[K] = 1 - sigma[K - 1];
      
//       // Baseline column of Jacobian
//       for (k in 1:K) J[k, 1] = 1;
      
//       // Diagonal entries of Jacobian
//       for (k in 2:K) {
//         // rho is the PDF of the latent distribution
//         real rho = exp(normal_lpdf(anchoredcutoffs[k - 1] | 0.0, 1.0)); ////  1.702 * sigma[k - 1] * (1 - sigma[k - 1]);
//         J[k, k] =  - rho;
//         J[k - 1, k] = + rho;
//       }
      
//       return dirichlet_lpdf(p | alpha) + log_determinant(J);
  
// }


////
//// Convert from cumul_probs -> ord_probs:
////
vector cumul_probs_to_ord_probs(vector cumul_probs) {
  
        int n_thr = num_elements(cumul_probs);
        int n_cat = n_thr + 1;
        vector[n_cat] ord_probs;
       
        ord_probs[1] = cumul_probs[1] - 0.0;
        for (k in 2:n_thr) {
             ord_probs[k] = cumul_probs[k] - cumul_probs[k - 1]; // since probs are INCREASING with k
        }
        ord_probs[n_cat] =  1.0 - cumul_probs[n_cat - 1];
        
        return ord_probs;
        
}
////
//// Convert from cumul_probs -> ord_probs (row-vector overload):
////
row_vector cumul_probs_to_ord_probs(row_vector cumul_probs) {
  
       return to_row_vector(cumul_probs_to_ord_probs(to_vector(cumul_probs)));
        
}
////
//// Convert from cumul_probs -> ord_probs (matrix overload)
//// NOTE: This will be more efficient when n_rows > n_cols (e.g. if each row is a study this will often be the case).
//// NOTE: assumed each col corresponds to a threshold/cutpoint.
////
matrix cumul_probs_to_ord_probs(matrix cumul_probs) {
  
        int n_cols = cols(cumul_probs);
        int n_rows = rows(cumul_probs);
        int n_thr = n_cols;
        int n_cat = n_thr + 1;
        matrix[n_rows, n_cat] ord_probs;
        
        // for (k in 1:n_cols) { //// col-major storage
           for (s in 1:n_rows) {
                 ord_probs[s, 1] = cumul_probs[s, 1] - 0.0;
                 for (k in 2:n_thr) {
                     ord_probs[s, k] = cumul_probs[s, k] - cumul_probs[s, k - 1]; // since probs are INCREASING with k
                 }
                 ord_probs[s, n_cat] =  1.0 - cumul_probs[s, n_cat - 1];
           }
        // }
        
        return ord_probs;
         
}



////
//// Convert from ord_probs -> cumul_probs:
////
vector ord_probs_to_cumul_probs( vector ord_probs) {
  
        int n_cat = num_elements(ord_probs);
        int n_thr = n_cat - 1;
        vector[n_thr] cumul_probs;
       
        cumul_probs[1] = ord_probs[1];
        for (k in 2:n_thr) {
             cumul_probs[k] = cumul_probs[k - 1] + ord_probs[k]; // since probs are INCREASING with k
        }
        
        return cumul_probs;
        
}
row_vector ord_probs_to_cumul_probs(row_vector ord_probs) {
  
       return to_row_vector(ord_probs_to_cumul_probs(to_vector(ord_probs)));
        
}




////
//// Convert from cumul_probs -> C (w/ induced dirichlet):
////
vector ID_cumul_probs_to_C( vector cumul_probs, 
                            int use_probit_link) {
  
        int n_thr = num_elements(cumul_probs);
        int n_cat = n_thr + 1;
        vector[n_thr] C;
        
        for (k in 1:n_thr) {
              if (cumul_probs[k] < 1e-300) {
                    C[k] = -37.5; ////  prob = 1e-38;
              } else if (cumul_probs[k] > 0.99999999999999999) {
                    C[k] = +8.20; ////  prob = 0.9999999999999;
              } else {
                    if (use_probit_link == 1) C[k] = nico_inv_Phi(cumul_probs[k]); 
                    else                      C[k] = nico_logit(cumul_probs[k]); 
              }
        }
        
        return C;
        
}
row_vector ID_cumul_probs_to_C( row_vector cumul_probs, 
                                int use_probit_link) {
  
       return to_row_vector(ID_cumul_probs_to_C(to_vector(cumul_probs), use_probit_link));
        
}



////
//// C_array -> C_mat function:
////
matrix convert_C_array_to_mat(array[] vector C_array) {
  
  
      int n_rows = size(C_array);
      int n_cols = num_elements(C_array[1]);
      matrix[n_rows, n_cols] C_mat;
  
      for (k in 1:n_cols) { //// col-major storage
         for (s in 1:n_rows) {
                  C_mat[s, k] = C_array[s][k];
             }
      }
  
      return C_mat;
 
}












 
 
////
////  Get the cutpoint index (k) to map "latent_surv_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_prob_to_fixed_C_HSROC( vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      int n_thr = num_elements(C); 
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = cutpoint_index[c, s, cut_i + 1];
                            if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = -(C[k] + locations[s, c])*scales[s, c];
                  }
          }
      }
      
      return latent_surv;
                                    
}
////
//// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_prob_to_random_C_HSROC( matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                                  int k = cutpoint_index[c, s, cut_i + 1];
                                  if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = -(C[s, k] + locations[s, c])*scales[s, c];
                     
                  }
          }
      }
      
      return latent_surv;
                                    
}
















// 
// ////
// ////  Get the cutpoint index (k) to map "latent_surv_prob[c][s, cut_i]" to correct cutpoint "C[k]":
// ////
// array[] matrix map_latent_surv_to_fixed_homog_C( vector C, 
//                                   matrix locations, 
//                                   matrix scales, 
//                                   data int n_studies,
//                                   data array[] int n_obs_cutpoints,
//                                   data array[,,] int cutpoint_index) {
//           
//       int n_thr = num_elements(C); 
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                             int k = cutpoint_index[c, s, cut_i + 1];
//                             if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = - (C[k] - locations[s, c])/scales[s, c];
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }
// array[] matrix map_latent_surv_to_fixed_homog_C( vector C, 
//                                   matrix locations, 
//                                   real scale, 
//                                   data int n_studies,
//                                   data array[] int n_obs_cutpoints,
//                                   data array[,,] int cutpoint_index) {
//           
//       int n_thr = num_elements(C); 
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                             int k = cutpoint_index[c, s, cut_i + 1];
//                             if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = - (C[k] - locations[s, c])/scale;
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }


////
////  Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_to_fixed_hetero_C(  array[] vector C, 
                                                matrix locations, 
                                                matrix scales, 
                                                data int n_studies,
                                                data array[] int n_obs_cutpoints,
                                                data array[,,] int cutpoint_index) {
                                    
      int n_thr = num_elements(C[1]); 
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
    
      for (c in 1:2) {
          for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                               int k = cutpoint_index[c, s, cut_i + 1];
                               if (k != -1) { 
                                   latent_surv[c][s, cut_i] = - (C[c][k] - locations[s, c])/scales[s, c];
                              }
                      }
          }
      }
      
      return latent_surv;
                                    
}
array[] matrix map_latent_surv_to_fixed_hetero_C(  array[] vector C, 
                                                matrix locations, 
                                                data int n_studies,
                                                data array[] int n_obs_cutpoints,
                                                data array[,,] int cutpoint_index) {
          
      int n_thr = num_elements(C[1]); 
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
 
      for (c in 1:2) {
          for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              int k = cutpoint_index[c, s, cut_i + 1];
                              if (k != -1)  latent_surv[c][s, cut_i] =   locations[s, c] - C[c][k];
                      }
          }
      }
      
      return latent_surv;
                                    
}


// ////
// //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
// ////
// array[] matrix map_latent_surv_to_random_homog_C( matrix C, 
//                                                   data int n_thr,
//                                                   matrix locations, 
//                                                   matrix scales, 
//                                                   data int n_studies,
//                                                   data array[] int n_obs_cutpoints,
//                                                   data array[,,] int cutpoint_index) {
//           
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                                   int k = cutpoint_index[c, s, cut_i + 1];
//                                   if (k != -1) latent_surv[c][s, cut_i] =  (locations[s, c] - C[s, k])/scales[s, c];
//                      
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }
// array[] matrix map_latent_surv_to_random_homog_C( matrix C,  
//                                                   data int n_thr,
//                                                   matrix locations, 
//                                                   data int n_studies,
//                                                   data array[] int n_obs_cutpoints,
//                                                   data array[,,] int cutpoint_index) {
//           
//       int n_cat = n_thr + 1;
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                                   int k = cutpoint_index[c, s, cut_i + 1];
//                                   if (k != -1)  latent_surv[c][s, cut_i] =   locations[s, c] - C[s, k];
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }
// 



////
//// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());

      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                                 int k = cutpoint_index[c, s, cut_i + 1];
                                 if (k != -1)  latent_surv[c][s, cut_i] = (locations[s, c] - C[c][s, k])/scales[s, c];
                  }
          }
      }
      
      return latent_surv;
                                    
}
array[] matrix map_latent_surv_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                               int k = cutpoint_index[c, s, cut_i + 1];
                               if (k != -1)  latent_surv[c][s, cut_i] = locations[s, c] - C[c][s, k];
                        
                  }
          }
      }
      
      return latent_surv;
                                    
}

























////
////  Get the cutpoint index (k) to map "latent_surv_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_fixed_C( vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_prob_to_fixed_C( vector C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}


////
////  Get the cutpoint index (k) to map "latent_cumul[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_fixed_hetero_C( array[] vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
                                    
      int n_thr = num_elements(C[1]); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_prob_to_fixed_hetero_C( array[] vector C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C[1]); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}


////
//// ---- Get the cutpoint index (k) to map "latent_cumul[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_random_C( matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[s, k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_prob_to_random_C( matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C); 
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[s, k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}




////
//// ---- Get the cutpoint index (k) to map "latent_cumul[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C[1]); 
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());

      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][s, k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C[1]); 
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][s, k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}










// ////
// ////
// ////
// vector convert_Dirichlet_alpha_to_SDs_jacobian(vector alpha) {
// 
//       ////
//       //// Pre-compute some quantities needed for Jacobian for the (alpha_k, alpha_0) -> dirichlet_cat_SDs_sigma transformation:
//       ////
//       int n_cat = num_elements(alpha);
//       real log_half = -0.6931471805599453;
//       real alpha_0 = sum(alpha);
//       real alpha_0_sq = square(alpha_0);
//       real alpha_0_cube = alpha_0_sq * alpha_0;
//       vector[n_cat] alpha_sq = square(alpha);
//       ////
//       //// NOTE: SD for each category probability in a Dirichlet is sqrt(α_k(α_0-α_k)/(α_0^2(α_0+1))) - 
//       //// where α_0 is the sum of all alphas:
//       ////
//       vector[n_cat] dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) ./ (alpha_0_sq * (alpha_0 + 1.0))); 
//       ////
//       //// Then compute the Jacobian adjustment:
//       ////
//       vector[n_cat] log_sigma = log(dirichlet_cat_SDs_sigma);
//       vector[n_cat] deriv_A = alpha_sq + alpha_0 - 2.0*alpha;
//       vector[n_cat] deriv_B = (1.0 ./ alpha) .* (3.0*alpha_0_sq*alpha + 2.0*alpha*alpha_0);
//       vector[n_cat] deriv_var  = deriv_A .* (alpha_0_cube + alpha_0_sq) + deriv_B .* (alpha*alpha_0 - alpha_sq);
//       vector[n_cat] log_abs_deriv_var_wrt_alpha = log(abs(deriv_var));
//       vector[n_cat] log_abs_deriv_SD_wrt_alpha = log_half - log_sigma + log_abs_deriv_var_wrt_alpha;
//       ////
//       jacobian += sum(log_abs_deriv_SD_wrt_alpha);
//       ////
//       //// Output:
//       ////
//       return(dirichlet_cat_SDs_sigma);
// 
// }
// 


        

// for (k in 1:n_cat) {
//         // // real deriv_A[k] = alpha_sq[k] + alpha_0 - 2.0*alpha[k];
//         // // real deriv_B[k] = (1.0 / alpha[k]) * (3.0*alpha_0_sq*alpha[k] + 2.0*alpha[k]*alpha_0);
//         // // //// To compute on log-scale (need to keep track of signs / use special log_abs_sum_exp() function):
//         // // //// Term 1:
//         // // real deriv_term_1 = deriv_A[k] * (alpha_0_cube + alpha_0_sq);
//         // // real log_abs_deriv_term_1 = log(abs(deriv_A[k])) + log(abs(alpha_0_cube + alpha_0_sq));
//         // // real sign_deriv_term_1 = sign(deriv_A[k]) * sign(alpha_0_cube + alpha_0_sq);
//         // // //// Term 2:
//         // // real deriv_term_2 =  deriv_B[k] * (alpha[k]*alpha_0 - alpha_sq[k]);
//         // // real log_abs_deriv_term_2 = log(abs( deriv_B[k])) + log(abs(alpha[k]*alpha_0 - alpha_sq[k]));
//         // // real sign_deriv_term_2 = sign( deriv_B[k]) * sign(alpha[k]*alpha_0 - alpha_sq[k]);
//         // // //// Full term for deriv_var_k:
//         // // real log_abs_deriv_var_k_wrt_alpha = log_abs_sum_exp(log_abs_deriv_term_1, sign_deriv_term_1, log_abs_deriv_term_2, sign_deriv_term_2);
//         // //// To compute on normal scale:
//         // //// Since we have: SD_k = sqrt(Var_k), this means:
//         // //// deriv_SD_k_wrt_alpha = 0.5 * pow(Var_k, -0.5) * deriv_var_k_wrt_alpha. Hence:
//         // //// log_abs_deriv_SD_k_wrt_alpha = log(0.5) -0.5*log(Var_k) + log_abs_deriv_var_k_wrt_alpha. Hence:
//         // //// log_abs_deriv_SD_k_wrt_alpha = log(0.5) - log_sigma[k] + log_abs_deriv_var_k_wrt_alpha.
//         // // real deriv_var_k  = deriv_A[k] * (alpha_0_cube + alpha_0_sq) + deriv_B[k] * (alpha[k]*alpha_0 - alpha_sq[k]);
//         // // real log_abs_deriv_var_k_wrt_alpha = log(abs(deriv_var_k));
//         // ////
//         // real log_abs_deriv_SD_k_wrt_alpha = log_half - log_sigma[k] + log_abs_deriv_var_k_wrt_alpha;
//         // Jacobian_for_alpha_k_to_category_SDs += log_abs_deriv_SD_k_wrt_alpha;
// }
        






















////
//// ---------------------------------------------------------------------------------------
//// The induced-Dirichlet log density of a set of ordered cutpoints C (at the latent location 0) with every category
//// probability on the LOG scale (Stan's std_normal_lcdf / std_normal_lccdf with their log-scale tails, or
//// log_inv_logit / log1m_inv_logit, and log_diff_exp), no clamps: the same values as induced_dirichlet_given_C_lpdf
//// of the category probabilities of C wherever those are representable, and a finite log density (or -inf for two
//// equal cutpoints) where Phi(C_k) - Phi(C_(k-1)) rounds to 0 in double precision, which made dirichlet_lpdf throw.
//// The same calculation as NicoSIMD's nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf, for plain Stan math
//// models (so that the plain and the AVX versions of a model compute the same thing):
////     log p_1 = log F(C_1),  log p_(K+1) = log(1 - F(C_K)),
////     log p_k = log_diff_exp(log F(C_k), log F(C_(k-1)))                 (C_(k-1) < 0),
////             = log_diff_exp(log(1 - F(C_(k-1))), log(1 - F(C_k)))       (C_(k-1) >= 0).
////
real induced_dirichlet_given_C_on_the_log_scale_lpdf(vector C,
                                                     vector alpha,
                                                     int use_probit_link) {

        int n_thr = num_elements(C);
        vector[n_thr] log_F;
        vector[n_thr] log_S;
        if (use_probit_link == 1) {
          for (k in 1:n_thr) {
            log_F[k] = std_normal_lcdf(C[k] |);
            log_S[k] = std_normal_lccdf(C[k] |);
          }
        } else {
          log_F = nico_log_inv_logit(C);
          log_S = nico_log1m_inv_logit(C);
        }
        vector[n_thr + 1] log_p;
        log_p[1] = log_F[1];
        log_p[n_thr + 1] = log_S[n_thr];
        for (k in 2:n_thr) {
          log_p[k] = (C[k - 1] >= 0) ? nico_log_diff_exp(log_S[k - 1], log_S[k]) : nico_log_diff_exp(log_F[k], log_F[k - 1]);
        }
        real log_prob = lgamma(sum(alpha)) - sum(nico_lgamma(alpha));
        for (j in 1:(n_thr + 1)) {
          if (alpha[j] != 1.0) log_prob += (alpha[j] - 1.0) * log_p[j];
        }
        if (use_probit_link == 1) log_prob += nico_normal_lupdf(C | 0, 1);
        else                      log_prob += sum(log_F + log_S);
        return log_prob;

}
////
//// Row-vector overload:
////
real induced_dirichlet_given_C_on_the_log_scale_lpdf(row_vector C,
                                                     vector alpha,
                                                     int use_probit_link) {

        return induced_dirichlet_given_C_on_the_log_scale_lupdf(to_vector(C) | alpha, use_probit_link);

}
////
//// ---- Probability-scale Dirichlet or Beta, or link-scale Normal, at chosen cutoffs:
////
real Se_Sp_prior_at_chosen_cutoffs_lpdf(vector C_baseline_at_chosen_cutoffs,
                                        int Se_Sp_prior_type,
                                        vector prior_Se_Sp_dirichlet_A,
                                        matrix prior_Se_Sp_beta,
                                        matrix prior_Se_Sp_link_mean_SD,
                                        int use_probit_link) {

        real log_prob = 0.0;
        //// Selected baseline values follow the ordered curve.
        if (Se_Sp_prior_type == 1 || Se_Sp_prior_type == 2) {
             vector[num_elements(C_baseline_at_chosen_cutoffs)] P_at_chosen_cutoffs = (use_probit_link == 1) ?
                   nico_Phi(C_baseline_at_chosen_cutoffs) : nico_inv_logit(C_baseline_at_chosen_cutoffs);
             if (Se_Sp_prior_type == 1) {
                  //// Dirichlet on the segment sums:
                  log_prob += nico_dirichlet_lupdf(cumul_probs_to_ord_probs(P_at_chosen_cutoffs) | prior_Se_Sp_dirichlet_A);
             } else {
                  log_prob += nico_beta_lupdf(P_at_chosen_cutoffs | prior_Se_Sp_beta[, 1], prior_Se_Sp_beta[, 2]);
             }
             if (use_probit_link == 1) {
                  log_prob += nico_normal_lupdf(C_baseline_at_chosen_cutoffs | 0, 1);
             } else {
                  log_prob += sum(nico_log_inv_logit(C_baseline_at_chosen_cutoffs) +
                                  nico_log1m_inv_logit(C_baseline_at_chosen_cutoffs));
             }
        } else {
             log_prob += nico_normal_lupdf(C_baseline_at_chosen_cutoffs | prior_Se_Sp_link_mean_SD[, 1],
                                     prior_Se_Sp_link_mean_SD[, 2]);
        }
        return log_prob;

}
////
//// ---- Induced Dirichlet with the Dirichlet of the segment sums (between chosen cutoffs) swapped for
//// ---- a prior on Sp (c = 1) or 1 - Se (c = 2) at the chosen cutoffs, all at the constant anchor:
////
real induced_dirichlet_with_Se_Sp_prior_lpdf(vector C_shifted,
                                             vector alpha,
                                             vector alpha_of_the_segments,
                                             array[] int chosen_cutoffs,
                                             int Se_Sp_prior_type,
                                             vector prior_Se_Sp_dirichlet_A,
                                             matrix prior_Se_Sp_beta,
                                             matrix prior_Se_Sp_link_mean_SD,
                                             int use_probit_link) {

        int n_thr = num_elements(C_shifted);
        int n_chosen_cutoffs = num_elements(chosen_cutoffs);
        vector[n_thr] P = (use_probit_link == 1) ? nico_Phi(C_shifted) : nico_inv_logit(C_shifted);
        vector[n_thr + 1] p_ord = cumul_probs_to_ord_probs(P);
        vector[n_chosen_cutoffs] P_at_chosen_cutoffs = P[chosen_cutoffs];
        vector[n_chosen_cutoffs + 1] segment_probs = cumul_probs_to_ord_probs(P_at_chosen_cutoffs);
        ////
        //// ---- induced Dirichlet minus the Dirichlet of the segment sums:
        ////
        real log_prob = induced_dirichlet_given_C_lupdf(p_ord | C_shifted, alpha, use_probit_link);
        log_prob += - nico_dirichlet_lupdf(segment_probs | alpha_of_the_segments);
        ////
        //// ---- prior on P at the chosen cutoffs:
        ////
        if (Se_Sp_prior_type == 1) {
             //// Dirichlet on the segment sums:
             log_prob += nico_dirichlet_lupdf(segment_probs | prior_Se_Sp_dirichlet_A);
        } else if (Se_Sp_prior_type == 2) {
             //// Beta on each P, restricted to rising P (normaliser constant):
             log_prob += nico_beta_lupdf(P_at_chosen_cutoffs | prior_Se_Sp_beta[, 1], prior_Se_Sp_beta[, 2]);
        } else if (Se_Sp_prior_type == 4) {
             //// Normal on the link scale, restricted to rising P (with the Jacobian from P to the link):
             log_prob += nico_normal_lupdf(C_shifted[chosen_cutoffs] | prior_Se_Sp_link_mean_SD[, 1],
                                     prior_Se_Sp_link_mean_SD[, 2]);
             if (use_probit_link == 1) {
                  log_prob += - nico_normal_lupdf(C_shifted[chosen_cutoffs] | 0, 1);
             } else {
                  log_prob += - sum(nico_log_inv_logit(C_shifted[chosen_cutoffs]) +
                                    nico_log1m_inv_logit(C_shifted[chosen_cutoffs]));
             }
        } else {
             //// Beta on each (P_j - P_{j-1}) / (1 - P_{j-1}):
             vector[n_chosen_cutoffs] P_previous = append_row(0.0, P_at_chosen_cutoffs[1:(n_chosen_cutoffs - 1)]);
             log_prob += nico_beta_lupdf(segment_probs[1:n_chosen_cutoffs] ./ (1.0 - P_previous) |
                                   prior_Se_Sp_beta[, 1], prior_Se_Sp_beta[, 2]);
             log_prob += - sum(nico_log1m(P_previous));
        }
        return log_prob;

}
////
//// ---- The median of the Beta of a cutpoint under the Dirichlet, on the link scale
////      (shapes kappa * P_lower and kappa * P_upper; the smaller shape first, then reflected):
////
real dirichlet_cutpoint_median(real P_lower,
                               real P_upper,
                               real kappa,
                               int use_probit_link) {

        int reflect = (P_lower > P_upper);
        real shape_a = kappa * ((reflect == 1) ? P_upper : P_lower);
        real shape_b = kappa * ((reflect == 1) ? P_lower : P_upper);
        real p = inv_inc_beta(shape_a, shape_b, 0.5);
        real cutpoint = (use_probit_link == 1) ? nico_inv_Phi(p) : nico_logit(p);
        return (reflect == 1) ? -cutpoint : cutpoint;

}
////
//// ---- Induced Dirichlet with the Dirichlet of the segment sums (between chosen cutoffs) swapped for
//// ---- a prior on baseline Sp (c = 1) or 1 - Se (c = 2) at the chosen cutoffs:
////
real induced_dirichlet_with_Se_Sp_prior_lpdf(vector C_shifted,
                                             vector C_baseline,
                                             vector alpha,
                                             vector alpha_of_the_segments,
                                             array[] int chosen_cutoffs,
                                             int Se_Sp_prior_type,
                                             vector prior_Se_Sp_dirichlet_A,
                                             matrix prior_Se_Sp_beta,
                                             matrix prior_Se_Sp_link_mean_SD,
                                             int use_probit_link) {

        int n_thr = num_elements(C_shifted);
        int n_chosen_cutoffs = num_elements(chosen_cutoffs);
        vector[n_thr] P = (use_probit_link == 1) ? nico_Phi(C_shifted) : nico_inv_logit(C_shifted);
        vector[n_thr + 1] p_ord = cumul_probs_to_ord_probs(P);
        vector[n_chosen_cutoffs] P_at_chosen_cutoffs = (use_probit_link == 1) ? nico_Phi(C_baseline[chosen_cutoffs]) :
                                                    nico_inv_logit(C_baseline[chosen_cutoffs]);
        vector[n_chosen_cutoffs + 1] segment_probs = cumul_probs_to_ord_probs(P_at_chosen_cutoffs);
        vector[n_chosen_cutoffs + 1] segment_probs_at_the_anchor = cumul_probs_to_ord_probs(P[chosen_cutoffs]);
        ////
        //// ---- induced Dirichlet minus the Dirichlet of the segment sums:
        ////
        real log_prob = induced_dirichlet_given_C_lupdf(p_ord | C_shifted, alpha, use_probit_link);
        log_prob += - nico_dirichlet_lupdf(segment_probs_at_the_anchor | alpha_of_the_segments);
        ////
        //// ---- Jacobian between selected anchor and baseline probabilities:
        ////
        for (j in 1:n_chosen_cutoffs) {
             int k = chosen_cutoffs[j];
             if (use_probit_link == 1) {
                  log_prob += - nico_normal_lupdf(C_shifted[k] | 0, 1);
                  if (Se_Sp_prior_type != 4) {
                       log_prob += nico_normal_lupdf(C_baseline[k] | 0, 1);
                  }
             } else {
                  log_prob += - nico_log_inv_logit(C_shifted[k]) - nico_log1m_inv_logit(C_shifted[k]);
                  if (Se_Sp_prior_type != 4) {
                       log_prob += nico_log_inv_logit(C_baseline[k]) + nico_log1m_inv_logit(C_baseline[k]);
                  }
             }
        }
        ////
        //// ---- prior on P at the chosen cutoffs:
        ////
        if (Se_Sp_prior_type == 1) {
             //// Dirichlet on the segment sums:
             log_prob += nico_dirichlet_lupdf(segment_probs | prior_Se_Sp_dirichlet_A);
        } else if (Se_Sp_prior_type == 2) {
             //// Beta on each P, restricted to rising P (normaliser constant):
             log_prob += nico_beta_lupdf(P_at_chosen_cutoffs | prior_Se_Sp_beta[, 1], prior_Se_Sp_beta[, 2]);
        } else if (Se_Sp_prior_type == 4) {
             //// Normal on the link scale, restricted to rising P:
             log_prob += nico_normal_lupdf(C_baseline[chosen_cutoffs] | prior_Se_Sp_link_mean_SD[, 1],
                                     prior_Se_Sp_link_mean_SD[, 2]);
        } else {
             //// Beta on each (P_j - P_{j-1}) / (1 - P_{j-1}):
             vector[n_chosen_cutoffs] P_previous = append_row(0.0, P_at_chosen_cutoffs[1:(n_chosen_cutoffs - 1)]);
             log_prob += nico_beta_lupdf(segment_probs[1:n_chosen_cutoffs] ./ (1.0 - P_previous) |
                                   prior_Se_Sp_beta[, 1], prior_Se_Sp_beta[, 2]);
             log_prob += - sum(nico_log1m(P_previous));
        }
        return log_prob;

}
