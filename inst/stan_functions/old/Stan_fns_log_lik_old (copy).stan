

real min_reals( real a, real b) { 
  vector[2] vec;
  vec[1] = a;
  vec[2] = b;
  return(min(vec));
}
  

real max_reals( real a, real b) { 
  vector[2] vec;
  vec[1] = a;
  vec[2] = b;
  return(max(vec));
}



real compute_log_lik_binomial_fact_lp(    array[] matrix latent_surv, 
                                          data int use_probit_link,
                                          data int n_thr,
                                          data int n_studies,
                                          data array[,,] int x_2, 
                                          data array[,,] int n, 
                                          data array[,] int N_total, 
                                          data array[] int n_obs_cutpoints
) {
          ////
          // vector[n_total_non_missing_obs] log_lik   = rep_vector(0.0, n_total_non_missing_obs);
          ////
          // array[2] matrix[n_studies, n_thr] surv_prob;
          // array[2] matrix[n_studies, n_thr] cond_prob;
          ////
          //// ---- Calculate surv  probabilities:
          ////
          ////
          //// ---- Factorised binomial Likelihood:
          ////
          real log_lik = 0.0;
          for (c in 1:2) {
                ////
                matrix[n_studies, n_thr] surv_prob;
                if (use_probit_link == 1) surv_prob = Phi(latent_surv[c]);
                else                      surv_prob = inv_logit(latent_surv[c]);
                ////
              
                ////
                //// ---- First compute conditional probs:
                ////
                for (s in 1:n_studies) {
                       vector[n_obs_cutpoints[s]] cond_prob;
                       ////
                       cond_prob[1] = surv_prob[s, 1];
                       ////
                       for (k in 2:n_obs_cutpoints[s]) {
                                  cond_prob[k] = surv_prob[s, k] / surv_prob[s, k - 1];
                       }
                       ////
                       //// ---- Likelihood:
                       ////
                       log_lik += binomial_lpmf(
                                  x_2[c, s, 1:n_obs_cutpoints[s]] | 
                                  n[c, s, 1:n_obs_cutpoints[s]], cond_prob);
                 }
          }
            
          return(log_lik);
           
}
array[,] matrix compute_log_lik_binomial_fact_data(   data array[] matrix latent_surv, 
                                                      data int use_probit_link,
                                                      data int n_thr,
                                                      data int n_studies,
                                                      data array[,,] int x_2, 
                                                      data array[,,] int n, 
                                                      data array[,] int N_total, 
                                                      data array[] int n_obs_cutpoints
) {
  
          ////
          array[2] matrix[n_studies, n_thr] surv_prob;
          array[2] matrix[n_studies, n_thr] cond_prob;
          array[2] matrix[n_studies, n_thr] log_lik;
          ////
          for (c in 1:2) { 
            surv_prob[c] = rep_matrix(1.0, n_studies, n_thr); 
            cond_prob[c] = rep_matrix(1.0, n_studies, n_thr);
            log_lik[c]   = rep_matrix(0.0, n_studies, n_thr);
          }
          ////
          //// ---- Calculate surv  probabilities:
          ////
          for (c in 1:2) {
              if (use_probit_link == 1) surv_prob[c] = Phi(latent_surv[c]);
              else                      surv_prob[c] = inv_logit(latent_surv[c]);
          }
          ////
          //// ---- Factorised binomial Likelihood:
          ////
          for (c in 1:2) {
                ////
                //// ---- First compute conditional probs:
                ////
                for (s in 1:n_studies) {
                     cond_prob[c][s, 1] = surv_prob[c][s, 1];
                     ////
                     for (k in 2:n_obs_cutpoints[s]) {
                                cond_prob[c][s, k] = surv_prob[c][s, k] / surv_prob[c][s, k - 1];
                     }
                }
                ////
                //// ---- Likelihood:
                ////
                for (s in 1:n_studies) {
                   for (k in 1:n_obs_cutpoints[s]) {
                      log_lik[c][s, k] = binomial_lpmf(x_2[c, s, k] | n[c, s, k], cond_prob[c][s, k]);
                   }
                }
          }
          
          array[3, 2] matrix[n_studies, n_thr] out;
          out[1] = log_lik;
          out[2] = cond_prob;
          out[3] = surv_prob;
          return(out);
          
}
    
    
    
    
    
real safe_binomial_lpmf( int x, 
                         int n,
                         real p) {
                           
      real p_safe = fmax(1e-10, fmin(1 - 1e-10, p));
      
      if (x == 0) {
        return n * log1m(p_safe);  // log(1-p)
      } else if (x == n) {
        return n * log(p_safe);
      } else {
        return binomial_lpmf(x | n, p_safe);
      }
  
}   
      
real compute_NMA_log_lik_binomial_fact_lp(    array[] matrix latent_surv, 
                                              data int use_probit_link,
                                              data int n_thr,
                                              data int n_studies,
                                              data array[,,] int x_2, 
                                              data array[,,] int n, 
                                              data array[,] int N_total, 
                                              data array[] int n_obs_cutpoints,
                                              data array[] int indicator_index_test_t_in_study,
                                              data array[] int K_fold_CV_indicator
) {
        
          ////
          array[2] matrix[n_studies, n_thr] surv_prob;
          array[2] matrix[n_studies, n_thr] cond_prob;
          ////
          //// ---- Calculate surv  probabilities:
          ////
          for (c in 1:2) {
              if (use_probit_link == 1) surv_prob[c] = Phi(latent_surv[c]);
              else                      surv_prob[c] = inv_logit(latent_surv[c]);
          }
          ////
          //// ---- Factorised binomial Likelihood:
          ////
          real log_lik = 0.0;
          for (c in 1:2) {
                ////
                //// ---- First compute conditional probs:
                ////
                for (s in 1:n_studies) {
                           cond_prob[c][s, 1] = surv_prob[c][s, 1];
                           ////
                           for (k in 2:n_obs_cutpoints[s]) {
                                      cond_prob[c][s, k] = surv_prob[c][s, k] / surv_prob[c][s, k - 1];
                           }
                }
                ////
                //// ---- Likelihood:
                ////
                for (s in 1:n_studies) {
                     // ////
                     // int dim = n_obs_cutpoints[s];
                     // array[dim] int x_2_int;
                     // array[dim] int n_int;
                     //  if (indicator_index_test_t_in_study[s] == 1) {
                     //     for (k in 1:n_obs_cutpoints[s]) {
                     //             x_2_int[k] = to_int(x_2[c][s, k]);
                     //             n_int[k]   = to_int(n[c][s, k]); //// [k]);
                     //             // if ((x_2_int[k] != -1) && (n_int[k] != -1)) {
                     //             //     if ( (!(is_nan(x_2_int[k]))) && (!(is_nan(n_int[k]))) &&
                     //             //          (!(is_nan(cond_prob[c][s, k]))) &&
                     //             //          (cond_prob[c][s, k] > 0.0)) {
                     //             //              real log_lik_i = binomial_lpmf(x_2_int[k] | n_int[k], 
                     //             //                               cond_prob[c][s, k]);
                     //             //           // if (!(is_inf(log_lik_i))) {
                     //             //              log_lik += log_lik_i;
                     //             //     }
                     //             // }
                     //         }
                     // }
                     if ((indicator_index_test_t_in_study[s] == 1) && (K_fold_CV_indicator[s] == 1)) {
                         // log_lik += binomial_lpmf(
                         //            x_2[c, s ,1:n_obs_cutpoints[s]] | 
                         //            n[c, s ,1:n_obs_cutpoints[s]], 
                         //            to_vector(cond_prob[c][s, 1:n_obs_cutpoints[s]]));
                                    for (k in 1:n_obs_cutpoints[s]) {
                                        log_lik += safe_binomial_lpmf( x_2[c, s, k] | n[c, s, k], cond_prob[c][s, k]);
                                    }
                     }
                 }
          }
          
          return(log_lik);
           
}
array[,] matrix compute_NMA_log_lik_binomial_fact_data(  array[] matrix latent_surv, 
                                                  data int use_probit_link,
                                                  data int n_thr,
                                                  data int n_studies,
                                                  data array[,,] int x_2, 
                                                  data array[,,] int n, 
                                                  data array[,] int N_total, 
                                                  data array[] int n_obs_cutpoints,
                                                  data array[] int indicator_index_test_t_in_study,
                                                  data array[] int K_fold_CV_indicator
) {
        
          ////
          array[2] matrix[n_studies, n_thr] surv_prob;
          array[2] matrix[n_studies, n_thr] cond_prob;
          array[2] matrix[n_studies, n_thr] log_lik;
          ////
          for (c in 1:2) { 
            surv_prob[c] = rep_matrix(1.0, n_studies, n_thr); 
            cond_prob[c] = rep_matrix(1.0, n_studies, n_thr);
            log_lik[c]   = rep_matrix(0.0, n_studies, n_thr);
          }
          ////;
          ////
          //// ---- Calculate surv  probabilities:
          ////
          for (c in 1:2) {
              if (use_probit_link == 1) surv_prob[c] = Phi(latent_surv[c]);
              else                      surv_prob[c] = inv_logit(latent_surv[c]);
          }
          ////
          //// ---- Factorised binomial Likelihood:
          ////
          for (c in 1:2) {
                ////
                //// ---- First compute conditional probs:
                ////
                for (s in 1:n_studies) {
                           cond_prob[c][s, 1] = surv_prob[c][s, 1];
                           ////
                           for (k in 2:n_obs_cutpoints[s]) {
                                      cond_prob[c][s, k] = surv_prob[c][s, k] / surv_prob[c][s, k - 1];
                           }
                }
                ////
                //// ---- Likelihood:
                ////
                for (s in 1:n_studies) {
                     // ////
                     // int dim = n_obs_cutpoints[s];
                     // array[dim] int x_2_int;
                     // array[dim] int n_int;
                     //  if (indicator_index_test_t_in_study[s] == 1) {
                     //     for (k in 1:n_obs_cutpoints[s]) {
                     //             x_2_int[k] = to_int(x_2[c][s, k]);
                     //             n_int[k]   = to_int(n[c][s, k]); //// [k]);
                     //             // if ((x_2_int[k] != -1) && (n_int[k] != -1)) {
                     //             //     if ( (!(is_nan(x_2_int[k]))) && (!(is_nan(n_int[k]))) &&
                     //             //          (!(is_nan(cond_prob[c][s, k]))) &&
                     //             //          (cond_prob[c][s, k] > 0.0)) {
                     //             //              real log_lik_i = binomial_lpmf(x_2_int[k] | n_int[k], 
                     //             //                               cond_prob[c][s, k]);
                     //             //           // if (!(is_inf(log_lik_i))) {
                     //             //              log_lik += log_lik_i;
                     //             //     }
                     //             // }
                     //         }
                     // }
                     if ((indicator_index_test_t_in_study[s] == 1) && (K_fold_CV_indicator[s] == 1)) {
                          for (k in 1:n_obs_cutpoints[s]) {
                            log_lik[c][s, k] = safe_binomial_lpmf(x_2[c, s, k] | n[c, s, k], cond_prob[c][s, k]);
                          }
                     }
                 }
          }
          
          array[3, 2] matrix[n_studies, n_thr] out;
          out[1] = log_lik;
          out[2] = cond_prob;
          out[3] = surv_prob;
          return(out);
           
}
    
    
        
    
         
           
//            
// real compute_log_lik_binomial_fact_NMA_lp(    array[] matrix latent_surv, 
//                                               data int use_probit_link,
//                                               data array[] matrix x_2, 
//                                               data array[] matrix n, 
//                                               data array[] vector N_total, 
//                                               data array[] int n_obs_cutpoints,
//                                               ata array[] int indicator_index_test_t_in_study
// ) {
//         
//           int n_studies = rows(x_2[1]);
//           int n_thr     = cols(x_2[1]);
//           ////
//           // vector[n_total_non_missing_obs] log_lik   = rep_vector(0.0, n_total_non_missing_obs);
//           ////
//           array[2] matrix[n_studies, n_thr] surv_prob;
//           array[2] matrix[n_studies, n_thr] cond_prob;
//           ////
//           //// ---- Calculate surv  probabilities:
//           ////
//           for (c in 1:2) {
//               if (use_probit_link == 1) surv_prob[c] = Phi(latent_surv[c]);
//               else                      surv_prob[c] = inv_logit(latent_surv[c]);
//           }
//           ////
//           //// ---- Factorised binomial Likelihood:
//           ////
//           real log_lik = 0.0;
//           for (c in 1:2) {
//                     ////
//                     //// ---- First compute conditional probs:
//                     ////
//                     for (s in 1:n_studies) {
//                                cond_prob[c][s, 1] = surv_prob[c][s, 1];
//                                ////
//                                for (k in 2:n_obs_cutpoints[s]) {
//                                           cond_prob[c][s, k] = surv_prob[c][s, k] / surv_prob[c][s, k - 1];
//                                }
//                     }
//                     ////
//                     //// ---- Likelihood:
//                     ////
//                     for (s in 1:n_studies) {
//                                ////
//                                int dim = n_obs_cutpoints[s];
//                                array[dim] int x_2_int;
//                                array[dim] int n_int;
//                                for (k in 1:n_obs_cutpoints[s]) {
//                                    x_2_int[k] = to_int(x_2[c][s, k]);
//                                    n_int[k]   = to_int(n[c][s, k]); //// [k]);
//                                }
//                                for (k in 1:n_obs_cutpoints[s]) {
//                                  if (indicator_index_test_t_in_study[s] == 1) {
//                                      log_lik += binomial_lpmf(x_2_int | n_int, to_vector(cond_prob[c][s, 1:n_obs_cutpoints[s]]));
//                                  }
//                                }
//                             
//       
//                      }
//                          
//           }
//             
//           return(log_lik);
//            
// }
//     
//     
//     
//     
//     
    




// array[,] matrix compute_log_lik_binomial_probit_fact_NMA (  array[] matrix latent_surv,
//                                                             data array[] matrix x, 
//                                                             data array[] matrix n,
//                                                             data array[] int n_obs_cutpoints,
//                                                             data array[] int indicator_index_test_t_in_study
//     ) {
//           
//             int n_studies = rows(x[1]);
//             int n_thr     = cols(x[1]);
//             array[2] matrix[n_studies, n_thr] surv_prob;
//             array[2] matrix[n_studies, n_thr] cond_prob;
//             ////
//             array[2] matrix[n_studies, n_thr] log_surv_prob;
//             array[2] matrix[n_studies, n_thr] log_cond_prob;
//             ////
//             array[2] matrix[n_studies, n_thr] log_lik;
//             for (c in 1:2) { 
//                 log_surv_prob[c] = rep_matrix(1.0, n_studies, n_thr);
//                 cond_prob[c]  = rep_matrix(1.0, n_studies, n_thr);
//                 log_lik[c]    = rep_matrix(0.0, n_studies, n_thr);
//             }
//             ////
//             //// ---- Calculate survATIVE probabilities:
//             ////
//             for (s in 1:n_studies) {
//                if (indicator_index_test_t_in_study[s] == 1) {
//                    for (c in 1:2) {
//                       log_surv_prob[c][s, 1:n_thr] = inv_logit(latent_surv[c][s, 1:n_thr]); //// INCREASING sequence (as C_k > C_{k - 1})
//                    }
//                }
//             }
//             ////
//             //// ---- Multinomial (factorised Binomial) likelihood:
//             ////
//             for (c in 1:2) {
//                 for (s in 1:n_studies) {
//                   
//                       int observed = to_int(indicator_index_test_t_in_study[s]);
//                       if (observed == 1)  {
//                               for (cut_i in 1:n_obs_cutpoints[s]) {
//                                       //// Current and next survative counts
//                                       int x_current = to_int(x[c][s, cut_i]);
//                                       int x_next    = to_int(n[c][s, cut_i]);
//                                       
//                                       //// Skip if the current count is zero (no observations to classify)
//                                       if (x_current != 0)  {
//                                       
//                                              //// Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint:
//                                              if (cut_i == n_obs_cutpoints[s]) { 
//                                                    vector[2] probs;
//                                                    probs[1] = 0.999999999999999;
//                                                    probs[2] = surv_prob[c][s, cut_i] / 1.0;
//                                                    cond_prob[c][s, cut_i] = min(probs);
//                                              } else {
//                                                   if (x_next > 0) { 
//                                                      cond_prob[c][s, cut_i] = surv_prob[c][s, cut_i] / surv_prob[c][s, cut_i + 1];
//                                                   } else { 
//                                                      cond_prob[c][s, cut_i] = 1.0;
//                                                   }
//                                              }
//                                             
//                                              //// Binomial for observations at or below current cutpoint out of those at or below next cutpoint
//                                              // log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c][s, cut_i]);
//                                              real binomial_n = x_current;
//                                              real binomial_N = x_next;
//                                              real binomial_N_minus_n = x_next - x_current;
//                                              ////
//                                              real log_cond_prob =    log(cond_prob[c][s, cut_i]);
//                                              if (is_inf(log_cond_prob) == 0) { 
//                                                  log_lik[c][s, cut_i] = binomial_n * log_cond_prob;
//                                              }
//                                              real log_1m_cond_prob = log1m(cond_prob[c][s, cut_i]);
//                                              if (is_inf(log_1m_cond_prob) == 0) { 
//                                                  log_lik[c][s, cut_i] += binomial_N_minus_n * log_1m_cond_prob;
//                                              }
//                                       } else { 
//                                             log_lik[c][s, cut_i] = 0.0;
//                                       }
//                               }
//                         } 
//                         
//                 }
//             }
//            
//             array[2, 3] matrix[n_studies, n_thr] out;
//             out[, 1] = log_lik;
//             out[, 2] = cond_prob;
//             out[, 3] = surv_prob;
//             return(out);
// 
// }
// 
// 



           