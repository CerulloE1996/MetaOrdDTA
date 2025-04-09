

 
  

array[,] matrix compute_log_lik_binomial_probit_fact (  array[] matrix latent_surv, 
                                                        int use_probit_link,
                                                        data array[] matrix x, 
                                                        // data array[] matrix n,
                                                        data array[] int n_obs_cutpoints
) {
        
          int n_studies = rows(x[1]);
          int n_thr     = cols(x[1]) - 1;
          ////
          array[2] matrix[n_studies, n_thr] surv_prob;
          array[2] matrix[n_studies, n_thr] cond_prob;
          array[2] matrix[n_studies, n_thr] log_lik;
          ////
          for (c in 1:2) { 
              surv_prob[c]  = rep_matrix(1.0, n_studies, n_thr);
              cond_prob[c]  = rep_matrix(1.0, n_studies, n_thr);
              log_lik[c]    = rep_matrix(0.0, n_studies, n_thr);
          }
          //// ---- Calculate surv  probabilities:
          for (c in 1:2) {
              if (use_probit_link == 1) surv_prob[c] = Phi(latent_surv[c]); //// INCREASING sequence (as C_k > C_{k - 1})
              else                      surv_prob[c] = inv_logit(latent_surv[c]); //// INCREASING sequence (as C_k > C_{k - 1})
          }
          //// ---- Multinomial (factorised Binomial) likelihood:
          // matrix[n_studies, n_thr] probs;
          // matrix[n_studies, n_thr] bin_n;
          // matrix[n_studies, n_thr] bin_N_m_n;
          for (c in 1:2) {
              for (s in 1:n_studies) {
                      for (i in 1:n_obs_cutpoints[s]) {
                         
                              //// ---- Current and PREVIOUS survative counts //// next
                              int x_current = to_int(x[c][s, i + 1]);
                              int x_prev    = to_int(x[c][s, i + 0]); //// to_int(n[c][s, i]); // x_prev > x_next (counts are DECREASING)
                              int x_diff    = x_prev - x_current;
                
                              //// ---- Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                              
                                       //// ---- Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint
                                       if (i == 1) {
                                              //// cond_prob[c][s, 1] = surv_prob[c][s, 1]; // / 1.0; = surv_prob[1] /  surv_prob[0], and surv_prob[0] = 1.0.
                                              cond_prob[c][s, i] = surv_prob[c][s, i]; //// surv_idx - 2]; // Adjust back by 2
                                       } else if (i > 1) { //// so i = {3, 4, .... K}
                                              if (x_prev > 0) {
                                                  // For thresholds > 1, maintain the ratio calculation but with corrected indices
                                                  real curr_surv = surv_prob[c][s, i]; ////  (surv_idx <= n_thr) ? surv_prob[c][s, surv_idx - 2] : 0.0;
                                                  real prev_surv = surv_prob[c][s, i - 1] ; //// (prev_surv_idx <= n_thr) ? surv_prob[c][s, prev_surv_idx - 2] : 1.0;
                                                  cond_prob[c][s, i] = curr_surv / prev_surv;
                                              } else {
                                                  cond_prob[c][s, i] = 0.0;
                                              }
                                       }
                                       // probs[s, i] = cond_prob[c][s, i];
                                       // bin_n[s, i]     = x_current;
                                       // bin_N_m_n[s, i] = x_diff;
                                       ////
                                       log_lik[c][s, i] = binomial_lpmf(x_current | x_prev, cond_prob[c][s, i]);
                                       
                              }
                              
                     }

                  } // end of n_studies loop
                   // log_lik[c] = bin_n .* log(probs) +  bin_N_m_n .* log1m(probs);
       
          }
        
            array[2, 3] matrix[n_studies, n_thr] out;
            out[, 1] = log_lik;
            out[, 2] = cond_prob;
            out[, 3] = surv_prob;
            return(out);
           
}
    
    
    
    
           
           
           
           
           





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



           