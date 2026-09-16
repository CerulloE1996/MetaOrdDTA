

 
  

array[,] matrix compute_log_lik_binomial_probit_fact (  array[] matrix latent_surv, 
                                                        int use_probit_link,
                                                        data array[] matrix x, 
                                                        // data array[] matrix n,
                                                        data array[] int n_obs_cutpoints
) {
        
          int n_studies = rows(x[1]);
          int n_thr     = cols(x[1]) - 1;
          // int n_cat = n_thr + 1;
          ////
          array[2] matrix[n_studies, n_thr] surv_prob;
          array[2] matrix[n_studies, n_thr] cond_prob;
          ////
          // array[2] matrix[n_studies, n_thr] log_surv_prob;
          // array[2] matrix[n_studies, n_thr] log_cond_prob;
          array[2] matrix[n_studies, n_thr] log_lik;
          for (c in 1:2) { 
              surv_prob[c]  = rep_matrix(1.0, n_studies, n_thr);
              cond_prob[c]  = rep_matrix(1.0, n_studies, n_thr);
              ////
              // log_surv_prob[c]  = rep_matrix(0.0, n_studies, n_thr);
              // log_cond_prob[c]  = rep_matrix(0.0, n_studies, n_thr);
              log_lik[c]        = rep_matrix(0.0, n_studies, n_thr);
          }
          ////
          //// ---- Calculate surv  probabilities:
          ////
          for (c in 1:2) {
              if (use_probit_link == 1) surv_prob[c][, ] = (Phi(latent_surv[c][, ])); //// INCREASING sequence (as C_k > C_{k - 1})
              else                      surv_prob[c][, ] = inv_logit(latent_surv[c][, ]); //// INCREASING sequence (as C_k > C_{k - 1})
          }
          // ////
          // //// ---- Calculate surv  probabilities:
          // ////
          // for (c in 1:2) {
          //     if (use_probit_link == 1) log_surv_prob[c][, ] = log(Phi(latent_surv[c][, ])); //// INCREASING sequence (as C_k > C_{k - 1})
          //     else                      log_surv_prob[c][, ] = log_inv_logit(latent_surv[c][, ]); //// INCREASING sequence (as C_k > C_{k - 1})
          // }
          ////
          //// ---- Multinomial (factorised Binomial) likelihood:
          ////
          for (c in 1:2) {
              for (s in 1:n_studies) {
                      {
                            //// Current and PREVIOUS survative counts //// next
                              int cut_i = 1;
                              int x_current = to_int(x[c][s, 2]);
                              int x_prev    = to_int(x[c][s, 1]); //// to_int(n[c][s, cut_i]); // x_prev > x_next (counts are DECREASING)
                              // int x_diff = x_prev - x_current;
                              ////
                              cond_prob[c][s, 1] = surv_prob[c][s, 1] - 0.0;
                              real prob = cond_prob[c][s, 1];
                              ////
                              log_lik[c][s, 1] = binomial_lpmf(x_current | x_prev, prob);
                              //// log_lik[c][s, 1] =  x_current * log_prob  +  x_diff * log1m_exp(log_prob);

                      }
                      for (cut_i in 2:(n_obs_cutpoints[s])) {
                        
                              //// Current and PREVIOUS survative counts //// next
                              int x_current = to_int(x[c][s, cut_i + 1]);
                              int x_prev    = to_int(x[c][s, cut_i + (-1)]); //// to_int(n[c][s, cut_i]); // x_prev > x_next (counts are DECREASING)
                              // int x_diff = x_prev - x_current;
                              
                              //// Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                              
                                       //// Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint
                                       if (x_prev > 0)   cond_prob[c][s, cut_i] = surv_prob[c][s, cut_i] / surv_prob[c][s, cut_i - 1];  
                                       ////  surv_prob[c][s, cut_i] / surv_prob[c][s, cut_i + 1];
                                       else  cond_prob[c][s, cut_i] = 0.0; //// 1.0;
                                       ////
                                       real prob = cond_prob[c][s, cut_i];
                                       ////
                                       // log_lik[c][s, cut_i]  = lchoose(binom_N, binom_n);  //// + log( p_pow_n * one_m_p_pow_N_m_n);
                                       // log_lik[c][s, cut_i] += binom_n * log_prob + binom_N_m_n * log1m_exp(log_prob); //  + binom_N_m_n * log1m(current_prob); 
                                       // log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_prev, exp(log_prob));
                                       log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_prev, prob);
                                       //// log_lik[c][s, 1] =  x_current * log_prob  +  x_diff * log1m_exp(log_prob);
                                       
                              }
                              
                     }

                  }
       
                 // log_lik[c] = binom_n .* log(cond_prob[c]) + binom_N_m_n .* log1m(cond_prob[c]);
                 // // log_lik[c] += lgamma(binom_N + 1) - lgamma(binom_n + 1) + lgamma(binom_N_m_n + 1);
          }
        
            array[2, 3] matrix[n_studies, n_thr] out;
            out[, 1] = log_lik;
            out[, 2] = cond_prob;
            out[, 3] = surv_prob;
            return(out);
           
}
    
    
    
    
           
           
           
           
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



           