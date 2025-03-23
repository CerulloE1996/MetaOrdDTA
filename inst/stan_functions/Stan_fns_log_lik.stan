

        
        
        


array[,] matrix compute_log_lik_binomial_fact (  array[] matrix cumul_prob,
                                                data array[] matrix x, 
                                                data array[] matrix n,
                                                data array[] int n_obs_cutpoints
    ) {
        
        int n_studies = rows(x[1]);
        int n_thr     = cols(x[1]);
        array[2] matrix[n_studies, n_thr] cond_prob; //// = init_array_of_matrices(n_studies, n_thr, 2, 1.0);
        array[2] matrix[n_studies, n_thr] log_lik; ////   = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        for (c in 1:2) { 
            cond_prob[c] = rep_matrix(1.0, n_studies, n_thr);
            log_lik[c]   = rep_matrix(0.0, n_studies, n_thr);
        }
        
        ////
        //// ---- Multinomial (factorised Binomial) likelihood:
        ////
        for (c in 1:2) {
            matrix[n_studies, n_thr] binomial_n;
            matrix[n_studies, n_thr] binomial_N;
            matrix[n_studies, n_thr] binomial_N_minus_n;
            for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              //// Current and next cumulative counts
                              int x_current = to_int(x[c][s, cut_i]);
                              int x_next    = to_int(n[c][s, cut_i]);
                              
                              //// Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                              
                                    //// Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint
                                    if (cut_i == n_obs_cutpoints[s]) { 
                                             cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / 1.0;
                                    } else {
                                          if (x_next > 0) { 
                                             cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / cumul_prob[c][s, cut_i + 1];
                                          } else { 
                                             cond_prob[c][s, cut_i] = 1.0;
                                          }
                                    }
                                    
                                    //// Binomial for observations at or below current cutpoint out of those at or below next cutpoint
                                    // log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c][s, cut_i]);
                                    binomial_n[s, cut_i] = x_current;
                                    binomial_N[s, cut_i] = x_next;
                                    binomial_N_minus_n[s, cut_i] = x_next - x_current;
                                    
                              }
                      }
                }

     
               log_lik[c] = binomial_n .* log(cond_prob[c]) + binomial_N_minus_n .* log1m(cond_prob[c]);
               log_lik[c] += lgamma(binomial_N + 1) - lgamma(binomial_n + 1) + lgamma(binomial_N_minus_n + 1);
        }
        
        
           for (c in 1:2) {
              for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              //// Current and next cumulative counts
                              int x_current = to_int(x[c][s, cut_i]);
                              int x_next    = to_int(n[c][s, cut_i]);
                              //// Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                                //// log_lik[c][s, cut_i] = 
                              } else { 
                                 log_lik[c][s, cut_i] = 0.0;
                              }
                      }
              }
           }
           
           
                    
           array[2, 2] matrix[n_studies, n_thr] out;
           out[1] = log_lik;
           out[2] = cond_prob;
           return(out);
           


}
           
           
           
           
           

array[,] matrix compute_log_lik_binomial_fact_NMA (  array[] matrix cumul_prob,
                                                    data array[] matrix x, 
                                                    data array[] matrix n,
                                                    data array[] int n_obs_cutpoints,
                                                    data array[] int indicator_index_test_t_in_study
    ) {
        
        int n_studies = rows(x[1]);
        int n_thr     = cols(x[1]);
        array[2] matrix[n_studies, n_thr] cond_prob; //// = init_array_of_matrices(n_studies, n_thr, 2, 1.0);
        array[2] matrix[n_studies, n_thr] log_lik; ////   = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        for (c in 1:2) { 
            cond_prob[c] = rep_matrix(1.0, n_studies, n_thr);
            log_lik[c]   = rep_matrix(0.0, n_studies, n_thr);
        }
        
        ////
        //// ---- Multinomial (factorised Binomial) likelihood:
        ////
        for (c in 1:2) {
            for (s in 1:n_studies) {
              
                  int observed = to_int(indicator_index_test_t_in_study[s]);
                  if (observed == 1)  {
                          for (cut_i in 1:n_obs_cutpoints[s]) {
                                  //// Current and next cumulative counts
                                  int x_current = to_int(x[c][s, cut_i]);
                                  int x_next    = to_int(n[c][s, cut_i]);
                                  
                                  //// Skip if the current count is zero (no observations to classify)
                                  if (x_current != 0)  {
                                  
                                         //// Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint:
                                         if (cut_i == n_obs_cutpoints[s]) { 
                                                 cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / 1.0;
                                         } else {
                                              if (x_next > 0) { 
                                                 cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / cumul_prob[c][s, cut_i + 1];
                                              } else { 
                                                 cond_prob[c][s, cut_i] = 1.0;
                                              }
                                         }
                                        
                                         //// Binomial for observations at or below current cutpoint out of those at or below next cutpoint
                                         // log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c][s, cut_i]);
                                         real binomial_n = x_current;
                                         real binomial_N = x_next;
                                         real binomial_N_minus_n = x_next - x_current;
                                         ////
                                         real log_cond_prob =    log(cond_prob[c][s, cut_i]);
                                         if (is_inf(log_cond_prob) == 0) { 
                                             log_lik[c][s, cut_i] = binomial_n * log_cond_prob;
                                         }
                                         real log_1m_cond_prob = log1m(cond_prob[c][s, cut_i]);
                                         if (is_inf(log_1m_cond_prob) == 0) { 
                                             log_lik[c][s, cut_i] += binomial_N_minus_n * log_1m_cond_prob;
                                         }
                                         real gamma_fn_stuff =   lgamma(binomial_N + 1);
                                         gamma_fn_stuff     += - lgamma(binomial_n + 1);
                                         gamma_fn_stuff     += + lgamma(binomial_N_minus_n + 1);
                                         if (is_inf(gamma_fn_stuff) == 0) { 
                                             log_lik[c][s, cut_i] += gamma_fn_stuff;
                                         }
                                         // // log_lik[c][s, cut_i] = binomial_n .* log(cond_prob[c][s, cut_i]);
                                         // // log_lik[c][s, cut_i] += binomial_N_minus_n .* log1m(cond_prob[c][s, cut_i]);
                                         // log_lik[c][s, cut_i] += lgamma(binomial_N + 1) - lgamma(binomial_n + 1) + lgamma(binomial_N_minus_n + 1);
                                  }
                          }
                    } 
                    
            }
     
               // log_lik[c] = binomial_n .* log(cond_prob[c]) + binomial_N_minus_n .* log1m(cond_prob[c]);
               // log_lik[c] += lgamma(binomial_N + 1) - lgamma(binomial_n + 1) + lgamma(binomial_N_minus_n + 1);
        }
        
          
           // for (c in 1:2) {
           //   for (s in 1:n_studies) {
           //     for (k in 1:n_thr[t]) {
           //         log_lik[c][s, k] = 0.0;
           //     }
           //   }
           // }
        
           for (c in 1:2) {
              for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              //// Current and next cumulative counts
                              int x_current = to_int(x[c][s, cut_i]);
                              int x_next    = to_int(n[c][s, cut_i]);
                              //// Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                                //// log_lik[c][s, cut_i] = 
                              } else { 
                                 log_lik[c][s, cut_i] = 0.0;
                              }
                      }
              }
           }
           
           array[2, 2] matrix[n_studies, n_thr] out;
           out[1] = log_lik;
           out[2] = cond_prob;
           return(out);
           


}





           