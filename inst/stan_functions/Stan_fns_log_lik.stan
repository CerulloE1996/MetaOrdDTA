

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


real compute_log_lik_binomial_fact_lp(    array[] matrix latent_surv, 
                                          data int use_probit_link,
                                          data int n_thr,
                                          data int n_studies,
                                          data array[,,] int x_2, 
                                          data array[,,] int n, 
                                          data array[,] int N_total, 
                                          data array[] int n_obs_cutpoints,
                                          data array[] int K_fold_CV_indicator
) {
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
                   if (K_fold_CV_indicator[s] == 1) {
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
                       for (k in 1:n_obs_cutpoints[s]) {
                          log_lik += safe_binomial_lpmf(x_2[c, s, k] | n[c, s, k], cond_prob[k]);
                       }
                   }
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
                                                      data array[] int n_obs_cutpoints,
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
                   if (K_fold_CV_indicator[s] == 1) {
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
                   if (indicator_index_test_t_in_study[s] == 1) {
                           cond_prob[c][s, 1] = surv_prob[c][s, 1];
                           ////
                           for (k in 2:n_obs_cutpoints[s]) {
                                      cond_prob[c][s, k] = surv_prob[c][s, k] / surv_prob[c][s, k - 1];
                           }
                   }
                }
                ////
                //// ---- Likelihood:
                ////
                for (s in 1:n_studies) {
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
    
    
        
    
          
          
          
          
