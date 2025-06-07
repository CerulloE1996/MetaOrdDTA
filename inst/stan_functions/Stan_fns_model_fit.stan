    
    
    
    
    
    
    
array[] matrix compute_deviance(    array[] matrix cond_prob, 
                                    data int n_thr,
                                    data array[] matrix x_2, 
                                    data array[] matrix n,
                                    data array[] int n_obs_cutpoints
) {
    
          int n_studies = rows(x_2[1]);
          // int n_thr = cols(x_2[1]);
          ////
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          {
              array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
              array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
              ////
              //// ---- Model-predicted ("re-constructed") data:
              ////
              for (s in 1:n_studies) {
                  for (c in 1:2) {
                     for (i in 1:to_int(n_obs_cutpoints[s])) {
                           // dev[i,j,t] <- 2*(x[i,j,t]*(log(x[i,j,t]) - log(xhat[i,j,t])) + (n[i,j,t] - x[i,j,t])*(log(n[i,j,t] - x[i,j,t])     
                           //                             - log(n[i,j,t] - xhat[i,j,t])))  # Residual deviance contribution
    
                          //// 
                          //// ---- Model-estimated data:
                          ////
                          x_hat[c][s, i] = cond_prob[c][s, i] * n[c][s, i];  	 // Fitted values
                          ////
                          //// ---- Compute residual deviance contribution:
                          ////
                          real n_i =  n[c][s, i];
                          real x_i =  x_2[c][s, i];
                          ////
                          real x_hat_i =  x_hat[c][s, i];
                          ////
                          real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                          real log_diff_n_minus_x = log(n_i - x_i);
                          real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
                          ////
                          dev[c][s, i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) );
    
                     }
                  }
              }
    
              x_hat_nd = x_hat[1];
              dev_nd = dev[1];
              x_hat_d = x_hat[2];
              dev_d =  dev[2];
          }
          
          array[4] matrix[n_studies, n_thr] outs;
          outs[1] = x_hat_nd;
          outs[2] = dev_nd;
          outs[3] = x_hat_d;
          outs[4] = dev_d;
          return(outs);
          
}
