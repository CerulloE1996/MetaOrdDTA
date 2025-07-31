    
    
    
    
    
    
array[] matrix compute_deviance(array[] matrix cond_prob,
                                 data int n_thr,
                                 data int n_studies,
                                 data array[,,] int x_2,
                                 data array[,,] int n,
                                 data array[] int n_obs_cutpoints) {
    
    matrix[n_studies, n_thr] x_hat_nd = rep_matrix(0.0, n_studies, n_thr);
    matrix[n_studies, n_thr] x_hat_d  = rep_matrix(0.0, n_studies, n_thr);
    matrix[n_studies, n_thr] dev_nd   = rep_matrix(0.0, n_studies, n_thr);
    matrix[n_studies, n_thr] dev_d    = rep_matrix(0.0, n_studies, n_thr);
    
    {
        array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        
        for (s in 1:n_studies) {
            for (c in 1:2) {
                for (i in 1:n_obs_cutpoints[s]) {
                    // Model-estimated data
                    x_hat[c][s, i] = cond_prob[c][s, i] * n[c, s, i];
                    
                    // Get data values
                    real n_i = n[c, s, i];
                    real x_i = x_2[c, s, i];
                    real x_hat_i = x_hat[c][s, i];
                    
                    // Skip if this is missing data
                    if (n_i == -1 || x_i == -1) {
                        dev[c][s, i] = 0.0;
                        continue;
                    }
                    
                    // Ensure x_hat_i is within valid bounds
                    x_hat_i = fmax(1e-10, fmin(n_i - 1e-10, x_hat_i));
                    
                    // Compute deviance components with safe guards
                    real dev_term1 = 0.0;
                    real dev_term2 = 0.0;
                    
                    // First term: x_i * (log(x_i) - log(x_hat_i))
                    if (x_i > 0) {
                        dev_term1 = x_i * (log(x_i) - log(x_hat_i));
                    } else {
                        // When x_i = 0, the limit of x*log(x) as x->0 is 0
                        // But we still have -x_i * log(x_hat_i) = 0
                        dev_term1 = 0.0;
                    }
                    
                    // Second term: (n_i - x_i) * (log(n_i - x_i) - log(n_i - x_hat_i))
                    real n_minus_x = n_i - x_i;
                    real n_minus_x_hat = n_i - x_hat_i;
                    
                    if (n_minus_x > 0) {
                        dev_term2 = n_minus_x * (log(n_minus_x) - log(n_minus_x_hat));
                    } else {
                        // When n_i = x_i (all successes), this term is 0
                        dev_term2 = 0.0;
                    }
                    
                    // Total deviance
                    dev[c][s, i] = 2.0 * (dev_term1 + dev_term2);
                    
                    // Additional check for NaN or Inf
                    if (is_nan(dev[c][s, i]) || is_inf(dev[c][s, i])) {
                        dev[c][s, i] = 0.0;
                    }
                }
            }
        }
        
        x_hat_nd = x_hat[1];
        dev_nd = dev[1];
        x_hat_d = x_hat[2];
        dev_d = dev[2];
    }
    
    array[4] matrix[n_studies, n_thr] outs;
    outs[1] = x_hat_nd;
    outs[2] = dev_nd;
    outs[3] = x_hat_d;
    outs[4] = dev_d;
    
    return(outs);
}