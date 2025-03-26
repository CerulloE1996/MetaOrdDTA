



array[,] matrix compute_RnG_NMA_summaries(
                                        int n_tests,
                                        real mu_beta,
                                        real mu_gamma,
                                        real sigma_beta,
                                        real sigma_gamma) {
    
            array[2] matrix[n_tests, n_tests] sigma_sq;
            array[2] matrix[n_tests, n_tests] rho;
            array[2] matrix[n_tests, n_tests] rho12;
            
            //// Function to compute f(x)
            real f(real x) {
                return exp(x); // Or use softplus if preferred
            }
            
            //// Function to compute f'(x)/f(x)
            real f_prime_over_f(real x) {
                return 1.0; // For exp function, derivative/original = 1
                // If using softplus, this would be exp(x)/(1+exp(x))
            }
            
            //// Compute variances for Se:
            real var_se = (0.25 * square(sigma_beta)) / square(f(0.5 * mu_gamma)) + (0.25 * square(sigma_gamma) * square(f_prime_over_f(0.5 * mu_gamma)));
            
            //// Compute variances for Sp:
            real var_sp = (0.25 * square(sigma_beta)) / square(f(-0.5 * mu_gamma)) + (0.25 * square(sigma_gamma) * square(f_prime_over_f(-0.5 * mu_gamma)));
            
            //// The implicit covariance between Se and Sp due to shared parameters:
            real cov_se_sp = -0.25 * square(sigma_beta) / (f(0.5 * mu_gamma) * f(-0.5 * mu_gamma));
            
            //// Fill the matrices:
            for (c in 1:2) {
                real var_c = (c == 1) ? var_se : var_sp;
                
                for (t1 in 1:n_tests) {
                    for (t2 in 1:n_tests) {
                        
                          //// For variances:
                          sigma_sq[c, t1, t2] = var_c * var_c;
                          
                          //// For within-outcome correlations (perfect for same test, 0 for different tests):
                          rho[c, t1, t2] = (t1 == t2) ? 1.0 : 0.0;
                          
                          //// For between-outcome correlations:
                          real implicit_corr = cov_se_sp / sqrt(var_se * var_sp); 
                          rho12[c, t1, t2] = (t1 == t2) ? implicit_corr : 0.0;
                        
                    }
                }
            }
            
            array[3, 2] matrix[n_tests, n_tests] outs;
            outs[1, ] = sigma_sq;
            outs[2, ] = rho;
            outs[3, ] = rho12;
            return(outs);
    
}



  
  
 



























