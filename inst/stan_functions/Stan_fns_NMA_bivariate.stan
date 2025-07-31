



////
//// ---- Function to compute DTA-NMA (Nyaga et al-based) summaries:
////
array[,] matrix compute_Nyaga_NMA_summaries(  int n_tests,
                                              matrix tau,
                                              matrix Sigma_bivariate,
                                              matrix Omega_bivariate) {
                    
            array[2] matrix[n_tests, n_tests] sigma_sq;
            array[2] matrix[n_tests, n_tests] rho;
            matrix[n_tests, n_tests] rho12;
            
            {
                array[n_tests] vector[2] sigma_b_sq;
                array[n_tests] vector[2] sigma_b;
                real rho_overall = Omega_bivariate[1, 2];
                matrix[n_tests, 2] tau_sq = square(tau);
                
                for (t in 1:n_tests) {
                    sigma_b_sq[t] = diagonal(Sigma_bivariate);
                    sigma_b[t] = sqrt(sigma_b_sq[t]);
                }
                
                //// Within-class correlations (Se with Se, Sp with Sp across tests)
                for (c in 1:2) {
                    for (t1 in 1:n_tests) {
                        for (t2 in 1:n_tests) {
                            //// Total variance for each test
                            real var_t1 = sigma_b_sq[t1][c] + tau_sq[t1, c];
                            real var_t2 = sigma_b_sq[t2][c] + tau_sq[t2, c];
                            
                            //// Correlation from shared eta scaled by test-specific sigmas
                            real covariance = sigma_b[t1][c] * sigma_b[t2][c];
                            
                            sigma_sq[c, t1, t2] = var_t1 * var_t2;
                            rho[c, t1, t2] = covariance / sqrt(sigma_sq[c, t1, t2]);
                        }
                    }
                }
                
                //// Between-class correlation (Se with Sp)
                for (t in 1:n_tests) {
                    real numerator = rho_overall * sigma_b[t][1] * sigma_b[t][2];
                    real denominator = sqrt((sigma_b_sq[t][1] + tau_sq[t, 1]) * 
                                            (sigma_b_sq[t][2] + tau_sq[t, 2]));
                    rho12[t, t] = numerator / denominator;
                }
            }
            
            array[3, 2] matrix[n_tests, n_tests] outs;
            outs[1, ] = sigma_sq;
            outs[2, ] = rho;
            outs[3, 1] = rho12;
            outs[3, 2] = rho12;
            return(outs);
    
    
}







////
//// ---- Function to compute DTA-NMA (Nyaga et al-based) summaries:
////
array[,] matrix compute_Nyaga_NMA_summaries_hetero_SDs(int n_tests,
                                                       matrix tau,
                                                       array[] matrix Sigma_bivariate,
                                                       matrix Omega_bivariate) {
    
            array[2] matrix[n_tests, n_tests] sigma_sq;
            array[2] matrix[n_tests, n_tests] rho;
            matrix[n_tests, n_tests] rho12;
            
            {
                array[n_tests] vector[2] sigma_b_sq;
                array[n_tests] vector[2] sigma_b;
                real rho_overall = Omega_bivariate[1, 2];
                matrix[n_tests, 2] tau_sq = square(tau);
                
                for (t in 1:n_tests) {
                    sigma_b_sq[t] = diagonal(Sigma_bivariate[t]);
                    sigma_b[t] = sqrt(sigma_b_sq[t]);
                }
                
                //// Within-class correlations (Se with Se, Sp with Sp across tests)
                for (c in 1:2) {
                    for (t1 in 1:n_tests) {
                        for (t2 in 1:n_tests) {
                            //// Total variance for each test
                            real var_t1 = sigma_b_sq[t1][c] + tau_sq[t1, c];
                            real var_t2 = sigma_b_sq[t2][c] + tau_sq[t2, c];
                            
                            //// Correlation from shared eta scaled by test-specific sigmas
                            real covariance = sigma_b[t1][c] * sigma_b[t2][c];
                            
                            sigma_sq[c, t1, t2] = var_t1 * var_t2;
                            rho[c, t1, t2] = covariance / sqrt(sigma_sq[c, t1, t2]);
                        }
                    }
                }
                
                //// Between-class correlation (Se with Sp)
                for (t in 1:n_tests) {
                    real numerator = rho_overall * sigma_b[t][1] * sigma_b[t][2];
                    real denominator = sqrt((sigma_b_sq[t][1] + tau_sq[t, 1]) * 
                                            (sigma_b_sq[t][2] + tau_sq[t, 2]));
                    rho12[t, t] = numerator / denominator;
                }
            }
            
            array[3, 2] matrix[n_tests, n_tests] outs;
            outs[1, ] = sigma_sq;
            outs[2, ] = rho;
            outs[3, 1] = rho12;
            outs[3, 2] = rho12;
            return(outs);
    
}



























