



////
//// ---- Function to compute DTA-NMA (Nyaga et al-based) summaries:
////
array[,] matrix compute_Nyaga_NMA_summaries(  int n_tests,
                                              matrix tau_sq,
                                              matrix Sigma_bivariate,
                                              matrix Omega_bivariate) {
 
            array[2] matrix[n_tests, n_tests] sigma_sq;
            array[2] matrix[n_tests, n_tests] rho;
            array[2] matrix[n_tests, n_tests] rho12;
            
            //// sigma_sq[c, t1, t2] - This calculates the product of the total variances for test t1 and test t2 for outcome c (i.e. Se and Sp). 
            //// Each total variance is the sum of the between-study variance (sigma_b_sq[c]) and the test-specific heterogeneity (tau_sq[t1, c] or tau_sq[t2, c]).
            ////
            //// rho[c, t1, t2] - This calculates the correlation induced by the common between-study variance component. 
            //// It's the ratio of the between-study variance to the geometric mean of the total variances. This represents how
            //// much of the total variance is explained by the common between-study heterogeneity.
            ////
            //// rho12[c, t1, t2] - This is calculating the correlation between tests t1 and t2 across different outcomes (1 and 2 - i.e. Se and Sp). 
            //// The numerator contains the overall correlation (rho_overall) multiplied by the between-study standard deviations, and the denominator
            //// is the geometric mean of the total variances.
            //// 
            //// The difference between rho and rho12:
            //// rho is the within-class correlation between tests (how tests correlate for the same class)
            //// rho12 is the between-class correlation (how Se of one test correlates with Se of another)
            
            {
              
              vector[2] sigma_b_sq = diagonal(Sigma_bivariate);
              vector[2] sigma_b = sqrt(sigma_b_sq);
              real rho_overall  = Omega_bivariate[1, 2];
                      
              for (c in 1:2) {
                  for (t1 in 1:n_tests) {
                      for (t2 in 1:n_tests) {
                          sigma_sq[c, t1, t2] = (sigma_b_sq[c] + tau_sq[t1, c]) * ((sigma_b_sq[c] + tau_sq[t2, c]));
                          rho[c, t1, t2]      = sigma_b_sq[c] / sqrt(sigma_sq[c, t1, t2]);
                          ////
                          //// rho12 is the correlation between the t1'th and t2'th test (t1=t2 and t1 =/=t2 both possible):
                          ////
                          real numerator = (rho_overall * sigma_b[1] * sigma_b[2]);
                          real denominator = sqrt( (sigma_b_sq[1] + tau_sq[t1, c]) * (sigma_b_sq[2] + tau_sq[t2, c]) );
                          rho12[c, t1, t2] = numerator / denominator;
                      }
                  }
              }
              
            }
          
            array[3, 2] matrix[n_tests, n_tests] outs;
            outs[1, ] = sigma_sq;
            outs[2, ] = rho;
            outs[3, ] = rho12;
            return(outs);
 
}

 


  
 

































