


////
//// Flexible Correlation matrix functions:: -----------------------------------------------------------------------------------
//// NOTE: adapted from Pinkney et al, 2024  -----------------------------------------------------------------------------------
////

      //// Bounds y between lb and ub (WITHOUT jacobian adjustment!):
      real lb_ub ( real y, 
                   real lb, 
                   real ub) {
          
            real neg_log_2 = -0.6931471805599453; // 64-bit (i.e. "double") precision
            real tanh_y = tanh(y);
            
            return lb +  (ub - lb) *  0.5 * (1.0 + tanh_y);
    
      }
      
      //// Overload "lb_ub" for vector input:
      vector lb_ub (  vector y, 
                      real lb, 
                      real ub) {
 
            real neg_log_2 = -0.6931471805599453; // 64-bit (i.e. "double") precision
            int N = num_elements(y); 
            vector[N] tanh_y = tanh(y);
            
            return lb +  (ub - lb) *  0.5 * (1.0 + tanh_y);
    
      }
            
      //// Overload "lb_ub" for vector input and vector of lower/upper bounds:
      vector lb_ub ( vector y, 
                     vector lb, 
                     vector ub) {
 
            real neg_log_2 = -0.6931471805599453; // 64-bit (i.e. "double") precision
            int N = num_elements(y); 
            vector[N] tanh_y = tanh(y);
            
            return lb + (ub - lb) .* (0.5 * (1.0 + tanh_y));
    
      }
      
      // ////
      // //// Same as "lb_ub" but with Jacobian adjustment increment ("jacobian +="):
      // ////
      // real lb_ub_jacobian (real y, 
      //                      real lb, 
      //                      real ub) {
      //     
      //       real neg_log_2 = -0.6931471805599453; // 64-bit (i.e. "double") precision
      //       real tanh_y = tanh(y);
      //       
      //       jacobian += neg_log_2  +  log( (ub - lb) * (1.0 - square(tanh_y)) ) ;
      //       return lb +  (ub - lb) *  0.5 * (1.0 + tanh_y);
      // 
      // }
      // 
      // //// Overload "lb_ub_jacobian" for vector input:
      // vector lb_ub_jacobian ( vector y, 
      //                         real lb, 
      //                         real ub) {
      // 
      //       real neg_log_2 = -0.6931471805599453; // 64-bit (i.e. "double") precision
      //       int N = num_elements(y); 
      //       vector[N] tanh_y = tanh(y);
      //       
      //       jacobian += neg_log_2  +  log( (ub - lb) * (1.0 - square(tanh_y)) ) ;
      //       return lb +  (ub - lb) *  0.5 * (1.0 + tanh_y);
      // 
      // }
      // 
      // //// Overload "lb_ub_jacobian" for vector input and vector of lower/upper bounds:
      // vector lb_ub_jacobian ( vector y, 
      //                         vector lb, 
      //                         vector ub) {
      // 
      //       real neg_log_2 = -0.6931471805599453; // 64-bit (i.e. "double") precision
      //       int N = num_elements(y); 
      //       vector[N] tanh_y = tanh(y);
      //       
      //       jacobian += sum(neg_log_2 + log( (ub - lb) .* (1.0 - square(tanh_y))));
      //       return lb + (ub - lb) .* (0.5 * (1.0 + tanh_y));
      // 
      // }
      // 
      
      

      // matrix corr_to_chol( real x, 
      //                      int J) {
      //                            
      //        matrix[J, J] Omega = add_diag(rep_matrix(x, J, J), 1 - x);
      //        return cholesky_decompose(Omega);
      //    
      // }
      
      
      ////
      //// Function to construct Omega (bivariate) from single correlation parameter:
      ////
      matrix make_bivariate_Omega(real corr) { 
        
              matrix[2, 2] Omega = diag_matrix(rep_vector(1.0, 2));
              Omega[1, 2] = corr;
              Omega[2, 1] = corr;
              return Omega;
        
      }
      
      
      ////
      //// Function to construct Omega (bivariate) from single correlation parameter:
      ////
      matrix make_bivariate_L_Omega(real corr) { 
  
              matrix[2, 2] L_Omega = rep_matrix(0.0, 2, 2);
              L_Omega[1, 1] = 1.0;
              L_Omega[2, 1] = corr;
              L_Omega[2, 2] = sqrt(1.0 - square(corr));
              return L_Omega;
        
      }


            
      // matrix make_restricted_bivariate_Omega_jacobian(  real corr, 
      //                                                   real lb_corr, 
      //                                                   real ub_corr) {
      //     
      //        real corr_restricted = lb_ub_jacobian(corr, lb_corr, ub_corr); //// this also increments log_det_J
      //        ////
      //        matrix[2, 2] Omega = diag_matrix(rep_vector(1.0, 2));
      //        Omega[1, 2] = corr_restricted;
      //        Omega[2, 1] = Omega[1, 2];
      //        ////
      //        return Omega;
      //                                         
      // }
            
     //  matrix make_restricted_bivariate_L_Omega_jacobian( real corr, 
     //                                                     real lb_corr, 
     //                                                     real ub_corr) {
     //   
     //        real corr_restricted = lb_ub_jacobian(corr, lb_corr, ub_corr); //// this also increments log_det_J
     //        ////
     //        matrix[2, 2] L_Omega = rep_matrix(0.0, 2, 2);
     //        L_Omega[1, 1] = 1.0;
     //        L_Omega[2, 1] = corr_restricted;
     //        L_Omega[2, 2] = sqrt(1.0 - square(corr_restricted));
     //        ////
     //        jacobian += lkj_corr_cholesky_lpdf(L_Omega | 1.0); //// Jacobian for transaformation L_Omega -> Omega
     //        return L_Omega;
     //                                          
     //  }
     // 
     //  // real lb_trans = 2.0*(lb + 1.0); //// Transform bounds to [0,1]
     //  // real ub_trans = 2.0*(ub + 1.0); //// Transform bounds to [0,1]
     //  // real mass_in_range = beta_cdf(ub_trans, eta, eta) - beta_cdf(lb_trans, eta, eta);
     //  // target += (eta - 1.0) * log1m(square(rho)) - log(mass_in_range);
     // 
     //  ////
     //  //// This Stan function was provided by Sean Pinkney (may be slightly adapted/different formatting):
     //  ////
     //  matrix Pinkney_LDL_bounds_opt_jacobian(  vector Omega_raw_vec,
     //                                           matrix lb_corr,
     //                                           matrix ub_corr,
     //                                           matrix known_values_indicator,
     //                                           matrix known_values) {
     //                                             
     //        ////
     //        //// First construct "Omega_raw_mat" matrix from "Omega_raw_vec" using row-by-row access (as MUST match the BayesMVP C++ functions/models!):
     //        ////
     //        int dim = rows(lb_corr);
     //        matrix[dim, dim] Omega_raw_mat = rep_matrix(0.0, dim, dim);
     //        int counter = 1;
     //        for (i in 2:dim) { 
     //          for (j in 1:(i - 1)) { //// j < i
     //              Omega_raw_mat[i, j] = Omega_raw_vec[counter];
     //              counter += 1;
     //          }
     //        }
     //        ////
     //        //// Then get 1st raw col excluding first element (0.0):
     //        ////
     //        vector[dim - 1] col_one_raw = segment(col(Omega_raw_mat, 1), 2, dim - 1); //// Extract first col. excluding the first element (0.0):
     //        ////
     //        vector[dim - 1] z = lb_ub_jacobian(col_one_raw, segment(col(lb_corr, 1), 2, dim - 1), segment(col(lb_corr, 1), 2, dim - 1));
     //        ////
     //        matrix[dim, dim] L = diag_matrix(rep_vector(1.0, dim));
     //        vector[dim] D;
     //        D[1] = 1;
     //        L[2:dim, 1] = z[1:dim - 1]; 
     //        D[2] = 1 - L[2, 1]^2;
     //        ////
     //        //// Fix known values for first col (if any known):
     //        ////
     //        for (i in 2:dim) { ////  skip first element since this is fixed to 1.0 (as corr matrix).
     //            if (known_values_indicator[i, 1] == 1)    L[i, 0] = known_values[i, 0];
     //        }
     //        ////
     //        //// Compute L_Omega:
     //        ////
     //        for (i in 3:dim) {
     //          
     //               D[i] = 1 - L[i, 1]^2; 
     //               L[i, 2:i - 1] = rep_row_vector(1 - L[i, 1]^2, i - 2);
     //               real L_ij_old = L[i, 2];
     //               
     //               for (j in 2:(i - 1))  {
     //                    
     //                          //  real L_ij_old = L[i, j];
     //                          real b1 = dot_product(L[j, 1:(j - 1)], D[1:j - 1]' .* L[i, 1:(j - 1)]);
     //                          
     //                          // how to derive the bounds
     //                          // we know that the correlation value C is bound by
     //                          // b1 - Ljj * Lij_old <= C <= b1 + Ljj * Lij_old
     //                          // Now we want our bounds to be enforced too so
     //                          // max(lb_corr, b1 - Ljj * Lij_old) <= C <= min(ub_corr, b1 + Ljj * Lij_old)
     //                          // We have the Lij_new = (C - b1) / Ljj
     //                          // To get the bounds on Lij_new is
     //                          // (bound - b1) / Ljj 
     //                          
     //                          if (known_values_indicator[i, j] == 1)  {
     //                            
     //                                L[i, j] =  known_values[i, j] / D[j];
     //                          
     //                          } else { 
     //                            
     //                                real sqrt_L_ij_old = sqrt(L_ij_old); 
     //                                real low = max({-sqrt_L_ij_old * D[j], lb_corr[i, j] - b1});
     //                                real up  = min({ sqrt_L_ij_old * D[j], ub_corr[i, j] - b1}); 
     //                     
     //                                real x = lb_ub_jacobian(Omega_raw_mat[i, j], low, up);   ////  lb_ub_jacobian(off_raw[counter], low, up);
     //                                L[i, j] = x / D[j]; 
     //                      
     //                                jacobian += -0.5 * log(D[j]);
     //                                
     //                          }
     //                          
     //                          L_ij_old *= 1.0 - (D[j] * L[i, j]^2) / L_ij_old;
     //                        
     //              }
     //            
     //              D[i] = L_ij_old;
     //            
     //        }
     //        
     //        return diag_post_multiply(L, sqrt(D)); //// = L_Omega
     //        
     // }
     // 
 
      // // need to add citation to this (slight modification from a HP. calculators forum post)
      // real inv_Phi_approx_from_prob(real p) { 
      //       return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 * logit(p)  )) ;
      // }
      // 
      // // need to add citation to this (slight modification from a HP. calculators forum post)
      // vector inv_Phi_approx_from_prob(vector p) { 
      //       return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 * logit(p)  )) ;  
      // }
      // 
      // // need to add citation to this  (slight modification from a HP. calculators forum post)
      // real inv_Phi_approx_from_logit_prob(real logit_p) { 
      //       return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 * logit_p  )) ; 
      // }
      // 
      // // need to add citation to this (slight modification from a HP. calculators forum post)
      // vector inv_Phi_approx_from_logit_prob(vector logit_p) { 
      //       return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 *logit_p  )) ; 
      // }
      // 
      // vector rowwise_sum(matrix M) {      // M is (N x T) matrix
      //       return M * rep_vector(1.0, cols(M));
      // }
      // 
      // vector rowwise_max(matrix M) {      // M is (N x T) matrix
      //       int N =  rows(M);
      //       vector[N] rowwise_maxes;
      //       for (n in 1:N) {
      //         rowwise_maxes[n] = max(M[n, ]);
      //       }
      //       return rowwise_maxes;
      // }
      // 
      // vector log_sum_exp_2d(matrix array_2d_to_lse) { 
      //       int N = rows(array_2d_to_lse);
      //       matrix[N, 2] rowwise_maxes_2d_array;
      //       rowwise_maxes_2d_array[, 1] =  rowwise_max(array_2d_to_lse);
      //       rowwise_maxes_2d_array[, 2] =  rowwise_maxes_2d_array[, 1];
      //       return  rowwise_maxes_2d_array[, 1] + log(rowwise_sum(exp((array_2d_to_lse  -  rowwise_maxes_2d_array))));
      // }
 







