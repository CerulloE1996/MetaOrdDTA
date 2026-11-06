// reduce_sum_chunked version of DTA_MA_ord_HSROC_RANDthr_kappa_NicoStan_fused_AVX.stan, written by NicoStan: its nuisance blocks (beta_z, raw_scale_z, C_raw) are held by NicoStan inside the chunks
//// DTA_MA_ord_HSROC_RANDthr_kappa_NicoStan_fused_AVX.stan
//// The model of MetaOrdDTA's DTA_MA_ord_HSROC_RANDthr_kappa.stan with NicoSIMD's AVX functions:
//// the likelihood and the induced Dirichlet of the cutpoints each in one autodiff node.
////

functions {
    // (NicoStan's readers of its nuisance block: reduce_sum_chunked)
    real nicostan_nuisance_value(int term, int values_before_the_block, real type);
    vector nicostan_nuisance_vector(int first_term, int last_term, int values_before_the_block, real type);
    vector nicostan_nuisance_vector_at(array[] int terms, int values_before_the_block, real type);
    array[] real nicostan_nuisance_array_of_reals(int first_term, int last_term, int values_before_the_block, real type);
    array[] real nicostan_nuisance_array_of_reals_at(array[] int terms, int values_before_the_block, real type);
    array[] vector nicostan_nuisance_array_of_vectors(int first_term, int last_term, int values_per_term, int values_before_the_block, real type);
    array[] row_vector nicostan_nuisance_array_of_row_vectors(int first_term, int last_term, int values_per_term, int values_before_the_block, real type);
    matrix nicostan_nuisance_matrix_of_terms(int first_term, int last_term, int values_per_term, int values_before_the_block, real type);
    // (the constraints of the nuisance blocks: Stan math's own transforms, NicoStan holding the unconstrained values)
    real nicostan_constrain(real x, int kind, real a, real b);
    vector nicostan_constrain(vector x, int kind, real a, real b);
    row_vector nicostan_constrain(row_vector x, int kind, real a, real b);
    matrix nicostan_constrain(matrix x, int kind, real a, real b);
    array[] real nicostan_constrain(array[] real x, int kind, real a, real b);
    array[] vector nicostan_constrain(array[] vector x, int kind, real a, real b);
    array[] row_vector nicostan_constrain(array[] row_vector x, int kind, real a, real b);
    array[] matrix nicostan_constrain(array[] matrix x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(real x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(vector x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(row_vector x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(matrix x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(array[] real x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(array[] vector x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(array[] row_vector x, int kind, real a, real b);
    real nicostan_log_jacobian_of_constraint(array[] matrix x, int kind, real a, real b);

        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////



// 
// real safe_Phi(real x) {
//   
//     if (is_nan(x)) {
//         return -1.0;  // or whatever default makes sense
//     } else if (x < -37.5) {
//         return 0.0;
//     } else if (x > 8.25) {
//         return 1.0;
//     } else {
//         return Phi(x);
//     }
//     
// }
// vector safe_Phi(vector x) {
//   
//       int N = num_elements(x);
//       vector[N] out;
//       
//       for (n in 1:N) { 
//         out[n] = safe_Phi(x[n]);
//       }
//       
//       return(out);
//       
// }


////
//// "Basic" custom Stan functions: ---------------------------------------------------------------------------------------
////
array[] matrix init_array_of_matrices( data int n_rows, 
                                       data int n_cols,
                                       data int n_arrays,
                                       real init_val) {
                                         
      //// Create a local copy that you can modify (In Stan, function parameters are passed by value (not by reference), 
      //// unlike C++ where you can use "&"):
      array[n_arrays] matrix[n_rows, n_cols] array_out;
    
      for (c in 1:n_arrays) {
          array_out[c]  = rep_matrix(init_val, n_rows, n_cols);
      }
      
      return array_out;
          
}






array[,] matrix init_nested_array_of_matrices(  data int n_rows, 
                                                data int n_cols,
                                                data int n_arrays_rows,
                                                data int n_arrays_cols,
                                                real init_val) { 
                                 
        array[n_arrays_rows, n_arrays_cols] matrix[n_rows, n_cols] array_of_arrays_out;
                                 
        for (i in 1:n_arrays_rows) {
              for (j in 1:n_arrays_cols) {
                       array_of_arrays_out[i, j]  = rep_matrix(init_val, n_rows, n_cols);
              }
        }
          
          return array_of_arrays_out;
                                 
                                 
}



array[] matrix compute_ss_accuracy_from_cumul_prob( array[] matrix cumul_prob, 
                                                    data int n_studies,
                                                    data int n_thr) {
  
          array[3] matrix[n_studies, n_thr]  out_array;
          ////
          out_array[1] = 1.0 - cumul_prob[1];  //  = fp
          out_array[2] = 1.0 - out_array[1];   //  = sp
          out_array[3] = 1.0 - cumul_prob[2];  //  = se
          ////
          return(out_array);
  
}







//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" 
//// around 1 (just like exp(x)!) rather than log(2)
//// using "log_2_recip = 1.4426950408889634" as this can be much more efficient 
//// than doing 1.0/log(2) every time. 
////
real softplus_scaled(real x) {
     real log_2_recip = 1.4426950408889634;  
     return log_2_recip*log1p_exp(x);
}
vector softplus_scaled(vector x) {
     real log_2_recip = 1.4426950408889634;
     return log_2_recip*log1p_exp(x);
}
row_vector softplus_scaled(row_vector x) {
     real log_2_recip = 1.4426950408889634;
     return log_2_recip*log1p_exp(x);
}
matrix softplus_scaled(matrix x) {
     real log_2_recip = 1.4426950408889634;
     return log_2_recip*log1p_exp(x);
}




 
 
 
// //// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" 
// //// around 1 (just like exp(x)!) rather than log(2)
// real softplus_scaled_jacobian(real x) {
//      real log_log_2 = -0.3665129205816643;
//      jacobian += log_log_2 + log_inv_logit(x);
//      return softplus_scaled(x);
// }
// vector softplus_scaled_jacobian(vector x) {
//      real log_log_2 = -0.3665129205816643;
//      jacobian += sum(log_log_2 + log_inv_logit(x));
//      return softplus_scaled(x);
// }
// row_vector softplus_scaled_jacobian(row_vector x) {
//      real log_log_2 = -0.3665129205816643;
//      jacobian += sum(log_log_2 + log_inv_logit(x));
//      return softplus_scaled(x);
// }
// matrix softplus_scaled_jacobian(matrix x) {
//      real log_log_2 = -0.3665129205816643;
//      jacobian += sum(log_log_2 + log_inv_logit(x));
//      return softplus_scaled(x);
// }
// 





// 
// real exp_jacobian(real x_raw) {
//      jacobian += x_raw;
//      return exp(x_raw);
// }
// vector exp_jacobian(vector x_raw) {
//      jacobian += sum(x_raw);
//      return exp(x_raw);
// }
// row_vector exp_jacobian(row_vector x_raw) {
//      jacobian += sum(x_raw);
//      return exp(x_raw);
// }
// matrix exp_jacobian(matrix x_raw) {
//      jacobian += sum(x_raw);
//      return exp(x_raw);
// }
// 
// 





//// Normal PDF for a real input x (overloaded fn):
real std_normal_pdf( real x ) {
                   
     real sqrt_2_pi = 2.5066282746310002;
     real sqrt_2_pi_recip =  0.3989422804014327;
     return sqrt_2_pi_recip * exp(-0.5 * x * x);
  
}

//// Normal PDF for a vector input x (overloaded fn):
vector std_normal_pdf( vector x ) {

     real sqrt_2_pi = 2.5066282746310002;
     real sqrt_2_pi_recip =  0.3989422804014327;
     return sqrt_2_pi_recip * exp(-0.5 .* x .* x);

}

//// Normal PDF for a row_vector input x (overloaded fn):
row_vector std_normal_pdf( row_vector x ) {

     real sqrt_2_pi = 2.5066282746310002;
     real sqrt_2_pi_recip =  0.3989422804014327;
     return sqrt_2_pi_recip * exp(-0.5 .* x .* x);

}







////
//// Median for vector (overloaded fn):
////
real median(vector x) {

      int n = num_elements(x);
      vector[n] sorted_x = sort_asc(x);
      
      if (n % 2 == 0) {   // For even number of elements, average the middle two
          int left_element = to_int(n/2.0);
          int right_element = to_int(n/2.0 + 1.0);
          return (sorted_x[left_element] + sorted_x[right_element]) / 2.0;
      } else {            // For odd number of elements, return the middle one
          int middle_element = to_int((n+1.0)/2.0);
          return sorted_x[middle_element];
      }
}

////
//// Median for a row_vector (overloaded fn):
////
real median(row_vector x) {

      int n = num_elements(x);
      row_vector[n] sorted_x = sort_asc(x);
      
      if (n % 2 == 0) {   // For even number of elements, average the middle two
          int left_element = to_int(n/2.0);
          int right_element = to_int(n/2.0 + 1.0);
          return (sorted_x[left_element] + sorted_x[right_element]) / 2.0;
      } else {    // For odd number of elements, return the middle one
          int middle_element = to_int((n+1.0)/2.0);
          return sorted_x[middle_element];
      }
}


////
//// rowMedians function:
////
vector rowMedians(matrix x) {
  
      int dim = cols(x);
      vector[dim] out_vec;
      
      for (k in 1:dim) {
            out_vec[k] = median(col(x, k));
      }
      
      return out_vec;
  
}





// 
// // Function to compute quantile of a vector
// real vector_quantile( vector v, 
//                       real p) {
//                         
//     int n = num_elements(v);
//     vector[n] v_sorted;
//     real result;
//     
//     // Sort the vector
//     v_sorted = sort_asc(v);
//     
//     // Calculate index h (with interpolation)
//     real h = (n - 1) * p + 1;
//     
//     // Get the integer and fractional parts of h
//     int h_floor = floor(h);
//     real h_frac = h - h_floor;
//     
//     // Handle edge cases
//     if (h_floor <= 0) return v_sorted[1];
//     if (h_floor >= n) return v_sorted[n];
//     
//     // Interpolate between adjacent values
//     result = (1 - h_frac) * v_sorted[h_floor] + h_frac * v_sorted[h_floor + 1];
//     
//     return result;
//   
// }
// 
// // Convenience function to compute quantiles for a matrix column
// real matrix_column_quantile( matrix m, 
//                              int c, 
//                              real p) {
//     n_rows = rows(m);                            
//     vector[n_rows] v;
//     
//     for (i in 1:n_rows) {
//       v[i] = m[i, c];
//     }
//     
//     return vector_quantile(v, p);
//   
// }
// 


 

 

real std_normal_approx_pdf( real x, 
                            real Phi_x) { 
  
      real a_times_3 = 3.0 * 0.07056;
      real b = 1.5976;
      real x_sq = x*x;
      real three_a_x_sq_plus_b = fma(a_times_3, x_sq, b); //  3*a*x^2 + b
      
      return three_a_x_sq_plus_b * Phi_x * (1.0 - Phi_x);
  
}
vector std_normal_approx_pdf( vector x, 
                              vector Phi_x) { 
  
      int dim = num_elements(x);
      real a_times_3 = 3.0 * 0.07056;
      real b = 1.5976;
      vector[dim] x_sq = x .* x;
      vector[dim] three_a_x_sq_plus_b = fma(a_times_3, x_sq, b); //  3*a*x^2 + b
      
      return three_a_x_sq_plus_b .* (Phi_x .* (1.0 - Phi_x));
  
}
row_vector std_normal_approx_pdf( row_vector x, 
                                  row_vector Phi_x) { 
                              
      int dim = num_elements(x);
      real a_times_3 = 3.0 * 0.07056;
      real b = 1.5976;
      row_vector[dim] x_sq = x .* x;
      row_vector[dim] three_a_x_sq_plus_b = fma(a_times_3, x_sq, b);
      
      return three_a_x_sq_plus_b .* (Phi_x .* (1.0 - Phi_x)); //  3*a*x^2 + b
  
}
matrix std_normal_approx_pdf( matrix x, 
                              matrix Phi_x) { 
  
      int n_rows = rows(x);
      int n_cols = cols(x);
      real a_times_3 = 3.0 * 0.07056;
      real b = 1.5976;
      matrix[n_rows, n_cols] x_sq = x .* x;
      matrix[n_rows, n_cols] three_a_x_sq_plus_b = fma(a_times_3, x_sq, b);
      
      return three_a_x_sq_plus_b .* (Phi_x .* (1.0 - Phi_x)); //  3*a*x^2 + b
  
}










real  Phi_approx_2(real x)  {
  
      real a = 0.07056;
      real b = 1.5976;
      real x_sq = x*x;
      real a_x_sq_plus_b = fma(a, x_sq, b);
      real stuff_to_inv_logit =  x*a_x_sq_plus_b;
      
      return inv_logit(stuff_to_inv_logit);
  
}
vector  Phi_approx_2(vector x)  {
  
      int dim = num_elements(x);
      real a = 0.07056;
      real b = 1.5976;
      vector[dim] x_sq = x .* x;
      vector[dim] a_x_sq_plus_b = fma(a, x_sq, b);
      vector[dim] stuff_to_inv_logit =  x .* a_x_sq_plus_b;
      
      return inv_logit(stuff_to_inv_logit);
  
}
row_vector  Phi_approx_2(row_vector x)  {
  
      int dim = num_elements(x);
      real a = 0.07056;
      real b = 1.5976;
      row_vector[dim] x_sq = x .* x;
      row_vector[dim] a_x_sq_plus_b = fma(a, x_sq, b);
      row_vector[dim] stuff_to_inv_logit =  x .* a_x_sq_plus_b;
      
      return inv_logit(stuff_to_inv_logit);
  
}
matrix  Phi_approx_2(matrix x)  {
  
      int n_rows = rows(x);
      int n_cols = cols(x);
      real a = 0.07056;
      real b = 1.5976;
      matrix[n_rows, n_cols] x_sq = x .* x;
      matrix[n_rows, n_cols] a_x_sq_plus_b = fma(a, x_sq, b);
      matrix[n_rows, n_cols] stuff_to_inv_logit =  x .* a_x_sq_plus_b;
      
      return inv_logit(stuff_to_inv_logit);
  
}








// need to add citation to this (slight modification from a HP. calculators forum post)
real inv_Phi_approx_from_prob(real p) { 
      return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 * logit(p)  )) ;
}

// need to add citation to this (slight modification from a HP. calculators forum post)
vector inv_Phi_approx_from_prob(vector p) { 
      return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 * logit(p)  )) ;  
}

// need to add citation to this  (slight modification from a HP. calculators forum post)
real inv_Phi_approx_from_logit_prob(real logit_p) { 
      return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 * logit_p  )) ; 
}

// need to add citation to this (slight modification from a HP. calculators forum post)
vector inv_Phi_approx_from_logit_prob(vector logit_p) { 
      return 5.494 *  sinh(0.33333333333333331483 * asinh( 0.3418 *logit_p  )) ; 
}





// vector rowwise_sum(matrix M) {      // M is (N x T) matrix
//       return M * rep_vector(1.0, cols(M));
// }
// 
// vector rowwise_max(matrix M) {      // M is (N x T) matrix
// 
//       int N =  rows(M);
//       vector[N] rowwise_maxes;
//       
//       for (n in 1:N) {
//         rowwise_maxes[n] = max(M[n, ]);
//       }
//       
//       return rowwise_maxes;
// }
// 
// vector log_sum_exp_2d(matrix array_2d_to_lse) {
//   
//       int N = rows(array_2d_to_lse);
//       matrix[N, 2] rowwise_maxes_2d_array;
//       rowwise_maxes_2d_array[, 1] =  rowwise_max(array_2d_to_lse);
//       rowwise_maxes_2d_array[, 2] =  rowwise_maxes_2d_array[, 1];
//       return  rowwise_maxes_2d_array[, 1] + log(rowwise_sum(exp((array_2d_to_lse  -  rowwise_maxes_2d_array))));
//       
// }
// 

 
 
 
 
 
 








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
 









////
//// Custom Stan functions - ordinal regression + cutpoint functions  --------------------------------------------------------------
////




real compute_AUC(  vector se_vec, 
                   vector sp_vec, 
                   int n_thr, 
                   real missing_value_marker) {
                       
      // Calculate AUC for a single test
      // Inputs:
      //   se_vec: vector of sensitivities (length n_thr)
      //   sp_vec: vector of specificities (length n_thr)
      //   n_thr: number of thresholds for this test
      //   invalid_threshold: value below which Se/Sp are considered invalid (typically -0.5)
      // Returns:
      //   AUC value
      
      // Count valid (non-missing) thresholds
      int n_valid = 0;
      for (k in 1:n_thr) {
          // Check if this threshold has actual data (not missing)
          if (se_vec[k] != missing_value_marker && sp_vec[k] != missing_value_marker) {
              n_valid += 1;
          }
      }
      ////
      // If no valid thresholds, return NaN or -1.0
      if (n_valid == 0) {
          return -1.0;
      }
      ////
      // We need to add points (0,0) and (1,1) to complete the ROC curve
      int n_points = n_valid + 2;
      vector[n_points] FPR;
      vector[n_points] TPR;
      ////
      // Add (0,0) point
      FPR[1] = 0.0;
      TPR[1] = 0.0;
      ////
      // Add only the VALID threshold points
      int idx = 2;
      for (k in 1:n_thr) {
          if (se_vec[k] != missing_value_marker && sp_vec[k] != missing_value_marker) {
              FPR[idx] = 1.0 - sp_vec[k];
              TPR[idx] = se_vec[k];
              idx += 1;
          }
      }
      ////
      // Add (1,1) point
      FPR[n_points] = 1.0;
      TPR[n_points] = 1.0;
      ////
      // Sort by FPR using bubble sort
      for (i in 1:(n_points-1)) {
          for (j in 1:(n_points-i)) {
              if (FPR[j] > FPR[j+1]) {
                  // Swap FPR
                  real temp_fpr = FPR[j];
                  FPR[j] = FPR[j+1];
                  FPR[j+1] = temp_fpr;
                  // Swap TPR
                  real temp_tpr = TPR[j];
                  TPR[j] = TPR[j+1];
                  TPR[j+1] = temp_tpr;
              }
          }
      }
      ////
      // Calculate AUC using trapezoidal rule
      real auc = 0.0;
      for (i in 1:(n_points-1)) {
          auc += 0.5 * (TPR[i] + TPR[i+1]) * (FPR[i+1] - FPR[i]);
      }
      ////
      return auc;
    
}


real SD_approx_ID_ord_prob_to_C_probit( real mean_C_cutpoint_scale, 
                                        real SD_prob_scale) {
                                 
    real SD_probit_scale = (1.0 / std_normal_pdf(mean_C_cutpoint_scale)) * SD_prob_scale;
    return SD_probit_scale;
                            
}
vector SD_approx_ID_ord_prob_to_C_probit( vector mean_C_cutpoint_scale, 
                                          vector SD_prob_scale) {
                                            
    int dim = num_elements(SD_prob_scale);
    vector[dim] SD_probit_scale = (1.0 / std_normal_pdf(mean_C_cutpoint_scale)) .* SD_prob_scale;
    return SD_probit_scale;
                            
}
////
//// ---- logit scale:
////
real SD_approx_ID_ord_prob_to_C_logit(  real mean_prob_scale, 
                                        real SD_prob_scale) {
                                 
    real SD_logit_scale = (1.0 / (mean_prob_scale * (1.0 - mean_prob_scale))) * SD_prob_scale;
    return SD_logit_scale;
                            
}
vector SD_approx_ID_ord_prob_to_C_logit( vector mean_prob_scale, 
                                         vector SD_prob_scale) {
                                           
    int dim = num_elements(SD_prob_scale);
    vector[dim] SD_logit_scale = (1.0 ./ (mean_prob_scale .* (1.0 - mean_prob_scale))) .* SD_prob_scale;
    return SD_logit_scale;
                            
}
 
  
  
////
//// Function to manually construct a cutpoint vector from raw unconstrained parameters - using exp() / log-differences. 
////
vector construct_C(  vector C_raw_vec, 
                     int softplus) {
  
  
          int n_total_cutpoints = num_elements(C_raw_vec);
          vector[n_total_cutpoints] C_vec;  
          //// ----  1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          C_vec[1] = C_raw_vec[1];
          //// ---- Rest of cutpoints are made using LOG-differences (OR softplus-differences):
          if (softplus == 1) {
                  vector[n_total_cutpoints - 1] softplus_C_vec = log1p_exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + softplus_C_vec[k - 1];
                             // jacobian += log_inv_logit(C_raw_vec[k]);   //// Jacobian for trans. C_raw -> C:
                  }
          } else { 
                  vector[n_total_cutpoints - 1] exp_C_vec = exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + exp_C_vec[k - 1];
                             // jacobian += C_raw_vec[k];  //// Jacobian for trans. C_raw -> C:
                  }
          }
          //// ---- Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}
row_vector construct_C(  row_vector C_raw_vec, 
                         int softplus) {
  
          int n_total_cutpoints = num_elements(C_raw_vec);
          row_vector[n_total_cutpoints] C_vec;  
          //// ----  1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          C_vec[1] = C_raw_vec[1];
          //// ---- Rest of cutpoints are made using LOG-differences (OR softplus-differences):
          if (softplus == 1) {
                  row_vector[n_total_cutpoints - 1] softplus_C_vec = log1p_exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + softplus_C_vec[k - 1];
                             // jacobian += log_inv_logit(C_raw_vec[k]);   //// Jacobian for trans. C_raw -> C:
                  }
          } else { 
                  row_vector[n_total_cutpoints - 1] exp_C_vec = exp(C_raw_vec[2:n_total_cutpoints]);
                  for (k in 2:n_total_cutpoints) {
                             C_vec[k] = C_vec[k - 1] + exp_C_vec[k - 1];
                             // jacobian += C_raw_vec[k];  //// Jacobian for trans. C_raw -> C:
                  }
          }
          //// ---- Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}




//// ---------------------------------------------------------------------------------------
//// Induced-Dirichlet ("ind_dir") log-density function:
//// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
//// NOTE: adapted from: Betancourt et al (see: https://betanalpha.github.io/assets/case_studies/ordinal_regression.html),
//// HOWEVER my version has a (much) more computationally efficient (lower-trianglar) Jacobian computation, which is 
//// mathematically still valid. 
////
real induced_dirichlet_given_C_lpdf(    vector p_ord,
                                        vector C,
                                        vector alpha,
                                        int use_probit_link) {
          
        int n_cat = num_elements(p_ord);
        int n_thr = n_cat - 1;
        
        real log_prob = 0.0;

        for (k in 1:n_thr) {
          if (use_probit_link == 1)  log_prob += normal_lpdf(C[k] | 0.0, 1.0);
          else                       log_prob += log_inv_logit(C[k]) + log1m_inv_logit(C[k]);
        }
        log_prob += dirichlet_lpdf(p_ord | alpha);
        
        return log_prob;
  
}
////
//// Row-vector overload for "induced_dirichlet_given_rho_lpdf":
////
real induced_dirichlet_given_C_lpdf(  row_vector p_ord,
                                      row_vector C,
                                      row_vector alpha,
                                      int use_probit_link) {
          

        int n_cat = num_elements(p_ord);
        int n_thr = n_cat - 1;
        vector[n_cat] p_ord_col_vec = to_vector(p_ord);
        vector[n_thr] C_col_vec     = to_vector(C);
        vector[n_cat] alpha_col_vec = to_vector(alpha);
        
        return induced_dirichlet_given_C_lpdf(p_ord_col_vec | C_col_vec, alpha_col_vec, use_probit_link);
  
}





//// ---------------------------------------------------------------------------------------
//// Induced-Dirichlet ("ind_dir") log-density function:
//// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
//// NOTE: adapted from: Betancourt et al (see: https://betanalpha.github.io/assets/case_studies/ordinal_regression.html)
////
// real induced_dirichlet_given_rho_lpdf(  vector p_ord,
//                                         vector rho,
//                                         vector alpha) {
          
//         int n_cat = num_elements(p_ord);
//         matrix[n_cat, n_cat] J = rep_matrix(0.0, n_cat, n_cat);
         
//         //// Jacobian computation:
//         for (k in 1:n_cat) {
//                 J[k, 1] = 1.0;
//         }
//         // // // for (k in 2:(n_cat - 1)) {
//         // // //    J[k, k] = +rho[k];
//         // // // }
//         // for (k in 2:n_cat) { // p_ord[k] = Phi(C_{k} - anchor) - Phi(C_{k - 1} - anchor)
//         //          J[k, k]     = - rho[k - 1];
//         //          J[k - 1, k] = + rho[k - 1];
//         // }
                
//         // First row (p_1 depends on c_1)
//         J[1, 2] = -rho[1];  // ∂p_1/∂c_1 = -ρ(φ-c_1)  [SIGN FLIPPED]
        
//         // Middle rows
//         for (k in 2:(n_cat-1)) {
//             J[k, k]     = + rho[k - 1];     // ∂p_k/∂c_{k-1} = ρ(φ-c_{k-1})  [SIGN FLIPPED] 
//             J[k, k - 1] = - rho[k];    // ∂p_k/∂c_k = -ρ(φ-c_k)  [SIGN FLIPPED]
//         }
        
//         // Last row (p_n_cat depends on c_{n_cat-1})
//         J[n_cat, n_cat] = rho[n_cat-1];  // ∂p_n_cat/∂c_{n_cat-1} = ρ(φ-c_{n_cat-1})  [SIGN FLIPPED]

//         return dirichlet_lpdf(p_ord | alpha) + log_determinant(J);
  
// }
////
//// Row-vector overload for "induced_dirichlet_given_rho_lpdf":
////
// real induced_dirichlet_given_rho_lpdf(  row_vector p_ord,
//                                  row_vector rho,
//                                  row_vector alpha) {
          

//         int n_cat = num_elements(p_ord);
//         int n_thr = n_cat - 1;
//         vector[n_cat] p_ord_col_vec = to_vector(p_ord);
//         vector[n_thr] rho_col_vec   = to_vector(rho);
//         vector[n_cat] alpha_col_vec = to_vector(alpha);
        
//         return induced_dirichlet_given_rho_lpdf(p_ord_col_vec | rho_col_vec, alpha_col_vec);
  
// }





//// Induced dirichlet from Betancourt, 2019. 
//// See https://betanalpha.github.io/assets/case_studies/ordinal_regression.html#3_cut_to_the_chase
// real induced_dirichlet_lpdf( vector c, 
//                              vector alpha, 
//                              real phi) {
                                   
//       int K = num_elements(c) + 1;
//       vector[K - 1] anchoredcutoffs = c - phi;
      
//       vector[K] sigma;
//       vector[K] p;
//       matrix[K, K] J = rep_matrix(0, K, K);
      
//       sigma[1:(K-1)] = Phi(anchoredcutoffs); 
//       sigma[K] = 1;
      
//       p[1] =  sigma[1];
//       for (k in 2:(K - 1))
//         p[k] = sigma[k] - sigma[k - 1];
//       p[K] = 1 - sigma[K - 1];
      
//       // Baseline column of Jacobian
//       for (k in 1:K) J[k, 1] = 1;
      
//       // Diagonal entries of Jacobian
//       for (k in 2:K) {
//         // rho is the PDF of the latent distribution
//         real rho = exp(normal_lpdf(anchoredcutoffs[k - 1] | 0.0, 1.0)); ////  1.702 * sigma[k - 1] * (1 - sigma[k - 1]);
//         J[k, k] =  - rho;
//         J[k - 1, k] = + rho;
//       }
      
//       return dirichlet_lpdf(p | alpha) + log_determinant(J);
  
// }


////
//// Convert from cumul_probs -> ord_probs:
////
vector cumul_probs_to_ord_probs(vector cumul_probs) {
  
        int n_thr = num_elements(cumul_probs);
        int n_cat = n_thr + 1;
        vector[n_cat] ord_probs;
       
        ord_probs[1] = cumul_probs[1] - 0.0;
        for (k in 2:n_thr) {
             ord_probs[k] = cumul_probs[k] - cumul_probs[k - 1]; // since probs are INCREASING with k
        }
        ord_probs[n_cat] =  1.0 - cumul_probs[n_cat - 1];
        
        return ord_probs;
        
}
////
//// Convert from cumul_probs -> ord_probs (row-vector overload):
////
row_vector cumul_probs_to_ord_probs(row_vector cumul_probs) {
  
       return to_row_vector(cumul_probs_to_ord_probs(to_vector(cumul_probs)));
        
}
////
//// Convert from cumul_probs -> ord_probs (matrix overload)
//// NOTE: This will be more efficient when n_rows > n_cols (e.g. if each row is a study this will often be the case).
//// NOTE: assumed each col corresponds to a threshold/cutpoint.
////
matrix cumul_probs_to_ord_probs(matrix cumul_probs) {
  
        int n_cols = cols(cumul_probs);
        int n_rows = rows(cumul_probs);
        int n_thr = n_cols;
        int n_cat = n_thr + 1;
        matrix[n_rows, n_cat] ord_probs;
        
        // for (k in 1:n_cols) { //// col-major storage
           for (s in 1:n_rows) {
                 ord_probs[s, 1] = cumul_probs[s, 1] - 0.0;
                 for (k in 2:n_thr) {
                     ord_probs[s, k] = cumul_probs[s, k] - cumul_probs[s, k - 1]; // since probs are INCREASING with k
                 }
                 ord_probs[s, n_cat] =  1.0 - cumul_probs[s, n_cat - 1];
           }
        // }
        
        return ord_probs;
         
}



////
//// Convert from ord_probs -> cumul_probs:
////
vector ord_probs_to_cumul_probs( vector ord_probs) {
  
        int n_cat = num_elements(ord_probs);
        int n_thr = n_cat - 1;
        vector[n_thr] cumul_probs;
       
        cumul_probs[1] = ord_probs[1];
        for (k in 2:n_thr) {
             cumul_probs[k] = cumul_probs[k - 1] + ord_probs[k]; // since probs are INCREASING with k
        }
        
        return cumul_probs;
        
}
row_vector ord_probs_to_cumul_probs(row_vector ord_probs) {
  
       return to_row_vector(ord_probs_to_cumul_probs(to_vector(ord_probs)));
        
}




////
//// Convert from cumul_probs -> C (w/ induced dirichlet):
////
vector ID_cumul_probs_to_C( vector cumul_probs, 
                            int use_probit_link) {
  
        int n_thr = num_elements(cumul_probs);
        int n_cat = n_thr + 1;
        vector[n_thr] C;
        
        for (k in 1:n_thr) {
              if (cumul_probs[k] < 1e-300) {
                    C[k] = -37.5; ////  prob = 1e-38;
              } else if (cumul_probs[k] > 0.99999999999999999) {
                    C[k] = +8.20; ////  prob = 0.9999999999999;
              } else {
                    if (use_probit_link == 1) C[k] = inv_Phi(cumul_probs[k]); 
                    else                      C[k] = logit(cumul_probs[k]); 
              }
        }
        
        return C;
        
}
row_vector ID_cumul_probs_to_C( row_vector cumul_probs, 
                                int use_probit_link) {
  
       return to_row_vector(ID_cumul_probs_to_C(to_vector(cumul_probs), use_probit_link));
        
}



////
//// C_array -> C_mat function:
////
matrix convert_C_array_to_mat(array[] vector C_array) {
  
  
      int n_rows = size(C_array);
      int n_cols = num_elements(C_array[1]);
      matrix[n_rows, n_cols] C_mat;
  
      for (k in 1:n_cols) { //// col-major storage
         for (s in 1:n_rows) {
                  C_mat[s, k] = C_array[s][k];
             }
      }
  
      return C_mat;
 
}












 
 
////
////  Get the cutpoint index (k) to map "latent_surv_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_prob_to_fixed_C_HSROC( vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      int n_thr = num_elements(C); 
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = cutpoint_index[c, s, cut_i + 1];
                            if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = -(C[k] + locations[s, c])*scales[s, c];
                  }
          }
      }
      
      return latent_surv;
                                    
}
////
//// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_prob_to_random_C_HSROC( matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                                  int k = cutpoint_index[c, s, cut_i + 1];
                                  if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = -(C[s, k] + locations[s, c])*scales[s, c];
                     
                  }
          }
      }
      
      return latent_surv;
                                    
}
















// 
// ////
// ////  Get the cutpoint index (k) to map "latent_surv_prob[c][s, cut_i]" to correct cutpoint "C[k]":
// ////
// array[] matrix map_latent_surv_to_fixed_homog_C( vector C, 
//                                   matrix locations, 
//                                   matrix scales, 
//                                   data int n_studies,
//                                   data array[] int n_obs_cutpoints,
//                                   data array[,,] int cutpoint_index) {
//           
//       int n_thr = num_elements(C); 
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                             int k = cutpoint_index[c, s, cut_i + 1];
//                             if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = - (C[k] - locations[s, c])/scales[s, c];
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }
// array[] matrix map_latent_surv_to_fixed_homog_C( vector C, 
//                                   matrix locations, 
//                                   real scale, 
//                                   data int n_studies,
//                                   data array[] int n_obs_cutpoints,
//                                   data array[,,] int cutpoint_index) {
//           
//       int n_thr = num_elements(C); 
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                             int k = cutpoint_index[c, s, cut_i + 1];
//                             if ((k != 0) && (k != -1)) latent_surv[c][s, cut_i] = - (C[k] - locations[s, c])/scale;
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }


////
////  Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_to_fixed_hetero_C(  array[] vector C, 
                                                matrix locations, 
                                                matrix scales, 
                                                data int n_studies,
                                                data array[] int n_obs_cutpoints,
                                                data array[,,] int cutpoint_index) {
                                    
      int n_thr = num_elements(C[1]); 
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
    
      for (c in 1:2) {
          for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                               int k = cutpoint_index[c, s, cut_i + 1];
                               if (k != -1) { 
                                   latent_surv[c][s, cut_i] = - (C[c][k] - locations[s, c])/scales[s, c];
                              }
                      }
          }
      }
      
      return latent_surv;
                                    
}
array[] matrix map_latent_surv_to_fixed_hetero_C(  array[] vector C, 
                                                matrix locations, 
                                                data int n_studies,
                                                data array[] int n_obs_cutpoints,
                                                data array[,,] int cutpoint_index) {
          
      int n_thr = num_elements(C[1]); 
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
 
      for (c in 1:2) {
          for (s in 1:n_studies) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              int k = cutpoint_index[c, s, cut_i + 1];
                              if (k != -1)  latent_surv[c][s, cut_i] =   locations[s, c] - C[c][k];
                      }
          }
      }
      
      return latent_surv;
                                    
}


// ////
// //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
// ////
// array[] matrix map_latent_surv_to_random_homog_C( matrix C, 
//                                                   data int n_thr,
//                                                   matrix locations, 
//                                                   matrix scales, 
//                                                   data int n_studies,
//                                                   data array[] int n_obs_cutpoints,
//                                                   data array[,,] int cutpoint_index) {
//           
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                                   int k = cutpoint_index[c, s, cut_i + 1];
//                                   if (k != -1) latent_surv[c][s, cut_i] =  (locations[s, c] - C[s, k])/scales[s, c];
//                      
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }
// array[] matrix map_latent_surv_to_random_homog_C( matrix C,  
//                                                   data int n_thr,
//                                                   matrix locations, 
//                                                   data int n_studies,
//                                                   data array[] int n_obs_cutpoints,
//                                                   data array[,,] int cutpoint_index) {
//           
//       int n_cat = n_thr + 1;
//       array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
//       
//       for (c in 1:2) {
//           for (s in 1:n_studies) {
//                   for (cut_i in 1:n_obs_cutpoints[s]) {
//                                   int k = cutpoint_index[c, s, cut_i + 1];
//                                   if (k != -1)  latent_surv[c][s, cut_i] =   locations[s, c] - C[s, k];
//                   }
//           }
//       }
//       
//       return latent_surv;
//                                     
// }
// 



////
//// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_surv_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());

      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                                 int k = cutpoint_index[c, s, cut_i + 1];
                                 if (k != -1)  latent_surv[c][s, cut_i] = (locations[s, c] - C[c][s, k])/scales[s, c];
                  }
          }
      }
      
      return latent_surv;
                                    
}
array[] matrix map_latent_surv_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[,,] int cutpoint_index) {
          
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_surv  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                               int k = cutpoint_index[c, s, cut_i + 1];
                               if (k != -1)  latent_surv[c][s, cut_i] = locations[s, c] - C[c][s, k];
                        
                  }
          }
      }
      
      return latent_surv;
                                    
}

























////
////  Get the cutpoint index (k) to map "latent_surv_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_fixed_C( vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_prob_to_fixed_C( vector C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}


////
////  Get the cutpoint index (k) to map "latent_cumul[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_fixed_hetero_C( array[] vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
                                    
      int n_thr = num_elements(C[1]); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_prob_to_fixed_hetero_C( array[] vector C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C[1]); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}


////
//// ---- Get the cutpoint index (k) to map "latent_cumul[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_random_C( matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C); 
      // int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[s, k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_prob_to_random_C( matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C); 
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[s, k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}




////
//// ---- Get the cutpoint index (k) to map "latent_cumul[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C[1]); 
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());

      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][s, k] - locations[s, c])/scales[s, c];
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}
array[] matrix map_latent_cumul_to_random_hetero_C( array[] matrix C, 
                                  data int n_thr,
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      // int n_thr = num_elements(C[1]); 
      int n_cat = n_thr + 1;
      array[2] matrix[n_studies, n_thr] latent_cumul  = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
      
      for (c in 1:2) {
          for (s in 1:n_studies) {
                  // latent_cumul[c][s, 1] = 1.0;
                  for (cut_i in 1:n_obs_cutpoints[s]) {
                            int k = to_int(cutpoint_index[1][s, cut_i]);
                           //  if (k < n_thr + 1) {
                                  latent_cumul[c][s, cut_i] =  (C[c][s, k] - locations[s, c])/scale;
                            // }
                  }
          }
      }
      
      return latent_cumul;
                                    
}










// ////
// ////
// ////
// vector convert_Dirichlet_alpha_to_SDs_jacobian(vector alpha) {
// 
//       ////
//       //// Pre-compute some quantities needed for Jacobian for the (alpha_k, alpha_0) -> dirichlet_cat_SDs_sigma transformation:
//       ////
//       int n_cat = num_elements(alpha);
//       real log_half = -0.6931471805599453;
//       real alpha_0 = sum(alpha);
//       real alpha_0_sq = square(alpha_0);
//       real alpha_0_cube = alpha_0_sq * alpha_0;
//       vector[n_cat] alpha_sq = square(alpha);
//       ////
//       //// NOTE: SD for each category probability in a Dirichlet is sqrt(α_k(α_0-α_k)/(α_0^2(α_0+1))) - 
//       //// where α_0 is the sum of all alphas:
//       ////
//       vector[n_cat] dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) ./ (alpha_0_sq * (alpha_0 + 1.0))); 
//       ////
//       //// Then compute the Jacobian adjustment:
//       ////
//       vector[n_cat] log_sigma = log(dirichlet_cat_SDs_sigma);
//       vector[n_cat] deriv_A = alpha_sq + alpha_0 - 2.0*alpha;
//       vector[n_cat] deriv_B = (1.0 ./ alpha) .* (3.0*alpha_0_sq*alpha + 2.0*alpha*alpha_0);
//       vector[n_cat] deriv_var  = deriv_A .* (alpha_0_cube + alpha_0_sq) + deriv_B .* (alpha*alpha_0 - alpha_sq);
//       vector[n_cat] log_abs_deriv_var_wrt_alpha = log(abs(deriv_var));
//       vector[n_cat] log_abs_deriv_SD_wrt_alpha = log_half - log_sigma + log_abs_deriv_var_wrt_alpha;
//       ////
//       jacobian += sum(log_abs_deriv_SD_wrt_alpha);
//       ////
//       //// Output:
//       ////
//       return(dirichlet_cat_SDs_sigma);
// 
// }
// 


        

// for (k in 1:n_cat) {
//         // // real deriv_A[k] = alpha_sq[k] + alpha_0 - 2.0*alpha[k];
//         // // real deriv_B[k] = (1.0 / alpha[k]) * (3.0*alpha_0_sq*alpha[k] + 2.0*alpha[k]*alpha_0);
//         // // //// To compute on log-scale (need to keep track of signs / use special log_abs_sum_exp() function):
//         // // //// Term 1:
//         // // real deriv_term_1 = deriv_A[k] * (alpha_0_cube + alpha_0_sq);
//         // // real log_abs_deriv_term_1 = log(abs(deriv_A[k])) + log(abs(alpha_0_cube + alpha_0_sq));
//         // // real sign_deriv_term_1 = sign(deriv_A[k]) * sign(alpha_0_cube + alpha_0_sq);
//         // // //// Term 2:
//         // // real deriv_term_2 =  deriv_B[k] * (alpha[k]*alpha_0 - alpha_sq[k]);
//         // // real log_abs_deriv_term_2 = log(abs( deriv_B[k])) + log(abs(alpha[k]*alpha_0 - alpha_sq[k]));
//         // // real sign_deriv_term_2 = sign( deriv_B[k]) * sign(alpha[k]*alpha_0 - alpha_sq[k]);
//         // // //// Full term for deriv_var_k:
//         // // real log_abs_deriv_var_k_wrt_alpha = log_abs_sum_exp(log_abs_deriv_term_1, sign_deriv_term_1, log_abs_deriv_term_2, sign_deriv_term_2);
//         // //// To compute on normal scale:
//         // //// Since we have: SD_k = sqrt(Var_k), this means:
//         // //// deriv_SD_k_wrt_alpha = 0.5 * pow(Var_k, -0.5) * deriv_var_k_wrt_alpha. Hence:
//         // //// log_abs_deriv_SD_k_wrt_alpha = log(0.5) -0.5*log(Var_k) + log_abs_deriv_var_k_wrt_alpha. Hence:
//         // //// log_abs_deriv_SD_k_wrt_alpha = log(0.5) - log_sigma[k] + log_abs_deriv_var_k_wrt_alpha.
//         // // real deriv_var_k  = deriv_A[k] * (alpha_0_cube + alpha_0_sq) + deriv_B[k] * (alpha[k]*alpha_0 - alpha_sq[k]);
//         // // real log_abs_deriv_var_k_wrt_alpha = log(abs(deriv_var_k));
//         // ////
//         // real log_abs_deriv_SD_k_wrt_alpha = log_half - log_sigma[k] + log_abs_deriv_var_k_wrt_alpha;
//         // Jacobian_for_alpha_k_to_category_SDs += log_abs_deriv_SD_k_wrt_alpha;
// }
        






















////
//// ---------------------------------------------------------------------------------------
//// The induced-Dirichlet log density of a set of ordered cutpoints C (at the latent location 0) with every category
//// probability on the LOG scale (Stan's std_normal_lcdf / std_normal_lccdf with their log-scale tails, or
//// log_inv_logit / log1m_inv_logit, and log_diff_exp), no clamps: the same values as induced_dirichlet_given_C_lpdf
//// of the category probabilities of C wherever those are representable, and a finite log density (or -inf for two
//// equal cutpoints) where Phi(C_k) - Phi(C_(k-1)) rounds to 0 in double precision, which made dirichlet_lpdf throw.
//// The same calculation as NicoSIMD's nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf, for plain Stan math
//// models (so that the plain and the AVX versions of a model compute the same thing):
////     log p_1 = log F(C_1),  log p_(K+1) = log(1 - F(C_K)),
////     log p_k = log_diff_exp(log F(C_k), log F(C_(k-1)))                 (C_(k-1) < 0),
////             = log_diff_exp(log(1 - F(C_(k-1))), log(1 - F(C_k)))       (C_(k-1) >= 0).
////
real induced_dirichlet_given_C_on_the_log_scale_lpdf(vector C,
                                                     vector alpha,
                                                     int use_probit_link) {

        int n_thr = num_elements(C);
        vector[n_thr] log_F;
        vector[n_thr] log_S;
        if (use_probit_link == 1) {
          for (k in 1:n_thr) {
            log_F[k] = std_normal_lcdf(C[k] |);
            log_S[k] = std_normal_lccdf(C[k] |);
          }
        } else {
          log_F = log_inv_logit(C);
          log_S = log1m_inv_logit(C);
        }
        vector[n_thr + 1] log_p;
        log_p[1] = log_F[1];
        log_p[n_thr + 1] = log_S[n_thr];
        for (k in 2:n_thr) {
          log_p[k] = (C[k - 1] >= 0) ? log_diff_exp(log_S[k - 1], log_S[k]) : log_diff_exp(log_F[k], log_F[k - 1]);
        }
        real log_prob = lgamma(sum(alpha)) - sum(lgamma(alpha));
        for (j in 1:(n_thr + 1)) {
          if (alpha[j] != 1.0) log_prob += (alpha[j] - 1.0) * log_p[j];
        }
        if (use_probit_link == 1) log_prob += std_normal_lpdf(C);
        else                      log_prob += sum(log_F + log_S);
        return log_prob;

}
////
//// Row-vector overload:
////
real induced_dirichlet_given_C_on_the_log_scale_lpdf(row_vector C,
                                                     vector alpha,
                                                     int use_probit_link) {

        return induced_dirichlet_given_C_on_the_log_scale_lpdf(to_vector(C) | alpha, use_probit_link);

}


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


real safe_Phi(real x) {
  
    if (is_nan(x)) {
        return -1.0;  // or whatever default makes sense
    } else if (x < -37.5) {
        return 0.0;
    } else if (x > 8.25) {
        return 1.0;
    } else {
        return Phi(x);
    }
    
}
vector safe_Phi(vector x) {
  
      int N = num_elements(x);
      vector[N] out;
      
      for (n in 1:N) { 
        out[n] = safe_Phi(x[n]);
      }
      
      return(out);
      
}
matrix safe_Phi(matrix x) { 
  
      int N = rows(x);
      int M = cols(x);
      matrix[N, M] out;
      
      for (m in 1:M) {
        for (n in 1:N) { 
          out[n, m] = safe_Phi(x[n, m]);
        }
      }
      
      return(out);
  
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
                if (use_probit_link == 1) surv_prob = safe_Phi(latent_surv[c]);
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
              if (use_probit_link == 1) surv_prob[c] = safe_Phi(latent_surv[c]);
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
              if (use_probit_link == 1) surv_prob[c] = safe_Phi(latent_surv[c]);
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
              if (use_probit_link == 1) surv_prob[c] = safe_Phi(latent_surv[c]);
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
    
    
        
    
          
          
          
          



////
//// Custom Stan functions - simplex functions  ---------------------------------------------------------------
//// NOTE: adapted from:
//// https://github.com/mjhajharia/transforms/blob/main/transforms/simplex/StickbreakingLogistic.stan
////

 
////
//// Function to manually construct a simplex from raw unconstrained parameters.
////
// vector stickbreaking_logistic_simplex_constrain_jacobian(vector y) {
//             
//             //// Other variables:
//             int N = rows(y) + 1;
//             vector[N] x;
//             real log_zi;
//             real log_xi;
//             real log_cum_prod = 0;
//             
//             //// Construct simplex:
//             for (i in 1:(N - 1)) {
//                 log_zi = log_inv_logit(y[i] - log(N - i)); // logistic_lcdf(y[i] | log(N - i), 1)
//                 log_xi = log_cum_prod + log_zi;
//                 x[i] = exp(log_xi);
//                 log_cum_prod += log1m_exp(log_zi);
//                 jacobian += log_xi;
//             }   
//             
//             x[N] = exp(log_cum_prod);
//             
//             //// Increment jacobian adjustment:
//             jacobian += log_cum_prod;
//             
//             //// Output simplex:
//             return x;
// 
// }
vector stickbreaking_logistic_simplex_constrain(vector y) {
            
            //// Other variables:
            int N = rows(y) + 1;
            vector[N] x;
            real log_zi;
            real log_xi;
            real log_cum_prod = 0;
            
            //// Construct simplex:
            for (i in 1:(N - 1)) {
                log_zi = log_inv_logit(y[i] - log(N - i)); // logistic_lcdf(y[i] | log(N - i), 1)
                log_xi = log_cum_prod + log_zi;
                x[i] = exp(log_xi);
                log_cum_prod += log1m_exp(log_zi);
                // jacobian += log_xi;
            }   
            
            x[N] = exp(log_cum_prod);
            
            // //// Increment jacobian adjustment:
            // jacobian += log_cum_prod;
            
            //// Output simplex:
            return x;

}
real stickbreaking_logistic_simplex_constrain_lp(vector y) {
            
            //// Other variables:
            int N = rows(y) + 1;
            vector[N] x;
            real log_zi;
            real log_xi;
            real log_cum_prod = 0;
            real out = 0.0;
            
            //// Construct simplex:
            for (i in 1:(N - 1)) {
                log_zi = log_inv_logit(y[i] - log(N - i)); // logistic_lcdf(y[i] | log(N - i), 1)
                log_xi = log_cum_prod + log_zi;
                x[i] = exp(log_xi);
                log_cum_prod += log1m_exp(log_zi);
                // jacobian += log_xi;
            }   
            
            x[N] = exp(log_cum_prod);
            
            //// Increment jacobian adjustment:
            out += log_cum_prod;
            return out;

}

vector stickbreaking_logistic_simplex_constrain_inc_J(vector y) {
            
            //// Other variables:
            int N = rows(y) + 1;
            vector[N] x;
           
            real log_zi;
            real log_xi;
            real log_cum_prod = 0;
            ////
            real log_det_J = 0.0;
            
            //// Construct simplex:
            for (i in 1:(N - 1)) {
                log_zi = log_inv_logit(y[i] - log(N - i)); // logistic_lcdf(y[i] | log(N - i), 1)
                log_xi = log_cum_prod + log_zi;
                x[i] = exp(log_xi);
                log_cum_prod += log1m_exp(log_zi);
                log_det_J += log_xi;
            }   
            
            x[N] = exp(log_cum_prod);
            
            //// Increment jacobian adjustment:
            log_det_J += log_cum_prod;
            
            //// Output simplex + log_det_J:
            vector[N + 1] out;
            out[1] = log_det_J;
            out[2:(N+1)] = x;
            return out;

}

// vector stickbricking_angular_simplex_constrain_jacobian(vector y) {
//   
//           int N = rows(y) + 1;
//           vector[N] x;
//           real log_phi, phi, u, s, c;
//           real s2_prod = 1;
//           real log_halfpi = log(pi()) - log2();
//           int rcounter = 2 * N - 3;
//           
//             for (i in 1:(N-1)) {
//               u = log_inv_logit(y[i]);
//               log_phi = u + log_halfpi;
//               phi = exp(log_phi);
//               s = sin(phi);
//               c = cos(phi);
//               x[i] = s2_prod * c^2;
//               s2_prod *= s^2;
//               jacobian += log_phi + log1m_exp(u) + rcounter * log(s) + log(c);
//               rcounter -= 2;
//           }
//           
//           x[N] = s2_prod;
//           jacobian += (N - 1) * log2();
//           
//           return x;
//   
// }
// 
// vector stickbreaking_power_logistic_simplex_constrain_jacobian(vector y) {
//     
//           int N = rows(y) + 1;
//           vector[N] x;
//           real log_u, log_w, log_z;
//           real log_cum_prod = 0;
//           
//           for (i in 1:(N-1)) {
//             
//             log_u = log_inv_logit(y[i]); // logistic_lcdf(y[i] | 0, 1);
//             log_w = log_u / (N - i);
//             log_z = log1m_exp(log_w);
//             x[i] = exp(log_cum_prod + log_z);
//             jacobian += 2 * log_u - y[i]; // logistic_lupdf(y[i] | 0, 1);
//             log_cum_prod += log1m_exp(log_z);
//             
//           }
//           
//           x[N] = exp(log_cum_prod);
//           jacobian += -lgamma(N);
//           
//           return x;
//   
// }


// vector stickbreaking_power_normal_simplex_constrain_jacobian(vector y) {
//   
//            int N = rows(y) + 1;
//            vector[N] x;
//            real log_u, log_w, log_z;
//            real log_cum_prod = 0;
//            
//            for (i in 1:(N-1)) {
//              log_u = std_normal_lcdf(y[i] |);
//              log_w = log_u / (N - i);
//              log_z = log1m_exp(log_w);
//              x[i] = exp(log_cum_prod + log_z);
//              jacobian += std_normal_lpdf(y[i] |);
//              log_cum_prod += log1m_exp(log_z);
//            }
//            
//            x[N] = exp(log_cum_prod);
//            jacobian += -lgamma(N);
//            
//            return x;
//    
//  }
//  
//  
//  
//  vector stickbreaking_normal_simplex_constrain_jacobian(vector y) {
//     
//            int N = rows(y) + 1;
//            vector[N] x;
//            real log_zi, log_xi, wi;
//            real log_cum_prod = 0;
//            
//            for (i in 1:N - 1) {
//              wi = y[i] - log(N - i) / 2;
//              log_zi = std_normal_lcdf(wi);
//              log_xi = log_cum_prod + log_zi;
//              x[i] = exp(log_xi);
//              jacobian += std_normal_lpdf(wi) + log_cum_prod;
//              log_cum_prod += log1m_exp(log_zi);
//            }
//            
//            x[N] = exp(log_cum_prod);
//            
//            return x;
//    
//  }
//  
//  
  
  
  






real raw_scale_to_scale_log_det_J_lp(   vector raw_scale,
                                        int softplus) {
                                      
       real log_det_J = 0.0;
       if (softplus == 1)  log_det_J += sum(log_inv_logit(raw_scale));
       else                log_det_J += sum(raw_scale);
       return log_det_J;
  
}
real raw_scale_to_scale_log_det_J_lp(   row_vector raw_scale,
                                        int softplus) {
                                          
                                      
       real log_det_J = 0.0;
       if (softplus == 1)  log_det_J += sum(log_inv_logit(raw_scale));
       else                log_det_J += sum(raw_scale);
       return log_det_J;
  
}
real raw_scale_to_scale_log_det_J_lp(   matrix raw_scale,
                                        int softplus) {
                                          
       real log_det_J = 0.0;                                
       if (softplus == 1)  log_det_J += sum(log_inv_logit(raw_scale));
       else                log_det_J += sum(raw_scale);
       return log_det_J;
  
}











real raw_C_to_C_log_det_J_lp(   vector raw_C,
                                int softplus) {
                                  
       int n_cutpoints = num_elements(raw_C);
       real log_det_J = 0.0;                              
       // if (softplus == 1)  log_det_J += sum(raw_C[2:n_cutpoints]);
       // else                log_det_J += sum(log_inv_logit(raw_C[2:n_cutpoints]));
       //// ---- the same Jacobians as construct_C(): softplus -> log_inv_logit(C_raw), exp -> C_raw:
       if (softplus == 1)  log_det_J += sum(log_inv_logit(raw_C[2:n_cutpoints]));
       else                log_det_J += sum(raw_C[2:n_cutpoints]);
       return log_det_J;
  
}
real raw_C_to_C_log_det_J_lp(   row_vector raw_C,
                                int softplus) {
                                  
       int n_cutpoints = num_elements(raw_C);
       real log_det_J = 0.0;                             
       // if (softplus == 1)  log_det_J += sum(raw_C[2:n_cutpoints]);
       // else                log_det_J += sum(log_inv_logit(raw_C[2:n_cutpoints]));
       //// ---- the same Jacobians as construct_C(): softplus -> log_inv_logit(C_raw), exp -> C_raw:
       if (softplus == 1)  log_det_J += sum(log_inv_logit(raw_C[2:n_cutpoints]));
       else                log_det_J += sum(raw_C[2:n_cutpoints]);
       return log_det_J;
  
}
real raw_C_to_C_log_det_J_lp(   matrix raw_C,
                                int softplus) {
                                  
       int n_cutpoints = cols(raw_C);
       // int n_studies   = rows(raw_C);
       real log_det_J = 0.0;                             
       // if (softplus == 1)  log_det_J += sum(raw_C[, 2:n_cutpoints]);
       // else                log_det_J += sum(log_inv_logit(raw_C[, 2:n_cutpoints]));
       //// ---- the same Jacobians as construct_C(): softplus -> log_inv_logit(C_raw), exp -> C_raw:
       if (softplus == 1)  log_det_J += sum(log_inv_logit(raw_C[, 2:n_cutpoints]));
       else                log_det_J += sum(raw_C[, 2:n_cutpoints]);
       return log_det_J;
  
}


















    
    
    
    
    
    
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
        ////
        //// ---- NicoSIMD's AVX functions (their bodies in NicoSIMD's user header):
        ////
        real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
              array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
              data real bound);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(vector C, vector alpha, data int link,
                                                                    data real location);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(matrix C, vector alpha, data int link,
                                                                    data real location);
        vector nico_exp(vector x);
        vector nico_log1p_exp(vector x);
        vector nico_log_inv_logit(vector x);
        vector nico_log1m_inv_logit(vector x);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(vector C_raw, data int softplus);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(matrix C_raw, data int softplus);
        array[] matrix nico_map_latent_surv_prob_to_random_C_HSROC(data matrix C, data int n_thr,
                                                                   data matrix locations, data matrix scales,
                                                                   data int n_studies,
                                                                   data array[] int n_obs_cutpoints,
                                                                   data array[,,] int cutpoint_index);
        array[,] matrix nico_compute_log_lik_binomial_fact_data(data array[] matrix latent_surv,
                                                                data int use_probit_link, data int n_thr,
                                                                data int n_studies, data array[,,] int x_2,
                                                                data array[,,] int n, data array[,] int N_total,
                                                                data array[] int n_obs_cutpoints,
                                                                data array[] int K_fold_CV_indicator);
        vector nico_Phi(vector x);
        vector nico_inv_logit(vector x);
        real nico_normal_lpdf(vector y, real mu, real sigma);
        real nico_normal_lpdf(vector y, vector mu, vector sigma);
}


data {
        ////
        //// ---- Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;
        ////
        //// ---- Covariates:
        ////
        int n_covariates; //// For R+G - cannot have SEPERATE covariates for D+ and D-, so only one n_covariates!
        matrix[n_studies, n_covariates] X;
        vector[n_covariates] baseline_case; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
        array[2, n_studies, n_thr + 1] int x;
        array[2, n_studies, n_thr + 1] int cutpoint_index;
        ////
        //// ---- Priors for location ("beta"):
        ////
        vector[n_covariates] prior_beta_mu_mean;
        vector[n_covariates] prior_beta_mu_SD;
        real prior_beta_SD_mean;
        real prior_beta_SD_SD;
        ////
        //// ---- Priors for raw_scale ("gamma"):
        ////
        vector[n_covariates] prior_raw_scale_mu_mean;
        vector[n_covariates] prior_raw_scale_mu_SD;
        real prior_raw_scale_SD_mean;
        real prior_raw_scale_SD_SD;
        ////
        //// ---- Induced-Dirichlet priors:
        ////
        vector[n_thr + 1] prior_dirichlet_alpha;
        real prior_kappa_mean;
        real<lower=0.0> prior_kappa_SD;
        real<lower=0.0> prior_kappa_dof;
        ////
        //// ---- Other:
        ////
        real<lower=0.0> kappa_lb;
        int<lower=0, upper=1> softplus;
        ////
        int use_probit_link;
        ////
        array[n_studies] int<lower=0, upper=1> K_fold_CV_indicator;
}


transformed data {  
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        ////
        real mult_nd = -1.0;
        real mult_d  = +1.0; //// hence: mult_d > mult_nd
        ////
        array[2, n_studies, n_thr] int x_2;
        array[2, n_studies] int N_total;
        array[2, n_studies, n_thr] int n;
        ////
        for (s in 1:n_studies) {
              for (c in 1:2) {
                  N_total[c, s] = x[c, s, 1];
                  for (k in 1:n_thr) { 
                     x_2[c, s, k] = x[c, s, k + 1];
                  }
              }
              for (c in 1:2) {
                   n[c, s, 1] = N_total[c, s];
                   for (k in 2:n_obs_cutpoints[s]) { 
                              n[c, s, k] = x_2[c, s, k - 1];
                   }
              }
       }
        ////
        //// ---- The terms of the likelihood, one per reported cutpoint of each study and class (its counts and
        ////      its places in the locations and the cutpoints):
        ////
        int n_terms = 0;
        for (c in 1:2) {
              for (s in 1:n_studies) {
                    if (K_fold_CV_indicator[s] == 1) n_terms += n_obs_cutpoints[s];
              }
        }
        array[n_terms] int count_above_the_cutpoint_of_each_term;
        array[n_terms] int count_above_the_previous_cutpoint_of_each_term;
        vector[n_terms] is_the_first_reported_cutpoint_of_each_term;
        array[n_terms] int index_of_the_location_of_each_term;
        array[n_terms] int index_of_the_cutpoint_of_each_term;
        {
              int t = 0;
              for (c in 1:2) {
                    for (s in 1:n_studies) {
                          if (K_fold_CV_indicator[s] == 1) {
                                for (cut_i in 1:n_obs_cutpoints[s]) {
                                      int k = cutpoint_index[c, s, cut_i + 1];
                                      if (k < 1) reject("a reported cutpoint without its index: class ", c,
                                                        ", study ", s, ", reported cutpoint ", cut_i);
                                      t += 1;
                                      count_above_the_cutpoint_of_each_term[t] = x_2[c, s, cut_i];
                                      count_above_the_previous_cutpoint_of_each_term[t] = n[c, s, cut_i];
                                      is_the_first_reported_cutpoint_of_each_term[t] = (cut_i == 1) ? 1.0 : 0.0;
                                      index_of_the_location_of_each_term[t] = (c - 1) * n_studies + s;
                                      index_of_the_cutpoint_of_each_term[t] = (k - 1) * n_studies + s;
                                }
                          }
                    }
              }
        }


    // (the nuisance blocks held by NicoStan: reduce_sum_chunked)
    int nicostan_terms_of_beta_z = n_studies;
    int nicostan_values_per_term_of_beta_z = 1;
    int nicostan_offset_of_beta_z = 0;
    array[nicostan_terms_of_beta_z] int nicostan_term_indices_of_beta_z = linspaced_int_array(nicostan_terms_of_beta_z, 1, max(1, nicostan_terms_of_beta_z));
    int nicostan_terms_of_raw_scale_z = n_studies;
    int nicostan_values_per_term_of_raw_scale_z = 1;
    int nicostan_offset_of_raw_scale_z = nicostan_offset_of_beta_z + nicostan_terms_of_beta_z * nicostan_values_per_term_of_beta_z;
    array[nicostan_terms_of_raw_scale_z] int nicostan_term_indices_of_raw_scale_z = linspaced_int_array(nicostan_terms_of_raw_scale_z, 1, max(1, nicostan_terms_of_raw_scale_z));
    int nicostan_terms_of_C_raw = n_thr;
    int nicostan_values_per_term_of_C_raw = n_studies;
    int nicostan_offset_of_C_raw = nicostan_offset_of_raw_scale_z + nicostan_terms_of_raw_scale_z * nicostan_values_per_term_of_raw_scale_z;
    array[nicostan_terms_of_C_raw] int nicostan_term_indices_of_C_raw = linspaced_int_array(nicostan_terms_of_C_raw, 1, max(1, nicostan_terms_of_C_raw));
    int nicostan_size_of_the_nuisance_block = nicostan_offset_of_C_raw + nicostan_terms_of_C_raw * nicostan_values_per_term_of_C_raw;
}
  

parameters {
        // (held by NicoStan, reduce_sum_chunked:) vector[n_studies] beta_z
    //;  //// the nuisance block, declared first
        vector[n_covariates] beta_mu;  
        real<lower=0.0> beta_SD;   
        // vector[n_studies] beta_z;
        ////
        vector[n_covariates] raw_scale_mu;    
        real<lower=0.0> raw_scale_SD;   
        // (held by NicoStan, reduce_sum_chunked:) vector[n_studies] raw_scale_z
    //;
        
        // (held by NicoStan, reduce_sum_chunked:) matrix<lower=-7.5, upper=2.5>[n_studies, n_thr] C_raw
    //; //// restricting solves divs/numerical stab. issues w/ lots of sparse categories!
        real<lower=1.0> kappa;
        simplex[n_cat] dirichlet_phi; // props/probs
}


transformed parameters {
    // (the log-Jacobian of the constraint of the nuisance block C_raw, whose unconstrained values NicoStan holds: once per value)
    jacobian += nicostan_log_jacobian_of_constraint(
                    nicostan_nuisance_array_of_reals(1, nicostan_terms_of_C_raw * nicostan_values_per_term_of_C_raw,
                                                     nicostan_offset_of_C_raw, beta_mu[1]),
                    3, -7.5, 2.5);
 
        real log_kappa = log(kappa);
        vector<lower=0.0>[n_cat] alpha = 0.01 + dirichlet_phi * kappa;
        vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
        {
           real alpha_0 = sum(alpha);
           dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) / 
                                     (square(alpha_0) * (alpha_0 + 1.0)));
        }
        ////
        //// ---- Construct study-specific cutpoints:
        ////
        // matrix[n_studies, n_thr] C;
        // for (s in 1:n_studies) {
        //              C[s, ] =  construct_C(C_raw[s, 1:n_thr], softplus);
        // }

        //// (as construct_C, for all studies at once)
        matrix[n_studies, n_thr] C;
        C[, 1] = (nicostan_constrain(to_matrix(nicostan_nuisance_vector(1, nicostan_terms_of_C_raw * nicostan_values_per_term_of_C_raw, nicostan_offset_of_C_raw, beta_mu[1]), nicostan_values_per_term_of_C_raw, nicostan_terms_of_C_raw), 3, -7.5, 2.5))[, 1];
        if (n_thr > 1) {
              vector[n_studies * (n_thr - 1)] increments_of_the_cutpoints =
                    (softplus == 1) ? nico_log1p_exp(to_vector((nicostan_constrain(to_matrix(nicostan_nuisance_vector(1, nicostan_terms_of_C_raw * nicostan_values_per_term_of_C_raw, nicostan_offset_of_C_raw, beta_mu[1]), nicostan_values_per_term_of_C_raw, nicostan_terms_of_C_raw), 3, -7.5, 2.5))[, 2:n_thr])) :
                                      nico_exp(to_vector((nicostan_constrain(to_matrix(nicostan_nuisance_vector(1, nicostan_terms_of_C_raw * nicostan_values_per_term_of_C_raw, nicostan_offset_of_C_raw, beta_mu[1]), nicostan_values_per_term_of_C_raw, nicostan_terms_of_C_raw), 3, -7.5, 2.5))[, 2:n_thr]));
              for (k in 2:n_thr) {
                    C[, k] = C[, k - 1] +
                              increments_of_the_cutpoints[((k - 2) * n_studies + 1):((k - 1) * n_studies)];
              }
        }
        ////
        //// ---- Between-study model for location and scale:
        ////
        vector[n_studies] beta_random = beta_SD * (nicostan_nuisance_vector(1, nicostan_terms_of_beta_z, nicostan_offset_of_beta_z, beta_mu[1]));
        vector[n_studies] raw_scale_random = raw_scale_SD * (nicostan_nuisance_vector(1, nicostan_terms_of_raw_scale_z, nicostan_offset_of_raw_scale_z, beta_mu[1]));
        ////
        matrix[n_studies, 2] Xlocations; // NOT local
        matrix[n_studies, 2] Xscale;     // NOT local
        ////
        {
            vector[n_studies] Xbeta_baseline = X[, 1:n_covariates] * to_vector(beta_mu[1:n_covariates]) + beta_random;
            Xlocations[, 1] = mult_nd*Xbeta_baseline; 
            Xlocations[, 2] = mult_d*Xbeta_baseline; 
            ////
            vector[n_studies] Xraw_scale_baseline = X[, 1:n_covariates] * to_vector(raw_scale_mu[1:n_covariates]) + raw_scale_random;
            // Xscale[, 1] =  ((softplus == 1) ? softplus_scaled( -mult_nd*Xraw_scale_baseline) : exp( -mult_nd*Xraw_scale_baseline));
            Xscale[, 1] =  ((softplus == 1) ? softplus_scaled( -mult_nd*Xraw_scale_baseline) :
                                                nico_exp( -mult_nd*Xraw_scale_baseline)); 
            // Xscale[, 2] =  ((softplus == 1) ? softplus_scaled( -mult_d*Xraw_scale_baseline)  : exp( -mult_d*Xraw_scale_baseline));
            Xscale[, 2] =  ((softplus == 1) ? softplus_scaled( -mult_d*Xraw_scale_baseline)  :
                                                nico_exp( -mult_d*Xraw_scale_baseline));
        }
}


model {
    // (the nuisance block C_raw, held by NicoStan, read and constrained once)
    matrix[nicostan_values_per_term_of_C_raw, nicostan_terms_of_C_raw] C_raw_nicostan_constrained = nicostan_constrain(to_matrix(nicostan_nuisance_vector(1, nicostan_terms_of_C_raw * nicostan_values_per_term_of_C_raw, nicostan_offset_of_C_raw, beta_mu[1]), nicostan_values_per_term_of_C_raw, nicostan_terms_of_C_raw), 3, -7.5, 2.5);

        ////
        //// ---- Priors:
        ////
        // beta_mu[1:n_covariates] ~ normal(
        //                           prior_beta_mu_mean[1:n_covariates], 
        //                           prior_beta_mu_SD[1:n_covariates]);
        beta_mu[1:n_covariates] ~ nico_normal(prior_beta_mu_mean[1:n_covariates],
                                              prior_beta_mu_SD[1:n_covariates]);
        beta_SD ~ normal(
                  prior_beta_SD_mean, 
                  prior_beta_SD_SD);
        ////
        // raw_scale_mu[1:n_covariates] ~ normal(
        //                                prior_raw_scale_mu_mean[1:n_covariates], 
        //                                prior_raw_scale_mu_SD[1:n_covariates]);
        raw_scale_mu[1:n_covariates] ~ nico_normal(prior_raw_scale_mu_mean[1:n_covariates],
                                                   prior_raw_scale_mu_SD[1:n_covariates]);
        raw_scale_SD ~ normal(
                       prior_raw_scale_SD_mean, 
                       prior_raw_scale_SD_SD);
        ////
        {
           log_kappa ~ student_t(prior_kappa_dof, prior_kappa_mean, prior_kappa_SD);
           target += -log_kappa;
           dirichlet_phi ~ dirichlet(prior_dirichlet_alpha);
        }
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        // for (s in 1:n_studies) {
        //      vector[n_cat] Ind_Dir_ord_prob;
        //      if (use_probit_link == 1) Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(Phi(C[s, ])));
        //      else                      Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(inv_logit(C[s, ])));
        //      ////
        //      Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
        //                         to_vector(C[s, ]), alpha, use_probit_link);
        // }

        //// (every study's cutpoints in one node, as induced_dirichlet_given_C on each row)
        //// (anchored at the prior mean of the mid-location)
        target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(C | alpha, use_probit_link,
                                                                          0.5 * (mult_nd + mult_d) *
                                                                          prior_beta_mu_mean[1]);
        ////
        //// ---- Likelihood / Model:
        ////
        // target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior

        target += nico_normal_lupdf(to_vector((nicostan_nuisance_vector(1, nicostan_terms_of_beta_z, nicostan_offset_of_beta_z, beta_mu[1]))) | 0, 1); // part of between-study model, NOT prior
        // target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior

        target += nico_normal_lupdf(to_vector((nicostan_nuisance_vector(1, nicostan_terms_of_raw_scale_z, nicostan_offset_of_raw_scale_z, beta_mu[1]))) | 0, 1); // part of between-study model, NOT prior
        // target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus); // no Jacobian needed for R&G models as it cancels out (since + then -).
        // target += raw_C_to_C_log_det_J_lp(C_raw, softplus);
        target += nico_log_det_jacobian_of_cutpoints_from_raw_differences(C_raw_nicostan_constrained, softplus);
        ////
        //// ---- Log-likelihood:
        ////
        // {
        //     ////
        //     //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
        //     ////
        //     array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_random_C_HSROC(
        //                                                     C, n_thr, Xlocations, Xscale, n_studies, 
        //                                                     n_obs_cutpoints, cutpoint_index);
        //     target += compute_log_lik_binomial_fact_lp(
        //               latent_surv, use_probit_link, n_thr, n_studies, 
        //               x_2, n, N_total, n_obs_cutpoints,
        //               K_fold_CV_indicator);
        // }

        {
            ////
            //// ---- every term in one autodiff node:
            ////      nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf
            ////
            vector[n_studies * n_thr] cutpoints_of_every_study = to_vector(C);
            vector[2 * n_studies] locations_of_both_classes = to_vector(Xlocations);
            vector[2 * n_studies] scales_of_both_classes = to_vector(Xscale);
            //// (as map_latent_surv_prob_to_random_C_HSROC: -(C + location) * scale)
            vector[n_terms] latent_of_each_term =
                  -(cutpoints_of_every_study[index_of_the_cutpoint_of_each_term] +
                    locations_of_both_classes[index_of_the_location_of_each_term]) .*
                  scales_of_both_classes[index_of_the_location_of_each_term];
            target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                            count_above_the_cutpoint_of_each_term |
                            count_above_the_previous_cutpoint_of_each_term, latent_of_each_term,
                            is_the_first_reported_cutpoint_of_each_term, use_probit_link, 1e-10);
        }
}


generated quantities {
      ////
      //// ---- SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      ////
      vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma_full;
      {
          real alpha_0 = sum(alpha);
          dirichlet_cat_SDs_sigma_full = sqrt((alpha .* (alpha_0 - alpha) ) / (square(alpha_0) * (alpha_0 + 1.0)));
      }
      ////
      //// ---- Calculate summary cutpoints:
      ////
      vector[n_thr] C_mu;
      ////
      // {
      //     for (k in 1:n_thr) {
      //         vector[n_studies] observed_cutpoints;
      //         int n_obs = 0;
      //         
      //         for (s in 1:n_studies) {
      //             // Search through the cutpoint_index for this study to see if cutpoint k was observed
      //             for (j in 2:(n_obs_cutpoints[s] + 1)) {  // Start at 2 because index 1 is 0
      //                 if (cutpoint_index[1, s, j] == k) {
      //                     n_obs += 1;
      //                     observed_cutpoints[n_obs] = C[s, k];
      //                     break;  // Found it, no need to continue searching
      //                 }
      //             }
      //         }
      //         
      //         if (n_obs > 0) {
      //             C_mu[k] = median(observed_cutpoints[1:n_obs]);
      //         }
      //     }
      // }
      for (k in 1:n_thr) {
        C_mu[k] = median(C[, k]);
      }
      ////
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      // vector[n_covariates] location_nd_mu = mult_nd*beta_mu;
      // vector[n_covariates] location_d_mu  = mult_d*beta_mu; 
      ////
      // vector[n_covariates] scale_nd_mu =  ((softplus == 1) ? softplus_scaled(-mult_nd*raw_scale_mu) : exp(-mult_nd*raw_scale_mu));
      // vector[n_covariates] scale_d_mu  =  ((softplus == 1) ? softplus_scaled(-mult_d*raw_scale_mu)  : exp(-mult_d*raw_scale_mu));
      ////
      real Xbeta_baseline = dot_product(baseline_case, beta_mu);
      real Xbeta_baseline_nd = mult_nd * Xbeta_baseline;
      real Xbeta_baseline_d  = mult_d  * Xbeta_baseline;
      ////
      real Xraw_scale_baseline = dot_product(baseline_case, raw_scale_mu);
      real scale_nd_baseline = ((softplus == 1) ? 
                               softplus_scaled(-mult_nd * Xraw_scale_baseline) : exp(-mult_nd * Xraw_scale_baseline));
      real scale_d_baseline  = ((softplus == 1) ? 
                               softplus_scaled(-mult_d  * Xraw_scale_baseline) : exp(-mult_d  * Xraw_scale_baseline));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? 
                                  nico_Phi(-(C_mu + Xbeta_baseline_nd)*scale_nd_baseline)  : 
                                  nico_inv_logit(-(C_mu + Xbeta_baseline_nd)*scale_nd_baseline);
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  nico_Phi(-(C_mu + Xbeta_baseline_d)*scale_d_baseline)    : 
                                  nico_inv_logit(-(C_mu + Xbeta_baseline_d)*scale_d_baseline);
      ////
      //// ---- Calculate predictive accuracy:
      ////
      vector[n_thr] C_pred;
      {
        for (k in 1:n_thr) {
            real sigma_C_k;
            real SD_prob_scale = dirichlet_cat_SDs_sigma_full[k];
            if (use_probit_link == 1) {
                real mean_C_cutpoint_scale = C_mu[k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale, SD_prob_scale);
            } else { 
                real mean_prob_scale = dirichlet_phi[k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_logit(mean_prob_scale, SD_prob_scale);
            }
            if (!is_nan(sigma_C_k)) {
               C_pred[k] = C_mu[k] + normal_rng(0, sigma_C_k);
            }
        }
      }
      ////
      real beta_random_pred      = normal_rng(0.0, beta_SD);
      real raw_scale_random_pred = normal_rng(0.0, raw_scale_SD);
      ////
      real Xbeta_baseline_pred = dot_product(baseline_case, beta_mu) + beta_random_pred;
      real Xbeta_baseline_pred_nd = mult_nd * Xbeta_baseline_pred;
      real Xbeta_baseline_pred_d  = mult_d  * Xbeta_baseline_pred;
      ////
      real Xraw_scale_baseline_pred = dot_product(baseline_case, raw_scale_mu) + raw_scale_random_pred;
      real scale_baseline_pred_nd = ((softplus == 1) ? 
                                    softplus_scaled(-mult_nd * Xraw_scale_baseline_pred) : exp(-mult_nd * Xraw_scale_baseline_pred));
      real scale_baseline_pred_d  = ((softplus == 1) ? 
                                    softplus_scaled(-mult_d  * Xraw_scale_baseline_pred) : exp(-mult_d  * Xraw_scale_baseline_pred));
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C_pred + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd) : 
                                       nico_inv_logit(-(C_pred + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd);
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C_pred + Xbeta_baseline_pred_d)*scale_baseline_pred_d)   : 
                                       nico_inv_logit(-(C_pred + Xbeta_baseline_pred_d)*scale_baseline_pred_d);
      ////
      //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
      ////
      vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
      ////
      // matrix[n_studies, n_thr] fp; // global
      // matrix[n_studies, n_thr] sp; // global
      // matrix[n_studies, n_thr] se; // global
      ////
      // vector[n_studies] deviance_nd; // global
      // vector[n_studies] deviance_d;  // global
      // vector[n_studies] deviance;    // global
      {
          array[2] matrix[n_studies, n_thr] latent_surv = nico_map_latent_surv_prob_to_random_C_HSROC(
                                                          C, n_thr, Xlocations, Xscale, n_studies, 
                                                          n_obs_cutpoints, cutpoint_index);
          ////
          array[n_studies] int dummy_CV_indicator;
          for (s in 1:n_studies) {
              dummy_CV_indicator[s] = 1;
          }
          array[3, 2] matrix[n_studies, n_thr] outs = nico_compute_log_lik_binomial_fact_data(
                                                      latent_surv, use_probit_link, n_thr, n_studies,
                                                      x_2, n, N_total, n_obs_cutpoints,
                                                      dummy_CV_indicator);
          ////
          //// ---- Sum log-lik across thresholds and disease status within each study
          ////
          for (s in 1:n_studies) {
              for (k in 1:n_obs_cutpoints[s]) {
                log_lik_study[s] += outs[1][1][s, k];  // Non-diseased
                log_lik_study[s] += outs[1][2][s, k];  // Diseased
              }
          }
          ////
          array[2] matrix[n_studies, n_thr] cond_prob = outs[2];
          array[2] matrix[n_studies, n_thr] surv_prob = outs[3];
          ////
          // fp = surv_prob[1];
          // sp = 1.0 - fp;
          // se = surv_prob[2];
          // ////
          // //// ---- Model fit (deviance):
          // ////
          // array[4] matrix[n_studies, n_thr] outs_model_fit = compute_deviance(
          //                                                    cond_prob, n_thr, n_studies, 
          //                                                    x_2, n, n_obs_cutpoints);
          // ////
          // matrix[n_studies, n_thr] x_hat_nd = outs_model_fit[1];
          // matrix[n_studies, n_thr] dev_nd   = outs_model_fit[2];
          // matrix[n_studies, n_thr] x_hat_d  = outs_model_fit[3];
          // matrix[n_studies, n_thr] dev_d    = outs_model_fit[4];
          // ////
          // for (s in 1:n_studies) {
          //     for (k in 1:to_int(n_obs_cutpoints[s])) {
          //        deviance_nd[s] += is_nan(dev_nd[s, k]) ? 0.0 : dev_nd[s, k];
          //        deviance_d[s]  += is_nan(dev_d[s, k]) ? 0.0 : dev_d[s, k];
          //     }
          //     deviance[s] = deviance_nd[s] + deviance_d[s];
          // }
      }
      // ////
      // //// ---- "Ordinal-bivariate equivalent" params:
      // ////
      // array[2] matrix[n_thr, n_covariates] biv_equiv_C;
      // for (x_i in 1:n_covariates) {
      //     biv_equiv_C[1][, x_i] = to_vector(C_mu ./ scale_nd_mu[x_i]);
      //     biv_equiv_C[2][, x_i] = to_vector(C_mu ./ scale_d_mu[x_i]);
      // }
    // (the size of the nuisance block held by NicoStan, read by NicoStan: reduce_sum_chunked)
    int nicostan_size_of_the_nuisance_block_held_by_NicoStan = nicostan_size_of_the_nuisance_block;
}

