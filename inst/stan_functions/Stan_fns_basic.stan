



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




 
 
 
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" 
//// around 1 (just like exp(x)!) rather than log(2)
real softplus_scaled_jacobian(real x) {
     real log_log_2 = -0.3665129205816643;
     jacobian += log_log_2 + log_inv_logit(x);
     return softplus_scaled(x);
}
vector softplus_scaled_jacobian(vector x) {
     real log_log_2 = -0.3665129205816643;
     jacobian += sum(log_log_2 + log_inv_logit(x));
     return softplus_scaled(x);
}
row_vector softplus_scaled_jacobian(row_vector x) {
     real log_log_2 = -0.3665129205816643;
     jacobian += sum(log_log_2 + log_inv_logit(x));
     return softplus_scaled(x);
}
matrix softplus_scaled_jacobian(matrix x) {
     real log_log_2 = -0.3665129205816643;
     jacobian += sum(log_log_2 + log_inv_logit(x));
     return softplus_scaled(x);
}







real exp_jacobian(real x_raw) {
     jacobian += x_raw;
     return exp(x_raw);
}
vector exp_jacobian(vector x_raw) {
     jacobian += sum(x_raw);
     return exp(x_raw);
}
row_vector exp_jacobian(row_vector x_raw) {
     jacobian += sum(x_raw);
     return exp(x_raw);
}
matrix exp_jacobian(matrix x_raw) {
     jacobian += sum(x_raw);
     return exp(x_raw);
}







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

 
 
 
 
 
 





