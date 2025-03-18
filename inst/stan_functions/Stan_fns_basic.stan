


////
//// "Basic" custom Stan functions: ----------------------------------------------------------------------------------------------------
////
array[] matrix init_array_of_matrices( int n_rows, 
                                       int n_cols,
                                       int n_arrays,
                                       real init_val) {
                                         
      //// Create a local copy that you can modify (In Stan, function parameters are passed by value (not by reference), unlike C++ where you can use "&"):
      array[n_arrays] matrix[n_rows, n_cols] array_out;
    
      for (c in 1:n_arrays) {
          array_out[c]  = rep_matrix(init_val, n_rows, n_cols);
      }
      
      return array_out;
          
}




array[] matrix compute_ss_accuracy_from_cumul_prob( array[] matrix cumul_prob, 
                                                    int n_studies,
                                                    int n_thr) {
  
          array[3] matrix[n_studies, n_thr]  out_array;
          ////
          out_array[1] = 1.0 - cumul_prob[1];  //  = fp
          out_array[2] = 1.0 - out_array[1];   //  = sp
          out_array[3] = 1.0 - cumul_prob[2];  //  = se
          ////
          return(out_array);
  
}



//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
real softplus_scaled(real x) {
     real log_2_recip = 1.4426950408889634;  // can be much more efficient than doing 1.0/log(2) every time. 
     return log_2_recip*log1p_exp(x);
}
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
vector softplus_scaled(vector x) {
     real log_2_recip = 1.4426950408889634;
     return log_2_recip*log1p_exp(x);
}
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
row_vector softplus_scaled(row_vector x) {
     real log_2_recip = 1.4426950408889634;
     return log_2_recip*log1p_exp(x);
}
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
matrix softplus_scaled(matrix x) {
     real log_2_recip = 1.4426950408889634;
     return log_2_recip*log1p_exp(x);
}

 





//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
real softplus_scaled_jacobian(real x) {
     real log_log_2 = -0.3665129205816643; // can be much more efficient than doing log(log(2)) every time.
     jacobian += log_log_2 + log_inv_logit(x);
     return softplus_scaled(x);
}
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
vector softplus_scaled_jacobian(vector x) {
     real log_log_2 = -0.3665129205816643; // can be much more efficient than doing log(log(2)) every time.
     jacobian += sum(log_log_2 + log_inv_logit(x));
     return softplus_scaled(x);
}
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
row_vector softplus_scaled_jacobian(row_vector x) {
     real log_log_2 = -0.3665129205816643; // can be much more efficient than doing log(log(2)) every time.
     jacobian += sum(log_log_2 + log_inv_logit(x));
     return softplus_scaled(x);
}
//// re-scaled softplus(x) = log(1.0 + exp(x)) function so that it's "centered" around 1 (just like exp(x)!) rather than log(2)
matrix softplus_scaled_jacobian(matrix x) {
     real log_log_2 = -0.3665129205816643; // can be much more efficient than doing log(log(2)) every time.
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
real normal_pdf( real x, 
                 real mu, 
                 real sigma) {
                   
     real sqrt_2_pi = sqrt(2 * pi());
     return (1.0 / (sigma * sqrt_2_pi)) * exp(-0.5 * ((x - mu) / sigma)^2);
  
}

//// Normal PDF for a vector input x (overloaded fn):
vector normal_pdf( vector x, 
                   real mu, 
                   real sigma) {

      real sqrt_2_pi = sqrt(2 * pi());
      return (1.0 / (sigma .* sqrt_2_pi)) * exp(-0.5 * ((x - mu) ./ sigma)^2);

}

//// Overload for vector x, vector mu, vector sigma (overloaded fn):
vector normal_pdf( vector x, 
                   vector mu, 
                   vector sigma) {
  
      real sqrt_2_pi = sqrt(2 * pi());
      return (1.0 / (sigma * sqrt_2_pi)) .* exp(-0.5 * ((x - mu) ./ sigma)^2);

}

//// Normal PDF for a row_vector input x (overloaded fn):
row_vector normal_pdf( row_vector x, 
                       real mu, 
                       real sigma) {

      real sqrt_2_pi = sqrt(2 * pi());
      return (1.0 / (sigma .* sqrt_2_pi)) * exp(-0.5 * ((x - mu) ./ sigma)^2);

}

//// Overload for row_vector x, row_vector mu, row_vector sigma (overloaded fn):
row_vector normal_pdf( row_vector x, 
                       row_vector mu, 
                       row_vector sigma) {
  
      real sqrt_2_pi = sqrt(2 * pi());
      return (1.0 / (sigma * sqrt_2_pi)) .* exp(-0.5 * ((x - mu) ./ sigma)^2);

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
// 
//   for (k in 1:n_thr) {
//         C_MU_empirical[k] = median(C[, k]);
//   }
//           
//           
          
          






