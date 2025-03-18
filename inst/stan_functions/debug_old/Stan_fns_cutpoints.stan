

////
//// Custom Stan functions - cutpoint functions  ---------------------------------------------------------------------------------------
////


////
//// Function to manually construct a cutpoint vector from raw unconstrained parameters.
//// Using exp() / log-differences. 
////
vector construct_cutpoints_using_exp_jacobian(vector C_raw_vec) {
  
          int n_total_cutpoints = num_elements(C_raw_vec);
          vector[n_total_cutpoints] C_vec;  
          ////
          //// 1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          ////
          C_vec[1] = C_raw_vec[1];
          ////
          //// Rest of cutpoints are made using LOG-differences:
          ////
          int counter = 2; //// already set first element so skip first index
          for (k in 2:n_total_cutpoints) {
                 if (counter <= num_elements(C_raw_vec)) {
                     C_vec[k] = C_vec[k - 1] + exp(C_raw_vec[counter]);
                     jacobian += C_raw_vec[k]; //// Jacobian for trans. C_raw_vec -> C_vec
                     counter += 1;
                 }
          }
       
          ////
          //// Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}


/////
//// Function to manually construct a cutpoint vector from raw unconstrained parameters.
//// Using softplus / log1p_exp().
////
vector construct_cutpoints_using_softplus_jacobian(vector C_raw_vec) {
      
          int n_total_cutpoints = num_elements(C_raw_vec);
          vector[n_total_cutpoints] C_vec;  
          ////
          //// 1st cutpoint (for each test) is unconstrained (no Jacobian needed):
          ////
          C_vec[1] = C_raw_vec[1];
          ////
          //// Rest of cutpoints are made using LOG-differences:
          ////
          int counter = 2; //// already set first element so skip first index
          for (k in 2:n_total_cutpoints) {
                 if (counter <= num_elements(C_raw_vec)) {
                     C_vec[k] = C_vec[k - 1] + log1p_exp(C_raw_vec[counter]);
                     jacobian += log_inv_logit(C_raw_vec[k]); //// Jacobian for trans. C_raw -> C. (deriv of SP(x) = log(1.0 + exp(x)) is inv_logit(x) so log(deriv) = log_inv_logit(x)). 
                     counter += 1;
                 }
          }
       
          ////
          //// Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}

