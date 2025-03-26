

////
//// Custom Stan functions - ordinal regression + cutpoint functions  -----------------------------------------------------------------------------------
////

////
//// Function to manually construct a cutpoint vector from raw unconstrained parameters - using exp() / log-differences. 
////
vector construct_C_using_exp_jacobian(vector C_raw_vec) {
  
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
                     //// Jacobian for trans. C_raw_vec -> C_vec:
                     jacobian += C_raw_vec[k];
                     counter += 1;
                 }
          }
       
          ////
          //// Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}


////
//// Function to manually construct a cutpoint vector from raw unconstrained parameters - using softplus / log1p_exp().
////
vector construct_C_using_SP_jacobian(vector C_raw_vec) {
      
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
                     //// Jacobian for trans. C_raw -> C. 
                     //// (deriv of SP(x) = log(1.0 + exp(x)) is inv_logit(x) so log(deriv) = log_inv_logit(x)):
                     jacobian += log_inv_logit(C_raw_vec[k]); 
                     counter += 1;
                 }
          }
       
          ////
          //// Output cutpoints + Jacobian (as 1st element in output vector):
          return(C_vec);
  
}







//// ---------------------------------------------------------------------------------------
//// Induced-Dirichlet ("ind_dir") log-density function:
//// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
//// NOTE: adapted from: Betancourt et al (see: https://betanalpha.github.io/assets/case_studies/ordinal_regression.html)
////
real induced_dirichlet_given_rho_lpdf(  vector p_ord,
                                        vector rho,
                                        vector alpha) {
          
        int n_cat = num_elements(p_ord);
        matrix[n_cat, n_cat] J = rep_matrix(0.0, n_cat, n_cat);
        
        //// Jacobian computation:
        for (k in 1:n_cat) {
                J[k, 1] = 1.0;
        }
        for (k in 2:n_cat) {
                // real rho =  normal_pdf(inv_Phi( p_cumul[k - 1])); //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
                J[k, k] = +rho[k - 1];
                J[k - 1, k] = -rho[k - 1];
        }
        
        return dirichlet_lpdf(p_ord | alpha) + log_determinant(J);
  
}
////
//// Row-vector overload for "induced_dirichlet_given_rho_lpdf":
////
real induced_dirichlet_given_rho_lpdf(  row_vector p_ord,
                                 row_vector rho,
                                 row_vector alpha) {
          

        int n_cat = num_elements(p_ord);
        int n_thr = n_cat - 1;
        vector[n_cat] p_ord_col_vec = to_vector(p_ord);
        vector[n_thr] rho_col_vec   = to_vector(rho);
        vector[n_cat] alpha_col_vec = to_vector(alpha);
        
        return induced_dirichlet_given_rho_lpdf(p_ord_col_vec | rho_col_vec, alpha_col_vec);
  
}





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
vector ord_probs_to_cumul_probs(vector ord_probs) {
  
        int n_cat = num_elements(ord_probs);
        int n_thr = n_cat - 1;
        vector[n_thr] cumul_probs;
       
        cumul_probs[1] = ord_probs[1];
        for (k in 2:n_thr) {
             cumul_probs[k] = cumul_probs[k - 1] + ord_probs[k]; // since probs are INCREASING with k
        }
        
        return cumul_probs;
        
}
////
//// Convert from ord_probs -> cumul_probs (row-vector overload):
////
row_vector ord_probs_to_cumul_probs(row_vector ord_probs) {
  
       return to_row_vector(ord_probs_to_cumul_probs(to_vector(ord_probs)));
        
}



////
//// Convert from cumul_probs -> C:
////
vector cumul_probs_to_C( vector cumul_probs,
                         vector location, 
                         vector scale) {
  
        int n_thr = num_elements(cumul_probs);
        int n_cat = n_thr + 1;
        vector[n_thr] C;
        
        for (k in 1:n_thr) {
              if (cumul_probs[k] < 1e-38) {
                    C[k] = -15.0; ////  prob = 1e-38;
              } else if (cumul_probs[k] > 0.9999999999999) {
                    C[k] = +15.0; ////  prob = 0.9999999999999;
              } else {
                     C[k] = location[k] + scale[k]*inv_Phi(cumul_probs[k]);
              }
        }
        
        return C;
        
}
vector cumul_probs_to_C( vector cumul_probs,
                         real location, 
                         real scale) {
  
        int n_thr = num_elements(cumul_probs);
        int n_cat = n_thr + 1;
        vector[n_thr] C;
        
        for (k in 1:n_thr) {
              if (cumul_probs[k] < 1e-38) {
                    C[k] = -15.0; ////  prob = 1e-38;
              } else if (cumul_probs[k] > 0.9999999999999) {
                    C[k] = +15.0; ////  prob = 0.9999999999999;
              } else {
                    C[k] = location + scale*inv_Phi(cumul_probs[k]);
              }
        }
        
        return C;
        
}
row_vector cumul_probs_to_C( row_vector cumul_probs,
                             real location, 
                             real scale) {
  
       return to_row_vector(cumul_probs_to_C(to_vector(cumul_probs), location, location));
        
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
////  Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_fixed_C( vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C);  
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[k] - locations[c, s])/scales[c, s];
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}
array[] matrix map_latent_cumul_prob_to_fixed_C( vector C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C);  
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[k] - locations[c, s])/scale;
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}


////
////  Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_fixed_hetero_C( array[] vector C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C[1]);  
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[c][k] - locations[c, s])/scales[c, s];
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}
array[] matrix map_latent_cumul_prob_to_fixed_hetero_C( array[] vector C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      int n_thr = num_elements(C[1]);  
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[c][k] - locations[c, s])/scale;
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}





////
////  Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_random_C( matrix C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data int  n_thr,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[s, k] - locations[c, s])/scales[c, s];
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}
array[] matrix map_latent_cumul_prob_to_random_C( matrix C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data int  n_thr,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[s, k] - locations[c, s])/scale;
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}




////
////  Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
////
array[] matrix map_latent_cumul_prob_to_random_hetero_C( array[] matrix C, 
                                  matrix locations, 
                                  matrix scales, 
                                  data int n_studies,
                                  data int  n_thr,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[c][s, k] - locations[c, s])/scales[c, s];
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}
array[] matrix map_latent_cumul_prob_to_random_hetero_C( array[] matrix C, 
                                  matrix locations, 
                                  real scale, 
                                  data int n_studies,
                                  data int  n_thr,
                                  data array[] int n_obs_cutpoints,
                                  data array[] matrix cutpoint_index) {
          
      array[2] matrix[n_studies, n_thr] latent_cumul_prob;
      for (c in 1:2) {
          latent_cumul_prob[c] = rep_matrix(positive_infinity(), n_studies, n_thr);
      }
      for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                              int k = to_int(cutpoint_index[c][s, cut_i]);
                              latent_cumul_prob[c][s, cut_i] = (C[c][s, k] - locations[c, s])/scale;
                      }
                }
      }
      
      return latent_cumul_prob;
                                    
}





////
////
////
vector convert_Dirichlet_alpha_to_SDs_jacobian(vector alpha) {

      ////
      //// Pre-compute some quantities needed for Jacobian for the (alpha_k, alpha_0) -> dirichlet_cat_SDs_sigma transformation:
      ////
      int n_cat = num_elements(alpha);
      real log_half = -0.6931471805599453;
      real alpha_0 = sum(alpha);
      real alpha_0_sq = square(alpha_0);
      real alpha_0_cube = alpha_0_sq * alpha_0;
      vector[n_cat] alpha_sq = square(alpha);
      ////
      //// NOTE: SD for each category probability in a Dirichlet is sqrt(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      ////
      vector[n_cat] dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) ./ (alpha_0_sq * (alpha_0 + 1.0))); // We are putting a prior on this!!
      ////
      //// Then compute the Jacobian adjustment:
      ////
      vector[n_cat] log_sigma = log(dirichlet_cat_SDs_sigma);
      vector[n_cat] deriv_A = alpha_sq + alpha_0 - 2.0*alpha;
      vector[n_cat] deriv_B = (1.0 ./ alpha) .* (3.0*alpha_0_sq*alpha + 2.0*alpha*alpha_0);
      vector[n_cat] deriv_var  = deriv_A .* (alpha_0_cube + alpha_0_sq) + deriv_B .* (alpha*alpha_0 - alpha_sq);
      vector[n_cat] log_abs_deriv_var_wrt_alpha = log(abs(deriv_var));
      vector[n_cat] log_abs_deriv_SD_wrt_alpha = log_half - log_sigma + log_abs_deriv_var_wrt_alpha;
      ////
      jacobian += sum(log_abs_deriv_SD_wrt_alpha);
      ////
      //// Output:
      ////
      return(dirichlet_cat_SDs_sigma);

}



        

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
        






















