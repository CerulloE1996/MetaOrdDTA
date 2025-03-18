
 
////
//// Custom Stan functions - simplex functions  ---------------------------------------------------------------------------------------
//// NOTE: adapted from: Betancourt et al (see: https://betanalpha.github.io/assets/case_studies/ordinal_regression.html)
////



//// ---------------------------------------------------------------------------------------
//// Induced-Dirichlet ("ind_dir") log-density function:
//// NOTE: You can use this for both ind_dir PRIORS and ind_dir MODELS:
////
real induced_dirichlet_v2_lpdf(  vector p_ord,
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
//// Row-vector overload for "induced_dirichlet_v2_lpdf":
real induced_dirichlet_v2_lpdf(  row_vector p_ord,
                                 row_vector rho,
                                 row_vector alpha) {
          

        int n_cat = num_elements(p_ord);
        vector[n_cat] p_ord_col_vec = to_vector(p_ord);
        vector[n_cat] rho_col_vec   = to_vector(rho);
        vector[n_cat] alpha_col_vec = to_vector(alpha);
        
        
        return induced_dirichlet_v2_lpdf(p_ord_col_vec | rho_col_vec, alpha_col_vec);
  
}



