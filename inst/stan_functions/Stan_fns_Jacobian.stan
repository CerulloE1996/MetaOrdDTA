





real raw_scale_to_scale_log_det_J_lp(   vector raw_scale,
                                        int softplus) {
                                      
       real log_det_J = 0.0;
       if (softplus == 1)  log_det_J += sum(raw_scale);
       else                log_det_J += sum(log_inv_logit(raw_scale));
       return log_det_J;
  
}
real raw_scale_to_scale_log_det_J_lp(   row_vector raw_scale,
                                        int softplus) {
                                          
                                      
       real log_det_J = 0.0;
       if (softplus == 1)  log_det_J += sum(raw_scale);
       else                log_det_J += sum(log_inv_logit(raw_scale));
       return log_det_J;
  
}
real raw_scale_to_scale_log_det_J_lp(   matrix raw_scale,
                                        int softplus) {
                                          
       real log_det_J = 0.0;                                
       if (softplus == 1)  log_det_J += sum(raw_scale);
       else                log_det_J += sum(log_inv_logit(raw_scale));
       return log_det_J;
  
}











real raw_C_to_C_log_det_J_lp(   vector raw_C,
                                int softplus) {
                                  
       int n_cutpoints = num_elements(raw_C);
       real log_det_J = 0.0;                              
       if (softplus == 1)  log_det_J += sum(raw_C[2:n_cutpoints]);
       else                log_det_J += sum(log_inv_logit(raw_C[2:n_cutpoints]));
       return log_det_J;
  
}
real raw_C_to_C_log_det_J_lp(   row_vector raw_C,
                                int softplus) {
                                  
       int n_cutpoints = num_elements(raw_C);
       real log_det_J = 0.0;                             
       if (softplus == 1)  log_det_J += sum(raw_C[2:n_cutpoints]);
       else                log_det_J += sum(log_inv_logit(raw_C[2:n_cutpoints]));
       return log_det_J;
  
}
real raw_C_to_C_log_det_J_lp(   matrix raw_C,
                                int softplus) {
                                  
       int n_cutpoints = cols(raw_C);
       // int n_studies   = rows(raw_C);
       real log_det_J = 0.0;                             
       if (softplus == 1)  log_det_J += sum(raw_C[, 2:n_cutpoints]);
       else                log_det_J += sum(log_inv_logit(raw_C[, 2:n_cutpoints]));
       return log_det_J;
  
}


















