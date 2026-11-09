// NicoSIMD version of Stan_fns_Jacobian.stan written by NicoStan (included file md5 a4bdee454322724717c6233773e519d6; NicoSIMD declarations md5 7f272e2ac00d6303123b13fb8f59653e; writer implementation md5 f6843a34264d79c6ed96effa2d04b0c4)






real raw_scale_to_scale_log_det_J_lp(   vector raw_scale,
                                        int softplus) {
                                      
       real log_det_J = 0.0;
       if (softplus == 1)  log_det_J += sum(nico_log_inv_logit(raw_scale));
       else                log_det_J += sum(raw_scale);
       return log_det_J;
  
}
real raw_scale_to_scale_log_det_J_lp(   row_vector raw_scale,
                                        int softplus) {
                                          
                                      
       real log_det_J = 0.0;
       if (softplus == 1)  log_det_J += sum(nico_log_inv_logit(raw_scale));
       else                log_det_J += sum(raw_scale);
       return log_det_J;
  
}
real raw_scale_to_scale_log_det_J_lp(   matrix raw_scale,
                                        int softplus) {
                                          
       real log_det_J = 0.0;                                
       if (softplus == 1)  log_det_J += sum(nico_log_inv_logit(raw_scale));
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
       if (softplus == 1)  log_det_J += sum(nico_log_inv_logit(raw_C[2:n_cutpoints]));
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
       if (softplus == 1)  log_det_J += sum(nico_log_inv_logit(raw_C[2:n_cutpoints]));
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
       if (softplus == 1)  log_det_J += sum(nico_log_inv_logit(raw_C[, 2:n_cutpoints]));
       else                log_det_J += sum(raw_C[, 2:n_cutpoints]);
       return log_det_J;
  
}

















