#' @useDynLib BayesMVP
NULL

#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_wrapper_EIGEN_double_mat <- function(x, 
                                          fn, 
                                          vect_type,
                                          skip_checks) {
  
  .Call(`_BayesMVP_Rcpp_wrapper_EIGEN_double_mat`,
        x, 
        fn,
        vect_type,
        skip_checks)
  
}




#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_wrapper_EIGEN_double_colvec <- function( x, 
                                              fn, 
                                              vect_type,
                                              skip_checks) {
  
  .Call(`_BayesMVP_Rcpp_wrapper_EIGEN_double_colvec`,
        x, 
        fn,
        vect_type,
        skip_checks)
  
}



#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_wrapper_fn_lp_grad <- function(Model_type, 
                                    force_autodiff, 
                                    force_PartialLog, 
                                    multi_attempts, 
                                    theta_main_vec, 
                                    theta_us_vec, 
                                    y, 
                                    grad_option, 
                                    Model_args_as_Rcpp_List) {
  
  .Call(`_BayesMVP_Rcpp_wrapper_fn_lp_grad`, 
        Model_type,
        force_autodiff, 
        force_PartialLog, 
        multi_attempts, 
        theta_main_vec,
        theta_us_vec, 
        y, 
        grad_option, 
        Model_args_as_Rcpp_List)
  
}

 

#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_compute_chain_stats <- function(mcmc_3D_array, 
                                     stat_type, 
                                          n_threads) {
  
  .Call(`_BayesMVP_Rcpp_compute_chain_stats`, 
        mcmc_3D_array, 
        stat_type, 
        n_threads)
  
}



#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_compute_MCMC_diagnostics <- function(mcmc_3D_array, 
                                          diagnostic, 
                                          n_threads) {
  
  .Call(`_BayesMVP_Rcpp_compute_MCMC_diagnostics`, 
        mcmc_3D_array, 
        diagnostic, 
        n_threads)
  
}

#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
detect_vectorization_support <- function() {
  .Call(`_BayesMVP_detect_vectorization_support`)
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_Rcpp_wrapper_update_M_dense_main_Hessian <- function(M_dense_main, 
                                                        M_inv_dense_main, 
                                                        M_inv_dense_main_chol, 
                                                        shrinkage_factor, 
                                                        ratio_Hess_main,
                                                        interval_width, 
                                                        num_diff_e, 
                                                        Model_type, 
                                                        force_autodiff, 
                                                        force_PartialLog,
                                                        multi_attempts, 
                                                        theta_main_vec, 
                                                        theta_us_vec, 
                                                        y, 
                                                        Model_args_as_Rcpp_List, 
                                                        ii, 
                                                        n_burnin, 
                                                        metric_type) {
  
  .Call(`_BayesMVP_fn_Rcpp_wrapper_update_M_dense_main_Hessian`, M_dense_main, M_inv_dense_main, 
        M_inv_dense_main_chol, shrinkage_factor, ratio_Hess_main, interval_width, num_diff_e, 
        Model_type, force_autodiff, force_PartialLog, multi_attempts, theta_main_vec, 
        theta_us_vec, y, Model_args_as_Rcpp_List, ii, n_burnin, metric_type)
  
}



#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_find_initial_eps_main_and_us <- function( theta_main_vec_initial_ref, 
                                                  theta_us_vec_initial_ref, 
                                                  seed, 
                                                  Model_type, 
                                                  force_autodiff,
                                                  force_PartialLog, 
                                                  multi_attempts, 
                                                  y_ref, 
                                                  Model_args_as_Rcpp_List, 
                                                  EHMC_args_as_Rcpp_List, 
                                                  EHMC_Metric_as_Rcpp_List) {
  
  .Call(`_BayesMVP_fn_find_initial_eps_main_and_us`,
        theta_main_vec_initial_ref, 
        theta_us_vec_initial_ref, 
        seed, Model_type, force_autodiff, force_PartialLog, multi_attempts, y_ref, 
        Model_args_as_Rcpp_List, EHMC_args_as_Rcpp_List, EHMC_Metric_as_Rcpp_List)
  
}


 





#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_Rcpp_wrapper_adapt_eps_ADAM <- function(eps,
                                           eps_m_adam, 
                                           eps_v_adam, 
                                           iter,
                                           n_burnin, 
                                           LR, 
                                           p_jump,
                                           adapt_delta,
                                           beta1_adam, 
                                           beta2_adam, 
                                           eps_adam) {
  
  .Call(`_BayesMVP_fn_Rcpp_wrapper_adapt_eps_ADAM`, eps, eps_m_adam, eps_v_adam, iter, 
        n_burnin, LR, p_jump, adapt_delta, beta1_adam, beta2_adam, eps_adam)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_update_snaper_m_and_s <- function(snaper_m,
                                     snaper_s_empirical, 
                                     theta_vec_mean, ii) {
  
  .Call(`_BayesMVP_fn_update_snaper_m_and_s`, snaper_m, snaper_s_empirical, 
        theta_vec_mean, ii)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_update_eigen_max_and_eigen_vec <- function(eigen_max,
                                              eigen_vector,
                                              snaper_w_vec) {
  
  .Call(`_BayesMVP_fn_update_eigen_max_and_eigen_vec`, eigen_max, eigen_vector, snaper_w_vec)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_update_snaper_w_dense_M <- function(snaper_w_vec,
                                       eigen_vector,
                                       eigen_max,
                                       theta_vec, 
                                       snaper_m_vec, 
                                       ii, 
                                       M_dense_sqrt) {
  
  .Call(`_BayesMVP_fn_update_snaper_w_dense_M`, snaper_w_vec, eigen_vector, eigen_max, 
        theta_vec, snaper_m_vec, ii, M_dense_sqrt)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_update_snaper_w_diag_M <- function(snaper_w_vec, 
                                      eigen_vector,
                                      eigen_max, 
                                      theta_vec, 
                                      snaper_m_vec, 
                                      ii, 
                                      sqrt_M_vec) {
  
  .Call(`_BayesMVP_fn_update_snaper_w_diag_M`, snaper_w_vec, eigen_vector, eigen_max, 
        theta_vec, snaper_m_vec, ii, sqrt_M_vec)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_Rcpp_wrapper_update_tau_w_diag_M_ADAM <- function(eigen_vector, 
                                                     eigen_max, 
                                                     theta_vec_initial, 
                                                     theta_vec_prop,
                                                     snaper_m_vec,
                                                     velocity_prop, 
                                                     velocity_0, 
                                                     tau, 
                                                     LR, 
                                                     ii, 
                                                     n_burnin, 
                                                     sqrt_M_vec, 
                                                     tau_m_adam, 
                                                     tau_v_adam,
                                                     tau_ii) {
  
  .Call(`_BayesMVP_fn_Rcpp_wrapper_update_tau_w_diag_M_ADAM`, eigen_vector, eigen_max, 
        theta_vec_initial, theta_vec_prop, snaper_m_vec, velocity_prop, velocity_0, 
        tau, LR, ii, n_burnin, sqrt_M_vec, tau_m_adam, tau_v_adam, tau_ii)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_Rcpp_wrapper_update_tau_w_dense_M_ADAM <- function(eigen_vector, 
                                                      eigen_max, 
                                                      theta_vec_initial, 
                                                      theta_vec_prop, 
                                                      snaper_m_vec, 
                                                      velocity_prop, 
                                                      velocity_0, 
                                                      tau, 
                                                      LR, 
                                                      ii, 
                                                      n_burnin, 
                                                      M_dense_sqrt, 
                                                      tau_m_adam,
                                                      tau_v_adam, 
                                                      tau_ii) {
  
  .Call(`_BayesMVP_fn_Rcpp_wrapper_update_tau_w_dense_M_ADAM`, eigen_vector, eigen_max, 
        theta_vec_initial, theta_vec_prop, snaper_m_vec, velocity_prop, velocity_0, 
        tau, LR, ii, n_burnin, M_dense_sqrt, tau_m_adam, tau_v_adam, tau_ii)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_det <- function(mat) {
  .Call(`_BayesMVP_Rcpp_det`, mat)
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_log_det <- function(mat) {
  .Call(`_BayesMVP_Rcpp_log_det`, mat)
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_solve <- function(mat) {
  .Call(`_BayesMVP_Rcpp_solve`, mat)
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_Chol <- function(mat) {
  .Call(`_BayesMVP_Rcpp_Chol`, mat)
}



#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_wrapper_fn_sample_HMC_multi_iter_single_thread <- function(chain_id, 
                                                                seed, 
                                                                n_iter, 
                                                                partitioned_HMC, 
                                                                Model_type, 
                                                                sample_nuisance, 
                                                                force_autodiff, 
                                                                force_PartialLog, 
                                                                multi_attempts,
                                                                n_nuisance_to_track,
                                                                theta_main_vector_from_single_chain_input_from_R,
                                                                theta_us_vector_from_single_chain_input_from_R,
                                                                y_Eigen_i,
                                                                Model_args_as_Rcpp_List,
                                                                EHMC_args_as_Rcpp_List,
                                                                EHMC_Metric_as_Rcpp_List) {
  
  .Call(`_BayesMVP_Rcpp_wrapper_fn_sample_HMC_multi_iter_single_thread`,
        chain_id, 
        seed, 
        n_iter, 
        partitioned_HMC, 
        Model_type, 
        sample_nuisance, 
        force_autodiff, 
        force_PartialLog, 
        multi_attempts, 
        n_nuisance_to_track, 
        theta_main_vector_from_single_chain_input_from_R, 
        theta_us_vector_from_single_chain_input_from_R, 
        y_Eigen_i, 
        Model_args_as_Rcpp_List, 
        EHMC_args_as_Rcpp_List,
        EHMC_Metric_as_Rcpp_List)
  
}



#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_compute_param_constrain_from_trace_parallel <- function(unc_params_trace_input_main, 
                                                           unc_params_trace_input_nuisance, 
                                                           pars_indicies_to_track, 
                                                           n_params_full, 
                                                           n_nuisance, 
                                                           n_params_main, 
                                                           include_nuisance, 
                                                           model_so_file, 
                                                           json_file_path) {
  
  .Call(`_BayesMVP_fn_compute_param_constrain_from_trace_parallel`, 
        unc_params_trace_input_main, 
        unc_params_trace_input_nuisance,
        pars_indicies_to_track, 
        n_params_full, 
        n_nuisance, 
        n_params_main,
        include_nuisance, 
        model_so_file, 
        json_file_path)
  
}


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
Rcpp_fn_RcppParallel_EHMC_sampling <- function(n_threads_R, 
                                               seed_R, 
                                               n_iter_R, 
                                               iter_one_by_one, 
                                               partitioned_HMC_R, 
                                               Model_type_R, 
                                               sample_nuisance_R, 
                                               force_autodiff_R, 
                                               force_PartialLog_R, 
                                               multi_attempts_R, 
                                               n_nuisance_to_track, 
                                               theta_main_vectors_all_chains_input_from_R, 
                                               theta_us_vectors_all_chains_input_from_R, 
                                               y_Eigen_R, 
                                               Model_args_as_Rcpp_List, 
                                               EHMC_args_as_Rcpp_List, 
                                               EHMC_Metric_as_Rcpp_List) {
  
  .Call(`_BayesMVP_Rcpp_fn_RcppParallel_EHMC_sampling`, n_threads_R, seed_R, n_iter_R, 
        iter_one_by_one, partitioned_HMC_R, Model_type_R, sample_nuisance_R, 
        force_autodiff_R, force_PartialLog_R, multi_attempts_R, n_nuisance_to_track, 
        theta_main_vectors_all_chains_input_from_R, theta_us_vectors_all_chains_input_from_R, 
        y_Eigen_R, Model_args_as_Rcpp_List, EHMC_args_as_Rcpp_List, EHMC_Metric_as_Rcpp_List)
  
}



 


#' @useDynLib BayesMVP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fn_R_RcppParallel_EHMC_single_iter_burnin <- function(n_threads_R,
                                                      seed_R, 
                                                      n_iter_R, 
                                                      current_iter_R,
                                                      n_adapt, 
                                                      burnin_indicator, 
                                                      Model_type_R, 
                                                      sample_nuisance_R, 
                                                      force_autodiff_R, 
                                                      force_PartialLog_R, 
                                                      multi_attempts_R, 
                                                      n_nuisance_to_track,
                                                      max_eps_main, 
                                                      max_eps_us, 
                                                      partitioned_HMC_R, 
                                                      metric_type_main, 
                                                      shrinkage_factor, 
                                                      metric_type_nuisance, 
                                                      tau_main_target, 
                                                      tau_us_target, 
                                                      clip_iter,
                                                      gap, 
                                                      main_L_manual,
                                                      us_L_manual, 
                                                      L_main_if_manual,
                                                      L_us_if_manual, 
                                                      max_L, 
                                                      tau_mult, 
                                                      ratio_M_us, 
                                                      ratio_Hess_main,
                                                      M_interval_width, 
                                                      theta_main_vectors_all_chains_input_from_R, 
                                                      theta_us_vectors_all_chains_input_from_R, 
                                                      y_Eigen_R, 
                                                      Model_args_as_Rcpp_List, 
                                                      EHMC_args_as_Rcpp_List, 
                                                      EHMC_Metric_as_Rcpp_List, 
                                                      EHMC_burnin_as_Rcpp_List) {
  
  .Call(`_BayesMVP_fn_R_RcppParallel_EHMC_single_iter_burnin`, n_threads_R, seed_R, 
        n_iter_R, current_iter_R, n_adapt, burnin_indicator, Model_type_R, sample_nuisance_R, 
        force_autodiff_R, force_PartialLog_R, multi_attempts_R, n_nuisance_to_track, 
        max_eps_main, max_eps_us, partitioned_HMC_R, metric_type_main, shrinkage_factor, 
        metric_type_nuisance, tau_main_target, tau_us_target, clip_iter, gap, 
        main_L_manual, us_L_manual, L_main_if_manual, L_us_if_manual, max_L, 
        tau_mult, ratio_M_us, ratio_Hess_main, M_interval_width, 
        theta_main_vectors_all_chains_input_from_R, theta_us_vectors_all_chains_input_from_R, 
        y_Eigen_R,
        Model_args_as_Rcpp_List,
        EHMC_args_as_Rcpp_List,
        EHMC_Metric_as_Rcpp_List,
        EHMC_burnin_as_Rcpp_List)

}




 









