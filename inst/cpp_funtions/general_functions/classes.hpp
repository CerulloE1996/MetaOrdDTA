#pragma once


class HMCResult {
           
         /////////// -------------- PRIVATE members (ONLY accessible WITHIN this class)
         private:
           //// Making states private prevents accidental misuse
           Eigen::Matrix<double, -1, 1> lp_and_grad_outs_;
           
           Eigen::Matrix<double, -1, 1> main_theta_vec_0_;
           Eigen::Matrix<double, -1, 1> main_theta_vec_;
           Eigen::Matrix<double, -1, 1> main_theta_vec_proposed_;
           Eigen::Matrix<double, -1, 1> main_velocity_0_vec_;
           Eigen::Matrix<double, -1, 1> main_velocity_vec_proposed_;
           Eigen::Matrix<double, -1, 1> main_velocity_vec_;
           
           double main_p_jump_;
           int main_div_;
           
           Eigen::Matrix<double, -1, 1> us_theta_vec_0_;
           Eigen::Matrix<double, -1, 1> us_theta_vec_;
           Eigen::Matrix<double, -1, 1> us_theta_vec_proposed_;
           Eigen::Matrix<double, -1, 1> us_velocity_0_vec_;
           Eigen::Matrix<double, -1, 1> us_velocity_vec_proposed_;
           Eigen::Matrix<double, -1, 1> us_velocity_vec_;
           
           double us_p_jump_;
           int us_div_;
           
           int N_;
           
         /////////// --------------- PUBLIC members (accessible from OUTSIDE this class [e.g. can do: "class.public_member"])
         public:
           //// Constructor 
           HMCResult(int n_params_main, 
                     int n_nuisance, 
                     int N)
           : lp_and_grad_outs_(Eigen::Matrix<double, -1, 1>::Zero((1 + N + n_params_main + n_nuisance)))
           //// main
           , main_theta_vec_0_(Eigen::Matrix<double, -1, 1>::Zero(n_params_main))
           , main_theta_vec_(Eigen::Matrix<double, -1, 1>::Zero(n_params_main))
           , main_theta_vec_proposed_(Eigen::Matrix<double, -1, 1>::Zero(n_params_main))
           , main_velocity_0_vec_(Eigen::Matrix<double, -1, 1>::Zero(n_params_main))
           , main_velocity_vec_proposed_(Eigen::Matrix<double, -1, 1>::Zero(n_params_main))
           , main_velocity_vec_(Eigen::Matrix<double, -1, 1>::Zero(n_params_main))
           , main_p_jump_(0.0)
           , main_div_(0)
           //// nuisance
           , us_theta_vec_0_(Eigen::Matrix<double, -1, 1>::Zero(n_nuisance))
           , us_theta_vec_(Eigen::Matrix<double, -1, 1>::Zero(n_nuisance))
           , us_theta_vec_proposed_(Eigen::Matrix<double, -1, 1>::Zero(n_nuisance))
           , us_velocity_0_vec_(Eigen::Matrix<double, -1, 1>::Zero(n_nuisance))
           , us_velocity_vec_proposed_(Eigen::Matrix<double, -1, 1>::Zero(n_nuisance))
           , us_velocity_vec_(Eigen::Matrix<double, -1, 1>::Zero(n_nuisance))
           , us_p_jump_(0.0)
           , us_div_(0),
           N_(N)
           {}
           
           //// "getters"/"setters" to access the private members
           Eigen::Matrix<double, -1, 1> &lp_and_grad_outs() { 
             return lp_and_grad_outs_; 
           }
           // const Eigen::Matrix<double, -1, 1> &lp_and_grad_outs() const { return lp_and_grad_outs_; }
           Eigen::Matrix<double, -1, 1> log_lik() { 
             return lp_and_grad_outs_.tail(N_);
           }
           ////
           //// main
           Eigen::Matrix<double, -1, 1> &main_theta_vec_0() { return main_theta_vec_0_; }
           //const Eigen::Matrix<double, -1, 1> &main_theta_vec_0() const { return main_theta_vec_0_; }
           Eigen::Matrix<double, -1, 1> &main_theta_vec() { return main_theta_vec_; }
           //const Eigen::Matrix<double, -1, 1> &main_theta_vec() const { return main_theta_vec_; }
           Eigen::Matrix<double, -1, 1> &main_theta_vec_proposed() { return main_theta_vec_proposed_; }
           //const Eigen::Matrix<double, -1, 1> &main_theta_vec_proposed() const { return main_theta_vec_proposed_; }
           Eigen::Matrix<double, -1, 1> &main_velocity_0_vec() { return main_velocity_0_vec_; }
           //const Eigen::Matrix<double, -1, 1> &main_velocity_0_vec() const { return main_velocity_0_vec_; }
           Eigen::Matrix<double, -1, 1> &main_velocity_vec_proposed() { return main_velocity_vec_proposed_; }
           //const Eigen::Matrix<double, -1, 1> &main_velocity_vec_proposed() const { return main_velocity_vec_proposed_; }
           Eigen::Matrix<double, -1, 1> &main_velocity_vec() { return main_velocity_vec_; }
           //const Eigen::Matrix<double, -1, 1> &main_velocity_vec() const { return main_velocity_vec_; }
           double &main_p_jump() { return main_p_jump_; }
           //const double &main_p_jump() const { return main_p_jump_; }
           int &main_div() { return main_div_; }
           //const int &main_div() const { return main_div_; }
           ////
           //// nuisance
           Eigen::Matrix<double, -1, 1> &us_theta_vec_0() { return us_theta_vec_0_; }
           //const Eigen::Matrix<double, -1, 1> &us_theta_vec_0() const { return us_theta_vec_0_; }
           Eigen::Matrix<double, -1, 1> &us_theta_vec() { return us_theta_vec_; }
           //const Eigen::Matrix<double, -1, 1> &us_theta_vec() const { return us_theta_vec_; }
           Eigen::Matrix<double, -1, 1> &us_theta_vec_proposed() { return us_theta_vec_proposed_; }
           //const Eigen::Matrix<double, -1, 1> &us_theta_vec_proposed() const { return us_theta_vec_proposed_; }
           Eigen::Matrix<double, -1, 1> &us_velocity_0_vec() { return us_velocity_0_vec_; }
           //const Eigen::Matrix<double, -1, 1> &us_velocity_0_vec() const { return us_velocity_0_vec_; }
           Eigen::Matrix<double, -1, 1> &us_velocity_vec_proposed() { return us_velocity_vec_proposed_; }
           //const Eigen::Matrix<double, -1, 1> &us_velocity_vec_proposed() const { return us_velocity_vec_proposed_; }
           Eigen::Matrix<double, -1, 1> &us_velocity_vec() { return us_velocity_vec_; }
           //const Eigen::Matrix<double, -1, 1> &us_velocity_vec() const { return us_velocity_vec_; }
           double &us_p_jump() { return us_p_jump_; }
           //const double &us_p_jump() const { return us_p_jump_; }
           int &us_div() { return us_div_; }
           //const int &us_div() const { return us_div_; }
           
           
           /////////// --------------  PUBLIC HELPER FNS - local to this class 
           //// HMC helper fns
           void store_current_state() {
             main_theta_vec_0_ = main_theta_vec_;
             us_theta_vec_0_ = us_theta_vec_;
             main_velocity_0_vec_ = main_velocity_vec_;
             us_velocity_0_vec_ = us_velocity_vec_;
           }  
           
           void accept_proposal_main() {
             main_theta_vec_ = main_theta_vec_proposed_;
             main_velocity_vec_ = main_velocity_vec_proposed_;
           }  
           
           void accept_proposal_us() {
             us_theta_vec_ = us_theta_vec_proposed_;
             us_velocity_vec_ = us_velocity_vec_proposed_;
           }  
           
           void reject_proposal_main() {
             main_theta_vec_ = main_theta_vec_0_;
             main_velocity_vec_ = main_velocity_0_vec_;
             //// lp_and_grad_outs_.setZero();
           } 
            
           void reject_proposal_us() {
             us_theta_vec_ = us_theta_vec_0_;
             us_velocity_vec_ = us_velocity_0_vec_;
             //// lp_and_grad_outs_.setZero();
           }  
           
           //// fn to check state validity
           bool check_state_valid() const {
               return  main_theta_vec_.allFinite() && 
                       main_velocity_vec_.allFinite() && 
                       us_theta_vec_.allFinite() && 
                       us_velocity_vec_.allFinite();
           }
            

   
};














class HMC_output_single_chain {
  
       /////////// -------------- PRIVATE members (ONLY accessible WITHIN this class)
       private:
          // Make result_input private since it's internal state
          HMCResult result_input_;
          
          // Internal trace storage structures
          struct TraceBuffers {
            Eigen::Matrix<double, -1, -1> main;
            Eigen::Matrix<double, -1, -1> div;
            Eigen::Matrix<double, -1, -1> nuisance;
            Eigen::Matrix<double, -1, -1> log_lik;
            
            TraceBuffers(int n_params_main, int n_iter, int n_nuisance_to_track, int N) 
              : main(Eigen::Matrix<double, -1, -1>::Zero(n_params_main, n_iter))
              , div(Eigen::Matrix<double, -1, -1>::Zero(1, n_iter))
              , nuisance(Eigen::Matrix<double, -1, -1>::Zero(n_nuisance_to_track, n_iter))
              , log_lik(Eigen::Matrix<double, -1, -1>::Zero(N, n_iter)) 
              {}
          };
          
          struct DiagnosticBuffers {
            Eigen::Matrix<int, -1, 1> div_us;
            Eigen::Matrix<int, -1, 1> div_main;
            Eigen::Matrix<double, -1, 1> p_jump_us;
            Eigen::Matrix<double, -1, 1> p_jump_main;
            
            DiagnosticBuffers(int n_iter)
              : div_us(Eigen::Matrix<int, -1, 1>::Zero(n_iter))
              , div_main(Eigen::Matrix<int, -1, 1>::Zero(n_iter))
              , p_jump_us(Eigen::Matrix<double, -1, 1>::Zero(n_iter))
              , p_jump_main(Eigen::Matrix<double, -1, 1>::Zero(n_iter)) 
              {}
          };
          
          TraceBuffers traces_;
          DiagnosticBuffers diagnostics_;
          
        /////////// --------------- PUBLIC members (accessible from OUTSIDE this class [e.g. can do: "class.public_member"])
        public:
          //// Constructor
          HMC_output_single_chain(int n_iter,
                                  int n_nuisance_to_track,
                                  int n_params_main,
                                  int n_nuisance,
                                  int N)
            : 
            result_input_(n_params_main, n_nuisance, N)
          , traces_(n_params_main, n_iter, n_nuisance_to_track, N)
          , diagnostics_(n_iter) 
          {}
          
          //// Getters for result_input
          HMCResult &result_input() { return result_input_; }
          //const HMCResult &result_input() const { return result_input_; }
          
          //// Getters for traces
          Eigen::Matrix<double, -1, -1> &trace_main() { return traces_.main; }
          //const Eigen::Matrix<double, -1, -1> &trace_main() const { return traces_.main; }
          Eigen::Matrix<double, -1, -1> &trace_div() { return traces_.div; }
          //const Eigen::Matrix<double, -1, -1> &trace_div() const { return traces_.div; }
          Eigen::Matrix<double, -1, -1> &trace_nuisance() { return traces_.nuisance; }
          //const Eigen::Matrix<double, -1, -1> &trace_nuisance() const { return traces_.nuisance; }
          Eigen::Matrix<double, -1, -1> &trace_log_lik() { return traces_.log_lik; }
          //const Eigen::Matrix<double, -1, -1> &trace_log_lik() const { return traces_.log_lik; }
          
          //// Getters for diagnostics
          Eigen::Matrix<int, -1, 1> &diagnostics_div_us() { return diagnostics_.div_us; }
          //const Eigen::Matrix<int, -1, 1> &diagnostics_div_us() const { return diagnostics_.div_us; }
          Eigen::Matrix<int, -1, 1> &diagnostics_div_main() { return diagnostics_.div_main; }
          //const Eigen::Matrix<int, -1, 1> &diagnostics_div_main() const { return diagnostics_.div_main; }
          Eigen::Matrix<double, -1, 1> &diagnostics_p_jump_us() { return diagnostics_.p_jump_us; }
          //const Eigen::Matrix<double, -1, 1> &diagnostics_p_jump_us() const { return diagnostics_.p_jump_us; }
          Eigen::Matrix<double, -1, 1> &diagnostics_p_jump_main() { return diagnostics_.p_jump_main; }
          //const Eigen::Matrix<double, -1, 1> &diagnostics_p_jump_main() const { return diagnostics_.p_jump_main; }
          
          
          /////////// --------------  PUBLIC HELPER FNS - local to this class 
          
          //// Helper function to store current iteration results
          void store_iteration(int ii, 
                               bool sample_nuisance) {
            
              trace_main().col(ii) = result_input().main_theta_vec();
              
              if (sample_nuisance == true) {
                trace_div()(0, ii) = 0.5 * (result_input().main_div() + result_input().us_div());
                trace_nuisance().col(ii) = result_input().us_theta_vec();
              } else {
                trace_div()(0, ii) = result_input().main_div();
              }
              
              //// Store diagnostics
              diagnostics_div_us()(ii) = result_input().us_div();
              diagnostics_div_main()(ii) = result_input().main_div();
              diagnostics_p_jump_us()(ii) = result_input().us_p_jump();
              diagnostics_p_jump_main()(ii) = result_input().main_p_jump();
              
              //// Store log-lik 
              trace_log_lik().col(ii) = result_input().log_lik();
            
          }
          
          
          //// Function to check if storage is valid
          bool check_storage_valid() const {
            
              return traces_.main.allFinite() &&
                traces_.div.allFinite() &&
                traces_.nuisance.allFinite() &&
                traces_.log_lik.allFinite() &&
                diagnostics_.p_jump_us.allFinite() &&
                diagnostics_.p_jump_main.allFinite();
            
          }
  
};








