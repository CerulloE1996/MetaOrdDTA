#ifndef STAN_MATH_PRIM_META_MVP_COLVEC_AND_ARRAY_FN_WRAPPERS
#define STAN_MATH_PRIM_META_MVP_COLVEC_AND_ARRAY_FN_WRAPPERS

 
 
 
 
 
// #include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/prob/std_normal_log_qf.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <stan/math/prim/fun/inv_Phi.hpp>
#include <stan/math/prim/fun/Phi_approx.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/fun/log_inv_logit.hpp>
 
#include <Eigen/Dense>
#include <Eigen/Core>
 
#include <immintrin.h>

 

 
using namespace Eigen;


 

 
// #if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2 && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2

//  //////////////////  for-loop and SIMD (for AVX2 and  AVX512) wrappers -  double -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//

typedef double (*FuncDouble)(double);
typedef double (*FuncDouble_wo_checks)(double);

#if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2
typedef __m256d (*FuncAVX2)(const __m256d);
typedef __m256d (*FuncAVX2_wo_checks)(const __m256d);
#endif


#if defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)
typedef __m512d (*FuncAVX512)(const __m512d);
typedef __m512d (*FuncAVX512_wo_checks)(const __m512d);
#endif




#if defined(__AVX2__) && !defined(__AVX512F__) // Only use if AVX2 is available and AVX512 is not   //// for AVX2


template <typename T, typename FuncAVX2, typename FuncDouble>
inline void fn_AVX2_row_or_col_vector(Eigen::Ref<T> x, 
                                     FuncAVX2 fn_AVX2, 
                                     FuncDouble fn_double) {
  
  const int N = x.size();
  const int N_divisible_by_4 = (N / 4) * 4;
  
  for (int i = 0; i + 4 <= N_divisible_by_4; i += 4) {
    __m256d AVX_array = _mm256_load_pd(&x(i));
    __m256d AVX_array_out = fn_AVX2(AVX_array);
    _mm256_store_pd(&x(i), AVX_array_out);
  }
  
  // Process remaining elements
  for (int i = N_divisible_by_4; i < N; ++i) {
    x(i) = fn_double(x(i));
  }
  
}


template <typename T, typename FuncAVX2, typename FuncDouble>
inline void fn_AVX2_matrix(Eigen::Ref<T> x, 
                           FuncAVX2 fn_AVX2,
                           FuncDouble fn_double) {
  
  const int rows = x.rows();
  const int cols = x.cols();
  
  if (cols > rows) {     // Wide format: Process row-by-row
    
    for (int i = 0; i < rows; ++i) {
      T row = x.row(i);
      fn_AVX2_row_vector(row, fn_AVX512, fn_double);
    }  
    
  } else { // Long format: Process column-by-column 
    
    for (int j = 0; j < cols; ++j) {
      T col = x.col(j);
      fn_AVX2_col_vector(col, fn_AVX512, fn_double);
    } 
    
  }
  
} 

 
 

template <typename T, typename FuncAVX2, typename FuncDouble>
inline void fn_AVX2_dbl_Eigen(Eigen::Ref<T> x, 
                              FuncAVX2 fn_AVX2, 
                              FuncDouble fn_double) {
  
     constexpr int n_rows = T::RowsAtCompileTime;
     constexpr int n_cols = T::ColsAtCompileTime;
     
     if constexpr (n_rows == 1 && n_cols == -1) {  // Row vector case
      
       fn_AVX2_row_or_col_vector(x, fn_AVX2, fn_double);
       
     } else if constexpr (n_rows == -1 && n_cols == 1) {    // Column vector case
    
       fn_AVX2_row_or_col_vector(x, fn_AVX2, fn_double);
       
     } else {        // General matrix case

       fn_AVX2_matrix(x, fn_AVX2, fn_double);
       
     }
     
}

 
#endif
 





#if defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__) //// for AVX-512
 

template <typename T, typename FuncAVX512, typename FuncDouble>
inline void fn_AVX512_row_or_col_vector(   Eigen::Ref<T>  x, 
                                           FuncAVX512 fn_AVX512, 
                                           FuncDouble fn_double) {
  
  
  const int N = x.size();
  const int vect_size = 8;
  const double vect_siz_dbl = 8.0;
  const int N_divisible_by_8 = std::floor( static_cast<double>(N) / vect_siz_dbl) * vect_size;
 
  T x_temp = x; // make a copy 
  
  if (N >= vect_size) {
    
          for (int i = 0; i + 8 <= N_divisible_by_8; i += vect_size) {
            const __m512d AVX_array = _mm512_loadu_pd(&x(i));
            const __m512d AVX_array_out = fn_AVX512(AVX_array);
            _mm512_storeu_pd(&x_temp(i), AVX_array_out);
          }
          
          // if (N_divisible_by_8 != N) {    // Handle remainder 
          //    const Eigen::Matrix<double, -1, 1> x_tail = x.tail(vect_size); // copy of last 8 elements 
             const int start_index = N - vect_size;
             const int end_index = N;
               for (int i = start_index; i < end_index; ++i) {
                      ///  x(i) = fn_double(x_tail(i - start_index));
                      x_temp(i) = fn_double(x(i));
               }
           // }
           
  }  else {   // If N < 8, handle everything with scalar operations
    
        for (int i = 0; i < N; ++i) {
          x_temp(i) = fn_double(x(i));
        }
        
  }
  
  x = x_temp;

  
}

 
 
 
template<typename T, typename FuncAVX512, typename FuncDouble>
inline void fn_AVX512_matrix(  Eigen::Ref<T> x,
                               FuncAVX512 fn_AVX512,
                               FuncDouble fn_double) {
  
   const int n_rows = x.rows();
   const int n_cols = x.cols();
   const int vect_size = 8;
   const double vect_siz_dbl = 8.0;
   const int rows_divisible_by_8 = std::floor( static_cast<double>(n_rows) / vect_siz_dbl) * vect_size;
 
   T x_temp = x; // make a copy 
   
   for (int j = 0; j < n_cols; ++j) { /// loop through cols first as col-major storage

        // Make sure we have at least 8 rows before trying AVX
        if (n_rows >= vect_size) {
          
              for (int i = 0; i < rows_divisible_by_8; i += vect_size) {
                const __m512d AVX_array = _mm512_loadu_pd(&x(i, j));
                const __m512d AVX_array_out = fn_AVX512(AVX_array);
                _mm512_storeu_pd(&x_temp(i, j), AVX_array_out);
              }
              
              // Handle remaining rows with double fns
            // if (rows_divisible_by_8 != n_rows) {
            //   const Eigen::Matrix<double, -1, 1> x_tail = x.col(j).tail(vect_size); // copy of last 8 elements 
              const int start_index = n_rows - vect_size;
              const int end_index = n_rows;
                for (int i = start_index; i < end_index; ++i) {
                     /// x(i, j) = fn_double(x_tail(i - start_index));
                    x_temp(i, j) = fn_double(x(i, j));
                }
            // }
              
        } else {    // If n_rows < 8, handle entire row with double operations
          for (int i = 0; i < n_rows; ++i) {
            x_temp(i, j) = fn_double(x(i, j));
          } 
        }

  }
   
    x = x_temp;
   

}

 


template <typename T, typename FuncAVX512, typename FuncDouble>
inline void fn_AVX512_dbl_Eigen(Eigen::Ref<T> x, 
                                FuncAVX512 fn_AVX512, 
                                FuncDouble fn_double) {
  
    constexpr int n_rows = T::RowsAtCompileTime;
    constexpr int n_cols = T::ColsAtCompileTime;
    
    if constexpr (n_rows == 1 && n_cols == -1) {   // Row vector case
    
      fn_AVX512_row_or_col_vector(x, fn_AVX512, fn_double);
      
    } else if constexpr (n_rows == -1 && n_cols == 1) {  // Column vector case
     
      fn_AVX512_row_or_col_vector(x, fn_AVX512, fn_double);
      
    } else {   // General matrix case
    
      fn_AVX512_matrix(x, fn_AVX512, fn_double);
      
    }
  
}



#endif
 

 
 
 
 
 
 
 
template <typename T,  typename FuncDouble>
inline void fn_Loop_row_or_col_vector(      Eigen::Ref<T>  x, 
                                            FuncDouble fn_double) {
  
   const int N = x.size();
   
     for (int i = 0; i < N; ++i) {
       x(i) = fn_double(x(i));
     }
   
}
 
 


 
 
template<typename T, typename FuncDouble>
inline void fn_Loop_matrix(  Eigen::Ref<T> x,
                             FuncDouble fn_double) {
     
     const int rows = x.rows(); 
     const int cols = x.cols();
      
     for (int j = 0; j < cols; ++j) {
           for (int i = 0; i < rows; ++i) {
             x(i, j) = fn_double(x(i, j));
           }
     }
 
}

  
  
  
  
template <typename T, typename FuncDouble>
inline void fn_Loop_dbl_Eigen( Eigen::Ref<T> x, 
                               FuncDouble fn_double) {
  
     constexpr int n_rows = T::RowsAtCompileTime;
     constexpr int n_cols = T::ColsAtCompileTime;
     
     if constexpr (n_rows == 1 && n_cols == -1) {   // Row vector case
       
       fn_Loop_row_or_col_vector(x, fn_double);
       
     } else if constexpr (n_rows == -1 && n_cols == 1) {  // Column vector case
       
       fn_Loop_row_or_col_vector(x, fn_double);
       
     } else {   // General matrix case
       
       fn_Loop_matrix(x, fn_double);
       
     }
 
}

 
  
 
 
 
 
 

///////////////// for-loop and SIMD (for AVX2 and  AVX512) wrappers - "master" functions -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


 
 
 
#if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2 
 
template<typename FuncAVX2, typename FuncDouble, typename FuncAVX2_wo_checks, typename FuncDouble_wo_checks, typename T>
inline void    fn_process_double_AVX2_sub_function(     Eigen::Ref<T> x,  
                                                        FuncAVX2 fn_fast_AVX2_function,
                                                        FuncDouble fn_fast_double_function,
                                                        FuncAVX2_wo_checks fn_fast_AVX2_function_wo_checks,
                                                        FuncDouble_wo_checks fn_fast_double_function_wo_checks, 
                                                        const bool skip_checks) {
  
  if (skip_checks == false) {
    
    fn_AVX2_dbl_Eigen(x, fn_fast_AVX2_function, fn_fast_double_function);
    
  }   else  {
    
    fn_AVX2_dbl_Eigen(x, fn_fast_AVX2_function_wo_checks, fn_fast_double_function_wo_checks);
    
  }
  
}

#endif





#if defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)  /// use AVX-512 

template<typename FuncAVX512, typename FuncDouble, typename FuncAVX512_wo_checks, typename FuncDouble_wo_checks, typename T>
inline void                   fn_process_double_AVX512_sub_function(     Eigen::Ref<T> x, // since this is helper function we call x by reference "&" not "&&"
                                                                         FuncAVX512 fn_fast_AVX512_function,
                                                                         FuncDouble fn_fast_double_function,
                                                                         FuncAVX512_wo_checks fn_fast_AVX512_function_wo_checks,
                                                                         FuncDouble_wo_checks fn_fast_double_function_wo_checks, 
                                                                         const bool skip_checks) {
  
  if (skip_checks == false) {
    
    fn_AVX512_dbl_Eigen(x, fn_fast_AVX512_function, fn_fast_double_function);
    
  }   else  {
    
    fn_AVX512_dbl_Eigen(x, fn_fast_AVX512_function_wo_checks, fn_fast_double_function_wo_checks);
 
  }
  
}


#endif



template<typename FuncDouble, typename FuncDouble_wo_checks, typename T>
inline void    fn_process_double_Loop_sub_function(     Eigen::Ref<T> x,  
                                                        FuncDouble fn_fast_double_function,
                                                        FuncDouble_wo_checks fn_fast_double_function_wo_checks, 
                                                        const bool skip_checks) {
  
  if (skip_checks == false) {
    
    fn_Loop_dbl_Eigen(x, fn_fast_double_function);
    
  }   else  {
    
    fn_Loop_dbl_Eigen(x, fn_fast_double_function_wo_checks);
    
  }
  
}






inline double mvp_std_exp(double x) { 
     return stan::math::exp(x);
}
inline double mvp_std_log(double x) { 
  return stan::math::log(x);
}
inline double mvp_std_log1p(double x) { 
  return stan::math::log1p(x);
}
inline double mvp_std_log1m(double x) { 
  return stan::math::log1m(x);
} 
inline double mvp_std_logit(double x) { 
  return stan::math::logit(x);
} 
inline double mvp_std_tanh(double x) { 
  return stan::math::tanh(x);
} 
inline double mvp_std_Phi_approx(double x) { 
  return stan::math::Phi_approx(x);
}
inline double mvp_std_Phi(double x) { 
  return stan::math::Phi(x);
}
inline double mvp_std_inv_Phi(double x) { 
  return stan::math::inv_Phi(x);
}
inline double mvp_std_inv_logit(double x) { 
  return stan::math::inv_logit(x);
}
inline double mvp_std_log_inv_logit(double x) { 
  return stan::math::log_inv_logit(x); 
}



 
#if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2
template <typename T>
inline  void       fn_return_Ref_double_AVX2( Eigen::Ref<T> x,
                                              const std::string &fn,
                                              const bool &skip_checks) {
 
 
    if        (fn == "exp") {                               fn_process_double_AVX2_sub_function(x, fast_exp_1_AVX2,  mvp_std_exp,   fast_exp_1_wo_checks_AVX2, mvp_std_exp, skip_checks) ;
    } else if (fn == "log") {                               fn_process_double_AVX2_sub_function(x, fast_log_1_AVX2, mvp_std_log, fast_log_1_wo_checks_AVX2, mvp_std_log, skip_checks) ;
    } else if (fn == "log1p") {                             fn_process_double_AVX2_sub_function(x, fast_log1p_1_AVX2, mvp_std_log1p, fast_log1p_1_wo_checks_AVX2, mvp_std_log1p, skip_checks) ;
    } else if (fn == "log1m") {                             fn_process_double_AVX2_sub_function(x, fast_log1m_1_AVX2, mvp_std_log1m, fast_log1m_1_wo_checks_AVX2, mvp_std_log1m, skip_checks) ;
    } else if (fn == "logit") {                             fn_process_double_AVX2_sub_function(x, fast_logit_AVX2, mvp_std_logit, fast_logit_wo_checks_AVX2, mvp_std_logit, skip_checks) ;
    } else if (fn == "tanh") {                              fn_process_double_AVX2_sub_function(x, fast_tanh_AVX2, mvp_std_tanh, fast_tanh_wo_checks_AVX2, mvp_std_tanh, skip_checks) ;
    } else if (fn == "Phi_approx") {                        fn_process_double_AVX2_sub_function(x, fast_Phi_approx_AVX2, mvp_std_Phi_approx, fast_Phi_approx_wo_checks_AVX2, mvp_std_Phi_approx, skip_checks) ;
    } else if (fn == "log_Phi_approx") {                    fn_process_double_AVX2_sub_function(x, fast_log_Phi_approx_AVX2, fast_log_Phi_approx, fast_log_Phi_approx_wo_checks_AVX2, fast_log_Phi_approx_wo_checks, skip_checks) ;
    } else if (fn == "inv_Phi_approx") {                    fn_process_double_AVX2_sub_function(x, fast_inv_Phi_approx_AVX2, fast_inv_Phi_approx, fast_inv_Phi_approx_wo_checks_AVX2, fast_inv_Phi_approx_wo_checks, skip_checks) ;
    } else if (fn == "inv_Phi_approx_from_logit_prob") {    fn_process_double_AVX2_sub_function(x, fast_inv_Phi_approx_from_logit_prob_AVX2, fast_inv_Phi_approx_from_logit_prob, fast_inv_Phi_approx_from_logit_prob_wo_checks_AVX2, fast_inv_Phi_approx_from_logit_prob_wo_checks, skip_checks) ;
    } else if (fn == "Phi") {                               fn_process_double_AVX2_sub_function(x, fast_Phi_AVX2, mvp_std_Phi, fast_Phi_wo_checks_AVX2, mvp_std_Phi, skip_checks) ;
    } else if (fn == "inv_Phi") {                           fn_process_double_AVX2_sub_function(x, fast_inv_Phi_wo_checks_AVX2, mvp_std_inv_Phi, fast_inv_Phi_wo_checks_AVX2, mvp_std_inv_Phi, skip_checks) ;
    } else if (fn == "inv_logit") {                         fn_process_double_AVX2_sub_function(x, fast_inv_logit_AVX2, mvp_std_inv_logit, fast_inv_logit_wo_checks_AVX2, mvp_std_inv_logit, skip_checks) ;
    } else if (fn == "log_inv_logit") {                     fn_process_double_AVX2_sub_function(x, fast_log_inv_logit_AVX2, fast_log_inv_logit, fast_log_inv_logit_wo_checks_AVX2, fast_log_inv_logit_wo_checks, skip_checks) ;
    }
 

}
#endif





#if  defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)
template <typename T>
inline   void        fn_return_Ref_double_AVX512(  Eigen::Ref<T> x,
                                                   const std::string &fn,
                                                   const bool &skip_checks) {
  
  if        (fn == "exp") {                              fn_process_double_AVX512_sub_function(x, fast_exp_1_AVX512, mvp_std_exp, fast_exp_1_wo_checks_AVX512, mvp_std_exp, skip_checks) ;
  } else if (fn == "log") {                              fn_process_double_AVX512_sub_function(x, fast_log_1_AVX512, mvp_std_log, fast_log_1_wo_checks_AVX512, mvp_std_log, skip_checks) ;
  } else if (fn == "log1p") {                            fn_process_double_AVX512_sub_function(x, fast_log1p_1_AVX512, mvp_std_log1p, fast_log1p_1_wo_checks_AVX512, mvp_std_log1p, skip_checks) ;
  } else if (fn == "log1m") {                            fn_process_double_AVX512_sub_function(x, fast_log1m_1_AVX512, mvp_std_log1m, fast_log1m_1_wo_checks_AVX512, mvp_std_log1m, skip_checks) ;
  } else if (fn == "logit") {                            fn_process_double_AVX512_sub_function(x, fast_logit_AVX512, mvp_std_logit, fast_logit_wo_checks_AVX512, mvp_std_logit, skip_checks) ;
  } else if (fn == "tanh") {                             fn_process_double_AVX512_sub_function(x, fast_tanh_AVX512, mvp_std_tanh, fast_tanh_wo_checks_AVX512, mvp_std_tanh, skip_checks) ;
  } else if (fn == "Phi_approx") {                       fn_process_double_AVX512_sub_function(x, fast_Phi_approx_AVX512, mvp_std_Phi_approx, fast_Phi_approx_wo_checks_AVX512, mvp_std_Phi_approx, skip_checks) ;
  } else if (fn == "log_Phi_approx") {                   fn_process_double_AVX512_sub_function(x, fast_log_Phi_approx_AVX512, fast_log_Phi_approx, fast_log_Phi_approx_wo_checks_AVX512, fast_log_Phi_approx_wo_checks, skip_checks) ;
  } else if (fn == "inv_Phi_approx") {                   fn_process_double_AVX512_sub_function(x, fast_inv_Phi_approx_AVX512, fast_inv_Phi_approx, fast_inv_Phi_approx_wo_checks_AVX512, fast_inv_Phi_approx_wo_checks, skip_checks) ;
  } else if (fn == "inv_Phi_approx_from_logit_prob") {   fn_process_double_AVX512_sub_function(x, fast_inv_Phi_approx_from_logit_prob_AVX512, fast_inv_Phi_approx_from_logit_prob, fast_inv_Phi_approx_from_logit_prob_wo_checks_AVX512, fast_inv_Phi_approx_from_logit_prob_wo_checks, skip_checks) ;
  } else if (fn == "Phi") {                              fn_process_double_AVX512_sub_function(x, fast_Phi_AVX512, mvp_std_Phi, fast_Phi_wo_checks_AVX512, mvp_std_Phi, skip_checks) ;
  } else if (fn == "inv_Phi") {                          fn_process_double_AVX512_sub_function(x, fast_inv_Phi_wo_checks_AVX512, mvp_std_inv_Phi, fast_inv_Phi_wo_checks_AVX512, mvp_std_inv_Phi, skip_checks) ;
  } else if (fn == "inv_logit") {                        fn_process_double_AVX512_sub_function(x, fast_inv_logit_AVX512, mvp_std_inv_logit, fast_inv_logit_wo_checks_AVX512, mvp_std_inv_logit, skip_checks) ;
  } else if (fn == "log_inv_logit") {                    fn_process_double_AVX512_sub_function(x, fast_log_inv_logit_AVX512, mvp_std_log_inv_logit, fast_log_inv_logit_wo_checks_AVX512, mvp_std_log_inv_logit, skip_checks) ;
  }
  
  
}
#endif




  
template <typename T>
inline      void          fn_void_Ref_double_Stan(    Eigen::Ref<T> x,
                                                      const std::string &fn,
                                                      const bool &skip_checks) {


 if        (fn == "exp")   {   x =   stan::math::exp(x);
 } else if (fn == "log")   {   x =   stan::math::log(x);
 } else if (fn == "log1p") {   x =   stan::math::log1p(x);
 } else if (fn == "log1m") {   x =   stan::math::log1m(x);
 } else if (fn == "logit") {   x =   stan::math::logit(x);
 } else if (fn == "tanh")  {   x =   stan::math::tanh(x);
 } else if (fn == "Phi_approx") {      x =    stan::math::Phi_approx(x);
 } else if (fn == "log_Phi_approx") {    x =   stan::math::log_inv_logit((0.07056*x.array().square()*x.array()  +  1.5976*x.array()).matrix());  
 } else if (fn == "inv_Phi_approx") {  x =    inv_Phi_approx_Stan(x);
 } else if (fn == "Phi") {             x =    stan::math::Phi(x);
 } else if (fn == "inv_Phi") {         x =    stan::math::inv_Phi(x);
 } else if (fn == "inv_Phi_approx_from_logit_prob") {  x =  inv_Phi_approx_from_logit_prob_Stan(x);
 } else if (fn == "inv_Phi_from_log_prob") {  x =   stan::math::std_normal_log_qf(x.matrix());
 } else if (fn == "inv_logit") {  x =   stan::math::inv_logit(x.matrix());
 } else if (fn == "log_inv_logit") {  x =   stan::math::log_inv_logit(x.matrix());
 } else {

 }


}



template <typename T>
inline   void        fn_return_Loop_AVX512(  Eigen::Ref<T> x,
                                             const std::string &fn,
                                             const bool &skip_checks) {
  
      if        (fn == "exp") {                              fn_process_double_Loop_sub_function(x, fast_exp_1, fast_exp_1_wo_checks, skip_checks) ;
      } else if (fn == "log") {                              fn_process_double_Loop_sub_function(x, fast_log_1, fast_log_1_wo_checks, skip_checks) ;
      } else if (fn == "log1p") {                            fn_process_double_Loop_sub_function(x, fast_log1p_1, fast_log1p_1_wo_checks, skip_checks) ;
      } else if (fn == "log1m") {                            fn_process_double_Loop_sub_function(x, fast_log1m_1, fast_log1m_1_wo_checks, skip_checks) ;
      } else if (fn == "tanh") {                             fn_process_double_Loop_sub_function(x, fast_tanh, fast_tanh_wo_checks, skip_checks) ;
      } else if (fn == "Phi_approx") {                       fn_process_double_Loop_sub_function(x, fast_Phi_approx, fast_Phi_approx_wo_checks, skip_checks) ;
      } else if (fn == "log_Phi_approx") {                   fn_process_double_Loop_sub_function(x, fast_log_Phi_approx, fast_log_Phi_approx_wo_checks, skip_checks) ;
      } else if (fn == "inv_Phi_approx") {                   fn_process_double_Loop_sub_function(x, fast_inv_Phi_approx, fast_inv_Phi_approx_wo_checks, skip_checks) ;
      } else if (fn == "inv_Phi_approx_from_logit_prob") {   fn_process_double_Loop_sub_function(x, fast_inv_Phi_approx_from_logit_prob, fast_inv_Phi_approx_from_logit_prob_wo_checks, skip_checks) ;
      } else if (fn == "Phi") {                              fn_process_double_Loop_sub_function(x, fast_Phi, fast_Phi_wo_checks, skip_checks) ;
      } else if (fn == "inv_Phi") {                          fn_process_double_Loop_sub_function(x, fast_inv_Phi_wo_checks, fast_inv_Phi_wo_checks, skip_checks) ;
      } else if (fn == "inv_logit") {                        fn_process_double_Loop_sub_function(x, fast_inv_logit, fast_inv_logit_wo_checks, skip_checks) ;
      } else if (fn == "log_inv_logit") {                    fn_process_double_Loop_sub_function(x, fast_log_inv_logit, fast_log_inv_logit_wo_checks, skip_checks) ;
      }
  
}





 
///// ----------------------------------------------------------------------------- Colvec function callers / wrappers

 


template <typename T>
inline  void               fn_EIGEN_Ref_double(         Eigen::Ref<T> x,
                                                        const std::string &fn,
                                                        const std::string &vect_type,
                                                        const bool &skip_checks) {
   
   // T  x_copy = x; /// make a copy (for debug)
  
 /////  stan::math::check_finite(fn.c_str(), "x", x);  // using c_str() to convert std::string to const char*

  if (fn == "inv_Phi_from_log_prob") {

      x = stan::math::std_normal_log_qf(x); 
      //return x;

  } else {

    if (vect_type == "Stan") {

         fn_void_Ref_double_Stan(Eigen::Ref<T>(x), fn, skip_checks);
        // return x;
 
    } else if (vect_type == "AVX2" ) { // use AVX-512 or AVX2 or loop (i.e., rely on automatic vectorisation)

#if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2
           fn_return_Ref_double_AVX2(Eigen::Ref<T>(x), fn, skip_checks);
          // return x;
#endif

    } else if (vect_type == "AVX512" ) { // use AVX-512 or AVX2 or loop (i.e., rely on automatic vectorisation)

#if defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)
          fn_return_Ref_double_AVX512(Eigen::Ref<T>(x), fn, skip_checks);
         // return x;
#endif

    } else {
      
          fn_return_Loop_AVX512(Eigen::Ref<T>(x), fn, skip_checks);
 
        // throw std::invalid_argument( os.str() ); /// note: std::invalid_argument doesnt seem to work w/ Stan math lib

       //return x;

    }

  }


  //  return x;

}


 

 
 
 ////// ------- "master" function  w/ return  ------------------------------------------------------------- 
 
 
 //// R-value
 template <typename T>
 inline  auto          fn_EIGEN_double(              T  &&x_R_val,
                                                     const std::string &fn,
                                                     const std::string &vect_type = "Stan",
                                                     const bool &skip_checks = false) {
   
   using matrix_type = Eigen::Matrix<double, -1, -1>;
   matrix_type x_matrix = x_R_val;   
   fn_EIGEN_Ref_double(Eigen::Ref<matrix_type>(x_matrix), fn, vect_type, skip_checks);
   return x_matrix;
   
 }
 
 //// Eigen Ref (this will also accept L_value [&T] as well as other types)
 template <typename T>
 inline  auto          fn_EIGEN_double(  Eigen::Ref<T> x_L_val,
                                         const std::string &fn,
                                         const std::string &vect_type = "Stan",
                                         const bool &skip_checks = false) {
   
   fn_EIGEN_Ref_double(x_L_val, fn, vect_type, skip_checks);
   return x_L_val;
   
 }
 
 //// const Eigen Ref (this will also accept L_value [&T] as well as other types)
 template <typename T>
 inline  auto          fn_EIGEN_double(  const Eigen::Ref<const T> x_L_val,
                                         const std::string &fn,
                                         const std::string &vect_type = "Stan",
                                         const bool &skip_checks = false) {
   
   T x_copy = x_L_val;
   fn_EIGEN_Ref_double(Eigen::Ref<T>(x_copy), fn, vect_type, skip_checks);
   return x_copy;
   
 }
 
 
 //// blocks
 template <typename T, int n_rows = Eigen::Dynamic, int n_cols = Eigen::Dynamic>
 inline auto  fn_EIGEN_double(                       Eigen::Ref<Eigen::Block<T, n_rows, n_cols>> x_Ref,
                                                     const std::string &fn,
                                                     const std::string &vect_type = "Stan",
                                                     const bool &skip_checks = false) {
   
   T x_matrix = x_Ref;
   fn_EIGEN_Ref_double(Eigen::Ref<T>(x_matrix), fn, vect_type, skip_checks);
   return x_matrix;
   
 }
 
 
 
 //// arrays
 template <int n_rows = Eigen::Dynamic, int n_cols = Eigen::Dynamic>
 inline auto  fn_EIGEN_double(               const Eigen::Array<double, n_rows, n_cols> &x,
                                             const std::string &fn, 
                                             const std::string &vect_type = "Stan",
                                             const bool &skip_checks = false) {
   
   using T = Eigen::Matrix<double, n_rows, n_cols>;
   T x_matrix = x.matrix();
   fn_EIGEN_Ref_double(Eigen::Ref<T>(x_matrix), fn, vect_type, skip_checks);
   return x_matrix;
   
   
 }
 
 
 // New overload for general expressions
 template <typename Derived>
 inline auto fn_EIGEN_double(const Eigen::EigenBase<Derived> &x,
                             const std::string &fn,
                             const std::string &vect_type = "Stan",
                             const bool &skip_checks = false) {
   
   using T = Eigen::Matrix<typename Derived::Scalar, 
                           Derived::RowsAtCompileTime, 
                           Derived::ColsAtCompileTime>;
   T x_copy = x;
   fn_EIGEN_Ref_double(Eigen::Ref<T>(x_copy), fn, vect_type, skip_checks);
   return x_copy;
 }
 
 
 
 // Additional Matrix expression overload
 template <typename Derived>
 inline auto fn_EIGEN_double(const Eigen::MatrixBase<Derived> &x,
                             const std::string &fn,
                             const std::string &vect_type = "Stan",
                             const bool &skip_checks = false) {
   
   using T = Eigen::Matrix<typename Derived::Scalar,
                           Derived::RowsAtCompileTime,
                           Derived::ColsAtCompileTime>;
   T x_copy = x;
   fn_EIGEN_Ref_double(Eigen::Ref<T>(x_copy), fn, vect_type, skip_checks);
   return x_copy;
 }
 
 
 
 // Additional overload for array expressions
 template <typename Derived>
 inline auto fn_EIGEN_double(const Eigen::ArrayBase<Derived> &x,
                             const std::string &fn,
                             const std::string &vect_type = "Stan", 
                             const bool &skip_checks = false) {
   
   using T = Eigen::Matrix<typename Derived::Scalar,
                           Derived::RowsAtCompileTime,
                           Derived::ColsAtCompileTime>;
   T x_copy = x.matrix();
   fn_EIGEN_Ref_double(Eigen::Ref<T>(x_copy), fn, vect_type, skip_checks);
   return x_copy;
 }
 
 
 
 
 
 


 
 
 

 
 
 
 
  

// 
// 
// ////// ------------------------------------------------------------ Array function callers / wrappers




  


// 
// void log_sum_exp_pair(const Eigen::Matrix<double, -1, 1>  &log_a,
//                       const Eigen::Matrix<double, -1, 1>  &log_b,
//                       const std::string &vect_type_exp,
//                       const std::string &vect_type_log,
//                       Eigen::Matrix<double, -1, 1> &log_sum_abs_result) {       // output parameter
//   
//   // for each element i, find max(log_a[i], log_b[i])
//   Eigen::Matrix<double, -1, 1> max_logs = log_a.array().max(log_b.array());
//   // for each element i, compute sign_a[i]*exp_a[i] + sign_b[i]*exp_b[i]
//   Eigen::Matrix<double, -1, 1> combined = (fn_EIGEN_double(log_a - max_logs, "exp", vect_type_exp).array()  + 
//     fn_EIGEN_double(log_b - max_logs, "exp", vect_type_exp).array()).matrix(); 
//   // fill both output vectors
//   log_sum_abs_result = max_logs + fn_EIGEN_double(combined.array().abs().matrix(), "log", vect_type_log);
//   
// }
//  




inline void log_sum_exp_general(     const Eigen::Ref<const Eigen::Matrix<double, -1, -1>> log_vals,  
                                     const std::string &vect_type_exp,
                                     const std::string &vect_type_log,
                                     Eigen::Ref<Eigen::Matrix<double, -1, 1>> log_sum_result,
                                     Eigen::Ref<Eigen::Matrix<double, -1, 1>> container_max_logs) {
  
  // find max for each row across all columns
  container_max_logs = log_vals.rowwise().maxCoeff();
  // sum across columns
 //Eigen::Matrix<double, -1, 1> sum_exp =  (fn_EIGEN_double( (log_vals.colwise() - container_max_logs) , "exp", vect_type_exp).array()).matrix().rowwise().sum();
  // compute results
  log_sum_result = container_max_logs + fn_EIGEN_double( (fn_EIGEN_double( (log_vals.colwise() - container_max_logs) , "exp", vect_type_exp).array()).matrix().rowwise().sum().array().abs(), "log", vect_type_log);
   
}



  
 


struct LogSumVecSingedResult { 
  
     double log_sum;
     double sign;
 
};
 
 
 
 
 
 inline LogSumVecSingedResult log_sum_vec_signed_v1(   const Eigen::Ref<const Eigen::Matrix<double, -1, 1>> log_abs_vec,
                                                       const Eigen::Ref<const Eigen::Matrix<double, -1, 1>> signs,
                                                       const std::string &vect_type) {
   
             // const double huge_neg = -700.0;
             double max_log_abs = stan::math::max(log_abs_vec);  // find max 
           
             const Eigen::Matrix<double, -1, 1> &shifted_logs = (log_abs_vec.array() - max_log_abs);   ///// Shift logs and clip
             // shifted_logs = (shifted_logs.array() < huge_neg).select(huge_neg, shifted_logs);   ///// additionally clip (can comment out for no clipping)
             
             // Compute sum with signs carefully
             const Eigen::Matrix<double, -1, 1> &exp_terms = fn_EIGEN_double((log_abs_vec.array() - max_log_abs), "exp", vect_type);
             double sum_exp = (signs.array() * exp_terms.array()).sum();
             
             // // Handle near-zero sums (optional)
             // if (stan::math::abs(sum_exp) < stan::math::exp(huge_neg)) {
             //   return {huge_neg, 0.0};
             // }
             
             double log_abs_sum = max_log_abs + stan::math::log(stan::math::abs(sum_exp));   
             
             // // Clip final result if too large (optional)
             // if (log_abs_sum > 10.0) {  // exp(10) ≈ 22026, reasonable bound
             //   log_abs_sum = 10.0;
             // }
             
             return {log_abs_sum, sum_exp > 0 ? 1.0 : -1.0};
   
 }

 
 
 
 
 
 
 
//// with optional additional underflow protection (can be commented out easily)
inline void log_abs_sum_exp_general_v2(  const Eigen::Ref<const Eigen::Matrix<double, -1, -1>> log_abs_vals,
                                         const Eigen::Ref<const Eigen::Matrix<double, -1, -1>> signs,
                                         const std::string &vect_type_exp,
                                         const std::string &vect_type_log,
                                         Eigen::Ref<Eigen::Matrix<double, -1, 1>> log_sum_abs_result,
                                         Eigen::Ref<Eigen::Matrix<double, -1, 1>> sign_result,
                                         Eigen::Ref<Eigen::Matrix<double, -1, 1>> container_max_logs,
                                         Eigen::Ref<Eigen::Matrix<double, -1, 1>> container_sum_exp_signed) {

  const double min_exp_neg = -700.0 ;
  const double max_exp_arg =  700.0;
  const double tiny = stan::math::exp(min_exp_neg);

  container_max_logs = log_abs_vals.rowwise().maxCoeff();    // Find max log_abs value for each row 


  const Eigen::Matrix<double, -1, -1>  &shifted_logs = (log_abs_vals.colwise() - container_max_logs);   //// Compute shifted logs with underflow protection
  
  //  Eigen::Matrix<double, -1, -1>  shifted_logs = (log_abs_vals.colwise() - container_max_logs);   //// Compute shifted logs with underflow protection
  //  shifted_logs = (shifted_logs.array() < -max_exp_arg).select( -max_exp_arg, shifted_logs );  ////// Clip very negative values to avoid unnecessary exp computations

  //// Compute exp terms and sum over columns with signs 
  container_sum_exp_signed = (fn_EIGEN_double((log_abs_vals.colwise() - container_max_logs).matrix(), "exp", vect_type_exp).array() *  signs.array()).matrix().rowwise().sum();

  //// Compute sign_result and log_sum_abs_result
  sign_result = container_sum_exp_signed.array().sign();
  log_sum_abs_result.array() = container_max_logs.array() + fn_EIGEN_double( container_sum_exp_signed.array().abs(), "log", vect_type_log).array();
 
   // sign_result(i) = std::copysign(1.0, sum_exp);
   // log_sum_abs_result(i) = container_max_logs(i) +  stan::math::log(stan::math::abs(sum_exp));
 
 // for (Eigen::Index i = 0; i < container_sum_exp_signed.rows(); ++i) {
 // 
 //       double sum_exp = container_sum_exp_signed(i);
 // 
 //       if (stan::math::abs(sum_exp) < tiny) {   //  if exp's cancel out or are too small
 // 
 //             sign_result(i) = 0.0;
 //             log_sum_abs_result(i) = min_exp_neg;
 // 
 //       } else {  // Normal case
 // 
 //             sign_result(i) = std::copysign(1.0, sum_exp);
 //             log_sum_abs_result(i) = container_max_logs(i) +  stan::math::log(stan::math::abs(sum_exp));
 //       }
 // 
 // }
  

} 


 
 
 
 
 
 
 

 
inline  void log_abs_matrix_vector_mult_v1(  const Eigen::Ref<const Eigen::Matrix<double, -1, -1>> log_abs_matrix,
                                             const Eigen::Ref<const Eigen::Matrix<double, -1, -1>> sign_matrix,
                                             const Eigen::Ref<const Eigen::Matrix<double, -1, 1>> log_abs_vector,
                                             const Eigen::Ref<const Eigen::Matrix<double, -1, 1>> sign_vector,
                                             const std::string &vect_type_exp,
                                             const std::string &vect_type_log,
                                             Eigen::Ref<Eigen::Matrix<double, -1, 1>> log_abs_result,
                                             Eigen::Ref<Eigen::Matrix<double, -1, 1>> sign_result) {

   int n_rows = log_abs_matrix.rows();
   int n_cols = log_abs_matrix.cols();

             // Initialize temp storage for max finding pass
             Eigen::Matrix<double, -1, 1> max_logs = Eigen::Matrix<double, -1, 1>::Constant(n_rows, -700.0);

             // First pass: find max_log for each row
             for (int j = 0; j < n_cols; j++) {
               double log_vec_j = log_abs_vector(j);
               for (int i = 0; i < n_rows; i++) {
                 max_logs(i) = std::max(max_logs(i), log_abs_matrix(i,j) + log_vec_j);
               }
             }

             // Second pass: compute sums using exp-trick
             Eigen::Matrix<double, -1, 1> sums = Eigen::Matrix<double, -1, 1>::Zero(n_rows);
             for (int j = 0; j < n_cols; j++) {
               double log_vec_j = log_abs_vector(j);
               double sign_vec_j = sign_vector(j);

               for (int i = 0; i < n_rows; i++) {
                 double term = std::exp(log_abs_matrix(i,j) + log_vec_j - max_logs(i)) *
                   sign_matrix(i,j) * sign_vec_j;
                 sums(i) += term;
               }
             }

             // Final pass: compute results
             for (int i = 0; i < n_rows; i++) {
               sign_result(i) = (sums(i) >= 0) ? 1.0 : -1.0;
               log_abs_result(i) = std::log(std::abs(sums(i))) + max_logs(i);
             }

 }


 
 
 
 
 
 
 
 
 
 
#if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2
inline Eigen::Matrix<double, -1, 1  >   fast_log_sum_exp_2d_AVX2_double(  Eigen::Ref<Eigen::Matrix<double, -1, -1>>   x) {
   
   const int N = x.rows();
   const std::string vect_type_exp = "AVX2";
   const std::string vect_type_log = "AVX2";
   
   Eigen::Matrix<double, -1, 1> log_sum_abs_result(N);
   Eigen::Matrix<double, -1, 1> container_max_logs(N);
   
   log_sum_exp_general(x, 
                       vect_type_exp,
                       vect_type_log,
                       log_sum_abs_result,
                       container_max_logs);
   
   return log_sum_abs_result;
   
 } 
#endif
 
 
 
 
 
 
 
 
 
 
#if defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)
inline Eigen::Matrix<double, -1, 1  >   fast_log_sum_exp_2d_AVX512_double(  Eigen::Ref<Eigen::Matrix<double, -1, -1>> x) {
   
   const int N = x.rows();
   const std::string vect_type_exp = "AVX512";
   const std::string vect_type_log = "AVX512";
   
   Eigen::Matrix<double, -1, 1> log_sum_abs_result(N);
   Eigen::Matrix<double, -1, 1> container_max_logs(N);
   
   log_sum_exp_general(x, 
                       vect_type_exp,
                       vect_type_log,
                       log_sum_abs_result,
                       container_max_logs);
   
   return log_sum_abs_result;
     
 }
#endif
 
  
 
 
 
 
 
 
 ///////////////////  fns - log_sum_exp_2d (for 2d-array, vectorised)   ----------------------------------------------------------------------------------------------------------------------------------------
 //-----------
 
 inline Eigen::Matrix<double, -1, 1  >   log_sum_exp_2d_Eigen_double( const Eigen::Ref<const Eigen::Matrix<double, -1, -1>>  x )  {
   
   int N = x.rows();
   Eigen::Matrix<double, -1, -1> rowwise_maxes_2d_array(N, 2);
   rowwise_maxes_2d_array.col(0) = x.array().rowwise().maxCoeff().matrix();
   rowwise_maxes_2d_array.col(1) = rowwise_maxes_2d_array.col(0);
   
   /// Eigen::Matrix<double, -1, 1>  rowwise_maxes_1d_vec = rowwise_maxes_2d_array.col(0);
   Eigen::Matrix<double, -1, 1>  sum_exp_vec =  (  (x.array()  -  rowwise_maxes_2d_array.array()).matrix() ).array().exp().matrix().rowwise().sum() ;
   
   return     ( rowwise_maxes_2d_array.col(0).array()    +    sum_exp_vec.array().log() ).matrix() ;
   
   
 }




 



inline Eigen::Matrix<double, -1, 1  >   log_sum_exp_2d_Stan_double( const Eigen::Ref<const Eigen::Matrix<double, -1, -1>>  x )  {
  
  using namespace stan::math;
  
  int N = x.rows();
  Eigen::Matrix<double, -1, -1>   rowwise_maxes_2d_array(N, 2);
  rowwise_maxes_2d_array.col(0) = x.array().rowwise().maxCoeff().matrix(); 
  rowwise_maxes_2d_array.col(1) = rowwise_maxes_2d_array.col(0);
  
  //// Eigen::Matrix<double, -1, 1>  rowwise_maxes_1d_vec = rowwise_maxes_2d_array.col(0);
  Eigen::Matrix<double, -1, 1>  sum_exp_vec =  stan::math::exp(  (x.array()  -  rowwise_maxes_2d_array.array()).matrix() ).rowwise().sum() ;
  
  return     ( rowwise_maxes_2d_array.col(0).array()    +   stan::math::log(sum_exp_vec).array() ).matrix() ;
  
}



 








inline Eigen::Matrix<double, -1, 1  >   fast_log_sum_exp_2d_double(   const Eigen::Ref<const Eigen::Matrix<double, -1, -1>>  x )  {
  
  int N = x.rows();
  Eigen::Matrix<double, -1, -1  > rowwise_maxes_2d_array(N, 2);
  rowwise_maxes_2d_array.col(0) = x.array().rowwise().maxCoeff().matrix();
  rowwise_maxes_2d_array.col(1) = rowwise_maxes_2d_array.col(0); 
  
  // std::function<double(double)>  exp_fn =  static_cast<double(*)(double const &)>(fast_exp_1);
  // std::function<double(double)>  log_fn =  static_cast<double(*)(double const &)>(fast_log_1);
  
  return      rowwise_maxes_2d_array.col(0) +
    fn_colvec_loop_dbl_Eigen(  fn_mat_loop_dbl_Eigen( (x  - rowwise_maxes_2d_array).matrix(),  static_cast<double(*)(double const)>(fast_exp_1)).array().rowwise().sum().abs().matrix(), static_cast<double(*)(double const)>(fast_log_1)   )   ;
  
  
}




 
 
 
  
 
 
inline Eigen::Matrix<double, -1, 1  > fn_log_sum_exp_2d_double(     Eigen::Ref<Eigen::Matrix<double, -1, -1>>  x,    // Eigen::Matrix<double, -1, 2> &x,
                                                                     const std::string &vect_type = "Stan",
                                                                     const bool &skip_checks = false) {
   
   {
     if (vect_type == "Eigen") {
       return  log_sum_exp_2d_Eigen_double(x);
     } else if (vect_type == "Stan") {
       return  log_sum_exp_2d_Stan_double(x);
     } else if (vect_type == "AVX2") {
#if defined(__AVX2__) && ( !(defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)) ) // use AVX2
       if (skip_checks == false)   return  fast_log_sum_exp_2d_AVX2_double(x);
       else                        return  fast_log_sum_exp_2d_AVX2_double(x);
#endif
     } else if (vect_type == "AVX512") {
#if defined(__AVX512VL__) && defined(__AVX512F__)  && defined(__AVX512DQ__)
       if (skip_checks == false)   return  fast_log_sum_exp_2d_AVX512_double(x);
       else                        return  fast_log_sum_exp_2d_AVX512_double(x);
#endif
     } else {
              return  log_sum_exp_2d_Stan_double(x);
     }
     
   }
   
   return  log_sum_exp_2d_Stan_double(x);
   
}
 
 
  






#endif
// 
// 
// 
// 
//  
//  
//  
//  
//  