#ifndef STAN_MATH_REV_FUN_MVP_EXP_APPROX_HPP
#define STAN_MATH_REV_FUN_MVP_EXP_APPROX_HPP
 
 

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
 
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
 
#include <typeinfo>
#include <type_traits>
#include <sstream>
#include <stdexcept>
#include <stan/math/prim/err/invalid_argument.hpp>
 
#include <stan/math/prim/meta/is_var_matrix.hpp>
 
#if defined(__AVX2__) || defined(__AVX512F__)
#include <immintrin.h>
#endif
  
  
#include "approx_cpp_fns_for_Stan_prim.hpp"
 
 
  
  
  
  
//////////////  --------- exp  --------------------------------------------------------------------------------------------------------------------




namespace stan {
namespace math {


 

// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_exp_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                    std::ostream* pstream__ = nullptr) {
  
  // All matrices used in reverse pass callback must be arena allocated
  // This creates a matrix of vars
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  // This part:
  // 1. value_of(x) -> gets Matrix<double>
  // 2. mvp_exp_approx(value_of(x)) -> returns Matrix<double>
  // 3. result = ... -> Stan automatically converts the doubles to vars
  /// stan::arena_t<Eigen::Matrix<double, Rows, Cols>> x_double =    value_of(x);
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> result_double =  mvp_exp_approx(value_of(x));
  result = result_double ; // mvp_exp_approx(value_of(x));
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative 
  ////stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv = result_double; // deriv of exp is exp!
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, result_double]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(result_double); 
  });
  
  // Make a hard copy to avoid overwriting values used in reverse pass
  return Eigen::Matrix<stan::math::var, Rows, Cols>(result);
   
}




 

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_exp_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                             std::ostream* pstream__ = nullptr) {
  
  return mvp_exp_approx_template<-1, 1>(a, pstream__);
  
}  

// stan::math::var overload - row_vector
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_exp_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                             std::ostream* pstream__ = nullptr) {
  
  return mvp_exp_approx_template<1, -1>(a, pstream__);
  
} 

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_exp_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                              std::ostream* pstream__ = nullptr) {
  
  return mvp_exp_approx_template<-1, -1>(a, pstream__); 
  
}




// 
// // Forward mode template
// template <int Rows, int Cols, typename T>
// Eigen::Matrix<fvar<T>, Rows, Cols> mvp_exp_approx_template_fwd(   const Eigen::Matrix<fvar<T>, Rows, Cols>& x,
//                                                                   std::ostream* pstream__ = nullptr) {
//   
//   Eigen::Matrix<fvar<T>, Rows, Cols> result(x.rows(), x.cols());
//   
//   // Get values using vectorized operations
//   auto values = x.unaryExpr([](const fvar<T>& x) { return x.val_; });
//   auto exp_values = mvp_exp_approx(values);
//   
//   // Set values and derivatives (d/dx exp(x) = exp(x))
//   result.unaryExpr([&](const fvar<T>& x) {
//     fvar<T> res;
//     res.val_ = exp_values(x);
//     res.d_ = exp_values(x) * x.d_;
//     return res;
//   });
//   
//   return result;
//   
// }
// 
// 
// // Forward mode overloads
// template <typename T>
// inline Eigen::Matrix<fvar<T>, -1, 1> mvp_exp_approx( const Eigen::Matrix<fvar<T>, -1, 1> &x,
//                                                      std::ostream* pstream__ = nullptr) {
//   
//   return mvp_exp_approx_template_fwd<-1, 1>(x, pstream__);
//   
// } 
// 
// template <typename T>
// inline Eigen::Matrix<fvar<T>, 1, -1> mvp_exp_approx( const Eigen::Matrix<fvar<T>, 1, -1> &x,
//                                                      std::ostream* pstream__ = nullptr) {
//   
//   return mvp_exp_approx_template_fwd<1, -1>(x, pstream__);
//   
// } 
// 
// template <typename T>
// inline Eigen::Matrix<fvar<T>, -1, -1> mvp_exp_approx( const Eigen::Matrix<fvar<T>, -1, -1> &x,
//                                                       std::ostream* pstream__ = nullptr) {
//   
//   return mvp_exp_approx_template_fwd<-1, -1>(x, pstream__);
//   
// } 
// 
// 




}
}






//////////////  --------- log  --------------------------------------------------------------------------------------------------------------------

namespace stan {
namespace math {


 



// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_log_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                    std::ostream* pstream__ = nullptr) {
  
   // Arena allocate input for reverse pass
   stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  

    // All matrices used in reverse pass callback must be arena allocated
   stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>>  result =  fn_EIGEN_double(value_of(x), "log", vect_type, false);

   //// // Pre-compute derivative
  ///  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv = x_double.cwiseInverse();

    // Use reverse_pass_callback since our return type is an Eigen matrix
    reverse_pass_callback([x_arena, result]() mutable {
      x_arena.adj().array() += result.adj().array() * (1.0 / value_of(x_arena).array()).array();
    });

     return Eigen::Matrix<stan::math::var, Rows, Cols>(result);

}




 

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_log_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                              std::ostream* pstream__ = nullptr) {

  return mvp_log_approx_template<-1, 1>(a, pstream__);

}

// stan::math::var overload - row_vector
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_log_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                             std::ostream* pstream__ = nullptr) {

  return mvp_log_approx_template<1, -1>(a, pstream__);

}

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_log_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                              std::ostream* pstream__ = nullptr) {

  return mvp_log_approx_template<-1, -1>(a, pstream__);

}




 
 
 
 
 // 
 // 
 // 
 // // Forward mode version
 // template <int Rows, int Cols>
 // inline Eigen::Matrix<fvar<double>, Rows, Cols> mvp_log_approx_template_fwd( const Eigen::Ref<const Eigen::Matrix<fvar<double>, Rows, Cols>> x,
 //                                                                             std::ostream* pstream__ = nullptr) {
 //   
 //   Eigen::Matrix<fvar<double>, Rows, Cols> result(x.rows(), x.cols());
 //   Eigen::Matrix<double, Rows, Cols> x_double = x.val();
 //   Eigen::Matrix<double, Rows, Cols> res_double =  mvp_log_approx(x_double);
 //   
 //   Eigen::Matrix<double, Rows, Cols> deriv_double = ( 1.0 / x_double.array() ).matrix();
 //   
 //   for (int j = 0; j < x.cols(); j++) { // loop through cols first as col-major storage 
 //     for (int i = 0; i < x.rows(); i++) {
 //       result(i, j).val_ = res_double(i, j);
 //       result(i, j).d_ = x(i, j).d_ * deriv_double(i, j); // Chain rule: d/dx(inv_logit(x)) = inv_logit(x)(1-inv_logit(x))
 //     }
 //   }
 //   
 //   return result;
 //   
 // }
 // 
 // 
 // 
 // 
 // 
 // // Forward mode overloads
 // inline Eigen::Matrix<fvar<double>, -1, 1> mvp_log_approx( const Eigen::Matrix<fvar<double>, -1, 1> &x,
 //                                                          std::ostream* pstream__ = nullptr) {
 //   
 //   return mvp_log_approx_template_fwd<-1, 1>(x, pstream__);
 //   
 // }
 // 
 // inline Eigen::Matrix<fvar<double>, 1, -1> mvp_log_approx( const Eigen::Matrix<fvar<double>, 1, -1> &x,
 //                                                          std::ostream* pstream__ = nullptr) {
 //   
 //   return mvp_log_approx_template_fwd<1, -1>(x, pstream__);
 //   
 // }
 // 
 // 
 // inline Eigen::Matrix<fvar<double>, -1, -1> mvp_log_approx( const Eigen::Matrix<fvar<double>, -1, -1> &x,
 //                                                           std::ostream* pstream__ = nullptr) {
 //   
 //   return mvp_log_approx_template_fwd<-1, -1>(x, pstream__);
 //   
 // }
 // 
 // 
 // 
 // 
 
 
 
 
 


}
}






//////////////  --------- log1p  --------------------------------------------------------------------------------------------------------------------

namespace stan {
namespace math {


 


 
 // Main template header
 template <int Rows, int Cols>
 Eigen::Matrix<stan::math::var, Rows, Cols> mvp_log1p_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                     std::ostream* pstream__ = nullptr) {
   
   // All matrices used in reverse pass callback must be arena allocated
   stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
   stan::arena_t<Eigen::Matrix<double, Rows, Cols>> x_double =    value_of(x);
   stan::arena_t<Eigen::Matrix<double, Rows, Cols>> result_double =  mvp_log1p_approx(x_double);
   result = result_double;
   
   // Arena allocate input for reverse pass
   stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
   
   // Pre-compute derivative  
   stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv =  ( 1.0 / (1.0 + x_double.array()).array() ).matrix();
   // auto deriv = (1.0 + value_of(x_arena)).cwiseInverse();
   
   // Use reverse_pass_callback since our return type is an Eigen matrix
   reverse_pass_callback([x_arena, result, deriv]() mutable {
     x_arena.adj() += result.adj().cwiseProduct(deriv); 
   });
   
 
 return Eigen::Matrix<stan::math::var, Rows, Cols>(result);
   
 }
 
 
 
 
 
 // Overloads with corrected types:
 
 // stan::math::var overload - vector
 inline Eigen::Matrix<stan::math::var, -1, 1> mvp_log1p_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                              std::ostream* pstream__ = nullptr) {
   
   return mvp_log1p_approx_template<-1, 1>(a, pstream__);
   
 } 
 
 // stan::math::var overload - row_vector 
 inline Eigen::Matrix<stan::math::var, 1, -1> mvp_log1p_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a, 
                                                              std::ostream* pstream__ = nullptr) {
   
   return mvp_log1p_approx_template<1, -1>(a, pstream__);
   
 }
 
 // stan::math::var overload - matrix
 inline Eigen::Matrix<stan::math::var, -1, -1> mvp_log1p_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                               std::ostream* pstream__ = nullptr) {
   
   return mvp_log1p_approx_template<-1, -1>(a, pstream__);
   
 }
 
 
 
 



}
}





//////////////  --------- log1m  --------------------------------------------------------------------------------------------------------------------

namespace stan {
namespace math {







// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_log1m_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                      std::ostream* pstream__ = nullptr) {
  
          // All matrices used in reverse pass callback must be arena allocated
          stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
          stan::arena_t<Eigen::Matrix<double, Rows, Cols>> x_double =    value_of(x);
          stan::arena_t<Eigen::Matrix<double, Rows, Cols>> result_double =  mvp_log1m_approx(x_double);
          result = result_double;
          
          // Arena allocate input for reverse pass
          stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
          
          // Pre-compute derivative  
          stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv =  ( - 1.0 / (1.0 - x_double.array()).array() ).matrix();
          
          // Use reverse_pass_callback since our return type is an Eigen matrix
          reverse_pass_callback([x_arena, result, deriv]() mutable {
            x_arena.adj() += result.adj().cwiseProduct(deriv); 
          });
          
         
         return Eigen::Matrix<stan::math::var, Rows, Cols>(result);
  
}





// Overloads with corrected types:

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_log1m_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                               std::ostream* pstream__ = nullptr) {
  
  return mvp_log1m_approx_template<-1, 1>(a, pstream__);
  
} 

// stan::math::var overload - row_vector 
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_log1m_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a, 
                                                               std::ostream* pstream__ = nullptr) {
  
  return mvp_log1m_approx_template<1, -1>(a, pstream__);
  
}

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_log1m_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                                std::ostream* pstream__ = nullptr) {
  
  return mvp_log1m_approx_template<-1, -1>(a, pstream__);
  
}

 



}
}









//////////////  --------- logit  --------------------------------------------------------------------------------------------------------------------

namespace stan {
namespace math {






// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_logit_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                      std::ostream* pstream__ = nullptr) {
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> x_double =    value_of(x);
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> result_double =  mvp_logit_approx(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative  
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv =  ( 1.0 / ( x_double.array() * (1.0 - x_double.array()).array() ).array() ).matrix();
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj().array() += result.adj().array() * deriv.array();
  }); 
  
 
 return Eigen::Matrix<stan::math::var, Rows, Cols>(result);
  
} 





// Overloads with corrected types:

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_logit_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                              std::ostream* pstream__ = nullptr) { 
  
  return mvp_logit_approx_template<-1, 1>(a, pstream__);
  
} 

// stan::math::var overload - row_vector 
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_logit_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a, 
                                                              std::ostream* pstream__ = nullptr) { 
  
  return mvp_logit_approx_template<1, -1>(a, pstream__);
  
} 

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_logit_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                               std::ostream* pstream__ = nullptr) { 
  
  return mvp_logit_approx_template<-1, -1>(a, pstream__);
  
} 










}
}














//////////////  --------- inv_logit  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {


 



 



// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_inv_logit_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                   std::ostream* pstream__ = nullptr) {

  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena =  (x);

  // All matrices used in reverse pass callback must be arena allocated
  /// stan::arena_t<Eigen::Matrix<double, Rows, Cols>> x_double =     value_of(x_arena);
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> result_double =   fn_EIGEN_double(value_of(x_arena), "inv_logit", vect_type, false);
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result =   result_double;

  // Pre-compute derivative
  ///  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv =   ((result_double.array() * (1.0 - result_double.array()).array() ).matrix());

  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, result_double]() mutable {
    x_arena.adj().array() += result.adj().array() * (result_double.array() * (1.0 - result_double.array()).array() ).array(); 
  });


  return Eigen::Matrix<stan::math::var, Rows, Cols>(result);

}





//
 // stan::math::var overload - vector
 inline Eigen::Matrix<stan::math::var, -1, 1> mvp_inv_logit( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                             std::ostream* pstream__ = nullptr) {

   return mvp_inv_logit_template<-1, 1>(a, pstream__);

 }

 // stan::math::var overload - row_vector
 inline Eigen::Matrix<stan::math::var, 1, -1> mvp_inv_logit( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                             std::ostream* pstream__ = nullptr) {

   return mvp_inv_logit_template<1, -1>(a, pstream__);

 }

 // stan::math::var overload - matrix
 inline Eigen::Matrix<stan::math::var, -1, -1> mvp_inv_logit( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                              std::ostream* pstream__ = nullptr) {

   return mvp_inv_logit_template<-1, -1>(a, pstream__);

 }

 



// 
// 
// // Forward mode version
// template <int Rows, int Cols>
// inline Eigen::Matrix<fvar<double>, Rows, Cols> mvp_inv_logit_template_fwd( const Eigen::Ref<const Eigen::Matrix<fvar<double>, Rows, Cols>> x,
//                                                                            std::ostream* pstream__ = nullptr) {
// 
//   Eigen::Matrix<fvar<double>, Rows, Cols> result(x.rows(), x.cols());
//   Eigen::Matrix<double, Rows, Cols> x_double = x.val();
//   Eigen::Matrix<double, Rows, Cols> res_double =  mvp_inv_logit(x_double);
//   Eigen::Matrix<double, Rows, Cols> deriv_double = res_double.array() * (1.0 - res_double.array());
// 
//   for (int j = 0; j < x.cols(); j++) { // loop through cols first as col-major storage 
//       for (int i = 0; i < x.rows(); i++) {
//         result(i, j).val_ = res_double(i, j);
//         result(i, j).d_ = x(i, j).d_ * deriv_double(i, j); // Chain rule: d/dx(inv_logit(x)) = inv_logit(x)(1-inv_logit(x))
//      }
//   }
// 
//   return result;
//   
// }
// 
// 
//  
// 
// 
// // Forward mode overloads
// inline Eigen::Matrix<fvar<double>, -1, 1> mvp_inv_logit( const Eigen::Matrix<fvar<double>, -1, 1> &x,
//                                                      std::ostream* pstream__ = nullptr) {
// 
//   return mvp_inv_logit_template_fwd<-1, 1>(x, pstream__);
// 
// }
// 
// inline Eigen::Matrix<fvar<double>, 1, -1> mvp_inv_logit( const Eigen::Matrix<fvar<double>, 1, -1> &x,
//                                                      std::ostream* pstream__ = nullptr) {
// 
//   return mvp_inv_logit_template_fwd<1, -1>(x, pstream__);
// 
// }
// 
// 
// inline Eigen::Matrix<fvar<double>, -1, -1> mvp_inv_logit( const Eigen::Matrix<fvar<double>, -1, -1> &x,
//                                                      std::ostream* pstream__ = nullptr) {
// 
//   return mvp_inv_logit_template_fwd<-1, -1>(x, pstream__);
// 
// }
// 
//  




}
}










//////////////  --------- tanh  --------------------------------------------------------------------------------------------------------------------

namespace stan {
namespace math {






// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_tanh_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                     std::ostream* pstream__ = nullptr) {
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> x_double =    value_of(x);
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> result_double =  mvp_tanh_approx(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative  
  stan::arena_t<Eigen::Matrix<double, Rows, Cols>> deriv =  (1.0 - result_double.array().square()).matrix();
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(deriv); 
  }); 
  
  
  return Eigen::Matrix<stan::math::var, Rows, Cols>(result); 
  
} 





// Overloads with corrected types:

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_tanh_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                              std::ostream* pstream__ = nullptr) {
   
  return mvp_tanh_approx_template<-1, 1>(a, pstream__);
  
}  

// stan::math::var overload - row_vector 
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_tanh_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a, 
                                                              std::ostream* pstream__ = nullptr) { 
  
  return mvp_tanh_approx_template<1, -1>(a, pstream__);
  
} 

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_tanh_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a, 
                                                               std::ostream* pstream__ = nullptr) {
  
  return mvp_tanh_approx_template<-1, -1>(a, pstream__);
  
} 










}
}











//////////////  --------- log_inv_logit  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {




// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_log_inv_logit_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                   std::ostream* pstream__ = nullptr) {
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
  Eigen::Matrix<double, Rows, Cols> result_double =  mvp_log_inv_logit(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative  
  Eigen::Matrix<double, Rows, Cols> deriv =  ( (1.0 - result_double.array()).array() ).matrix();
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(deriv); 
  });  
  
 
  return result;
  
} 
 




// Overloads with corrected types:

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_log_inv_logit( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                            std::ostream* pstream__ = nullptr) { 
  
  return mvp_log_inv_logit_template<-1, 1>(a, pstream__);
   
} 

// stan::math::var overload - row_vector 
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_log_inv_logit( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a, 
                                                            std::ostream* pstream__ = nullptr) {  
  
  return mvp_log_inv_logit_template<1, -1>(a, pstream__);
  
} 
 
// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_log_inv_logit( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                             std::ostream* pstream__ = nullptr) {  
  
  return mvp_log_inv_logit_template<-1, -1>(a, pstream__);
  
}  











}
}
 





//////////////  --------- Phi  --------------------------------------------------------------------------------------------------------------------


namespace stan {
namespace math {






// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_Phi_template(  const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                              std::ostream* pstream__ = nullptr) {
  
    // All matrices used in reverse pass callback must be arena allocated
    stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
    Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
    Eigen::Matrix<double, Rows, Cols> result_double =  mvp_Phi(x_double);
    result = result_double;
    
    // Arena allocate input for reverse pass
    stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
    
    // Pre-compute derivative
    const double sqrt_2_pi_recip = 1.0 / sqrt(2.0 * M_PI);
    Eigen::Matrix<double, Rows, Cols> squared_x = x_double.array().square().matrix();
    Eigen::Matrix<double, Rows, Cols> neg_half_squared = (-0.5 * squared_x);
    ///// Eigen::Matrix<double, Rows, Cols> exp_term = mvp_exp_approx(neg_half_squared);
    Eigen::Matrix<double, Rows, Cols> deriv = sqrt_2_pi_recip * (  mvp_exp_approx(neg_half_squared).array() ).matrix() ;
    
    //  auto deriv = sqrt_2_pi_recip * mvp_exp_approx(  (-0.5 * x_double.array().square() ).matrix() );
    
    // Use reverse_pass_callback since our return type is an Eigen matrix
    reverse_pass_callback([x_arena, result, deriv]() mutable {
      x_arena.adj() += result.adj().cwiseProduct(deriv);  
    });
    
 
    return result;
  
}


 

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_Phi( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a,
                                                      std::ostream* pstream__ = nullptr) {
  
  return mvp_Phi_template<-1, 1>(a, pstream__);
  
} 

// stan::math::var overload - row_vector
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_Phi(  const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                       std::ostream* pstream__ = nullptr) {
  
  return mvp_Phi_template<1, -1>(a, pstream__);
  
}

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_Phi(  const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                        std::ostream* pstream__ = nullptr) {
  
    return mvp_Phi_template<-1, -1>(a, pstream__);
  
}







}
}








//////////////  --------- Phi_approx  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {







// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_Phi_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                    std::ostream* pstream__ = nullptr) {
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
  Eigen::Matrix<double, Rows, Cols> result_double =  mvp_Phi_approx(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative
  const double    aa =       0.07056;
  const double    bb =       1.5976;
  Eigen::Matrix<double, Rows, Cols> inv_logit_x = mvp_inv_logit(x_double);
  Eigen::Matrix<double, Rows, Cols> deriv =  (3.0*aa*aa + bb) *  ( inv_logit_x.array() * (1.0 - inv_logit_x.array()).array() ).matrix() ; 
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(deriv);  
  });
  
 
  return result;
  
} 





 
 

// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a, 
                                                             std::ostream* pstream__ = nullptr) {
  
  return mvp_Phi_approx_template<-1, 1>(a, pstream__);
  
}  

// stan::math::var overload - row_vector
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                             std::ostream* pstream__ = nullptr) {
  
  return mvp_Phi_approx_template<1, -1>(a, pstream__); 
  
} 

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1>& a,
                                                              std::ostream* pstream__ = nullptr) {
  return mvp_Phi_approx_template<-1, -1>(a, pstream__);
   
}





 
 


// also overload for Eigen expressions 
template <typename Derived>
inline auto mvp_Phi_approx(const Eigen::EigenBase<Derived> &x,
                                               std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, var_value<double>>::value, Eigen::Matrix<var_value<double>, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<var_value<double>, -1, 1> temp = x;
      return mvp_Phi_approx_template<-1, 1>(temp, pstream__);
    } else if (Derived::RowsAtCompileTime == 1) {
      Eigen::Matrix<var_value<double>, 1, -1> temp = x;
      return mvp_Phi_approx_template<1, -1>(temp, pstream__);
    } else { 
      Eigen::Matrix<var_value<double>, -1, -1> temp = x;
      return mvp_Phi_approx_template<-1, -1>(temp, pstream__);
    } 
    
     
  } 






}
}






//////////////  --------- inv_Phi  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {









// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_inv_Phi_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                    std::ostream* pstream__ = nullptr) {
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
  Eigen::Matrix<double, Rows, Cols> result_double =  mvp_inv_Phi(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute to avoid array operations in callback
  const double sqrt_2_pi_recip = INV_SQRT_TWO_PI;  // using Stan's built-in constant
  Eigen::Matrix<double, Rows, Cols> x_squared = x_double.array().square().matrix();
  Eigen::Matrix<double, Rows, Cols> stuff_to_exp =   ( - 0.5 * x_squared.array() ).matrix() ; 
  Eigen::Matrix<double, Rows, Cols> deriv = ( 1.0  / ( sqrt_2_pi_recip * mvp_exp_approx(stuff_to_exp) ).array() ).matrix();
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(deriv);   
  });
  
 
  return result;
  
} 
 


  
// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_inv_Phi( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a, 
                                                          std::ostream* pstream__ = nullptr) {
  
  return mvp_inv_Phi_template<-1, 1>(a, pstream__);
  
}

// stan::math::var overload - row_vector
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_inv_Phi( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                          std::ostream* pstream__ = nullptr) {
  
  return mvp_inv_Phi_template<1, -1>(a, pstream__);
  
}

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_inv_Phi( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                           std::ostream* pstream__ = nullptr) {
  return mvp_inv_Phi_template<-1, -1>(a, pstream__);
   
}




 
  



}
}








//////////////  --------- inv_Phi_approx_from_logit_prob  --------------------------------------------------------------------------------------------------------------------


////// ----- first make helper fns for "mvp_inv_Phi_approx_from_logit_prob" - can also make these Stan fns in their own right (do in future).
double mvp_DOUBLE_asinh_approx(const double x) {
  return stan::math::mvp_log_approx(x + sqrt(x*x + 1.0));
}

template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols> mvp_EIGEN_asinh_approx(const Eigen::Matrix<double, Rows, Cols> &x) {
  const Eigen::Matrix<double, Rows, Cols> stuff_to_log = (x.array() + (x.array().square() + 1.0).sqrt()).matrix();
  return stan::math::mvp_log_approx(stuff_to_log);
}
 


double mvp_DOBULE_cosh_approx(const double x) {
  return stan::math::mvp_log_approx(x + sqrt(x*x - 1.0));
}

template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols> mvp_EIGEN_cosh_approx(const Eigen::Matrix<double, Rows, Cols> &x) {
  const Eigen::Matrix<double, Rows, Cols> stuff_to_log = (x.array() + (x.array().square() - 1.0).sqrt()).matrix();
  return stan::math::mvp_log_approx(stuff_to_log);
}




double d_inv_Phi_approx_from_logit_prob(const double logit_p) {
  const double inner_part = 0.3418 * logit_p;
  const double middle_part = 0.33333333333333331483 * mvp_DOUBLE_asinh_approx(inner_part);
  return 5.494 * mvp_DOBULE_cosh_approx(middle_part) * 0.33333333333333331483 * (1.0/sqrt(1.0 + inner_part*inner_part)) * 0.3418;
}

template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols> d_inv_Phi_approx_from_logit_prob(const Eigen::Matrix<double, Rows, Cols> &logit_p) {
  const Eigen::Matrix<double, Rows, Cols> inner_part = 0.3418 * logit_p;
  const Eigen::Matrix<double, Rows, Cols> middle_part = 0.33333333333333331483 * mvp_EIGEN_asinh_approx(inner_part);
  return 5.494 * mvp_EIGEN_cosh_approx(middle_part).array() * 0.33333333333333331483 * (1.0/(1.0 + inner_part.array().square()).array().sqrt()) * 0.3418;
}







namespace stan {
namespace math {



 
 
 
 // Main template header
 template <int Rows, int Cols>
 Eigen::Matrix<stan::math::var, Rows, Cols> mvp_inv_Phi_approx_from_logit_prob_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                  std::ostream* pstream__ = nullptr) {
   
   check_finite("mvp_inv_Phi_approx_from_logit_prob", "x", x);  // using Stan's input validation check
   
   // All matrices used in reverse pass callback must be arena allocated
   stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
   Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
   Eigen::Matrix<double, Rows, Cols> result_double =  mvp_inv_Phi_approx_from_logit_prob(x_double);
   result = result_double;
   
   // Arena allocate input for reverse pass
   stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
   
   // Pre-compute derivative
   Eigen::Matrix<double, Rows, Cols> deriv =  d_inv_Phi_approx_from_logit_prob(x_double);
   
   // Use reverse_pass_callback since our return type is an Eigen matrix
   reverse_pass_callback([x_arena, result, deriv]() mutable {
     x_arena.adj() += result.adj().cwiseProduct(deriv);   
   });
    
   return result;
   
 } 
 
 
  
   
 
 
 
 
 // stan::math::var overload - vector
 inline Eigen::Matrix<stan::math::var, -1, 1> mvp_inv_Phi_approx_from_logit_prob( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a, 
                                                           std::ostream* pstream__ = nullptr) {
   
   return mvp_inv_Phi_approx_from_logit_prob_template<-1, 1>(a, pstream__);
   
 }
 
 // stan::math::var overload - row_vector 
 inline Eigen::Matrix<stan::math::var, 1, -1> mvp_inv_Phi_approx_from_logit_prob( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                           std::ostream* pstream__ = nullptr) {
   
   return mvp_inv_Phi_approx_from_logit_prob_template<1, -1>(a, pstream__);
   
 }
 
 // stan::math::var overload - matrix
 inline Eigen::Matrix<stan::math::var, -1, -1> mvp_inv_Phi_approx_from_logit_prob( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                            std::ostream* pstream__ = nullptr) {
   return mvp_inv_Phi_approx_from_logit_prob_template<-1, -1>(a, pstream__);
   
 }
 
 
 
 


 
  

 
 


} // namespace math
} // namespace stan




 
 
 
//////////////  --------- inv_Phi_approx  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {



 
 double d_inv_Phi_approx(const double p) {
   
   const double deriv_wrt_logit_p = d_inv_Phi_approx_from_logit_prob(mvp_logit_approx(p));
   const double deriv_logit_p =    1.0 / (p * (1.0 - p) ) ;    
   const double deriv =   deriv_wrt_logit_p * deriv_logit_p;
   
   return deriv;
   
 }

template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols> d_inv_Phi_approx(const Eigen::Matrix<double, Rows, Cols> &p) {
  
  /// Eigen::Matrix<double, Rows, Cols> logit_p = mvp_logit_approx(p);
  Eigen::Matrix<double, Rows, Cols> deriv_wrt_logit_p = d_inv_Phi_approx_from_logit_prob(mvp_logit_approx(p));
  Eigen::Matrix<double, Rows, Cols> deriv_logit_p =  ( ( 1.0 / (p.array() * (1.0 - p.array())) ).array() ).matrix();
  Eigen::Matrix<double, Rows, Cols> deriv =   (  deriv_wrt_logit_p.array() * deriv_logit_p.array() ).matrix();
  
  return deriv;
  
}
 



// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_inv_Phi_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                        std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_inv_Phi_approx", "x", x);  // using Stan's input validation check
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
  Eigen::Matrix<double, Rows, Cols> result_double =  mvp_inv_Phi_approx(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative
  Eigen::Matrix<double, Rows, Cols> deriv =  d_inv_Phi_approx(x_double);
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(deriv);   
  });
  
  return result;
  
} 








// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_inv_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a, 
                                                                                 std::ostream* pstream__ = nullptr) {
  
  return mvp_inv_Phi_approx_template<-1, 1>(a, pstream__);
  
}

// stan::math::var overload - row_vector 
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_inv_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                                                 std::ostream* pstream__ = nullptr) {
  
  return mvp_inv_Phi_approx_template<1, -1>(a, pstream__);
  
}

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_inv_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                                                  std::ostream* pstream__ = nullptr) {
  return mvp_inv_Phi_approx_template<-1, -1>(a, pstream__);
  
}






 




 
}
}











//////////////  --------- log_Phi_approx  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {




 




// Main template header
template <int Rows, int Cols>
Eigen::Matrix<stan::math::var, Rows, Cols> mvp_log_Phi_approx_template( const Eigen::Ref<const Eigen::Matrix<var_value<double>, Rows, Cols>> x,
                                                                        std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_log_Phi_approx", "x", x);  // using Stan's input validation check
  
  // All matrices used in reverse pass callback must be arena allocated
  stan::arena_t<Eigen::Matrix<stan::math::var, Rows, Cols>> result(x.rows(), x.cols());
  Eigen::Matrix<double, Rows, Cols> x_double =    value_of(x);
  Eigen::Matrix<double, Rows, Cols> result_double =  mvp_log_Phi_approx(x_double);
  result = result_double;
  
  // Arena allocate input for reverse pass
  stan::arena_t<Eigen::Matrix<var_value<double>, Rows, Cols>> x_arena = x;
  
  // Pre-compute derivative
  const double    aa =       0.07056;
  const double    bb =       1.5976; 
  Eigen::Matrix<double, Rows, Cols> deriv_wrt_Phi_approx = ( 1.0 / mvp_Phi_approx(x_double).array() ).matrix();
  Eigen::Matrix<double, Rows, Cols> deriv_of_Phi_approx = d_inv_Phi_approx(x_double);
  Eigen::Matrix<double, Rows, Cols> deriv = ( deriv_wrt_Phi_approx.array() * deriv_of_Phi_approx.array() ).matrix() ; 
 
  
  // Use reverse_pass_callback since our return type is an Eigen matrix
  reverse_pass_callback([x_arena, result, deriv]() mutable {
    x_arena.adj() += result.adj().cwiseProduct(deriv);   
  }); 
  
  return result;
  
}  








// stan::math::var overload - vector
inline Eigen::Matrix<stan::math::var, -1, 1> mvp_log_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, 1> &a, 
                                                                 std::ostream* pstream__ = nullptr) {
  
  return mvp_log_Phi_approx_template<-1, 1>(a, pstream__);
  
} 

// stan::math::var overload - row_vector 
inline Eigen::Matrix<stan::math::var, 1, -1> mvp_log_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, 1, -1> &a,
                                                                 std::ostream* pstream__ = nullptr) {
  
  return mvp_log_Phi_approx_template<1, -1>(a, pstream__);
  
} 

// stan::math::var overload - matrix
inline Eigen::Matrix<stan::math::var, -1, -1> mvp_log_Phi_approx( const Eigen::Matrix<stan::math::var_value<double>, -1, -1> &a,
                                                                  std::ostream* pstream__ = nullptr) {
  return mvp_log_Phi_approx_template<-1, -1>(a, pstream__);
  
} 







 

 



}
}




//
//
// 




  
  
  
#endif
 
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
