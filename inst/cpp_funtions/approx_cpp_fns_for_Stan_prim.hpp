#ifndef STAN_MATH_PRIM_FUN_MVP_EXP_APPROX_HPP
#define STAN_MATH_PRIM_FUN_MVP_EXP_APPROX_HPP

 
#include <stan/math/prim.hpp>
 

#include <Eigen/Dense>
#include <Eigen/Core>

#include <vector>

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

 
#include <typeinfo>
#include <sstream>
#include <stdexcept>
#include <stan/math/prim/err/invalid_argument.hpp>
 
#include <stan/math/prim/meta/is_matrix.hpp>
#include <stan/math/prim/meta/is_matrix_cl.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
 
 
 
#include <stan/math/prim/fun/Eigen.hpp>
 
 
#if defined(__AVX2__) || defined(__AVX512F__)
#include <immintrin.h>
#endif
 

 
#include "general_functions/double_fns.hpp"
 
 
//// --------- BayesMVP includes - General functions
#include "general_functions/fns_SIMD_and_wrappers/fn_wrappers_Stan.hpp"
#include "general_functions/fns_SIMD_and_wrappers/fn_wrappers_Loop.hpp"
////
#include "general_functions/fns_SIMD_and_wrappers/fast_and_approx_AVX2_fns.hpp" // will only compile if  AVX2 is available
#include "general_functions/fns_SIMD_and_wrappers/fast_and_approx_AVX512_fns.hpp" // will only compile if  AVX-512 is available
#include "general_functions/fns_SIMD_and_wrappers/fn_wrappers_SIMD_AVX_general.hpp" // will only compile if AVX-512 (1st choice) or AVX2 available
////
#include "general_functions/fns_SIMD_and_wrappers/fn_wrappers_overall.hpp"
#include "general_functions/fns_SIMD_and_wrappers/fn_wrappers_log_sum_exp_dbl.hpp"
#include "general_functions/fns_SIMD_and_wrappers/fn_wrappers_log_sum_exp_SIMD.hpp"
 
 
 
 
 
////////////////-------- first determine SIMD (i.e. vectorisation) type to use as global static variable - will add more in future -------------------------------------------
static const std::string vect_type = [] {
      #ifdef __AVX512F__
       return "AVX512"; 
      #elif defined(__AVX2__) && (!(defined(__AVX512F__)))
       return "AVX2";
      #else 
       return "Stan";
      #endif
}();
 

 

//////////////  --------- exp  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {
    
    
    
    

///// DOUBLE version
inline double mvp_exp_approx(double x,
                             std::ostream* pstream__ = nullptr) {

     return fast_exp_1(x);

}

 



// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_exp_approx(const Eigen::Matrix<double, -1, 1> &x,
                                                   std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_exp_approx", "x", x);
  
  return fn_EIGEN_double(x, "exp", vect_type, false);
  
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_exp_approx(const Eigen::Matrix<double, 1, -1> &x,
                                                   std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_exp_approx", "x", x);
  
  return fn_EIGEN_double(x, "exp", vect_type, false);
  
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_exp_approx(const Eigen::Matrix<double, -1, -1> &x,
                                                    std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_exp_approx", "x", x);
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_exp_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher." );
  }
  
  return fn_EIGEN_double(x, "exp", vect_type, false);
  
} 


// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_exp_approx(const Eigen::EigenBase<Derived> &x,
                           std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "exp", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) { 
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "exp", vect_type, false);
    } else {  
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "exp", vect_type, false);
    }
     
 
    
    
} 
 
     
    
}  // namespace math
}  // namespace stan

 
 
 

//////////////  --------- log  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {

 


///// DOUBLE version
inline double mvp_log_approx(double x,
                             std::ostream* pstream__ = nullptr) {
  
  return fast_log_1(x);
  
}



 

// 
// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_log_approx(const Eigen::Matrix<double, -1, 1> &x,
                                                    std::ostream* pstream__ = nullptr) {

  //check_finite("mvp_log_approx", "x", x);
  return fn_EIGEN_double(x, "log", vect_type, false);

}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_log_approx(const Eigen::Matrix<double, 1, -1> &x,
                                                   std::ostream* pstream__ = nullptr) {

 // check_finite("mvp_log_approx", "x", x);
  return fn_EIGEN_double(x, "log", vect_type, false);

}

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1>  mvp_log_approx(const Eigen::Matrix<double, -1, -1> &x,
                                                     std::ostream* pstream__ = nullptr) {

  //check_finite("mvp_log_approx", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_log_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }

  return fn_EIGEN_double(x, "log", vect_type, false);

}



 // also overload for Eigen expressions
 template <typename Derived>
 inline auto mvp_log_approx(const Eigen::EigenBase<Derived> &x,
                            std::ostream* pstream__ = nullptr)
   -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {

     /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
     if (Derived::ColsAtCompileTime == 1) {
       Eigen::Matrix<double, -1, 1> temp = x;
       return fn_EIGEN_double(temp, "log", vect_type, false);
     } else if (Derived::RowsAtCompileTime == 1) { 
       Eigen::Matrix<double, 1, -1> temp = x;
       return fn_EIGEN_double(temp, "log", vect_type, false);
     } else {  
       Eigen::Matrix<double, -1, -1> temp = x;
       return fn_EIGEN_double(temp, "log", vect_type, false);
     }
     
      
     
   }
 
 
 

// 


}  // namespace math
}  // namespace stan




//////////////  --------- log1p  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {
        
        
///// DOUBLE version
inline double mvp_log1p_approx(double x,
                               std::ostream* pstream__ = nullptr) {
  
  return fast_log1p_1(x);
  
}



 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_log1p_approx(const Eigen::Matrix<double, -1, 1> &x,
                                                     std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_log1p_approx", "x", x);
  
  return fn_EIGEN_double(x, "log1p", vect_type, false);
  
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_log1p_approx(const Eigen::Matrix<double, 1, -1> &x,
                                                     std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_log1p_approx", "x", x);
  
  return fn_EIGEN_double(x, "log1p", vect_type, false);
  
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_log1p_approx(const Eigen::Matrix<double, -1, -1> &x,
                                                      std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_log1p_approx", "x", x);
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_log1p_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  
  return fn_EIGEN_double(x, "log1p", vect_type, false); 
  
}


 
 // also overload for Eigen expressions  
 template <typename Derived>
 inline auto mvp_log1p_approx(const Eigen::EigenBase<Derived> &x,
                                                std::ostream* pstream__ = nullptr)
   -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
     
     /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
     if (Derived::ColsAtCompileTime == 1) {
       Eigen::Matrix<double, -1, 1> temp = x;
       return fn_EIGEN_double(temp, "log1p", vect_type, false);
     } else if (Derived::RowsAtCompileTime == 1) { 
       Eigen::Matrix<double, 1, -1> temp = x;
       return fn_EIGEN_double(temp, "log1p", vect_type, false);
     } else {  
       Eigen::Matrix<double, -1, -1> temp = x;
       return fn_EIGEN_double(temp, "log1p", vect_type, false);
     }
     
     
     
}
  


}  // namespace math
}  // namespace stan




//////////////  --------- log1m  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {

        

///// DOUBLE version
inline double mvp_log1m_approx(double x,
                               std::ostream* pstream__ = nullptr) {
  
  return fast_log1m_1(x);
   
}
 

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1>
mvp_log1m_approx(const Eigen::Matrix<double, -1, 1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_log1m_approx", "x", x);
  return fn_EIGEN_double(x, "log1m", vect_type, false);
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1>
mvp_log1m_approx(const Eigen::Matrix<double, 1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_log1m_approx", "x", x);
  return fn_EIGEN_double(x, "log1m", vect_type, false);
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1>
mvp_log1m_approx(const Eigen::Matrix<double, -1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_log1m_approx", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_log1m_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  return fn_EIGEN_double(x, "log1m", vect_type, false); 
}



        
// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_log1m_approx(const Eigen::EigenBase<Derived> &x,
                                               std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "log1m", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {  
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "log1m", vect_type, false);
    } else {   
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "log1m", vect_type, false);
    }
    
     
    
    
  } 
        


}  // namespace math
}  // namespace stan




//////////////  --------- logit  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {



///// DOUBLE version
inline double mvp_logit_approx(double x,
                              std::ostream* pstream__ = nullptr) {
  
  return fast_tanh(x);
  
}

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_logit_approx(const Eigen::Matrix<double, -1, 1> &x,
                                                     std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_logit_approx", "x", x); 
  
  return fn_EIGEN_double(x, "logit", vect_type, false);
  
}
 
// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_logit_approx(const Eigen::Matrix<double, 1, -1> &x,
                                                     std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_logit_approx", "x", x); 
  
  return fn_EIGEN_double(x, "logit", vect_type, false);
  
} 
 
// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_logit_approx(const Eigen::Matrix<double, -1, -1> &x,
                                                      std::ostream* pstream__ = nullptr) { 
  
  check_finite("mvp_logit_approx", "x", x);
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_logit_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  } 
  
  return fn_EIGEN_double(x, "logit", vect_type, false); 
  
} 


// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_logit_approx(const Eigen::EigenBase<Derived> &x,
                             std::ostream* pstream__ = nullptr) 
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {

    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "logit", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {  
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "logit", vect_type, false);
    } else {   
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "logit", vect_type, false);
    }
     
    
    
}





}  // namespace math
}  // namespace stan










//////////////  --------- tanh  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {



///// DOUBLE version
inline double mvp_tanh_approx(double x,
                       std::ostream* pstream__ = nullptr) {

       return fast_tanh(x);

}

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_tanh_approx(const Eigen::Matrix<double, -1, 1> &x,
                                                    std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_tanh_approx", "x", x);
  
  return fn_EIGEN_double(x, "tanh", vect_type, false);
  
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_tanh_approx(const Eigen::Matrix<double, 1, -1> &x,
                                                    std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_tanh_approx", "x", x);
  
  return fn_EIGEN_double(x, "tanh", vect_type, false);
  
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_tanh_approx(const Eigen::Matrix<double, -1, -1> &x,
                                                     std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_tanh_approx", "x", x);
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_tanh_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  
  return fn_EIGEN_double(x, "tanh", vect_type, false); 
  
}


// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_tanh_approx(const Eigen::EigenBase<Derived> &x,
                            std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "tanh", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {  
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "tanh", vect_type, false);
    } else {   
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "tanh", vect_type, false);
    }
     
    
    
    
  } 
 




}  // namespace math
}  // namespace stan







//////////////  --------- inv_logit  --------------------------------------------------------------------------------------------------------------------


namespace stan {
namespace math {




///// DOUBLE version
inline double mvp_inv_logit(double x,
                            std::ostream* pstream__ = nullptr) {
  
  return fast_inv_logit(x);
  
}
 

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_inv_logit(const Eigen::Matrix<double, -1, 1> &x,
                                                  std::ostream* pstream__ = nullptr) { 
  
  check_finite("mvp_inv_logit", "x", x);
  
  return fn_EIGEN_double(x, "inv_logit", vect_type, false);
  
} 

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_inv_logit(const Eigen::Matrix<double, 1, -1> &x,
                                                  std::ostream* pstream__ = nullptr) { 
  
  check_finite("mvp_inv_logit", "x", x);
  
  return fn_EIGEN_double(x, "inv_logit", vect_type, false);
  
} 
 
// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_inv_logit(const Eigen::Matrix<double, -1, -1> &x,
                                                   std::ostream* pstream__ = nullptr) { 
  
  check_finite("mvp_inv_logit", "x", x);
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_inv_logit"
                                  " called with vector input. This may be a bug in Stan's dispatcher."); 
  } 
  
  return fn_EIGEN_double(x, "inv_logit", vect_type, false);
  
}  



// also overload for Eigen expressions
template <typename Derived>
inline auto mvp_inv_logit(const Eigen::EigenBase<Derived> &x,
                          std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {

    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "inv_logit", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_logit", vect_type, false);
    } else {
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_logit", vect_type, false);
    }




  }

 


}  // namespace math
}  // namespace stan











//////////////  --------- log_inv_logit  --------------------------------------------------------------------------------------------------------------------


namespace stan {
namespace math {





///// DOUBLE version
inline double mvp_log_inv_logit(double x,
                                std::ostream* pstream__ = nullptr) {
   
  return fast_log_inv_logit(x);
  
} 


 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_log_inv_logit(const Eigen::Matrix<double, -1, 1> &x,
                                                      std::ostream* pstream__ = nullptr) { 
  
  check_finite("mvp_log_inv_logit", "x", x);
  
  return fn_EIGEN_double(x, "log_inv_logit", vect_type, false); 
  
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_log_inv_logit(const Eigen::Matrix<double, 1, -1> &x,
                                                      std::ostream* pstream__ = nullptr) { 
  
  check_finite("mvp_log_inv_logit", "x", x);
  
  return fn_EIGEN_double(x, "log_inv_logit", vect_type, false);
  
}
 
// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_log_inv_logit(const Eigen::Matrix<double, -1, -1> &x,
                                                       std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_log_inv_logit", "x", x); 
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_log_inv_logit"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  } 
  
  return fn_EIGEN_double(x, "log_inv_logit", vect_type, false); 
  
} 



// also overload for Eigen expressions  
template <typename Derived> 
inline auto mvp_log_inv_logit(const Eigen::EigenBase<Derived> &x,
                              std::ostream* pstream__ = nullptr) 
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "log_inv_logit", vect_type, false); 
    } else if (Derived::RowsAtCompileTime == 1) { 
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "log_inv_logit", vect_type, false); 
    } else {  
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "log_inv_logit", vect_type, false);
    }
     
    
    
    
  } 
 




}  // namespace math
}  // namespace stan








//////////////  --------- Phi  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {


///// DOUBLE version
inline double mvp_Phi(double x,
                      std::ostream* pstream__ = nullptr) {

  return fast_Phi_wo_checks(x);

}



 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1> mvp_Phi(const Eigen::Matrix<double, -1, 1> &x,
                                            std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_Phi", "x", x);
  
  return fn_EIGEN_double(x, "Phi", vect_type, false);
  
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1> mvp_Phi(const Eigen::Matrix<double, 1, -1> &x,
                                            std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_Phi", "x", x);
  
  return fn_EIGEN_double(x, "Phi", vect_type, false);
  
}
 
// Matrix - must come last
inline Eigen::Matrix<double, -1, -1> mvp_Phi(const Eigen::Matrix<double, -1, -1> &x,
                                             std::ostream* pstream__ = nullptr) {
  
  check_finite("mvp_Phi", "x", x);
  
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_Phi"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  
  return fn_EIGEN_double(x, "Phi", vect_type, false); 
  
}


// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_Phi(const Eigen::EigenBase<Derived> &x,
                    std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "Phi", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {  
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "Phi", vect_type, false);
    } else {   
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "Phi", vect_type, false);
    }
     
    
    
    
  } 
 



}  // namespace math
}  // namespace stan



//////////////  --------- inv_Phi  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {



///// DOUBLE version
inline double mvp_inv_Phi(double x,
                          std::ostream* pstream__ = nullptr) {

  return fast_inv_Phi_wo_checks(x);

}

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1>
mvp_inv_Phi(const Eigen::Matrix<double, -1, 1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi", "x", x);
  return fn_EIGEN_double(x, "inv_Phi", vect_type, false);
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1>
mvp_inv_Phi(const Eigen::Matrix<double, 1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi", "x", x);
  return fn_EIGEN_double(x, "inv_Phi", vect_type, false);
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1>
mvp_inv_Phi(const Eigen::Matrix<double, -1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_inv_Phi"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  return fn_EIGEN_double(x, "inv_Phi", vect_type, false); 
}


// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_inv_Phi(const Eigen::EigenBase<Derived> &x,
                                               std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) { 
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi", vect_type, false);
    } else {   
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi", vect_type, false);
    }
     
    
    
    
  } 
 


}  // namespace math
}  // namespace stan





//////////////  --------- Phi_approx  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {




///// DOUBLE version
inline double mvp_Phi_approx(double x,
                             std::ostream* pstream__ = nullptr) {

  return fast_Phi_approx(x);

}

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1>
mvp_Phi_approx(const Eigen::Matrix<double, -1, 1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_Phi_approx", "x", x);
  return fn_EIGEN_double(x, "Phi_approx", vect_type, false);
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1>
mvp_Phi_approx(const Eigen::Matrix<double, 1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_Phi_approx", "x", x);
  return fn_EIGEN_double(x, "Phi_approx", vect_type, false);
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1>
mvp_Phi_approx(const Eigen::Matrix<double, -1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_Phi_approx", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_Phi_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  return fn_EIGEN_double(x, "Phi_approx", vect_type, false); 
}



// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_Phi_approx(const Eigen::EigenBase<Derived> &x,
                                               std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "Phi_approx", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) { 
      Eigen::Matrix<double, 1, -1> temp = x; 
      return fn_EIGEN_double(temp, "Phi_approx", vect_type, false);
    } else {  
      Eigen::Matrix<double, -1, -1> temp = x; 
      return fn_EIGEN_double(temp, "Phi_approx", vect_type, false);
    }
     
    
    
  } 
 




}  // namespace math
}  // namespace stan








//////////////  --------- inv_Phi_approx  --------------------------------------------------------------------------------------------------------------------


namespace stan {
namespace math {







///// DOUBLE version
inline double mvp_inv_Phi_approx(double x,
                                 std::ostream* pstream__ = nullptr) {

  return fast_inv_Phi_approx(x);

}

 


// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1>
mvp_inv_Phi_approx(const Eigen::Matrix<double, -1, 1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi_approx", "x", x);
  return fn_EIGEN_double(x, "inv_Phi_approx", vect_type, false);
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1>
mvp_inv_Phi_approx(const Eigen::Matrix<double, 1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi_approx", "x", x);
  return fn_EIGEN_double(x, "inv_Phi_approx", vect_type, false);
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1>
mvp_inv_Phi_approx(const Eigen::Matrix<double, -1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi_approx", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_inv_Phi_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  return fn_EIGEN_double(x, "inv_Phi_approx", vect_type, false); 
}


// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_inv_Phi_approx(const Eigen::EigenBase<Derived> &x,
                                               std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi_approx", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {  
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi_approx", vect_type, false);
    } else {  
      Eigen::Matrix<double, -1, -1> temp = x; 
      return fn_EIGEN_double(temp, "inv_Phi_approx", vect_type, false);
    }
    
     
    
  } 
 



}  // namespace math
}  // namespace stan








//////////////  --------- inv_Phi_approx_from_logit_prob  --------------------------------------------------------------------------------------------------------------------
namespace stan {
namespace math {




 ///// DOUBLE version
 inline double mvp_inv_Phi_approx_from_logit_prob(double x,
                                  std::ostream* pstream__ = nullptr) {

   return fast_inv_Phi_approx_from_logit_prob(x);

 }




// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1>
mvp_inv_Phi_approx_from_logit_prob(const Eigen::Matrix<double, -1, 1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi_approx_from_logit_prob", "x", x);
  return fn_EIGEN_double(x, "inv_Phi_approx_from_logit_prob", vect_type, false);
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1>
mvp_inv_Phi_approx_from_logit_prob(const Eigen::Matrix<double, 1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi_approx_from_logit_prob", "x", x);
  return fn_EIGEN_double(x, "inv_Phi_approx_from_logit_prob", vect_type, false);
} 

// Matrix - must come last 
inline Eigen::Matrix<double, -1, -1>
mvp_inv_Phi_approx_from_logit_prob(const Eigen::Matrix<double, -1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_inv_Phi_approx_from_logit_prob", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_inv_Phi_approx_from_logit_prob"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  }
  return fn_EIGEN_double(x, "inv_Phi_approx_from_logit_prob", vect_type, false); 
}



// also overload for Eigen expressions  
template <typename Derived>
inline auto mvp_inv_Phi_approx_from_logit_prob(const Eigen::EigenBase<Derived> &x,
                                               std::ostream* pstream__ = nullptr)
  -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
    
    /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
    if (Derived::ColsAtCompileTime == 1) {
      Eigen::Matrix<double, -1, 1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi_approx_from_logit_prob", vect_type, false);
    } else if (Derived::RowsAtCompileTime == 1) {  
      Eigen::Matrix<double, 1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi_approx_from_logit_prob", vect_type, false);
    } else {   
      Eigen::Matrix<double, -1, -1> temp = x;
      return fn_EIGEN_double(temp, "inv_Phi_approx_from_logit_prob", vect_type, false);
    } 
    
    
    

} 
 



}  // namespace math
}  // namespace stan









//////////////  --------- log_Phi_approx  --------------------------------------------------------------------------------------------------------------------


namespace stan {
namespace math {




///// DOUBLE version
inline double mvp_log_Phi_approx(double x,
                                std::ostream* pstream__ = nullptr) {

  return fast_log_Phi_approx(x);

}
 



// Column vector - must come before matrix (note: using Eigen::Ref here for overloads creattes abiguuity errors so avoid it)
inline Eigen::Matrix<double, -1, 1>
mvp_log_Phi_approx(const Eigen::Matrix<double, -1, 1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_log_Phi_approx", "x", x);
  return fn_EIGEN_double(x, "log_Phi_approx", vect_type, false);
}

// Row vector - must come before matrix
inline Eigen::Matrix<double, 1, -1>
mvp_log_Phi_approx(const Eigen::Matrix<double, 1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_log_Phi_approx", "x", x);
  return fn_EIGEN_double(x, "log_Phi_approx", vect_type, false);
} 

// Matrix - must come last
inline Eigen::Matrix<double, -1, -1>
mvp_log_Phi_approx(const Eigen::Matrix<double, -1, -1> &x,
               std::ostream* pstream__ = nullptr) {
  check_finite("mvp_log_Phi_approx", "x", x);
  // Add runtime check to prevent vectors being passed here
  if (x.cols() == 1 || x.rows() == 1) {
    throw std::invalid_argument("Matrix version of mvp_log_Phi_approx"
                                  " called with vector input. This may be a bug in Stan's dispatcher.");
  } 
  return fn_EIGEN_double(x, "log_Phi_approx", vect_type, false);
}
 
 

 
 // also overload for Eigen expressions  
 template <typename Derived>
 inline auto mvp_log_Phi_approx(const Eigen::EigenBase<Derived> &x,
                                                std::ostream* pstream__ = nullptr)
   -> std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value, Eigen::Matrix<double, -1, -1>> {
     
     /// now look at rows/ccols at copmpile time to apply correct fn (as we have seperate overloads depending on rows/cols @ compile time)
     if (Derived::ColsAtCompileTime == 1) {
       Eigen::Matrix<double, -1, 1> temp = x;
       return fn_EIGEN_double(temp, "log_Phi_approx", vect_type, false);
     } else if (Derived::RowsAtCompileTime == 1) { 
       Eigen::Matrix<double, 1, -1> temp = x;
       return fn_EIGEN_double(temp, "log_Phi_approx", vect_type, false);
     } else {  
       Eigen::Matrix<double, -1, -1> temp = x;
       return fn_EIGEN_double(temp, "log_Phi_approx", vect_type, false);
     }
     
      
     
     
   } 
  




}  // namespace math
}  // namespace stan


 








 
#endif
    
     

 
