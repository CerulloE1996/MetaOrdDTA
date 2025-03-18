


////
//// ---- Box-Cox custom Stan functions: ----------------------------------------------------------------------------------------------
////
////    
//// ---- Box-Cox transform function:
////
real fn_Stan_box_cox( real x, 
                      real lambda) {
  
  if (lambda == 0.0) {
    if (x != 0.0) {
      return log(x);
    } else {
      return -700.0;
    }
  } else {
    return (pow(x, lambda) - 1.0) / lambda;
  }
  
}
////
//// ---- Vectorized Box-Cox transform function (overloaded function):
////
vector fn_Stan_box_cox( vector x, 
                        real lambda) {
      
      int N = num_elements(x);
      vector[N] result;
      
      for (n in 1:N) {
        if (lambda == 0.0) {
            if (x[n] != 0.0) {
              result[n] = log(x[n]);
            } else {
              result[n] = -700.0;
            }
        } else {
            result[n] = (pow(x[n], lambda) - 1.0) / lambda;
        }
      }
      
      return result;
      
}
////
//// ---- Row-vector overload:
////
row_vector fn_Stan_box_cox( row_vector x, 
                            real lambda) { 
    
      int N = num_elements(x);
      vector[N] x_col_vec = to_vector(x);
      vector[N] x_col_vec_trans = fn_Stan_box_cox(x_col_vec, lambda);
      row_vector[N] x_row_vec_trans = to_row_vector(x_col_vec_trans);
      
      return(x_row_vec_trans);
  
}
