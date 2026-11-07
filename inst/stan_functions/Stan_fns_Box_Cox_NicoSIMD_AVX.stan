// NicoSIMD version of Stan_fns_Box_Cox.stan written by NicoStan (included file md5 d1475e1c5a808be4881b335f219dfe24; NicoSIMD declarations md5 e3c567fd115d2ccf49ec504688009e6b; writer implementation md5 aef1406a74d6694a1d260c4cd31c7a35)



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
      return nico_log(x);
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
              result[n] = nico_log(x[n]);
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


