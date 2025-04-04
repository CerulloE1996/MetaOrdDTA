


////
//// Custom Stan functions - simplex functions  ---------------------------------------------------------------
//// NOTE: adapted from:
//// https://github.com/mjhajharia/transforms/blob/main/transforms/simplex/StickbreakingLogistic.stan
////

 
////
//// Function to manually construct a simplex from raw unconstrained parameters.
////
vector stickbreaking_logistic_simplex_constrain_jacobian(vector y) {
            
            //// Other variables:
            int N = rows(y) + 1;
            vector[N] x;
            real log_zi;
            real log_xi;
            real log_cum_prod = 0;
            
            //// Construct simplex:
            for (i in 1:(N - 1)) {
                log_zi = log_inv_logit(y[i] - log(N - i)); // logistic_lcdf(y[i] | log(N - i), 1)
                log_xi = log_cum_prod + log_zi;
                x[i] = exp(log_xi);
                log_cum_prod += log1m_exp(log_zi);
                jacobian += log_xi;
            }   
            
            x[N] = exp(log_cum_prod);
            
            //// Increment jacobian adjustment:
            jacobian += log_cum_prod;
            
            //// Output simplex:
            return x;

}



// 
//   vector stickbricking_angular_simplex_constrain_lp(vector y) {
//     int N = rows(y) + 1;
//     vector[N] x;
//     real log_phi, phi, u, s, c;
//     real s2_prod = 1;
//     real log_halfpi = log(pi()) - log2();
//     int rcounter = 2 * N - 3;
//     for (i in 1:(N-1)) {
//       u = log_inv_logit(y[i]);
//       log_phi = u + log_halfpi;
//       phi = exp(log_phi);
//       s = sin(phi);
//       c = cos(phi);
//       x[i] = s2_prod * c^2;
//       s2_prod *= s^2;
//       target += log_phi + log1m_exp(u) + rcounter * log(s) + log(c);
//       rcounter -= 2;
//     }
//     
//     x[N] = s2_prod;
//     target += (N - 1) * log2();
//     return x;
//     
//   }