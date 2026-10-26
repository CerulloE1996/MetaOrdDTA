////
//// ---- Sample an induced-Dirichlet cutpoint prior from a target distribution:
////
//// Each row contains one complete ordered cutpoint draw from the target distribution.
//// The anchor and link are supplied explicitly. The likelihood is the joint induced-Dirichlet
//// density of the whole cutpoint vector, including the change-of-variables Jacobian.
//// There are no hyperpriors, interval penalties or restrictions on the concentration beyond positivity.
//// HMC includes the parameter-transform Jacobian and uses the flat measure on positive concentrations.
//// Posterior propriety and sampling diagnostics must be checked for the supplied target and free inputs.
//// This model derives the cutpoint distribution itself; it does not integrate another model's
//// random locations, random scales or between-study cutpoint hierarchy.
//// Concentrations belong to the intervals separated by the supplied cutpoints. Omitted raw-score
//// cutpoints are not reconstructed, and their separate category concentrations are not identified.
////


data {

      int<lower=2> n_draws;
      int<lower=1> n_thr;
      array[n_draws] ordered[n_thr] cutpoint_draws;
      real anchor;
      int<lower=0, upper=1> use_probit_link;

}


transformed data {

      int n_cat = n_thr + 1;
      array[n_draws] vector[n_cat] log_category_probability;
      vector[n_draws] log_cutpoint_Jacobian;
      int target_varies = 0;

      if (is_nan(anchor) || is_inf(anchor)) {
            reject("The cutpoint anchor must be finite.");
      }

      for (draw_i in 1:n_draws) {

            vector[n_thr] C = cutpoint_draws[draw_i] - anchor;
            vector[n_thr] log_F;
            vector[n_thr] log_S;

            for (k in 1:n_thr) {

                  if (is_nan(C[k]) || is_inf(C[k])) {
                        reject("The anchored target cutpoints must be finite: draw ", draw_i, ", cutoff ", k);
                  }

                  if (use_probit_link == 1) {
                        log_F[k] = std_normal_lcdf(C[k] |);
                        log_S[k] = std_normal_lccdf(C[k] |);
                  } else {
                        log_F[k] = log_inv_logit(C[k]);
                        log_S[k] = log1m_inv_logit(C[k]);
                  }

            }

            log_category_probability[draw_i, 1] = log_F[1];
            log_category_probability[draw_i, n_cat] = log_S[n_thr];

            if (n_thr > 1) {
                  for (k in 2:n_thr) {
                        log_category_probability[draw_i, k] = (C[k - 1] >= 0)
                              ? log_diff_exp(log_S[k - 1], log_S[k])
                              : log_diff_exp(log_F[k], log_F[k - 1]);
                  }
            }

            log_cutpoint_Jacobian[draw_i] = (use_probit_link == 1)
                  ? std_normal_lpdf(C)
                  : sum(log_F + log_S);

            for (k in 1:n_cat) {
                  if (is_nan(log_category_probability[draw_i, k])
                        || is_inf(log_category_probability[draw_i, k])) {
                        reject("The target category log probability is not finite: draw ", draw_i,
                               ", category ", k, ". No probabilities are clipped or replaced.");
                  }

                  if (draw_i > 1) {
                        if (log_category_probability[draw_i, k] != log_category_probability[1, k]) {
                              target_varies = 1;
                        }
                  }

            }

            if (is_nan(log_cutpoint_Jacobian[draw_i]) || is_inf(log_cutpoint_Jacobian[draw_i])) {
                  reject("The target cutpoint Jacobian is not finite: draw ", draw_i);
            }

      }

      if (target_varies == 0) {
            reject("Identical target draws have no finite maximum-likelihood Dirichlet concentration.");
      }

}


parameters {

      vector<lower=0>[n_cat] alpha;

}


model {

      real log_normalising_constant = lgamma(sum(alpha)) - sum(lgamma(alpha));
      vector[n_cat] alpha_minus_one = alpha - rep_vector(1.0, n_cat);

      for (draw_i in 1:n_draws) {
            target += log_normalising_constant
                  + dot_product(alpha_minus_one, log_category_probability[draw_i])
                  + log_cutpoint_Jacobian[draw_i];
      }

}






















