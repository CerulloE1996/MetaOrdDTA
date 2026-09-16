




data {
  int<lower=0> N; // Number of observations
  vector[N] y; // Observations
  real<lower=0> sigma; // Measurement variability
  int<lower=0> K;
  array[N] int<lower=1, upper=K> indiv_idx;
  vector<lower=0, upper=1>[K] w; // weights 
}


parameters {
  real mu;            // Population location
  real<lower=0> tau;  // Population scale
  vector[K] z;      // Non-centered individual parameters
}

transformed parameters {
  
 vector[K] tau_wt = tau - w * (tau - 1);
 vector[K] theta = (z * tau + mu * w) ./ tau_wt;  // Recentered individual parameters
}


model {
  mu ~ normal(0, 5);                   // Prior model
  tau ~ normal(0, 5);                  // Prior model
  
  z ~ normal(mu * (1 - w), tau_wt);
  
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}










// The code is copyrighted by Michael Betancourt and licensed under the new BSD (3-clause) license

data {
  int<lower=0> N; // Number of observations
  vector[N] y; // Observations
  real<lower=0> sigma; // Measurement variability
  int<lower=0> K;   // Number of individual contexts in hierarchy
  array[N] int<lower=1, upper=K> indiv_idx;  // Individual context from which each observation is generated
}

parameters {
  real mu; // Population location
  real<lower=0> tau; // Population scale
  vector[K] z; // Non-centered individual parameters
}

transformed parameters {
 
  vector[K] theta = mu + tau * z;  // Recentered individual parameters
}

model {
  mu ~ normal(0, 5); // Prior model
  tau ~ normal(0, 5); // Prior model
  
  z ~ std_normal(); // Non-centered hierarchical model
  
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}







