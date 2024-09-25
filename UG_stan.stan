functions {
  // Log-verossimilhança da distribuição Gama Unitária
  real unit_gamma_lpdf(real y, real mu, real phi) {
    return phi * log(mu^(1/phi) / (1 - mu^(1/phi))) - lgamma(phi) + 
           ((mu^(1/phi) / (1 - mu^(1/phi))) - 1) * log(y) + 
           (phi - 1) * log(-log(y));
  }
}

data {
  int<lower=1> N;                    // Tamanho da amostra
  vector<lower=0, upper=1>[N] y;      // Variável resposta (entre 0 e 1)
}

parameters {
  real<lower=0, upper=1> mu;           // Parâmetro de média
  real<lower=0> phi;                   // Parâmetro de dispersão
}

model {
  
  // Priori para mu
  mu ~ uniform(0, 1);
  
  // Priori para o parâmetro de dispersão phi (gamma com baixa informação)
  phi ~ gamma(0.001, 0.001);
  
  // Verossimilhança para cada observação
  for (i in 1:N) {
    target += unit_gamma_lpdf(y[i] | mu, phi); // Distribuição Gama Unitária
  }
}
