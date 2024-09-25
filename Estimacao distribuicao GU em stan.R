# Carregando pacotes necessários
library(rstan)       # Interface com o Stan para modelos Bayesianos
library(tidyverse)   # Pacote de manipulação de dados e visualização

# Configurações do Stan para escrever modelos compilados no disco
rstan_options(auto_write = TRUE)

# Detectando o número de núcleos de CPU disponíveis para paralelização
options(mc.cores = parallel::detectCores())

# Definindo o diretório de trabalho e verificando a localização atual
setwd("C:/Users/Lenovo/Desktop/Mestrado/2023.1/Recursos Computacionais/V0/Stan")
getwd()

# Carregando pacotes adicionais
library(stats)

# Função para calcular a densidade da distribuição unitária gama (dugamma)
dugamma = function(x, mu, sigma, log.p = FALSE){
  
  # Verificação se o parâmetro mu está dentro do intervalo permitido (0, 1)
  if (any(mu <= 0 | mu >= 1))
    warning("The mean parameter must be on (0, 1).")
  
  # Transformando vetor x em matriz
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  
  n <- dim(x)[1]  # Número de linhas (observações)
  d <- dim(x)[2]  # Número de colunas (variáveis)
  
  # Ajustando mu e sigma para o mesmo tamanho de x
  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  
  # Verificando conformidade das dimensões
  if(n %% dim(mu)[1] != 0)
    stop("x and mean have non-conforming size")
  if(n %% dim(sigma)[1] != 0)
    stop("x and sigma have non-conforming size")
  
  # Replicando mu e sigma para igualar o tamanho de x
  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))
  
  # Inicializando matriz de probabilidades com -Inf
  pmf <- matrix(-Inf, n, d)
  
  # Índices de valores inválidos para mu
  NaNindex <- rbind(which(mu <= 0 | mu >= 1, arr.ind = TRUE))
  pmf[NaNindex] <- NaN
  
  # Índices de valores válidos para o cálculo da densidade
  index <- which(0 < x & x < 1 & !is.nan(pmf), arr.ind = TRUE)
  
  # Calculando densidade para os valores válidos
  mu.star <- mu[index]^(1 / sigma[index])
  dp <- mu.star / (1 - mu.star)
  
  # Fórmula da densidade
  pmf[index] <- (sigma[index]) * log(dp) +
    (dp - 1) * log(x[index]) +
    (sigma[index] - 1) * log(log(1 / x[index])) -
    lgamma(sigma[index])
  
  # Retornando a densidade, no logaritmo se solicitado
  if (n == 1){
    pmf <- as.numeric(pmf)
  }
  
  if (log.p == TRUE){
    pmf
  } else {
    exp(pmf)
  }
}

# Função para calcular a quantil da distribuição unitária gama (qugamma)
qugamma = function(p, mu, sigma, lower.tail = TRUE){
  
  # Verificação se mu está dentro do intervalo permitido
  if (any((mu <= 0) | (mu >= 1)))
    warning("The mean parameter must be on the unit interval (0, 1).")
  
  # Transformando vetor p em matriz
  if (is.vector(p))
    p <- matrix(p, ncol = length(p))
  
  n <- dim(p)[1]  # Número de linhas (observações)
  d <- dim(p)[2]  # Número de colunas (variáveis)
  
  # Ajustando mu e sigma para o mesmo tamanho de p
  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  
  # Verificando conformidade das dimensões
  if(n %% dim(mu)[1] != 0)
    stop("p and mean have non-conforming size")
  if(n %% dim(sigma)[1] != 0)
    stop("p and sigma have non-conforming size")
  
  # Replicando mu e sigma para igualar o tamanho de p
  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))
  
  quanti <- matrix(0, n, d)  # Inicializando a matriz de quantis
  
  # Índices de valores inválidos
  NaNindex <- rbind(which(mu <= 0 | mu >= 1, arr.ind = TRUE))
  quanti[NaNindex] <- NaN
  
  # Índices de valores válidos para o cálculo do quantil
  index <- which(0 < p & p < 1 & !is.nan(quanti), arr.ind = TRUE)
  
  # Ajuste para a cauda superior, se necessário
  if(lower.tail == FALSE)
    p[index] <- 1 - p[index]
  
  # Calculando o quantil
  alpha <- mu[index]^(1 / sigma[index])
  dp <- alpha / (1 - alpha)
  quanti[index] <- exp(-stats::qgamma(1 - p[index], sigma[index], dp))
  
  # Retornando o resultado como vetor, se for o caso
  if (n == 1)
    quanti <- as.numeric(quanti)
  
  quanti
}

# Função para gerar valores aleatórios da distribuição unitária gama (rugamma)
rugamma = function(n, mu, sigma, rep = 1L){
  ifelse(rep == 1L, u <- stats::runif(n),
         u <- matrix(stats::runif(n * rep), ncol = n))
  
  qugamma(u, mu, sigma)
}

# Gerando dados simulados da distribuição unitária gama
set.seed(20222)
N = 200
dat = rugamma(N, 0.2, 5)

# Visualizando os dados simulados com um gráfico de densidade
qplot(dat, geom = 'density')

###############################################################################

# Ajustando o modelo Stan
fit = stan(file = 'unit_gamma2.stan', data = list(Y = dat, n = N), 
           warmup = 1000, iter = 2500, chains = 3, verbose = FALSE)

# Exibindo os resultados do ajuste
print(fit)

###############################################################################

# Carregando pacote para visualização interativa de modelos Stan
library(shinystan)

# Lançando o visualizador ShinyStan
launch_shinystan(fit)

###############################################################################

# Extraindo amostras da distribuição posterior
fit_value = rstan::extract(fit)
mu_p = fit_value$mu
phi_p = fit_value$phi

# Convertendo para o formato coda para análise das cadeias MCMC
library(coda)
line.coda = as.mcmc(as.matrix(fit))
dim(line.coda)

# Calculando o desvio padrão ajustado pelo tamanho efetivo da amostra
round(apply(line.coda, 2, sd) / (sqrt(coda::effectiveSize(line.coda))), 6)

# Gerando novos dados a partir das estimativas posteriores
y_rep = matrix(0, N, N)
for (i in 1:N) {
  y_rep[i,] <- rugamma(N, mu_p[i], phi_p[i])
}

# Visualizando sobreposição de densidades entre dados observados e gerados
library(bayesplot)
ppc_dens_overlay(dat, y_rep[1:100, ], color_scheme_set("darkgray")) + xlim(0, 1)

# Alternativa de visualização com verificação posterior (pp_check)
pp_check(stan_data$y, fun = 'dens_overlay')
