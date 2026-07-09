rm(list = ls())

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

B <- 5000
K <- 5
T_total <- 1000
sigma <- 1
prior_mean <- 0
prior_var <- 100
c_ucb <- 2

pasta_saida <- "resultados_parte1"
if (!dir.exists(pasta_saida)) dir.create(pasta_saida)

politicas <- c("Uniforme", "Greedy", "UCB", "Thompson Sampling")

calc_post <- function(soma, n) {
  v <- 1 / (1 / prior_var + n / sigma^2)
  m <- v * (prior_mean / prior_var + soma / sigma^2)
  list(m = m, v = v)
}

simula <- function(politica, seed) {
  set.seed(seed)
  soma <- matrix(rnorm(B * K, 0, sigma), B, K)
  n <- matrix(1, B, K)
  id <- 1:B
  
  for (t in (K + 1):T_total) {
    if (politica == "Uniforme") {
      a <- sample(1:K, B, replace = TRUE)
    }
    
    if (politica == "Greedy") {
      a <- max.col(soma / n, ties.method = "first")
    }
    
    if (politica == "UCB") {
      p <- calc_post(soma, n)
      ucb <- p$m + c_ucb * sqrt(p$v)
      a <- max.col(ucb, ties.method = "first")
    }
    
    if (politica == "Thompson Sampling") {
      p <- calc_post(soma, n)
      amostra <- matrix(rnorm(B * K, as.vector(p$m), sqrt(as.vector(p$v))), B, K)
      a <- max.col(amostra, ties.method = "first")
    }
    
    r <- rnorm(B, 0, sigma)
    pos <- cbind(id, a)
    soma[pos] <- soma[pos] + r
    n[pos] <- n[pos] + 1
  }
  
  medias <- soma / n
  braco <- max.col(medias, ties.method = "first")
  n_alpha <- n[cbind(id, braco)]
  muhat_alpha <- medias[cbind(id, braco)]
  li <- muhat_alpha - 1.96 / sqrt(n_alpha)
  ls <- muhat_alpha + 1.96 / sqrt(n_alpha)
  cobriu <- as.integer(li <= 0 & 0 <= ls)
  
  data.frame(
    politica = politica,
    braco_selecionado = braco,
    muhat_alpha = muhat_alpha,
    n_alpha = n_alpha,
    limite_inferior = li,
    limite_superior = ls,
    cobriu = cobriu
  )
}

resultados <- data.frame()

for (i in 1:length(politicas)) {
  cat("Rodando", politicas[i], "\n")
  resultados <- rbind(resultados, simula(politicas[i], 12345 + 1000 * i))
}

resumo <- data.frame()

for (pol in politicas) {
  d <- resultados[resultados$politica == pol, ]
  resumo <- rbind(resumo, data.frame(
    politica = pol,
    cobertura = mean(d$cobriu),
    media_muhat_alpha = mean(d$muhat_alpha),
    dp_muhat_alpha = sd(d$muhat_alpha),
    q025_muhat_alpha = as.numeric(quantile(d$muhat_alpha, 0.025)),
    mediana_muhat_alpha = median(d$muhat_alpha),
    q975_muhat_alpha = as.numeric(quantile(d$muhat_alpha, 0.975)),
    media_n_alpha = mean(d$n_alpha),
    dp_n_alpha = sd(d$n_alpha),
    q025_n_alpha = as.numeric(quantile(d$n_alpha, 0.025)),
    mediana_n_alpha = median(d$n_alpha),
    q975_n_alpha = as.numeric(quantile(d$n_alpha, 0.975))
  ))
}

freq <- prop.table(table(resultados$politica, resultados$braco_selecionado), 1)
freq_df <- as.data.frame.matrix(freq)
freq_df$politica <- rownames(freq_df)
rownames(freq_df) <- NULL
freq_df <- freq_df[, c("politica", setdiff(names(freq_df), "politica"))]

write.csv(resultados, file.path(pasta_saida, "resultados_completos.csv"), row.names = FALSE)
write.csv(resumo, file.path(pasta_saida, "resumo_parte1.csv"), row.names = FALSE)
write.csv(freq_df, file.path(pasta_saida, "frequencia_braco_selecionado.csv"), row.names = FALSE)

g1 <- ggplot(resumo, aes(x = politica, y = cobertura)) +
  geom_col() +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(ylim = c(0.80, 1.00)) +
  labs(title = "Cobertura empírica do intervalo ingênuo de 95%", x = "Política de alocação", y = "Cobertura empírica") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(pasta_saida, "fig_cobertura.png"), g1, width = 9, height = 6, dpi = 300)

g2 <- ggplot(resultados, aes(x = muhat_alpha, fill = politica)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.35, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  labs(title = expression(paste("Distribuição de ", hat(mu)[hat(alpha)])), x = expression(hat(mu)[hat(alpha)]), y = "Densidade", fill = "Política") +
  theme_minimal(base_size = 13)

ggsave(file.path(pasta_saida, "fig_muhat_alpha.png"), g2, width = 9, height = 6, dpi = 300)

g3 <- ggplot(resultados, aes(x = politica, y = n_alpha)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression(paste("Tamanho amostral do braço selecionado: ", n[hat(alpha)])), x = "Política de alocação", y = expression(n[hat(alpha)])) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(pasta_saida, "fig_n_alpha.png"), g3, width = 9, height = 6, dpi = 300)

resultados$braco_selecionado <- factor(resultados$braco_selecionado)

g4 <- ggplot(resultados, aes(x = politica, fill = braco_selecionado)) +
  geom_bar(position = "fill") +
  labs(title = "Frequência de seleção final de cada braço", x = "Política de alocação", y = "Frequência relativa", fill = "Braço") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(pasta_saida, "fig_freq_braco_selecionado.png"), g4, width = 9, height = 6, dpi = 300)

cat("\nResumo da Parte 1:\n")
print(resumo)
cat("\nFrequência dos braços:\n")
print(freq_df)
