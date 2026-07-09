rm(list = ls())

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

B <- 1000
K <- 5
T_total <- 1000
lambda <- c(18, 20, 21, 23, 27)
alpha0 <- 2
beta0 <- 0.1
c_ucb <- 2
epsilon <- 0.10

pasta_saida <- "resultados_parte2"
if (!dir.exists(pasta_saida)) dir.create(pasta_saida)

politicas <- c("Uniforme", "Greedy", "Thompson Sampling", "UCB", "Epsilon-greedy")

resultados <- data.frame()
curvas <- data.frame(tempo = 1:T_total)

for (i in 1:length(politicas)) {
  politica <- politicas[i]
  set.seed(12345 + 1000 * i)
  cat("Rodando", politica, "\n")
  curva_politica <- matrix(0, B, T_total)
  
  for (b in 1:B) {
    n <- rep(1, K)
    soma <- rpois(K, lambda)
    recompensa_acumulada <- numeric(T_total)
    recompensa_acumulada[1:K] <- cumsum(soma)
    
    for (t in (K + 1):T_total) {
      if (politica == "Uniforme") {
        a <- sample(1:K, 1)
      }
      
      if (politica == "Greedy") {
        a <- which.max(soma / n)
      }
      
      if (politica == "Thompson Sampling") {
        alpha_post <- alpha0 + soma
        beta_post <- beta0 + n
        amostra <- rgamma(K, shape = alpha_post, rate = beta_post)
        a <- which.max(amostra)
      }
      
      if (politica == "UCB") {
        alpha_post <- alpha0 + soma
        beta_post <- beta0 + n
        media_post <- alpha_post / beta_post
        dp_post <- sqrt(alpha_post) / beta_post
        ucb <- media_post + c_ucb * dp_post
        a <- which.max(ucb)
      }
      
      if (politica == "Epsilon-greedy") {
        if (runif(1) < epsilon) {
          a <- sample(1:K, 1)
        } else {
          a <- which.max(soma / n)
        }
      }
      
      r <- rpois(1, lambda[a])
      soma[a] <- soma[a] + r
      n[a] <- n[a] + 1
      recompensa_acumulada[t] <- recompensa_acumulada[t - 1] + r
    }
    
    lambda_hat <- soma / n
    curva_politica[b, ] <- recompensa_acumulada
    
    resultados <- rbind(resultados, data.frame(
      politica = politica,
      simulacao = b,
      recompensa_total = recompensa_acumulada[T_total],
      melhor_braco_estimado = which.max(lambda_hat),
      n1 = n[1],
      n2 = n[2],
      n3 = n[3],
      n4 = n[4],
      n5 = n[5],
      lambda_hat_1 = lambda_hat[1],
      lambda_hat_2 = lambda_hat[2],
      lambda_hat_3 = lambda_hat[3],
      lambda_hat_4 = lambda_hat[4],
      lambda_hat_5 = lambda_hat[5]
    ))
  }
  
  curvas[[politica]] <- colMeans(curva_politica)
}

resumo <- data.frame()

for (pol in politicas) {
  d <- resultados[resultados$politica == pol, ]
  resumo <- rbind(resumo, data.frame(
    politica = pol,
    media_recompensa = mean(d$recompensa_total),
    dp_recompensa = sd(d$recompensa_total),
    q025_recompensa = as.numeric(quantile(d$recompensa_total, 0.025)),
    mediana_recompensa = median(d$recompensa_total),
    q975_recompensa = as.numeric(quantile(d$recompensa_total, 0.975))
  ))
}

alocacoes <- data.frame(politica = politicas)
for (j in 1:K) {
  alocacoes[[paste0("braco_", j)]] <- sapply(politicas, function(pol) mean(resultados[resultados$politica == pol, paste0("n", j)]))
}

prop <- alocacoes
prop[, -1] <- prop[, -1] / T_total

lambda_hat <- data.frame(politica = politicas)
for (j in 1:K) {
  lambda_hat[[paste0("lambda_hat_", j)]] <- sapply(politicas, function(pol) mean(resultados[resultados$politica == pol, paste0("lambda_hat_", j)]))
}

eqm <- data.frame(politica = politicas)
for (j in 1:K) {
  eqm[[paste0("EQM_lambda_", j)]] <- sapply(politicas, function(pol) mean((resultados[resultados$politica == pol, paste0("lambda_hat_", j)] - lambda[j])^2))
}
eqm$EQM_medio <- rowMeans(eqm[, -1])

freq <- prop.table(table(resultados$politica, resultados$melhor_braco_estimado), 1)
freq_df <- as.data.frame.matrix(freq)
freq_df$politica <- rownames(freq_df)
rownames(freq_df) <- NULL
freq_df <- freq_df[, c("politica", setdiff(names(freq_df), "politica"))]

write.csv(resultados, file.path(pasta_saida, "resultados_completos.csv"), row.names = FALSE)
write.csv(resumo, file.path(pasta_saida, "resumo_recompensa.csv"), row.names = FALSE)
write.csv(alocacoes, file.path(pasta_saida, "alocacoes_medias.csv"), row.names = FALSE)
write.csv(prop, file.path(pasta_saida, "proporcao_alocacoes.csv"), row.names = FALSE)
write.csv(lambda_hat, file.path(pasta_saida, "lambda_hat_medias.csv"), row.names = FALSE)
write.csv(eqm, file.path(pasta_saida, "eqm_lambda.csv"), row.names = FALSE)
write.csv(curvas, file.path(pasta_saida, "curvas_recompensa_media.csv"), row.names = FALSE)
write.csv(freq_df, file.path(pasta_saida, "frequencia_melhor_braco_estimado.csv"), row.names = FALSE)

curvas_longas <- data.frame()
for (pol in politicas) {
  curvas_longas <- rbind(curvas_longas, data.frame(tempo = curvas$tempo, politica = pol, recompensa_acumulada_media = curvas[[pol]]))
}

prop_longa <- data.frame()
for (i in 1:nrow(prop)) {
  for (j in 1:K) {
    prop_longa <- rbind(prop_longa, data.frame(politica = prop$politica[i], braco = paste0("Braço ", j), proporcao = prop[i, paste0("braco_", j)]))
  }
}

g1 <- ggplot(resumo, aes(x = politica, y = media_recompensa)) +
  geom_col() +
  labs(title = "Recompensa acumulada média por política", x = "Política de alocação", y = "Recompensa acumulada média") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(pasta_saida, "fig_recompensa_media.png"), g1, width = 9, height = 6, dpi = 300)

g2 <- ggplot(curvas_longas, aes(x = tempo, y = recompensa_acumulada_media, color = politica)) +
  geom_line(linewidth = 1) +
  labs(title = "Recompensa acumulada média ao longo do tempo", x = "Tempo", y = "Recompensa acumulada média", color = "Política") +
  theme_minimal(base_size = 13)

ggsave(file.path(pasta_saida, "fig_recompensa_tempo.png"), g2, width = 9, height = 6, dpi = 300)

g3 <- ggplot(prop_longa, aes(x = politica, y = proporcao, fill = braco)) +
  geom_col(position = "stack") +
  labs(title = "Proporção média de alocação em cada braço", x = "Política de alocação", y = "Proporção média", fill = "Braço") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(pasta_saida, "fig_proporcao_alocacoes.png"), g3, width = 9, height = 6, dpi = 300)

g4 <- ggplot(eqm, aes(x = politica, y = EQM_medio)) +
  geom_col() +
  labs(title = "Erro quadrático médio das estimativas de lambda", x = "Política de alocação", y = "EQM médio") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(pasta_saida, "fig_eqm_lambda.png"), g4, width = 9, height = 6, dpi = 300)

cat("\nResumo da recompensa acumulada:\n")
print(resumo)
cat("\nProporção média de alocações:\n")
print(prop)
cat("\nEstimativas médias finais de lambda:\n")
print(lambda_hat)
cat("\nEQM:\n")
print(eqm)
cat("\nFrequência do melhor braço estimado:\n")
print(freq_df)
