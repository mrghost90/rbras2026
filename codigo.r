rm(list = ls())
library(lme4)
library(ggplot2)
library(dplyr)

source("Constr_Lat_Primo.R")
source("SBGD.R")
planilha_base <- read.table("planilha_base.csv", header = TRUE)
colnames(planilha_base) <- c("R", "B", "Tr")

set.seed(210826)
p <- 16 # Número de genitores
a <- 2 * qnorm((1:16) / 17)    # Efeitos Aditivos
d <- qnorm((1:120) / 121)      # Efeitos de Dominância

g <- numeric(120)
k <- 0
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    k <- k + 1
    g[k] <- (a[i] + a[j]) / 2 + d[k]
  }
}

sigma_g <- var(g)
v <- 121 # 120 cruzamentos + 1 testemunha

et <- numeric(v)
et[1:120] <- g
et[121] <- qnorm(0.9, mean = 0, sd = sqrt(sigma_g))


X_A <- matrix(0, nrow = v, ncol = p)
k <- 0
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    k <- k + 1
    X_A[k, i] <- 0.5
    X_A[k, j] <- 0.5
  }
}

planilha <- planilha_base
planilha$Local <- ceiling(planilha$R / 2)
n <- nrow(planilha) # 1452 parcelas

X_T <- design(planilha$Tr, v)$design
X_E <- design(planilha$Local, 6)$design
X_R <- design(planilha$R, 12)$design
X_B <- design(planilha$B, 132)$design

eL  <- qnorm((1:6) / 7)
eBL <- qnorm((1:132) / 133)
eR  <- runif(12, min = -0.5, max = 0.5)

yp <- (X_E %*% eL) + (X_R %*% eR) + (X_B %*% eBL) + (X_T %*% et)
y  <- yp + rnorm(n)

planilha$Y <- as.vector(y)
write.csv(planilha, file = "Delineamento_Simulado.csv", row.names = FALSE)

planilha$Tr    <- as.factor(planilha$Tr)
planilha$Local <- as.factor(planilha$Local)
planilha$R     <- as.factor(planilha$R)
planilha$B     <- as.factor(planilha$B)

modelo <- lmer(Y ~ (1 | Local) + (1 | Local:R) + (1 | Local:R:B) + (1 | Tr), data = planilha)

blups_bib <- ranef(modelo)$Tr
ghat_bib  <- blups_bib[order(as.numeric(rownames(blups_bib))), 1]

a_hat_bib  <- as.vector(solve(t(X_A) %*% X_A) %*% t(X_A) %*% ghat_bib)
gA_hat_bib <- as.vector(X_A %*% a_hat_bib)
d_hat_bib  <- ghat_bib - gA_hat_bib

# Vetores de Erros Quadráticos Individuais (BIB)
eq_vec_bib_a  <- (a - a_hat_bib)^2
eq_vec_bib_d  <- (d - d_hat_bib[1:120])^2
eq_vec_bib_et <- (et - ghat_bib)^2

BIB <- data.frame(
  Delineamento = "Látice Balanceado (BIB)",
  rho_a  = cor(a, a_hat_bib),
  rho_d  = cor(d, d_hat_bib[1:120]),
  rho_et = cor(et, ghat_bib),
  EQM_a  = mean(eq_vec_bib_a),
  EQM_d  = mean(eq_vec_bib_d),
  EQM_et = mean(eq_vec_bib_et)
)

latice_simples_base <- subset(planilha_base, R %in% c(1, 2))
planilha_simples <- data.frame()


for (loc in 1:6) {
  temp_planilha <- latice_simples_base
  temp_planilha$Local <- loc
  temp_planilha$R <- temp_planilha$R + (loc - 1) * 2
  temp_planilha$B <- temp_planilha$B + (loc - 1) * 22
  
  sorteio_tratamentos <- sample(1:v, v, replace = FALSE)
  temp_planilha$Tr <- sorteio_tratamentos[temp_planilha$Tr]
  
  planilha_simples <- rbind(planilha_simples, temp_planilha)
}


a_check <- 2 * qnorm((1:16) / 17)
d_check <- qnorm((1:120) / 121)

X_T <- design(planilha_simples$Tr, v)$design
X_E <- design(planilha_simples$Local, 6)$design
X_R <- design(planilha_simples$R, 12)$design
X_B <- design(planilha_simples$B, 132)$design

eL  <- qnorm((1:6) / 7)
eBL <- qnorm((1:132) / 133)
eR  <- runif(12, min = -0.5, max = 0.5)

yp <- (X_E %*% eL) + (X_R %*% eR) + (X_B %*% eBL) + (X_T %*% et)
y  <- yp + rnorm(n)

planilha_simples$Y <- as.vector(y)
write.csv(planilha_simples, file = "Delineamento_Simulado_LaticeSimples.csv", row.names = FALSE)

planilha_simples$Tr    <- as.factor(planilha_simples$Tr)
planilha_simples$Local <- as.factor(planilha_simples$Local)
planilha_simples$R     <- as.factor(planilha_simples$R)
planilha_simples$B     <- as.factor(planilha_simples$B)

modelo_simples <- lmer(Y ~ (1 | Local) + (1 | Local:R) + (1 | Local:R:B) + (1 | Tr), data = planilha_simples)

blups_bid <- ranef(modelo_simples)$Tr
ghat_bid  <- blups_bid[order(as.numeric(rownames(blups_bid))), 1]

a_hat_bid  <- as.vector(solve(t(X_A) %*% X_A) %*% t(X_A) %*% ghat_bid)
gA_hat_bid <- as.vector(X_A %*% a_hat_bid)
d_hat_bid  <- ghat_bid - gA_hat_bid

# Vetores de Erros Quadráticos Individuais (BID / Simples Independente)
eq_vec_bid_a  <- (a - a_hat_bid)^2
eq_vec_bid_d  <- (d - d_hat_bid[1:120])^2
eq_vec_bid_et <- (et - ghat_bid)^2

BID <- data.frame(
  Delineamento = "Látice Simples Independente",
  rho_a  = cor(a, a_hat_bid),
  rho_d  = cor(d, d_hat_bid[1:120]),
  rho_et = cor(et, ghat_bid),
  EQM_a  = mean(eq_vec_bid_a),
  EQM_d  = mean(eq_vec_bid_d),
  EQM_et = mean(eq_vec_bid_et)
)

# Tabela comparativa global
print(BIB)
print(BID)



#GRAFICOS
df_distribuicao <- bind_rows(
  data.frame(Delineamento = "Látice Balanceado (BIB)",     Parametro = "Genotípico Total (n = 121)", EQ = eq_vec_bib_et),
  data.frame(Delineamento = "Látice Balanceado (BIB)",     Parametro = "Aditivo (n = 16)",           EQ = eq_vec_bib_a),
  data.frame(Delineamento = "Látice Balanceado (BIB)",     Parametro = "Dominância (n = 120)",       EQ = eq_vec_bib_d),
  data.frame(Delineamento = "Látice Simples Independente", Parametro = "Genotípico Total (n = 121)", EQ = eq_vec_bid_et),
  data.frame(Delineamento = "Látice Simples Independente", Parametro = "Aditivo (n = 16)",           EQ = eq_vec_bid_a),
  data.frame(Delineamento = "Látice Simples Independente", Parametro = "Dominância (n = 120)",       EQ = eq_vec_bid_d)
)


resumo_dispersao <- df_distribuicao %>%
  group_by(Parametro, Delineamento) %>%
  summarise(
    Media_EQM = mean(EQ),
    Mediana   = median(EQ),
    IIQ       = IQR(EQ),
    Maximo    = max(EQ),
    .groups   = "drop"
  )
print(resumo_dispersao)




df_limpo <- df_distribuicao %>%
  mutate(
    Delineamento = factor(
      ifelse(grepl("Simples|PBIB", Delineamento), "PBIB", "BIB"),
      levels = c("BIB", "PBIB")
    ),
    Parametro = factor(
      case_when(
        grepl("Aditivo", Parametro) ~ "Aditivo",
        grepl("Dominância", Parametro) ~ "Dominância",
        TRUE ~ "Genotípico Total"
      ),
      levels = c("Aditivo", "Dominância", "Genotípico Total")
    )
  )



df_limpo <- df_distribuicao %>%
  mutate(
    Delineamento = factor(
      ifelse(grepl("Simples|PBIB", Delineamento), "PBIB", "BIB"),
      levels = c("BIB", "PBIB")
    ),
    Parametro = factor(
      case_when(
        grepl("Aditivo", Parametro) ~ "Aditivo",
        grepl("Dominância", Parametro) ~ "Dominância",
        TRUE ~ "Genotípico Total"
      ),
      levels = c("Aditivo", "Dominância", "Genotípico Total")
    )
  )

grafico_limpo <- ggplot(df_limpo, aes(x = Delineamento, y = EQ, fill = Delineamento)) +
  geom_boxplot(width = 0.45, alpha = 0.82, color = "gray15", outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.45, color = "gray15") +
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 3.5,
    fill = "white", color = "black", stroke = 1.2
  ) +
  facet_wrap(~Parametro, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("BIB" = "#073A70", "PBIB" = "#F9BD38")) +
  labs(
    x = NULL,
    y = "Erro Quadrático"
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill = "gray95", color = "gray75"),
    strip.text         = element_text(face = "bold", size = 13, color = "gray10"),
    axis.text.x        = element_text(face = "bold", size = 13, color = "gray10"),
    axis.text.y        = element_text(size = 11, color = "gray20"),
    axis.title.y       = element_text(face = "bold", size = 13, margin = margin(r = 10))
  )

print(grafico_limpo)

ggsave("Distribuicao_EQM_Limpo.png", plot = grafico_limpo, width = 11, height = 4.5, dpi = 300)




library(ggplot2)
library(dplyr)


planos_bib <- list(
  # Ambiente 1 (Plano 1 - Linhas; Plano 2 - Colunas)
  c(1, 2, 3), c(4, 5, 6), c(7, 8, 9),
  c(1, 4, 7), c(2, 5, 8), c(3, 6, 9),
  # Ambiente 2 (Planos 3 e 4 - Quadrados Latinos Ortogonais)
  c(1, 5, 9), c(2, 6, 7), c(3, 4, 8),
  c(1, 6, 8), c(2, 4, 9), c(3, 5, 7)
)

latice_simples_base <- planos_bib[1:6] # Apenas Planos 1 e 2


perm_e2 <- sample(1:9)

planos_e1 <- latice_simples_base
planos_e2 <- lapply(latice_simples_base, function(bloco) perm_e2[bloco])


gerar_grafo <- function(lista_blocos, cx, cy, raio, grupo, subtitulo = NA) {
  Lambda <- matrix(0, nrow = 9, ncol = 9)
  for (b in lista_blocos) {
    Lambda[b, b] <- Lambda[b, b] + 1
  }
  
  angulos <- pi/2 - 2 * pi * (0:8) / 9
  nos <- data.frame(
    id        = 1:9,
    x         = cx + raio * cos(angulos),
    y         = cy + raio * sin(angulos),
    grupo     = grupo,
    subtitulo = subtitulo,
    cx        = cx,
    cy        = cy,
    raio      = raio
  )
  
  arestas <- expand.grid(i = 1:9, j = 1:9) %>%
    filter(i < j) %>%
    mutate(lambda = mapply(function(r, c) Lambda[r, c], i, j)) %>%
    filter(lambda > 0) %>%
    left_join(nos %>% select(id, x1 = x, y1 = y), by = c("i" = "id")) %>%
    left_join(nos %>% select(id, x2 = x, y2 = y), by = c("j" = "id")) %>%
    mutate(grupo = grupo)
  
  list(nos = nos, arestas = arestas)
}

# BIB unificado (E1 + E2 juntos -> 36 arestas, lambda_ij = 1 para todos os pares)
g_bib <- gerar_grafo(planos_bib, cx = -1.25, cy = -0.05, raio = 0.92, grupo = "BIB")

# PBIBs locais (E1 e E2 independentes -> 18 arestas cada, metade dos pares sem contato)
g_e1 <- gerar_grafo(planos_e1, cx = 0.78, cy = -0.05, raio = 0.48, grupo = "PBIB", subtitulo = "E1")
g_e2 <- gerar_grafo(planos_e2, cx = 2.02, cy = -0.05, raio = 0.48, grupo = "PBIB", subtitulo = "E2")

nos_bib     <- g_bib$nos
nos_pbib    <- bind_rows(g_e1$nos, g_e2$nos)
arestas_all <- bind_rows(g_bib$arestas, g_e1$arestas, g_e2$arestas)

rotulos_amb <- nos_pbib %>%
  distinct(subtitulo, cx, cy, raio)

titulos <- data.frame(
  x     = c(-1.25, 1.40),
  y     = c(1.12, 1.12),
  label = c("BIB", "PBIB")
)

grafo_intro <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.6) +
  
 
  geom_segment(
    data = arestas_all,
    aes(x = x1, y = y1, xend = x2, yend = y2, color = grupo, linewidth = grupo),
    alpha = 0.88
  ) +
  
 
  geom_point(
    data = nos_bib,
    aes(x = x, y = y),
    size = 8.5, shape = 21, fill = "#073A70", color = "white", stroke = 1.2
  ) +
  geom_text(
    data = nos_bib,
    aes(x = x, y = y, label = id),
    color = "white", fontface = "bold", size = 4
  ) +
  
  
  geom_point(
    data = nos_pbib,
    aes(x = x, y = y),
    size = 5.2, shape = 21, fill = "#073A70", color = "white", stroke = 0.9
  ) +
  geom_text(
    data = nos_pbib,
    aes(x = x, y = y, label = id),
    color = "white", fontface = "bold", size = 2.6
  ) +
  
  
  geom_text(
    data = rotulos_amb,
    aes(x = cx, y = cy + raio + 0.13, label = subtitulo),
    color = "gray30", fontface = "bold", size = 4
  ) +
  

  geom_text(
    data = titulos,
    aes(x = x, y = y, label = label),
    color = "gray10", fontface = "bold", size = 5.5
  ) +
  
  scale_color_manual(values = c("BIB" = "#073A70", "PBIB" = "#F9BD38")) +
  scale_linewidth_manual(values = c("BIB" = 0.85, "PBIB" = 1.05)) +
  coord_equal(xlim = c(-2.35, 2.65), ylim = c(-1.10, 1.25)) +
  theme_void() +
  theme(legend.position = "none")

print(grafo_intro)
ggsave("Grafo_Concorrencia_Intro.png", plot = grafo_intro, width = 9, height = 4.3, dpi = 300)
