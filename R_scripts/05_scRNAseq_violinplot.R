
load("data16genesNew.rda")

library(ggplot2)
library(dplyr)
library(pals)

df <- all %>% 
  mutate(conca = factor(conca))

palette25 <- cols25()

ggplot(df, aes(x = conca, y = CORO2A, fill = conca, color = conca)) +
  geom_violin(trim = FALSE, alpha = 1) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 1) +
  scale_fill_manual(values = palette25) +
  scale_color_manual(values = palette25) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(
    x = "",
    y = "CORO2A expression",
    title = ""
  )
