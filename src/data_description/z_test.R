library(tidyverse)

dist1 <- tibble(x = seq(-7, 10, by = 0.01), y = dnorm(x, mean = 0, sd = 2), dist = "Sample background")
dist2 <- tibble(x = seq(-7, 10, by = 0.01), y = dnorm(x, mean = 3, sd = 2), dist = "Kinase targets")
tib <- bind_rows(dist1, dist2) %>%
  mutate(dist = fct_relevel(dist, "Sample background"))

plot <- tib %>%
  ggplot(mapping=aes(x = x, y = y, fill = dist, color = dist)) +
  geom_line(lwd = 1.5) +
  geom_linerange(data = filter(tib, x == 0, dist == "Sample background"), mapping = aes(x = x, ymin = 0, ymax = y), lwd = 1.5) +
  geom_ribbon(mapping = aes(x = x, ymin = 0, ymax = y), alpha = 0.5) +
  geom_linerange(data = filter(tib, x == 3, dist == "Kinase targets"), mapping = aes(x = x, ymin = 0, ymax = y), lwd = 1.5) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    #axis.text = element_text(size = 24, colour = "black"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 35, colour = "black")) +
  scale_fill_manual(values = c("#636363", "#feb24c")) +
  scale_color_manual(values = c("#636363", "#feb24c"), guide = F)
  #scale_x_continuous(breaks = c(0, 3), labels = c(expression(mu), "X"))

ggsave(filename="z_test.png", plot = plot, path = "./output/plots/data_description/", width = 9, height = 3)
ggsave(filename="z_test.pdf", plot = plot, path = "./output/plots/data_description/", width = 9, height = 3)

