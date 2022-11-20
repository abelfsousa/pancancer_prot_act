library(tidyverse)


nci60 <- read_tsv(file = "./data/nci60_crc65/metadata/nci60_cell_lines.txt") %>%
  select(cell=Name, panel=`Cell line panel`, tissue=`Tissue of origin`) %>%
  mutate(tissue = str_replace(tissue, "Colon", "Colorectal"))
crc65 <- read_tsv(file = "./data/nci60_crc65/metadata/crc65_cell_lines.txt") %>%
  select(cell=Name, panel=`Cell line panel`) %>%
  mutate(tissue = "Colorectal")

cells <- bind_rows(nci60, crc65)


panel_barplot <- cells %>%
  mutate(panel = fct_rev(fct_infreq(panel))) %>%
  ggplot(mapping = aes(x=panel, fill = panel)) +
  geom_bar() +
  theme() +
  coord_flip() +
  labs(x = "Panel", y = "Count", fill = "")

ggsave(filename = "nci60crc65_cells.png", plot = panel_barplot, path = "./output/plots/nci60_crc65/", width = 6, height = 2)
unlink("nci60crc65_cells.png")


tissue_barplot <- cells %>%
  mutate(tissue = fct_rev(fct_infreq(tissue)), panel = fct_relevel(panel, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x=tissue, fill = panel)) +
  geom_bar(position = "dodge") +
  theme() +
  coord_flip() +
  labs(x = "Tissue", y = "Count", fill = "Panel")

ggsave(filename = "nci60crc65_tissue.png", plot = tissue_barplot, path = "./output/plots/nci60_crc65/", width = 6, height = 3)
unlink("nci60crc65_tissue.png")

