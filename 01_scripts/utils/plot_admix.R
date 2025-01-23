# Plot NGSadmix results

library(tidyr)
library(ggplot2)

# 1. Access files ---------------------------------------------------------
argv <- commandArgs(T)
ADMIX <- argv[1]
ID_POP <- argv[2]

## Order of samples in ID_POP should be the same as in bamlist
ID_POP <- read.delim(ID_POP, header=FALSE)
admix <- read.table(ADMIX, 
                    quote="\"", comment.char="")


## Add IP and POP info 
admix$ID <-ID_POP$V1
admix$POP <-ID_POP$V2

# Convert to long format
admix_long <- pivot_longer(data = admix,
             cols = colnames(admix)[grepl(colnames(admix), pattern = 'V')],
             names_to = 'anc',
             values_to = 'admix')


# 2. Plot -----------------------------------------------------------------
admix_plot <- 
ggplot(data = admix_long) + 
  facet_grid(~POP, scales = 'free_x') +
  geom_bar(aes(fill = anc, x = ID, y = admix), position = 'fill', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, size = 2.6, hjust = 1),
        panel.spacing.x = unit(1, 'points')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = 'D')

# Export
ggsave(admix_plot, 
       filename = paste0(unlist(strsplit(ADMIX, split = '.qopt')), ".pdf"),
       device = "pdf", dpi = 320,
       height = 6, width = 5, units = "in")

