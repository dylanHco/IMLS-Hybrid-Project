#### Load libraries to make plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
# Set working directory
setwd("~/Desktop/Oaks_working_folder/R_Oaks")

#Input the data
admix_data <- read.table("Oak_AD_6.txt", header=FALSE)

#Name the Columns for adding labels later
colnames(admix_data) <- c("Sample", "Cluster1", "Cluster2", "Cluster3","Cluster4", "Cluster5","Cluster6")

# Cluster Names
cluster_labels <- c("Cluster1" = "Montana", 
                    "Cluster2" = "Alba", 
                    "Cluster3" = "Margarettae", 
                    "Cluster4" = "Boyntonii_Stellata", 
                    "Cluster5" = "Muehlenbergii",
                    "Cluster6" = "Boyntonii_Marg_Stell")

###### Order samples within each cluster
library(tidyverse)

# Convert to long format
admix_long <- pivot_longer(admix_data, cols = starts_with("Cluster"), 
                           names_to = "Cluster", values_to = "Ancestry")

# Determine the primary cluster for each sample (where they have the highest ancestry proportion)
admix_summary <- admix_long %>%
  group_by(Sample) %>%
  summarize(MainCluster = Cluster[which.max(Ancestry)],  # Find the dominant cluster
            MaxAncestry = max(Ancestry))                # Store the highest ancestry proportion

# Order samples within each main cluster from least to most admixed (low MaxAncestry first)
admix_summary <- admix_summary %>%
  arrange(MainCluster, desc(MaxAncestry))  # Reverse the order using `desc()`

# Update factor levels for Sample to reflect new order
admix_long$Sample <- factor(admix_long$Sample, levels = admix_summary$Sample)



### Oak plot
ggplot(admix_long, aes(x = Sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),        # Increase legend text size
    legend.title = element_text(size = 16),       # Increase legend title size
    legend.key.size = unit(.2, "cm")               # Increase size of color boxes
  ) +
  scale_fill_manual(
    values = c(
      "#E69F00",  # orange
      "#56B4E9",  # sky blue
      "#009E73",  # bluish green
      "#F0E442",  # yellow
      "#0072B2",  # blue
      "#D55E00",  # vermillion
      "#CC79A7",  # reddish purple
      "#999999",  # gray
      "#332288",  # dark blue
      "#88CCEE",  # light blue
      "#117733",  # dark green
      "#DDCC77",  # sand
      "#CC6677",  # red-pink
      "#AA4499",  # purple
      "#44AA99",  # teal
      "#882255",  # wine
      "#661100",  # brown
      "#6699CC"   # steel blue
    ),
    labels = cluster_labels,
    name = "Species"
  ) +
  labs(title = "ADMIXTURE Plot (K=6_Oaks)", x = "Samples", y = "Ancestry Proportion")

#Adjust width heigh and dpi for plat
ggsave("Oak6b.png", width = 24, height = 12, dpi = 500)
dev.off()
