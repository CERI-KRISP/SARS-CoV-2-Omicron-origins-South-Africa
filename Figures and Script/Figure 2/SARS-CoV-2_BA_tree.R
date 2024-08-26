# Required libraries
library(ggplot2)
library("readxl")
library(ggtree)
library(viridis)

setwd("~/Downloads/Omicron analysis/Full_analysis/ggtree")

# Load tree (from treefile output)
treefile <- read.newick("")
treefile$tip.label <- gsub("'", "", treefile$tip.label) # Remove single quotation marks from the tree labels for joining to metadata file

# Load metadata
metadata <- read.delim("", header = TRUE, sep = "\t")

## Alternative method for joining the tree and metadata
#tree_tbl <- as_tibble(tree) # Convert the tree object to a tibble with ggtree (to be able to join the two)
#tree_metadata <- left_join(tree_tbl, metadata, by = c("label" = "lineage")) # Join tree and metadata
#tree_annotated <- as.phylo(tree_metadata) # Convert back to a tree object and plot

# Customize and plot the tree
tree <- ggtree(treefile, mrsd = "", as.Date = TRUE) +
  #geom_tiplab() +  # Add tip labels
  #geom_tippoint(aes(color = metadata$clade)) +   # Add tip points
  theme_tree2() +   # Use a cleaner theme
  labs(color = "Sub-lineage") +   # Legend title
  theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +  # Legend position
  scale_x_date(date_labels = "%Y")  # Format the x-axis as dates

TreeAnnotated <- tree %<+% metadata +
  #geom_tippoint(aes(color = lineage))+ # Tip colour
  geom_tree(aes(color = lineage, group = lineage)) + # Line colour
  scale_color_viridis(discrete = TRUE, option="D")

# Plot the tree
print(TreeAnnotated) # Use 'dev.off()' if plot is not showing

# Save the tree to a file
ggsave("", TreeAnnotated, width = 11, height = 8)
