
title: "biometrics_ESD_plots"
author: "MMRocker"
date: "06 June 2025"
output: html_document




#set up working file & R

getwd() #check working directory


library(ggplot2)
library(readxl)
library(lubridate)
require(grid)
library(gridExtra)
library(reshape2)
library(dplyr) #NO plyr
library(patchwork)
library(ggtext)




st<-read_excel("ESDExpt_biometrics.xlsx",sheet="biometrics")

st<-as.data.frame(st)

colnames(st)

st1 <- as.data.frame(st[,c(1:5,7,9,8)])

hist(st1$initialcoverage)
hist(st1$CarapaceLength1)
hist(st1$CF1)


st2 <- as.data.frame(st[,c(1:5,6)])
st2$Source_Pop <- st2$Source_Pop
st2$Source_Pop <- as.character(st2$Source_Pop)

barplot(table(st2$MoltStage1))




# Reshape to long format for ggplot - disease coverage, CL, & K
long_df <- st1 %>%
  tidyr::pivot_longer(cols = c(CarapaceLength1, CF1, initialcoverage), names_to = "variable", values_to = "value")


# For bar plot: summarize category counts and label it as variable = "MoltStage1)
bar_df <- st2 %>%
  count(Source_Pop, MoltStage1) %>%
  rename(group = MoltStage1, subgroup = Source_Pop, value = n) %>%
  mutate(variable = "MoltStage1")

# Rename for consistency with long_df
colnames(bar_df) <- c("Source_Pop","group", "value", "variable")  # trick to match names

# Combine boxplot data and barplot data
plot_df <- bind_rows(long_df, bar_df)

plot_df <- plot_df %>%
  mutate(
    xvar = ifelse(variable == "MoltStage1", group, Source_Pop),
    xvar = as.factor(xvar)
  )




# Subset for boxplots
box_df <- subset(plot_df, variable != "MoltStage1")
box_vars <- unique(box_df$variable)

# Subset for barplot
bar_df <- subset(plot_df, variable == "MoltStage1")



#rename & reorder states
box_df$Source_Pop <- plyr::revalue(box_df$Source_Pop, c("ME" = "Casco Bay",
                                                  "MA" = "Buzzards Bay",
                                                  "CT" = "Niantic Bay"))
box_df$Source_Pop <- factor(box_df$Source_Pop, levels = c("Casco Bay",
                                                          "Buzzards Bay",
                                                          "Niantic Bay"))
#rename & reorder states
bar_df$Source_Pop <- plyr::revalue(bar_df$Source_Pop, c("ME" = "Casco Bay, ME",
                                                  "MA" = "Buzzards Bay, MA",
                                                  "CT" = "Niantic Bay, CT"))
bar_df$Source_Pop <- factor(bar_df$Source_Pop, levels = c("Casco Bay, ME",
                                                          "Buzzards Bay, MA",
                                                          "Niantic Bay, CT"))


# Generate individual plots
p1 <- ggplot(subset(box_df, variable == box_vars[1]), aes(x = Source_Pop, y = value, fill = Source_Pop)) +
  geom_boxplot(alpha =0.8) +
  geom_point(aes(fill = Source_Pop), shape = 21, size = 0.75) +
  theme_bw() + 
  labs( x = "", y = "") +
  scale_fill_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  scale_color_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  theme(legend.position = "none",
        #plot.title = element_text(hjust = 0.5, size = 10),
        plot.title = element_text("none"),
        plot.margin = margin(t = 5))  # avoid overlap
p1

# Create a separate "title" plot with a gray background
title_1 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Carapace length (mm)", size = 3) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

# Combine using patchwork
final1 <- title_1 / p1 +
  plot_layout(heights = c(0.15, 1))  # Adjust height ratio of title to plot
final1


p2 <- ggplot(subset(box_df, variable == box_vars[2]), aes(x = Source_Pop, y = value, fill = Source_Pop)) +
  geom_boxplot(alpha =0.8) +
  geom_point(aes(fill = Source_Pop), shape = 21, size = 0.75) +
  theme_bw() + labs( x = "", y = "") +
  scale_fill_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  scale_color_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  theme(legend.position = "none",
        #plot.title = element_text(hjust = 0.5, size = 10),
        plot.title = element_text("none"),
        plot.margin = margin(t = 5))  # avoid overlap

# Create a separate "title" plot with a gray background
title_2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Condition factor (K)", size = 3) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

# Combine using patchwork
final2 <- title_2 / p2 +
  plot_layout(heights = c(0.15, 1))  # Adjust height ratio of title to plot
final2


p3 <- ggplot(subset(box_df, variable == box_vars[3]), aes(x = Source_Pop, y = value, fill = Source_Pop)) +
  geom_boxplot(alpha =0.8) +
  geom_point(aes(fill = Source_Pop), shape = 21, size = 0.75) +
  theme_bw() + labs( x = "", y = "") +
  scale_fill_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  scale_color_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  theme(legend.position = "none",
        #plot.title = element_text(hjust = 0.5, size = 10),
        plot.title = element_text("none"),
        plot.margin = margin(t = 5))  # avoid overlap

# Create a separate "title" plot with a gray background
title_3 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "ESD severity (%)", size = 3) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

# Combine using patchwork
final3 <- title_3 / p3 +
  plot_layout(heights = c(0.15, 1))  # Adjust height ratio of title to plot
final3


p4 <- ggplot(bar_df, aes(x = group, y = value, fill = Source_Pop)) +
  geom_col(position = "stack", alpha = 0.85, color = "black") +
  theme_bw() + labs( x = "", y = "") +
  scale_fill_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  theme(legend.position = "right",legend.direction = "vertical",
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs( fill = "Source Population")

p4_noleg <- ggplot(bar_df, aes(x = group, y = value, fill = Source_Pop)) +
  geom_col(position = "stack", alpha = 0.85, color = "black") +
  theme_bw() + labs( x = "", y = "") +
  scale_fill_manual(
    values = viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
  ) +
  theme(legend.position = "none",
        #plot.title = element_text(hjust = 0.5, size = 10),
        plot.title = element_text("none"),
        plot.margin = margin(t = 5))  # avoid overlap

p4_noleg

# Create a separate "title" plot with a gray background
title_4 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Molt stage (n)", size = 3) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

# Combine using patchwork
final4 <- title_4 / p4_noleg +
  plot_layout(heights = c(0.15, 1))  # Adjust height ratio of title to plot
final4

# Combine in 2x2 layout
((p1 + p2) / (p3 + p4)) #+ theme(legend.position = "bottom") # & plot_layout(guides = "collect") &






# Load libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)

# Get US states and world map
states <- ne_states(country = "united states of america", returnclass = "sf")
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define main map limits
main_xlim <- c(-74, -66)
main_ylim <- c(40.5, 45)

# Define bay locations
stars <- data.frame(
  bay = c("Casco Bay", "Buzzards Bay", "Niantic Bay"),
  lon = c(-70.05, -70.80, -72.26),
  lat = c(43.65, 41.45, 41.31)
)
stars_sf <- st_as_sf(stars, coords = c("lon", "lat"), crs = 4326)
stars_sf$bay <- factor(stars_sf$bay, levels = c("Casco Bay", "Buzzards Bay", "Niantic Bay"))

# Define colors
star_colors <- viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
names(star_colors) <- c("Casco Bay", "Buzzards Bay", "Niantic Bay")


# Main map
main_map <- ggplot() +
  geom_sf(data = states, fill = "gray90", color = "black", size = 0.2) +
  geom_sf(data = stars_sf, aes(fill = bay), shape = 23, size = 4) +  # no fixed color here
  coord_sf(xlim = main_xlim, ylim = main_ylim, expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3, bar_cols = c("black", "white"), text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  scale_fill_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT")) +  # apply your viridis colors here
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white"),
    legend.position = c(0.95, 0.5),
    legend.justification = c("right", "center"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.box.background = element_blank()) +
  labs( fill = "Source Population")

main_map



#plotting all together:
bottom_grid <- plot_grid(final1, final2, final3, final4, 
                         ncol = 2#,  # 2 x 2 layout
                         #labels = c("B", "C", "D", "E")
                         )
final_plot <- plot_grid(
  main_map,         # your map
  bottom_grid,      # 4 plots underneath
  ncol = 1,         # stack vertically
  rel_heights = c(1, 1.25)#,  # control height ratio if needed
  #labels = c("A", "")     # optional: label map as "A"
)

final_plot

# Wrap with ggdraw() and set white background
final_plot_white <- ggdraw(final_plot) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

final_plot_white

ggsave("BiometricsAndMap_Fig1.png",
       device = "png",
       height = 6, # in inches, but you can set the units
       width = 6,
       dpi = 300)



