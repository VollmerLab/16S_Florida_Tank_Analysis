florida_plot <- ggplot(data = earth) +
  geom_sf(fill = "#D8CEBC") +
  geom_rect(aes(xmin = -82, xmax = -79.5, ymin = 24, ymax = 26.5), col = "red", fill = NA, linewidth = 0.75) +
  coord_sf(xlim = c(-88, -78), ylim = c(23.5, 33), expand = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "#B1B1B1", linetype = "dashed", 
                                        size = 0.2), panel.background = element_rect(fill = "#B9ECFF"),
        plot.background = element_rect(fill='transparent', color = NA),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#only the FL Keys, colored by region
fl_keys_region_plot <- ggplot(data = earth) +
  geom_sf(fill = "#D8CEBC") +
  coord_sf(xlim = c(-82, -79.5), ylim = c(24, 26.5), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "#B1B1B1", linetype = "dashed", 
                                        size = 0.2), 
        panel.background = element_rect(fill = "#B9ECFF"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(color = "Region")

# Final Region Plot (combine Keys and FL)
hmm <- ggdraw() +
  draw_plot(fl_keys_region_plot) +
  draw_plot(florida_plot, x = .47, y = 0.075, width = .3, height = .3) +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

ggsave("Documents/florida_envdata_transparent.png", hmm, bg = "transparent", 
       width = 17.5, height = 7, units = "in")




