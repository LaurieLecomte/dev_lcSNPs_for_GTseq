ggplot(data = maf0.05_pctind0.65_maxdepth10_combined_0.1_scored_100.good.5) +
+     facet_wrap( ~V1) +
+     geom_point(aes(x = V2, y = V2), alpha = 0.5, size = 0.5) +
+     theme(panel.spacing = unit(0.1, 'points'),
+           strip.text.x = element_text(size = 6),
+           axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
+           panel.background = element_rect(color = "gray60"),
+           strip.placement = "inside",
+           strip.background = element_rect(colour = 'gray60')
+     ) + 
+     # scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
+     scale_x_continuous(
+         labels = function(x) {
+             round(x/10^8, 1)
+         }
+     )