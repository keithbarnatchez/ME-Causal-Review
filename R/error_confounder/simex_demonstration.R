set_plot_theme <- function() {
  # theme_set(theme_bw() +
  #             theme(plot.title = element_text(hjust = 0, size = 16),
  #                   plot.subtitle = element_text(hjust = 0, size = 12),
  #                   axis.title = element_text(size = 12),
  #                   strip.text = element_text(size = 12),
  #                   legend.position = "bottom"))
  ggplot2::theme_set(
    ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.caption = ggplot2::element_text(hjust = 0),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        plot.subtitle = ggplot2::element_text(
          margin = ggplot2::margin(b = 10)
        ),
        strip.text = ggplot2::element_text(
          face = "bold",
          hjust = 0,
          margin = ggplot2::margin(b = 5),
          size = 11
        )
      )
  )
}

set_plot_theme()


# Create a data frame with x values
df <- data.frame(x = seq(-1, 2, by = 0.01))

# Calculate y values based on the function
df$y <- ((df$x - 3)^2) / 9

# Create a data frame for the jittered points
points_df <- data.frame(x = c(0, 0.4, 0.8, 1.2, 1.6, 2.0))

# Calculate y values for the jittered points
points_df$y <- ((points_df$x - 3)^2) / 9

point_lab <- TeX('$\\beta$')

extr_point <- data.frame(x=-1, y=(16)/9)

# Create the plot
simexplot <- ggplot(df, aes(x, y)) +
  geom_line(aes(linetype = ifelse(x < 0, "dashed", "solid"))) +
  geom_point(data = points_df, aes(x, y), position = position_jitter(width = 0.05), color = "black") +
  geom_text(data = points_df[1,], aes(x, y, label = 'Naive estimate', vjust = -1, hjust=-.1)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_point(data=extr_point, aes(x,y),color='red') +
  geom_text(data=extr_point, aes(x,y,label='SIMEX estimate', hjust=-.2)) +
  theme(legend.position = 'none') +
  labs(title='SIMEX bias correction',
       y='Coefficient estimate',
       x=expression(lambda)) ; simexplot

ggsave('../../output/figures/simex_demonstration.pdf',
       width=5,height=3,units='in')

