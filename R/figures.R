library("ggplot2")
plot_data <- performance_all %>%
  dplyr::select(
    model,
    share_RMSFE,
    share_RMSE,
    share_Correct_Sparsity,
    share_Model_Inlcuded,
    share_Relevant
  ) %>%
  pivot_longer(
    cols = -model,
    names_to = "metric",
    values_to = "share"
  ) %>%
  mutate(
    share = 100 * share,
    metric = recode(
      metric,
      share_RMSFE = "RMSFE",
      share_RMSE = "RMSE",
      share_Correct_Sparsity = "Correct sparsity",
      share_Model_Inlcuded = "True model included",
      share_Relevant = "Relevant coefficients\nselected"
    ),
    model = recode(
      model,
      "aLASSO_ridge" = "aLASSO(ridge)",
      "aLASSO_lasso" = "aLASSO(lasso)",
      "LASSO" = "Lasso"
    ),
    model = factor(
      model,
      levels = c(
        "BAR",
        "aLASSO(ridge)",
        "aLASSO(lasso)",
        "Lasso",
        "SCAD"
      )
    )
  )

ggplot(
  plot_data,
  aes(
    x = share,
    y = factor(
      metric,
      levels = c(
        "Relevant coefficients\nselected",
        "True model included",
        "Correct sparsity",
        "RMSE",
        "RMSFE"
      )
    ),
    fill = model
  )
) +
  geom_col(
    position = position_dodge(width = 0.75,reverse = TRUE),
    width = 0.5,
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = sprintf("%.1f", share)),
    position = position_dodge(width = 0.75,reverse=TRUE),
    hjust = -0.15,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "BAR" = "royalblue",
      "aLASSO(ridge)" = "orange",
      "aLASSO(lasso)" = "forestgreen",
      "Lasso" = "firebrick",
      "SCAD" = "purple"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    #title = "Random sign eigenvalues",
    x = "Average win rate (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 10)
  )