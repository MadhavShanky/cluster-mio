# theme_scs.R — publication ggplot theme + palette + figure helpers for the SCS paper.
# Aesthetic: clean, minimal, colorblind-aware. Hero method (L0-MIO) = forest; emphasis = coral.
suppressMessages({ library(ggplot2); library(scales) })

# ---- palette ---------------------------------------------------------------
scs_cols <- c(
  "L0-MIO"        = "#2D5A3D",  # forest — the hero method
  "L0-MIO+stab"   = "#1F7A4D",  # lighter forest (stability variant)
  "LMM-Gaussian"  = "#4C6E9C",  # muted blue
  "LMM-Laplace"   = "#8A6BBF",  # muted purple
  "OLS"           = "#9A9A9A",  # neutral grey
  "Truth"         = "#222222"
)
scs_emphasis <- "#C2453E"       # coral — at most one emphasis per figure

scale_color_scs <- function(...) scale_color_manual(values = scs_cols, ...)
scale_fill_scs  <- function(...) scale_fill_manual(values = scs_cols, ...)

# diverging fill for "advantage" phase diagrams: negative (shrinkage wins) coral -> 0 grey -> positive
# (L0 wins) forest. `limit` symmetrizes the scale around 0.
scale_fill_advantage <- function(limit = NULL, name = "L0 advantage\n(log2 ratio)", ...) {
  scale_fill_gradient2(low = scs_emphasis, mid = "grey92", high = "#2D5A3D",
                       midpoint = 0, limits = limit, oob = scales::squish, name = name, ...)
}

# ---- theme -----------------------------------------------------------------
theme_scs <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
      axis.title       = element_text(face = "plain"),
      axis.text        = element_text(colour = "grey25"),
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(colour = "grey35", size = rel(0.9)),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      strip.background = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = rel(0.85)),
      legend.key.height = unit(0.8, "lines"),
      plot.caption     = element_text(colour = "grey55", size = rel(0.75), hjust = 0)
    )
}

# ---- figure helpers --------------------------------------------------------
# F1-style phase diagram: tile of `value` over x (sparsity/pi) by y (K), optionally faceted.
scs_phase_tile <- function(df, x, y, value, diverging = TRUE, facet = NULL,
                           title = NULL, subtitle = NULL, xlab = NULL, ylab = "clusters K",
                           value_name = NULL, limit = NULL) {
  p <- ggplot(df, aes(x = factor(.data[[x]]), y = factor(.data[[y]]), fill = .data[[value]])) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", .data[[value]])), size = 3, colour = "grey15") +
    labs(title = title, subtitle = subtitle, x = xlab %||% x, y = ylab) +
    theme_scs() + theme(panel.grid.major = element_blank())
  p <- p + if (diverging) scale_fill_advantage(limit = limit, name = value_name %||% "advantage")
           else scale_fill_gradient(low = "grey92", high = "#2D5A3D", name = value_name %||% "value")
  if (!is.null(facet)) p <- p + facet_wrap(facet)
  p
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# consistent device saver (PNG + optionally PDF) at paper-ready size/res.
scs_save <- function(plot, file, w = 6.5, h = 4.2, dpi = 300, pdf = FALSE) {
  ggsave(file, plot, width = w, height = h, dpi = dpi)
  if (pdf) ggsave(sub("\\.png$", ".pdf", file), plot, width = w, height = h)
  invisible(file)
}
