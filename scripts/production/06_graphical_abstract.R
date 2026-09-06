# Reproducible conceptual graphical abstract for the MethodsX MOB workflow.
#
# The layout is designed at the manuscript display size (6.8 x 3.4 inches),
# not only at the larger export size. Short, explicitly wrapped labels prevent
# headings and body text from crossing panel boundaries in Word.

build_methodsx_graphical_abstract <- function(colors) {
  required_colors <- c(
    "blue", "gold", "orange", "olive", "ink", "mid_grey", "light_grey"
  )
  if (!all(required_colors %in% names(colors))) {
    stop(
      "The graphical abstract requires these named colors: ",
      paste(required_colors, collapse = ", "),
      call. = FALSE
    )
  }

  alpha_fill <- function(color, alpha = 0.12) {
    grDevices::adjustcolor(unname(color), alpha.f = alpha)
  }

  max_line_length <- function(x) {
    max(vapply(
      strsplit(x, "\n", fixed = TRUE),
      function(lines) max(nchar(lines, type = "width")),
      numeric(1)
    ))
  }

  title_text <- "HOW DOES NETWORK STRUCTURE CHANGE?"
  subtitle_text <- paste(
    "Compare two complementary evidence families",
    "at each annual transition"
  )
  footer_text <- paste(
    "Agreement strengthens the claim;",
    "disagreement identifies the change."
  )

  panel_titles <- c(
    "ANNUAL\nNETWORKS",
    "NETWORK\nINDICATORS",
    "GRAPHON\nDISTANCES",
    "CALIBRATE\n& CHECK",
    "INTERPRET"
  )
  panel_bodies <- c(
    paste(
      "19 annual snapshots",
      "2006-2024",
      "175 aligned nodes",
      sep = "\n"
    ),
    paste(
      "Activity  |  Density",
      "Reciprocity  |  Concentration",
      sep = "\n"
    ),
    paste(
      "Raw scale  |  Shape",
      "Frobenius  |  Cut bound",
      "Spectra  |  Roles",
      sep = "\n"
    ),
    paste(
      "Known-truth",
      "simulations",
      "",
      "Channel-specific scores",
      "",
      "Estimator + resolution",
      "checks",
      sep = "\n"
    ),
    paste(
      "AGREEMENT",
      "Broad structural change",
      "",
      "DISAGREEMENT",
      "Type + scale of change",
      sep = "\n"
    )
  )

  if (
    max_line_length(panel_titles) > 18L ||
      max_line_length(panel_bodies) > 29L ||
      nchar(title_text, type = "width") > 40L ||
      nchar(subtitle_text, type = "width") > 70L ||
      nchar(footer_text, type = "width") > 75L
  ) {
    stop(
      "A graphical-abstract label exceeds the manuscript-size layout limit.",
      call. = FALSE
    )
  }

  box_data <- data.frame(
    xmin = c(0.40, 3.70, 3.70, 8.00, 11.90),
    xmax = c(3.00, 7.30, 7.30, 11.20, 15.60),
    ymin = c(1.35, 4.65, 1.35, 1.35, 1.35),
    ymax = c(7.55, 7.55, 4.15, 7.55, 7.55),
    fill = c(
      alpha_fill(colors["ink"], 0.055),
      alpha_fill(colors["blue"], 0.13),
      alpha_fill(colors["orange"], 0.13),
      alpha_fill(colors["gold"], 0.15),
      alpha_fill(colors["olive"], 0.15)
    ),
    border = unname(colors[c("mid_grey", "blue", "orange", "gold", "olive")]),
    stringsAsFactors = FALSE
  )

  label_data <- data.frame(
    x = c(1.70, 5.50, 5.50, 9.60, 13.75),
    y = c(6.75, 6.82, 3.53, 6.75, 6.82),
    label = panel_titles,
    color = unname(colors[c("ink", "blue", "orange", "gold", "olive")]),
    stringsAsFactors = FALSE
  )

  body_data <- data.frame(
    x = c(1.70, 5.50, 5.50, 9.60, 13.75),
    y = c(2.75, 5.25, 1.92, 3.72, 3.82),
    label = panel_bodies,
    stringsAsFactors = FALSE
  )

  question_data <- data.frame(
    x = 1.70,
    y = 1.82,
    label = "WHEN?  WHAT?\nHOW ROBUST?",
    stringsAsFactors = FALSE
  )

  arrows <- data.frame(
    x = c(3.05, 3.05, 7.35, 7.35, 11.25),
    y = c(5.70, 2.75, 5.70, 2.75, 4.45),
    xend = c(3.62, 3.62, 7.92, 7.92, 11.82),
    yend = c(5.70, 2.75, 5.70, 2.75, 4.45),
    stringsAsFactors = FALSE
  )

  network_nodes <- data.frame(
    x = c(0.80, 1.18, 1.62, 2.02, 2.45, 2.78),
    y = c(5.02, 5.50, 4.95, 5.62, 5.02, 5.50),
    size = c(3.5, 4.8, 3.8, 5.1, 3.9, 3.5),
    stringsAsFactors = FALSE
  )
  network_edges <- data.frame(
    x = c(0.80, 0.80, 1.18, 1.18, 1.62, 2.02, 2.02, 2.45),
    y = c(5.02, 5.02, 5.50, 5.50, 4.95, 5.62, 5.62, 5.02),
    xend = c(1.18, 1.62, 1.62, 2.02, 2.02, 2.45, 2.78, 2.78),
    yend = c(5.50, 4.95, 4.95, 5.62, 5.62, 5.02, 5.50, 5.50),
    stringsAsFactors = FALSE
  )

  indicator_series <- data.frame(
    x = rep(seq(4.20, 6.80, length.out = 6), 3),
    y = c(
      5.80, 5.92, 5.86, 6.12, 6.02, 6.18,
      5.72, 5.78, 5.90, 5.83, 6.02, 5.98,
      5.66, 5.75, 5.72, 5.86, 5.82, 5.92
    ),
    series = factor(rep(1:3, each = 6)),
    stringsAsFactors = FALSE
  )
  indicator_colors <- c(
    "1" = unname(colors["blue"]),
    "2" = unname(colors["olive"]),
    "3" = unname(colors["gold"])
  )
  indicator_series$color <- unname(indicator_colors[indicator_series$series])
  anomaly_point <- data.frame(x = 5.76, y = 6.12)

  graphon_surface <- data.frame(
    group = rep(1:5, each = 4),
    x = c(
      4.45, 4.82, 5.10, 4.73,
      4.82, 5.19, 5.47, 5.10,
      5.19, 5.56, 5.84, 5.47,
      5.56, 5.93, 6.21, 5.84,
      5.93, 6.30, 6.58, 6.21
    ),
    y = c(
      2.84, 2.66, 2.79, 3.01,
      2.66, 2.54, 2.86, 2.79,
      2.54, 2.70, 3.12, 2.86,
      2.70, 2.96, 3.22, 3.12,
      2.96, 3.05, 3.10, 3.22
    ),
    z = rep(c(1, 2, 4, 3, 2), each = 4),
    stringsAsFactors = FALSE
  )
  surface_colors <- grDevices::colorRampPalette(c(
    alpha_fill(colors["orange"], 0.20),
    alpha_fill(colors["orange"], 0.75)
  ))(5)
  graphon_surface$fill <- surface_colors[graphon_surface$z]

  calibration_rings <- data.frame(
    size = c(9.0, 6.0, 3.0),
    color = c(
      alpha_fill(colors["gold"], 0.45),
      alpha_fill(colors["gold"], 0.72),
      unname(colors["gold"])
    ),
    stringsAsFactors = FALSE
  )

  interpretation_lines <- data.frame(
    x = c(12.75, 12.75, 13.55, 13.55),
    y = c(5.80, 5.42, 5.61, 5.61),
    xend = c(13.55, 13.55, 14.48, 14.48),
    yend = c(5.61, 5.61, 5.88, 5.34),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = box_data,
      ggplot2::aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill, color = border
      ),
      linewidth = 0.70, show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = unname(colors["mid_grey"]), linewidth = 0.75,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.075, "in"))
    ) +
    ggplot2::geom_segment(
      data = network_edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = alpha_fill(colors["blue"], 0.65), linewidth = 0.60,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.040, "in"))
    ) +
    ggplot2::geom_point(
      data = network_nodes,
      ggplot2::aes(x = x, y = y, size = size),
      shape = 21, fill = "white", color = unname(colors["blue"]),
      stroke = 0.70, show.legend = FALSE
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::geom_segment(
      ggplot2::aes(x = 4.15, y = 5.62, xend = 4.15, yend = 6.28),
      color = unname(colors["mid_grey"]), linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 4.15, y = 5.62, xend = 6.85, yend = 5.62),
      color = unname(colors["mid_grey"]), linewidth = 0.35
    ) +
    ggplot2::geom_line(
      data = indicator_series,
      ggplot2::aes(x = x, y = y, group = series, color = color),
      linewidth = 0.55, show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = anomaly_point,
      ggplot2::aes(x = x, y = 5.62, xend = x, yend = 6.17),
      color = "#C63C3C", linewidth = 0.45, linetype = 2,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = anomaly_point,
      ggplot2::aes(x = x, y = y),
      shape = 21, size = 2.2, stroke = 0.75,
      fill = "white", color = "#C63C3C"
    ) +
    ggplot2::geom_polygon(
      data = graphon_surface,
      ggplot2::aes(x = x, y = y, group = group, fill = fill),
      color = unname(colors["orange"]), linewidth = 0.30,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = calibration_rings,
      ggplot2::aes(x = 9.60, y = 5.55, size = size, color = color),
      shape = 1, stroke = 0.85, show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 9.40, y = 5.54, xend = 9.55, yend = 5.39),
      color = unname(colors["olive"]), linewidth = 0.75
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 9.55, y = 5.39, xend = 9.86, yend = 5.78),
      color = unname(colors["olive"]), linewidth = 0.75
    ) +
    ggplot2::geom_segment(
      data = interpretation_lines,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = unname(colors["olive"]), linewidth = 0.65,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.045, "in"))
    ) +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = x, y = y, label = label, color = color),
      fontface = "bold", size = 2.85, lineheight = 0.92,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = body_data,
      ggplot2::aes(x = x, y = y, label = label),
      color = unname(colors["ink"]), size = 2.50, lineheight = 1.04
    ) +
    ggplot2::geom_text(
      data = question_data,
      ggplot2::aes(x = x, y = y, label = label),
      color = unname(colors["ink"]), fontface = "bold",
      size = 2.22, lineheight = 0.98
    ) +
    ggplot2::annotate(
      "text", x = 8.0, y = 8.55,
      label = title_text,
      color = unname(colors["ink"]), fontface = "bold", size = 3.65
    ) +
    ggplot2::annotate(
      "text", x = 8.0, y = 8.05,
      label = subtitle_text,
      color = unname(colors["mid_grey"]), size = 2.55
    ) +
    ggplot2::annotate(
      "rect", xmin = 0.40, xmax = 15.60, ymin = 0.25, ymax = 1.05,
      fill = unname(colors["ink"]), color = NA
    ) +
    ggplot2::annotate(
      "text", x = 8.0, y = 0.65,
      label = footer_text,
      color = "white", fontface = "bold", size = 2.75
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 16), ylim = c(0, 9), expand = FALSE, clip = "on"
    ) +
    ggplot2::theme_void(base_family = "sans") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(4, 4, 4, 4),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}
