geoplot <- function(g,
                    mode = c("basemap", "rworldmap", "custom_shape"),
                    basemap = "OpenStreetMap",
                    shape = NULL,
                    show_tmap_on_custom = FALSE,
                    node_col = "black",
                    node_alpha = 0.8,
                    node_size = 0.8,
                    node_size_scale = c(0.5, 2.5),
                    colormap = grDevices::rainbow,
                    edge_col = "gray",
                    edge_alpha = 0.5,
                    edge_width = 1,
                    edge_curved = FALSE,
                    curvature = 0.15,
                    arrow_size = 0.0,
                    shape_bg = "#f0f0f0",
                    shape_border_col = "#d0d0d0",
                    # [NEW] Label parameters
                    show_labels = FALSE,
                    label_field = "nodeLabel",
                    label_col = "black",
                    label_size_mult = 1,
                    top_n_labels = NULL) {

  library(igraph); library(sf); library(tmap); library(scales)
  #tmap_options(facet.max = 2000)
  mode <- match.arg(mode)

  v_count <- igraph::vcount(g)
  e_count <- igraph::ecount(g)
  nodes_df <- igraph::as_data_frame(g, what = "vertices")

  # [NEW] Detect whether the graph is directed (controls arrowhead rendering)
  is_dir <- igraph::is_directed(g)

  # Detect pre-existing edge attributes on the igraph object
  has_edge_color_attr <- e_count > 0 && !is.null(igraph::edge_attr(g, "color"))
  has_edge_width_attr <- e_count > 0 && !is.null(igraph::edge_attr(g, "width"))

  # --- 1. NODE DATA (color: handles functions and fixed palettes) ---
  get_colors <- function(n_levels) {
    if (is.function(colormap)) {
      return(colormap(n_levels))
    } else if (is.character(colormap) && length(colormap) >= n_levels) {
      return(colormap[1:n_levels])
    } else {
      return(grDevices::rainbow(n_levels))
    }
  }

  if (length(node_col) == 1 && node_col %in% names(nodes_df)) {
    vals <- as.factor(nodes_df[[node_col]])
    node_colors_raw <- get_colors(length(levels(vals)))[vals]
  } else if (length(node_col) == v_count && !is.character(node_col)) {
    vals <- as.factor(node_col)
    node_colors_raw <- get_colors(length(levels(vals)))[vals]
  } else if (length(node_col) == v_count) {
    node_colors_raw <- node_col
  } else {
    node_colors_raw <- rep(node_col[1], v_count)
  }
  final_node_colors <- scales::alpha(node_colors_raw, node_alpha)

  # --- NODE SIZE ---
  rescale_vec <- function(x, range_to) {
    if (length(x) == 0) return(numeric(0))
    if (all(x == x[1])) return(rep(range_to[1], length(x)))
    x_num <- as.numeric(x)
    (x_num - min(x_num, na.rm = TRUE)) /
      (max(x_num, na.rm = TRUE) - min(x_num, na.rm = TRUE)) *
      (range_to[2] - range_to[1]) + range_to[1]
  }

  # [CHANGED] Guard: coerce list to atomic vector (igraph may store attrs as lists)
  if (is.list(node_size) && !is.data.frame(node_size)) {
    node_size <- unlist(node_size)
  }

  if (is.character(node_size) && length(node_size) == 1 && node_size %in% names(nodes_df)) {
    col_data <- nodes_df[[node_size]]
    if (is.list(col_data)) col_data <- unlist(col_data)
    final_node_size <- rescale_vec(as.numeric(col_data), node_size_scale)
  } else if (length(node_size) == v_count) {
    final_node_size <- rescale_vec(as.numeric(node_size), node_size_scale)
  } else {
    size_val <- suppressWarnings(as.numeric(node_size[1]))
    if (is.na(size_val)) {
      warning("node_size could not be interpreted as numeric. Falling back to 0.8.")
      size_val <- 0.8
    }
    final_node_size <- rep(size_val, v_count)
  }

  nodes_sf <- sf::st_as_sf(nodes_df, coords = c("lon", "lat"), crs = 4326)
  nodes_sf$node_color_val <- final_node_colors
  nodes_sf$node_size_val  <- final_node_size * 0.05

  # --- [NEW] NODE LABELS ---
  # Precompute label text, font sizes, and the display mask so that both
  # the tmap and the static rendering paths can use them.
  label_text <- NULL
  label_cex  <- NULL
  label_mask <- NULL

  if (show_labels) {
    # Resolve the label text source
    if (label_field %in% names(nodes_df)) {
      label_text <- as.character(nodes_df[[label_field]])
    } else if ("name" %in% names(nodes_df)) {
      warning(paste0("Label field '", label_field, "' not found in vertex ",
                     "attributes. Falling back to 'name'."))
      label_text <- as.character(nodes_df[["name"]])
    } else {
      warning(paste0("Label field '", label_field, "' not found and no 'name' ",
                     "attribute available. Using node indices."))
      label_text <- as.character(seq_len(v_count))
    }

    # Font size proportional to node size, scaled by user multiplier
    label_cex <- final_node_size * label_size_mult

    # Optionally restrict to the top-N largest nodes
    if (!is.null(top_n_labels) && is.numeric(top_n_labels) &&
        top_n_labels > 0 && top_n_labels < v_count) {
      top_idx    <- order(final_node_size, decreasing = TRUE)[1:top_n_labels]
      label_mask <- seq_len(v_count) %in% top_idx
    } else {
      label_mask <- rep(TRUE, v_count)
    }

    # Store in nodes_sf for tmap modes
    nodes_sf$node_label    <- ifelse(label_mask, label_text, NA_character_)
    nodes_sf$label_size_val <- ifelse(label_mask, label_cex * 0.4, 0)
  }

  # --- 2. EDGE DATA ---
  edges_sf  <- NULL
  arrows_sf <- NULL          # [NEW] Arrowhead polygons for directed graphs

  if (e_count > 0) {
    edge_list <- igraph::as_edgelist(g, names = TRUE)

    # -- [CHANGED] Edge colors: prefer E(g)$color when available --
    if (has_edge_color_attr) {
      e_cols_raw <- igraph::E(g)$color
      if (is.list(e_cols_raw)) e_cols_raw <- unlist(e_cols_raw)
    } else if (edge_col == "mixed") {
      from_idx   <- match(edge_list[, 1], nodes_df$name)
      e_cols_raw <- node_colors_raw[from_idx]
    } else {
      e_cols_raw <- rep(edge_col[1], e_count)
    }
    final_edge_colors <- scales::alpha(e_cols_raw, edge_alpha)

    # -- [NEW] Edge widths: prefer E(g)$width when available --
    if (has_edge_width_attr) {
      final_edge_widths <- igraph::E(g)$width
      if (is.list(final_edge_widths)) final_edge_widths <- unlist(final_edge_widths)
      final_edge_widths <- as.numeric(final_edge_widths)
      final_edge_widths[is.na(final_edge_widths)] <- edge_width
    } else {
      final_edge_widths <- rep(edge_width, e_count)
    }

    # -- Edge geometries (straight or Bézier) --
    calculate_bezier <- function(p1, p2, v, n = 30) {
      mid <- (p1 + p2) / 2
      diff <- p2 - p1
      norm <- c(-diff[2], diff[1])
      control <- mid + norm * v
      t <- seq(0, 1, length.out = n)
      coords <- matrix(NA, n, 2)
      for (i in 1:n)
        coords[i, ] <- (1 - t[i])^2 * p1 +
        2 * (1 - t[i]) * t[i] * control +
        t[i]^2 * p2
      return(coords)
    }

    edge_geoms <- lapply(1:e_count, function(i) {
      p1 <- as.numeric(nodes_df[nodes_df$name == edge_list[i, 1], c("lon", "lat")])
      p2 <- as.numeric(nodes_df[nodes_df$name == edge_list[i, 2], c("lon", "lat")])
      sf::st_linestring(
        if (edge_curved) calculate_bezier(p1, p2, curvature) else rbind(p1, p2)
      )
    })

    # [CHANGED] edges_sf now carries per-edge width as well
    edges_sf <- sf::st_sf(
      geometry       = sf::st_sfc(edge_geoms, crs = 4326),
      edge_color_val = as.character(final_edge_colors),
      edge_width_val = final_edge_widths
    )

    # ------------------------------------------------------------------
    # [NEW] Build arrowhead polygons for directed graphs when arrow_size > 0
    # ------------------------------------------------------------------
    # Arrowheads are drawn as filled triangle polygons at the target end of
    # each edge.  The size is relative to the spatial extent of the node
    # bounding box so that arrowheads look proportional regardless of zoom.
    # arrow_size acts as a user-controlled multiplier (default 0 = no arrows).
    # ------------------------------------------------------------------
    if (is_dir && arrow_size > 0) {

      # Compute a base arrowhead length from the spatial extent of all nodes
      bbox <- sf::st_bbox(nodes_sf)
      extent_diag <- sqrt((bbox["xmax"] - bbox["xmin"])^2 +
                            (bbox["ymax"] - bbox["ymin"])^2)
      # Fallback for degenerate cases (all nodes at the same location)
      if (!is.finite(extent_diag) || extent_diag == 0) extent_diag <- 1

      # Base length = 1.5 % of the diagonal, scaled by the user's arrow_size
      base_arrow_len <- extent_diag * 0.015 * arrow_size

      # Helper: create a single arrowhead triangle polygon
      # tip       – numeric(2), the point of the arrowhead (target end)
      # direction – numeric(2), vector pointing from penultimate to tip
      # size      – scalar, length of the arrowhead along the direction axis
      make_arrowhead <- function(tip, direction, size) {
        d_len <- sqrt(sum(direction^2))
        if (d_len == 0) return(NULL)
        d    <- direction / d_len          # unit vector along edge
        perp <- c(-d[2], d[1])             # perpendicular unit vector
        base_mid <- tip - d * size         # midpoint of the triangle base
        p1 <- base_mid + perp * size * 0.45
        p2 <- base_mid - perp * size * 0.45
        sf::st_polygon(list(rbind(tip, p1, p2, tip)))
      }

      arrow_polys  <- vector("list", e_count)
      arrow_colors <- character(e_count)
      keep         <- logical(e_count)

      for (i in seq_len(e_count)) {
        coords <- sf::st_coordinates(edges_sf$geometry[i])
        n_pts  <- nrow(coords)
        # Tip = last point (target node), prev = penultimate point
        tip  <- coords[n_pts, 1:2]
        prev <- coords[max(n_pts - 1, 1), 1:2]
        direction <- tip - prev

        poly <- make_arrowhead(as.numeric(tip),
                               as.numeric(direction),
                               base_arrow_len)
        if (!is.null(poly)) {
          arrow_polys[[i]] <- poly
          arrow_colors[i]  <- edges_sf$edge_color_val[i]
          keep[i]          <- TRUE
        }
      }

      if (any(keep)) {
        arrows_sf <- sf::st_sf(
          geometry    = sf::st_sfc(arrow_polys[keep], crs = 4326),
          arrow_color = arrow_colors[keep]
        )
      }
    }
  }

  # --- 3. RENDERING ---
  if (mode == "basemap" || (mode == "custom_shape" && show_tmap_on_custom)) {
    # ---- Interactive (tmap) mode ----
    tmap_mode("view")

    final_bmap <- if (basemap == "osm") "OpenStreetMap"
    else if (basemap == "dark") "CartoDB.DarkMatter"
    else if (basemap == "light") "CartoDB.Positron"
    else basemap

    mapa <- if (mode == "basemap") {
      tm_basemap(final_bmap)
    } else {
      tm_shape(st_transform(shape, 4326)) +
        tm_polygons(col = shape_bg, border.col = shape_border_col)
    }

    # [CHANGED] Edges: variable width from edge_width_val column
    if (!is.null(edges_sf)) {
      mapa <- mapa +
        tm_shape(edges_sf) +
        tm_lines(col     = "edge_color_val",
                 lwd     = "edge_width_val",
                 col.legend = tm_legend_hide(),
                 lwd.legend = tm_legend_hide())
    }

    # [NEW] Arrowheads (tmap): rendered as filled triangle polygons on top
    # of the edge lines so that the direction of each edge is visible.
    if (!is.null(arrows_sf)) {
      mapa <- mapa +
        tm_shape(arrows_sf) +
        tm_polygons(col        = "arrow_color",
                    border.col = "arrow_color",
                    col.legend = tm_legend_hide())
    }

    # Nodes
    mapa <- mapa +
      tm_shape(nodes_sf) +
      tm_symbols(col    = "node_color_val",
                 size   = "node_size_val",
                 border.lwd = 0,
                 col.legend  = tm_legend_hide(),
                 size.legend = tm_legend_hide())

    # [NEW] Labels (tmap)
    if (show_labels) {
      # Subset to labeled nodes only (non-NA label)
      labeled_sf <- nodes_sf[!is.na(nodes_sf$node_label), ]
      if (nrow(labeled_sf) > 0) {
        mapa <- mapa +
          tm_shape(labeled_sf) +
          tm_text("node_label",
                  size = "label_size_val",
                  col  = label_col,
                  size.legend = tm_legend_hide())
      }
    }

    return(mapa)

  } else {
    # ---- Static (base R) mode ----
    if (mode == "rworldmap") {
      world     <- rworldmap::getMap(resolution = "high")
      shape_gps <- sf::st_as_sf(world)
      bbox      <- st_bbox(nodes_sf)
      lims      <- list(x = c(bbox["xmin"], bbox["xmax"]),
                        y = c(bbox["ymin"], bbox["ymax"]))
    } else {
      shape_gps <- st_transform(shape, 4326)
      lims      <- NULL
    }

    plot(st_geometry(shape_gps),
         col = shape_bg, border = shape_border_col,
         xlim = lims$x, ylim = lims$y)

    # [CHANGED] Edges: per-edge color AND per-edge width
    if (!is.null(edges_sf)) {
      plot(st_geometry(edges_sf), add = TRUE,
           col = edges_sf$edge_color_val,
           lwd = edges_sf$edge_width_val)
    }

    # [NEW] Arrowheads (base R): filled triangle polygons overlaid on edges
    if (!is.null(arrows_sf)) {
      plot(st_geometry(arrows_sf), add = TRUE,
           col    = arrows_sf$arrow_color,
           border = arrows_sf$arrow_color)
    }

    # Nodes
    plot(st_geometry(nodes_sf), add = TRUE,
         col = nodes_sf$node_color_val,
         pch = 16, cex = final_node_size)

    # [NEW] Labels (base R)
    if (show_labels) {
      idx_to_label <- which(label_mask)
      if (length(idx_to_label) > 0) {
        coords <- sf::st_coordinates(nodes_sf)
        text(coords[idx_to_label, 1],
             coords[idx_to_label, 2],
             labels = label_text[idx_to_label],
             cex    = label_cex[idx_to_label],
             col    = label_col,
             pos    = 3,
             offset = 0.5)
      }
    }

    return(recordPlot())
  }
}




