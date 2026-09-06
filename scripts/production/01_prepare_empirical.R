# Prepare the node-aligned annual empirical networks.

mx_log("Stage 01: preparing empirical networks.")
mx_require_packages(c("igraph"))

nodes_path <- file.path(MX_PATHS$raw_empirical, "nodes.txt")
mx_assert_files(nodes_path)

edge_files <- list.files(
  MX_PATHS$raw_empirical,
  pattern = "^[0-9]{4}_from_to_kist_weight[.]txt$",
  full.names = TRUE
)
edge_years <- suppressWarnings(as.integer(substr(basename(edge_files), 1L, 4L)))
edge_files <- edge_files[order(edge_years)]
edge_years <- edge_years[order(edge_years)]

if (!identical(edge_years, MX_CONFIG$years)) {
  stop(
    "Empirical edge files must cover 2006-2024 exactly. Found: ",
    paste(edge_years, collapse = ", "), call. = FALSE
  )
}

nodes <- utils::read.table(
  nodes_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_node_columns <- c("nodeID", "nodeLabel", "nodeLat", "nodeLong")
if (!all(required_node_columns %in% names(nodes))) {
  stop(
    "nodes.txt must contain: ", paste(required_node_columns, collapse = ", "),
    call. = FALSE
  )
}
nodes$nodeID <- as.character(nodes$nodeID)
if (anyDuplicated(nodes$nodeID)) stop("nodes.txt contains duplicate nodeID values.")
if (nrow(nodes) != MX_CONFIG$expected_nodes) {
  stop(
    "Expected ", MX_CONFIG$expected_nodes, " nodes; found ", nrow(nodes), ".",
    call. = FALSE
  )
}

read_annual_network <- function(path) {
  edges <- utils::read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_edge_columns <- c("from", "to", "weight")
  if (!all(required_edge_columns %in% names(edges))) {
    stop(basename(path), " must contain from, to, and weight.", call. = FALSE)
  }
  edges <- edges[, required_edge_columns, drop = FALSE]
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$weight <- suppressWarnings(as.numeric(edges$weight))
  if (any(!is.finite(edges$weight)) || any(edges$weight < 0)) {
    stop(basename(path), " contains invalid weights.", call. = FALSE)
  }
  unknown <- setdiff(unique(c(edges$from, edges$to)), nodes$nodeID)
  if (length(unknown)) {
    stop(
      basename(path), " contains unknown node IDs: ",
      paste(head(unknown, 10L), collapse = ", "), call. = FALSE
    )
  }
  if (anyDuplicated(edges[c("from", "to")])) {
    edges <- stats::aggregate(weight ~ from + to, data = edges, FUN = sum)
  }
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  if (!identical(as.character(igraph::V(graph)$name), nodes$nodeID)) {
    stop("Vertex order changed while reading ", basename(path), ".", call. = FALSE)
  }
  graph
}

networks <- lapply(edge_files, read_annual_network)
names(networks) <- as.character(edge_years)

annual_checks <- do.call(rbind, lapply(names(networks), function(year) {
  graph <- networks[[year]]
  weights <- as.numeric(igraph::E(graph)$weight)
  data.frame(
    year = as.integer(year),
    nodes = igraph::vcount(graph),
    directed_edges = igraph::ecount(graph),
    total_weight = sum(weights),
    finite_weights = all(is.finite(weights)),
    nonnegative_weights = all(weights >= 0),
    stringsAsFactors = FALSE
  )
}))

prepared <- list(
  nodes = nodes,
  networks = networks,
  years = edge_years,
  source_files = edge_files,
  annual_checks = annual_checks
)

prepared_path <- file.path(
  MX_PATHS$derived_empirical, "EMPIRICAL_node_aligned_networks.rds"
)
saveRDS(prepared, prepared_path, compress = "gzip")
mx_write_csv(
  annual_checks,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_annual_input_checks.csv")
)

input_paths <- c(nodes_path, edge_files)
input_manifest <- data.frame(
  file = gsub("\\\\", "/", input_paths),
  bytes = unname(file.info(input_paths)$size),
  md5 = unname(tools::md5sum(input_paths)),
  stringsAsFactors = FALSE
)
mx_write_csv(
  input_manifest,
  file.path(MX_PATHS$diagnostics, "EMPIRICAL_input_manifest.csv")
)

mx_log("Prepared ", length(networks), " annual networks at ", prepared_path, ".")
