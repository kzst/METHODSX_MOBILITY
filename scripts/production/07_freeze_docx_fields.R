# Prepare and verify Word OOXML packages used by the manuscript pipeline.
#
# The MethodsX reference document contains a valid header, but its relationship
# identifier is named `rIdMethodsXHeader`. officer/officedown expects numeric
# `rId<number>` identifiers and consequently emits repeated coercion warnings.
# The source template is preserved: mx_prepare_reference_docx() makes a clean
# runtime copy with numeric relationship identifiers before rendering.
#
# officedown also writes figure numbers as dirty SEQ fields. Microsoft Word can
# then show a generic security prompt claiming that fields may refer to other
# files, even though all figures are embedded and the fields only perform local
# caption numbering. mx_freeze_docx_fields() replaces those SEQ fields with
# their visible integer values and verifies the final archive.

mx_count_regex <- function(pattern, text) {
  matches <- gregexpr(pattern, text, perl = TRUE)[[1L]]
  if (length(matches) == 1L && matches[1L] < 0L) 0L else length(matches)
}

mx_require_docx_packages <- function() {
  missing <- c("xml2", "zip")[!vapply(
    c("xml2", "zip"), requireNamespace, logical(1L), quietly = TRUE
  )]
  if (length(missing)) {
    stop(
      "Required DOCX package(s) not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

mx_extract_docx <- function(docx_path, exdir) {
  mx_require_docx_packages()
  if (!file.exists(docx_path)) {
    stop("DOCX file not found: ", docx_path, call. = FALSE)
  }

  archive_path <- normalizePath(
    docx_path, winslash = "/", mustWork = TRUE
  )
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  extraction_error <- NULL
  tryCatch(
    zip::unzip(archive_path, exdir = exdir),
    error = function(error) extraction_error <<- conditionMessage(error)
  )

  document_path <- file.path(exdir, "word", "document.xml")
  if (!is.null(extraction_error) || !file.exists(document_path)) {
    members <- tryCatch(
      zip::zip_list(archive_path)$filename,
      error = function(error) character()
    )
    archive_has_document <- "word/document.xml" %in% members
    detail <- if (is.null(extraction_error)) {
      "the extractor returned without creating word/document.xml"
    } else {
      extraction_error
    }
    stop(
      "Could not extract the DOCX archive at its absolute path: ",
      archive_path, ". Archive index contains word/document.xml: ",
      archive_has_document, ". Detail: ", detail, ".",
      call. = FALSE
    )
  }

  invisible(archive_path)
}

mx_pack_docx <- function(work_dir, output_path) {
  mx_require_docx_packages()
  work_dir <- normalizePath(work_dir, winslash = "/", mustWork = TRUE)
  output_parent <- dirname(output_path)
  dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  output_parent <- normalizePath(
    output_parent, winslash = "/", mustWork = TRUE
  )
  output_path <- file.path(output_parent, basename(output_path))

  archive_members <- list.files(
    work_dir, recursive = TRUE, all.files = TRUE,
    no.. = TRUE, include.dirs = FALSE
  )
  packed_docx <- tempfile(
    pattern = "methodsx_docx_pack_", tmpdir = output_parent,
    fileext = ".docx"
  )
  on.exit(unlink(packed_docx, force = TRUE), add = TRUE)
  zip::zipr(
    zipfile = packed_docx,
    files = archive_members,
    root = work_dir,
    mode = "mirror",
    include_directories = FALSE,
    compression_level = 9
  )

  packed_members <- zip::zip_list(packed_docx)$filename
  if (!"word/document.xml" %in% packed_members) {
    stop("Repacked DOCX has no word/document.xml part.", call. = FALSE)
  }
  if (!file.copy(packed_docx, output_path, overwrite = TRUE)) {
    stop("Could not write the DOCX archive: ", output_path, call. = FALSE)
  }
  invisible(output_path)
}

mx_relationship_audit <- function(extracted_dir) {
  relationship_files <- list.files(
    extracted_dir, pattern = "[.]rels$", recursive = TRUE, full.names = TRUE
  )
  external_file_targets <- character()
  external_hyperlinks <- character()
  invalid_relationship_ids <- character()

  for (relationship_file in relationship_files) {
    relationship_doc <- xml2::read_xml(relationship_file)
    relationship_nodes <- xml2::xml_find_all(
      relationship_doc, ".//*[local-name()='Relationship']"
    )
    if (!length(relationship_nodes)) next

    ids <- xml2::xml_attr(relationship_nodes, "Id")
    invalid <- is.na(ids) | !grepl("^rId[0-9]+$", ids)
    if (any(invalid)) {
      relative_path <- substring(
        normalizePath(relationship_file, winslash = "/", mustWork = TRUE),
        nchar(normalizePath(extracted_dir, winslash = "/", mustWork = TRUE)) + 2L
      )
      invalid_relationship_ids <- c(
        invalid_relationship_ids,
        paste0(relative_path, "::", ids[invalid])
      )
    }

    is_external <- xml2::xml_attr(relationship_nodes, "TargetMode") == "External"
    is_external[is.na(is_external)] <- FALSE
    if (!any(is_external)) next
    targets <- xml2::xml_attr(relationship_nodes[is_external], "Target")
    types <- xml2::xml_attr(relationship_nodes[is_external], "Type")
    is_hyperlink <- grepl("/hyperlink$", types)
    is_file_target <- grepl(
      "^(?:file:|[A-Za-z]:[/\\\\]|\\\\\\\\|/)", targets,
      ignore.case = TRUE, perl = TRUE
    )
    external_file_targets <- c(
      external_file_targets, targets[is_file_target | !is_hyperlink]
    )
    external_hyperlinks <- c(
      external_hyperlinks, targets[is_hyperlink & !is_file_target]
    )
  }

  list(
    invalid_relationship_ids = unique(invalid_relationship_ids),
    external_file_targets = unique(external_file_targets),
    external_hyperlinks = unique(external_hyperlinks)
  )
}

mx_audit_docx_fields <- function(docx_path) {
  mx_require_docx_packages()
  audit_dir <- tempfile("methodsx_docx_audit_")
  dir.create(audit_dir, recursive = TRUE)
  on.exit(unlink(audit_dir, recursive = TRUE, force = TRUE), add = TRUE)
  mx_extract_docx(docx_path, audit_dir)

  word_xml_paths <- c(
    file.path(audit_dir, "word", "document.xml"),
    list.files(
      file.path(audit_dir, "word"),
      pattern = "^(header|footer)[0-9]+[.]xml$", full.names = TRUE
    )
  )
  word_xml <- paste(vapply(
    word_xml_paths,
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  settings_path <- file.path(audit_dir, "word", "settings.xml")
  settings_xml <- if (file.exists(settings_path)) {
    paste(readLines(settings_path, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  relationship_audit <- mx_relationship_audit(audit_dir)

  c(
    list(
      seq_fields = mx_count_regex("<w:instrText[^>]*>\\s*SEQ fig", word_xml),
      field_markers = mx_count_regex("<w:fldChar\\b", word_xml),
      dirty_fields = mx_count_regex("w:dirty=", word_xml),
      update_fields = mx_count_regex("<w:updateFields\\b", settings_xml)
    ),
    relationship_audit
  )
}

mx_prepare_reference_docx <- function(reference_path, output_path) {
  mx_require_docx_packages()
  source_path <- normalizePath(
    reference_path, winslash = "/", mustWork = TRUE
  )
  output_parent <- dirname(output_path)
  dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  target_path <- file.path(
    normalizePath(output_parent, winslash = "/", mustWork = TRUE),
    basename(output_path)
  )
  if (identical(source_path, target_path)) {
    stop(
      "The cleaned reference DOCX must not overwrite the source template.",
      call. = FALSE
    )
  }

  work_dir <- tempfile("methodsx_reference_clean_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  mx_extract_docx(source_path, work_dir)

  relationships_path <- file.path(
    work_dir, "word", "_rels", "document.xml.rels"
  )
  document_path <- file.path(work_dir, "word", "document.xml")
  relationship_doc <- xml2::read_xml(relationships_path)
  relationship_nodes <- xml2::xml_find_all(
    relationship_doc, ".//*[local-name()='Relationship']"
  )
  ids <- xml2::xml_attr(relationship_nodes, "Id")
  types <- xml2::xml_attr(relationship_nodes, "Type")
  numeric_ids <- suppressWarnings(as.integer(sub("^rId", "", ids)))
  next_id <- if (all(is.na(numeric_ids))) 1L else max(numeric_ids, na.rm = TRUE) + 1L
  repair <- grepl("/header$", types) & !grepl("^rId[0-9]+$", ids)

  replacements <- character()
  if (any(repair)) {
    document_xml <- paste(readLines(document_path, warn = FALSE), collapse = "\n")
    for (node_index in which(repair)) {
      old_id <- ids[node_index]
      new_id <- paste0("rId", next_id)
      next_id <- next_id + 1L
      old_reference <- paste0('r:id="', old_id, '"')
      new_reference <- paste0('r:id="', new_id, '"')
      if (!grepl(old_reference, document_xml, fixed = TRUE)) {
        stop(
          "Header relationship ", old_id,
          " is not referenced from word/document.xml.",
          call. = FALSE
        )
      }
      xml2::xml_set_attr(relationship_nodes[[node_index]], "Id", new_id)
      document_xml <- gsub(
        old_reference, new_reference, document_xml, fixed = TRUE
      )
      replacements <- c(replacements, paste0(old_id, " -> ", new_id))
    }
    xml2::write_xml(relationship_doc, relationships_path, options = "format")
    writeLines(document_xml, document_path, useBytes = TRUE)
  }

  mx_pack_docx(work_dir, target_path)
  audit <- mx_audit_docx_fields(target_path)
  if (length(audit$invalid_relationship_ids)) {
    stop(
      "The cleaned Word reference still contains invalid relationship IDs: ",
      paste(audit$invalid_relationship_ids, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(replacements)
}

mx_freeze_docx_fields <- function(docx_path) {
  mx_require_docx_packages()
  work_dir <- tempfile("methodsx_docx_freeze_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  mx_extract_docx(docx_path, work_dir)

  document_path <- file.path(work_dir, "word", "document.xml")
  document_doc <- xml2::read_xml(document_path)
  namespaces <- xml2::xml_ns(document_doc)
  all_instructions <- xml2::xml_find_all(document_doc, ".//w:instrText", namespaces)
  sequence_instructions <- xml2::xml_find_all(
    document_doc,
    ".//w:instrText[normalize-space(.)='SEQ fig \\* Arabic']",
    namespaces
  )

  if (length(all_instructions) != length(sequence_instructions)) {
    stop(
      "The Word document contains field instructions other than figure SEQ fields; ",
      "they were not flattened automatically.",
      call. = FALSE
    )
  }

  if (length(sequence_instructions)) {
    for (index in seq_along(sequence_instructions)) {
      instruction <- sequence_instructions[[index]]
      instruction_run <- xml2::xml_parent(instruction)
      begin_run <- xml2::xml_find_first(
        instruction_run, "preceding-sibling::w:r[1]", namespaces
      )
      end_run <- xml2::xml_find_first(
        instruction_run, "following-sibling::w:r[1]", namespaces
      )
      begin_marker <- xml2::xml_find_first(
        begin_run, ".//w:fldChar[@w:fldCharType='begin']", namespaces
      )
      end_marker <- xml2::xml_find_first(
        end_run, ".//w:fldChar[@w:fldCharType='end']", namespaces
      )
      if (inherits(begin_marker, "xml_missing") || inherits(end_marker, "xml_missing")) {
        stop("Unexpected Word field structure around a figure number.", call. = FALSE)
      }

      xml2::xml_name(instruction) <- "t"
      xml2::xml_text(instruction) <- as.character(index)
      xml2::xml_remove(begin_run)
      xml2::xml_remove(end_run)
    }
  }

  xml2::write_xml(document_doc, document_path, options = "format")
  document_xml <- readLines(document_path, warn = FALSE)
  document_xml <- gsub(' w:dirty="true"', "", document_xml, fixed = TRUE)
  writeLines(document_xml, document_path, useBytes = TRUE)

  settings_path <- file.path(work_dir, "word", "settings.xml")
  if (file.exists(settings_path)) {
    settings_doc <- xml2::read_xml(settings_path)
    settings_ns <- xml2::xml_ns(settings_doc)
    update_nodes <- xml2::xml_find_all(settings_doc, ".//w:updateFields", settings_ns)
    if (length(update_nodes)) xml2::xml_remove(update_nodes)
    xml2::write_xml(settings_doc, settings_path, options = "format")
  }

  frozen_docx <- tempfile(fileext = ".docx")
  on.exit(unlink(frozen_docx, force = TRUE), add = TRUE)
  mx_pack_docx(work_dir, frozen_docx)

  audit <- mx_audit_docx_fields(frozen_docx)
  if (
    audit$seq_fields != 0L || audit$field_markers != 0L ||
      audit$dirty_fields != 0L || audit$update_fields != 0L ||
      length(audit$external_file_targets) != 0L ||
      length(audit$invalid_relationship_ids) != 0L
  ) {
    stop("The field-free DOCX verification did not pass.", call. = FALSE)
  }

  if (!file.copy(frozen_docx, docx_path, overwrite = TRUE)) {
    stop("Could not replace the rendered DOCX with its field-free copy.", call. = FALSE)
  }
  invisible(length(sequence_instructions))
}
