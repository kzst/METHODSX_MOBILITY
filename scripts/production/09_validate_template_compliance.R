# Validate the rendered manuscript against the mandatory MethodsX method-article
# form. This is a structural pre-submission gate; it complements, but does not
# replace, visual inspection in Microsoft Word.

mx_word_count <- function(text) {
  cleaned <- trimws(gsub("[^[:alnum:]'’-]+", " ", text, perl = TRUE))
  if (!nzchar(cleaned)) return(0L)
  length(strsplit(cleaned, "[[:space:]]+", perl = TRUE)[[1L]])
}

mx_validate_methodsx_template <- function(
  docx_path = file.path("output", "MethodsX_MOB.docx"),
  report_path = file.path(
    "output", "diagnostics", "METHODSX_TEMPLATE_compliance.csv"
  )
) {
  if (!exists("mx_extract_docx", mode = "function")) {
    source(file.path(
      "scripts", "production", "07_freeze_docx_fields.R"
    ), local = parent.frame())
  }
  mx_require_docx_packages()

  audit_dir <- tempfile("methodsx_template_audit_")
  dir.create(audit_dir, recursive = TRUE)
  on.exit(unlink(audit_dir, recursive = TRUE, force = TRUE), add = TRUE)
  mx_extract_docx(docx_path, audit_dir)

  document_path <- file.path(audit_dir, "word", "document.xml")
  document_doc <- xml2::read_xml(document_path)
  namespaces <- xml2::xml_ns(document_doc)
  paragraph_nodes <- xml2::xml_find_all(
    document_doc, ".//w:body/w:p", namespaces
  )
  paragraph_text <- trimws(gsub(
    "[[:space:]]+", " ",
    vapply(paragraph_nodes, xml2::xml_text, character(1L)),
    perl = TRUE
  ))

  required_headings <- c(
    "Article information",
    "Abstract",
    "Graphical abstract",
    "Specifications table",
    "Background",
    "Method details",
    "Method validation",
    "Limitations",
    "Ethics statements",
    "CRediT author statement",
    "Acknowledgments",
    "Declaration of interests",
    "Supplementary material and/or additional information [OPTIONAL]",
    "References"
  )
  heading_positions <- match(required_headings, paragraph_text)
  headings_present <- all(!is.na(heading_positions))
  headings_in_order <- headings_present && all(diff(heading_positions) > 0L)

  expected_title <- paste(
    "Comparing Network Indicators and Graphon Distances for Structural",
    "Anomaly Detection: A Reproducible Workflow for Directed Weighted",
    "Longitudinal Networks"
  )
  title_position <- match("Article title", paragraph_text)
  author_position <- match("Authors", paragraph_text)
  title_candidates <- if (
    !is.na(title_position) && !is.na(author_position) &&
      author_position > title_position + 1L
  ) {
    paragraph_text[(title_position + 1L):(author_position - 1L)]
  } else {
    character()
  }
  title_candidates <- title_candidates[nzchar(title_candidates)]
  article_title <- paste(title_candidates, collapse = " ")
  title_words <- mx_word_count(article_title)

  abstract_position <- match("Abstract", paragraph_text)
  graphical_position <- match("Graphical abstract", paragraph_text)
  abstract_indices <- if (
    !is.na(abstract_position) && !is.na(graphical_position) &&
      graphical_position > abstract_position + 1L
  ) {
    (abstract_position + 1L):(graphical_position - 1L)
  } else {
    integer()
  }
  abstract_text <- paste(paragraph_text[abstract_indices], collapse = " ")
  abstract_words <- mx_word_count(abstract_text)
  abstract_bullets <- if (length(abstract_indices)) {
    sum(vapply(
      paragraph_nodes[abstract_indices],
      function(node) length(xml2::xml_find_all(
        node, ".//w:numPr", namespaces
      )) > 0L || grepl("^[[:space:]]*•", xml2::xml_text(node)),
      logical(1L)
    ))
  } else {
    0L
  }

  specification_fields <- c(
    "Subject area", "More specific subject area", "Name of your method",
    "Name and reference of original method", "Resource availability"
  )
  table_cell_text <- trimws(gsub(
    "[[:space:]]+", " ",
    vapply(
      xml2::xml_find_all(document_doc, ".//w:tbl//w:tc", namespaces),
      xml2::xml_text, character(1L)
    ),
    perl = TRUE
  ))

  graphical_indices <- if (
    !is.na(graphical_position) && !is.na(match("Specifications table", paragraph_text))
  ) {
    start <- graphical_position + 1L
    end <- match("Specifications table", paragraph_text) - 1L
    if (end >= start) start:end else integer()
  } else {
    integer()
  }
  graphical_drawings <- if (length(graphical_indices)) {
    sum(vapply(
      paragraph_nodes[graphical_indices],
      function(node) length(xml2::xml_find_all(
        node, ".//w:drawing | .//w:pict", namespaces
      )),
      integer(1L)
    ))
  } else {
    0L
  }

  full_text <- paste(paragraph_text, collapse = "\n")
  forbidden_instruction_fragments <- c(
    "AUTHOR INSTRUCTIONS",
    "Max. 20 words. The title should be unique.",
    "Please Select Subject Area from dropdown list",
    "Please delete this line and everything above it",
    "Reminder: Before you submit, please delete all the instructional text",
    "Our goal is to publish methods that can be replicated"
  )
  forbidden_present <- vapply(
    forbidden_instruction_fragments,
    grepl, logical(1L), x = full_text, fixed = TRUE
  )
  comment_anchors <- length(xml2::xml_find_all(
    document_doc,
    ".//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference",
    namespaces
  ))

  checks <- data.frame(
    check = c(
      "first_visible_section_is_article_information",
      "required_template_headings_present",
      "required_template_headings_in_order",
      "article_title_exact",
      "article_title_max_20_words",
      "no_duplicate_front_matter_title",
      "both_author_names_present",
      "corresponding_author_marked",
      "both_orcids_present",
      "both_institutional_emails_present",
      "abstract_max_200_words",
      "abstract_has_1_to_3_bullets",
      "graphical_abstract_declared_as_separate_file",
      "graphical_abstract_not_embedded_in_manuscript",
      "all_specifications_fields_present",
      "both_competing_interest_options_retained",
      "template_instruction_text_removed",
      "template_comment_anchors_removed",
      "obsolete_conclusion_heading_removed"
    ),
    observed = c(
      paragraph_text[which(nzchar(paragraph_text))[1L]],
      paste(sum(!is.na(heading_positions)), "of", length(required_headings)),
      paste(heading_positions, collapse = ","),
      article_title,
      as.character(title_words),
      as.character(sum(paragraph_text == expected_title)),
      paste(vapply(
        c("Zsolt Tibor Kosztyán", "Kornél Dénes"),
        grepl, logical(1L), x = full_text, fixed = TRUE
      ), collapse = ","),
      as.character(grepl("Kosztyán.*\\*", full_text, perl = TRUE)),
      paste(vapply(
        c("0000-0001-7345-8336", "0009-0002-9527-6124"),
        grepl, logical(1L), x = full_text, fixed = TRUE
      ), collapse = ","),
      paste(vapply(
        c("kosztyan.zsolt@gtk.uni-pannon.hu", "morgosz@student.elte.hu"),
        grepl, logical(1L), x = full_text, fixed = TRUE
      ), collapse = ","),
      as.character(abstract_words),
      as.character(abstract_bullets),
      as.character(grepl("Submitted as separate files", full_text, fixed = TRUE)),
      as.character(graphical_drawings),
      paste(sum(specification_fields %in% table_cell_text), "of", length(specification_fields)),
      paste(vapply(
        c(
          "The authors declare that they have no known competing financial interests",
          "The authors declare the following financial interests/personal relationships"
        ), grepl, logical(1L), x = full_text, fixed = TRUE
      ), collapse = ","),
      paste(forbidden_instruction_fragments[forbidden_present], collapse = "; "),
      as.character(comment_anchors),
      as.character(grepl(
        "Practical interpretation and conclusions", full_text, fixed = TRUE
      ))
    ),
    passed = c(
      identical(paragraph_text[which(nzchar(paragraph_text))[1L]], "Article information"),
      headings_present,
      headings_in_order,
      identical(article_title, expected_title),
      title_words <= 20L,
      sum(paragraph_text == expected_title) == 1L,
      all(vapply(
        c("Zsolt Tibor Kosztyán", "Kornél Dénes"),
        grepl, logical(1L), x = full_text, fixed = TRUE
      )),
      grepl("Kosztyán.*\\*", full_text, perl = TRUE),
      all(vapply(
        c("0000-0001-7345-8336", "0009-0002-9527-6124"),
        grepl, logical(1L), x = full_text, fixed = TRUE
      )),
      all(vapply(
        c("kosztyan.zsolt@gtk.uni-pannon.hu", "morgosz@student.elte.hu"),
        grepl, logical(1L), x = full_text, fixed = TRUE
      )),
      abstract_words <= 200L,
      abstract_bullets >= 1L && abstract_bullets <= 3L,
      grepl("Submitted as separate files", full_text, fixed = TRUE),
      graphical_drawings == 0L,
      all(specification_fields %in% table_cell_text),
      all(vapply(c(
        "The authors declare that they have no known competing financial interests",
        "The authors declare the following financial interests/personal relationships"
      ), grepl, logical(1L), x = full_text, fixed = TRUE)),
      !any(forbidden_present),
      comment_anchors == 0L,
      !grepl("Practical interpretation and conclusions", full_text, fixed = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(checks, report_path, row.names = FALSE, na = "")
  if (!all(checks$passed)) {
    stop(
      "MethodsX template-compliance checks failed: ",
      paste(checks$check[!checks$passed], collapse = ", "),
      ". See ", report_path, ".",
      call. = FALSE
    )
  }

  invisible(list(
    title_words = title_words,
    abstract_words = abstract_words,
    abstract_bullets = abstract_bullets,
    report = report_path,
    checks = checks
  ))
}
