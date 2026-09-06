# Reproducible one-page cover letter for the MethodsX submission.

mx_build_cover_letter <- function(
  output_path = file.path("output", "Cover_Letter_MethodsX.docx")
) {
  mx_require_packages("officer")
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  submission_date <- Sys.getenv("METHODSX_SUBMISSION_DATE", unset = "")
  if (!nzchar(submission_date)) submission_date <- as.character(Sys.Date())
  submission_date <- as.Date(submission_date)
  if (is.na(submission_date)) {
    stop("METHODSX_SUBMISSION_DATE must use YYYY-MM-DD format.", call. = FALSE)
  }
  english_months <- c(
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
  )
  date_label <- sprintf(
    "%d %s %d",
    as.integer(format(submission_date, "%d")),
    english_months[as.integer(format(submission_date, "%m"))],
    as.integer(format(submission_date, "%Y"))
  )

  ink <- "#25313A"
  blue <- "#2E74B5"
  grey <- "#657078"
  body_text <- officer::fp_text(
    font.family = "Calibri", font.size = 11, color = ink
  )
  body_bold <- officer::fp_text(
    font.family = "Calibri", font.size = 11, bold = TRUE, color = ink
  )
  small_text <- officer::fp_text(
    font.family = "Calibri", font.size = 9.5, color = grey
  )
  name_text <- officer::fp_text(
    font.family = "Calibri", font.size = 15, bold = TRUE, color = blue
  )
  subject_text <- officer::fp_text(
    font.family = "Calibri", font.size = 11, bold = TRUE, color = blue
  )
  paragraph_format <- officer::fp_par(
    text.align = "left", padding.bottom = 6, line_spacing = 1.10
  )
  compact_format <- officer::fp_par(
    text.align = "left", padding.bottom = 1, line_spacing = 1.00
  )

  add_text_paragraph <- function(doc, text, style = body_text, fp_p = paragraph_format) {
    officer::body_add_fpar(
      doc,
      officer::fpar(officer::ftext(text, style), fp_p = fp_p)
    )
  }

  doc <- officer::read_docx()
  doc <- officer::body_add_fpar(
    doc,
    officer::fpar(
      officer::ftext("Zsolt T. Kosztyán", name_text),
      fp_p = compact_format
    )
  )
  doc <- add_text_paragraph(
    doc,
    "BRIDGE and Department of Quantitative Methods, University of Pannonia",
    small_text, compact_format
  )
  doc <- add_text_paragraph(
    doc,
    "Egyetem Street 10, 8200 Veszprém, Hungary | kosztyan.zsolt@gtk.uni-pannon.hu",
    small_text, compact_format
  )
  doc <- officer::body_add_fpar(
    doc,
    officer::fpar(
      officer::ftext(date_label, body_text),
      fp_p = officer::fp_par(text.align = "right", padding.bottom = 10)
    )
  )
  doc <- add_text_paragraph(doc, "Dear Editor,", body_text)
  doc <- officer::body_add_fpar(
    doc,
    officer::fpar(
      officer::ftext("Re: Submission to MethodsX — ", body_bold),
      officer::ftext(
        paste0(
          "Comparing Network Indicators and Graphon Distances for Structural ",
          "Anomaly Detection: A Reproducible Workflow for Directed Weighted ",
          "Longitudinal Networks"
        ),
        subject_text
      ),
      fp_p = paragraph_format
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "Please consider our manuscript as a Method Article for publication in ",
      "MethodsX. It addresses a practical problem in longitudinal network ",
      "analysis: conventional annual indicators and graphon-inspired distances ",
      "capture different forms of structural change and should be compared ",
      "without treating either evidence family as a universal substitute for ",
      "the other."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "The proposed workflow retains interpretable indicators of activity, ",
      "connectivity, reciprocity, and concentration; adds node-aligned fitted-",
      "kernel, spectral, and structural-role distances; and integrates them ",
      "through channel-specific anomaly scoring, known-truth calibration, and ",
      "transparent consensus reporting. Agreement supports a broad structural-",
      "anomaly claim, while disagreement identifies the type, scale, or ",
      "specification sensitivity of the change."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "Validation combines controlled simulations with 19 annual Hungarian ",
      "higher-education application networks covering 2006–2024 and 175 aligned ",
      "micro-regions. The calibrated maximum null false-positive rate is 5%, and ",
      "the empirical analysis identifies 2017–2018 as the strongest six-",
      "perspective consensus. The complete R pipeline, aggregate inputs, ",
      "validation tables, and manuscript-generation sources accompany the ",
      "submission."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "We believe the manuscript fits MethodsX because its primary contribution ",
      "is an auditable, reusable workflow with explicit assumptions, validation ",
      "gates, robustness checks, and publication-ready outputs. The manuscript is ",
      "original and is not under consideration elsewhere. Both authors have ",
      "approved the submitted version and declare no known competing financial ",
      "interests or personal relationships that could have influenced the work."
    )
  )
  doc <- add_text_paragraph(
    doc,
    "Thank you for considering our manuscript.",
    body_text,
    officer::fp_par(text.align = "left", padding.bottom = 10)
  )
  doc <- add_text_paragraph(doc, "Sincerely,", body_text, compact_format)
  doc <- add_text_paragraph(doc, "Zsolt T. Kosztyán", body_bold, compact_format)
  doc <- add_text_paragraph(
    doc, "Corresponding author, on behalf of both authors", small_text, compact_format
  )

  print(doc, target = output_path)
  invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE))
}

# PRE-SUBMISSION GATE: before sending the letter, confirm that both authors have
# approved the exact submitted manuscript and that the originality/exclusivity
# and conflict-of-interest statements remain accurate.
