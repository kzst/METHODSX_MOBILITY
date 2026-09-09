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
      officer::ftext("Zsolt Tibor Kosztyán", name_text),
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
    paste0(
      "Egyetem Street 10, 8200 Veszprém, Hungary | ",
      "kosztyan.zsolt@gtk.uni-pannon.hu | ORCID 0000-0001-7345-8336"
    ),
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
      officer::ftext("Re: Revised submission MEX-D-26-02151 — ", body_bold),
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
      "Thank you for inviting us to amend and resubmit this MethodsX Method ",
      "Article. We have addressed the two pre-evaluation requirements in the ",
      "editorial letter."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "First, the manuscript has been rebuilt using the exact MethodsX Method ",
      "Article Template supplied through the link in the editorial decision. ",
      "Every mandatory field has been completed, the prescribed section ",
      "headings and order have been retained, the title contains 18 words, the ",
      "abstract contains three method bullets and remains below 200 words, and ",
      "the Background remains below 500 words. All instructional text and ",
      "comments have been removed. The graphical abstract is supplied ",
      "separately in PDF and high-resolution PNG formats."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "Second, the manuscript and submission metadata identify both authors and ",
      "their ORCIDs: Zsolt Tibor Kosztyán (0000-0001-7345-8336) and Kornél ",
      "Dénes (0009-0002-9527-6124). Both Editorial Manager author records will ",
      "use kosztyan.zsolt@gtk.uni-pannon.hu and morgosz@student.elte.hu, ",
      "respectively, and both authors will complete the ",
      "requested authorship verification before editorial evaluation."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "The methodological presentation now states the precise adjustment to ",
      "existing approaches: conventional network indicators are retained as ",
      "interpretable evidence, directed weighted block-kernel distances are ",
      "compared on raw and unit-mean scales, and both families are integrated ",
      "through within-channel anomaly profiles and known-truth calibration. The ",
      "complete R pipeline, aggregate inputs, validation tables, graphical-",
      "abstract source, and template-compliance report accompany the revision."
    )
  )
  doc <- add_text_paragraph(
    doc,
    paste0(
      "The manuscript is original and is not under consideration elsewhere. The authors declare no known ",
      "competing financial interests or personal relationships that could have ",
      "influenced the work."
    )
  )
  doc <- add_text_paragraph(
    doc,
    "Thank you for reconsidering the revised submission.",
    body_text,
    officer::fp_par(text.align = "left", padding.bottom = 10)
  )
  doc <- add_text_paragraph(doc, "Sincerely,", body_text, compact_format)
  doc <- add_text_paragraph(doc, "Zsolt Tibor Kosztyán", body_bold, compact_format)
  doc <- add_text_paragraph(
    doc, "Corresponding author, on behalf of both authors", small_text, compact_format
  )

  print(doc, target = output_path)
  invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE))
}

# PRE-SUBMISSION GATE: before sending the letter, confirm that both authors have
# approved the exact submitted manuscript and that the originality/exclusivity
# and conflict-of-interest statements remain accurate.
