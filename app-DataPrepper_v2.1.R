library(shiny)
library(dplyr)
library(stringr)
library(DT)
library(rentrez)
library(xml2)

# Increase the maximum request size to 100MB
options(shiny.maxRequestSize = 100 * 1024^2)  # 100MB

# Increase the timeout for Shiny app to 5 minutes (300 seconds)
options(shiny.timeout = 300)


ui <- fluidPage(
  titlePanel("DataPrepper"),
  sidebarLayout(
    sidebarPanel(
      div(style = "font-size: 11px; color: #bbb; margin-bottom: 8px;",
          "Session ID: ", textOutput("session_id_display", inline = TRUE)),
      fileInput("file1", "Upload OTU.txt File", accept = c(".txt")),
      radioButtons("datatype", "Select Data Type:", choices = c("16S", "ITS")),
      hr(),
      h5("Taxonomy Correction", style = "margin: 5px 0;"),
      radioButtons("tax_method", NULL,
        choices = c(
          "Skip (no correction)" = "none",
          "NCBI Entrez"          = "ncbi"
        ),
        selected = "ncbi"
      ),
      hr(),
      actionButton("process", "Process OTU File Data"),
      downloadButton("downloadProcessed1", "Download Processed OTU File"),
      downloadButton("downloadMetadata", "Download Metadata File"),
      uiOutput("log_download_ui"),
      hr(),
      fileInput("file2", "Upload % of Library OTU File", accept = c(".txt")),
      actionButton("process2", "Reformat for PAST"),
      downloadButton("downloadProcessed2", "Download File for PAST")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Preview Data",
          verbatimTextOutput("progressText"),
          h3("Preview Data"),
          helpText("Preview shows up to the first 100 rows. Downloaded files contain the full table."),
          tableOutput("outputTable")
        ),
        tabPanel("Verify Metadata",
          uiOutput("metadata_verification_ui")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  processed_data1 <- reactiveVal(NULL)
  processed_data2 <- reactiveVal(NULL)
  metadata_val <- reactiveVal(NULL)
  correction_log_val <- reactiveVal(NULL)
  datatype_val <- reactiveVal(NULL)

  session_dir <- file.path(tempdir(), session$token)
  dir.create(session_dir, showWarnings = FALSE)
  session$onSessionEnded(function() unlink(session_dir, recursive = TRUE))

  output$session_id_display <- renderText(substr(session$token, 1, 8))

  observeEvent(input$process, {
    req(input$file1)
    req(input$datatype)
    datatype_val(input$datatype)

    temp_input1 <- input$file1$datapath

    # Show progress
    withProgress(message = "Processing OTU counts file...", value = 0, {
      incProgress(0.5, detail = "Running transformation pipeline...")
      
      # Step 2: Read the input file
      input_file1 <- read.delim(temp_input1, header = TRUE, sep = "\t")
      
      # Step 3: Separate the taxonomy column and the rest of the columns
      first_column <- input_file1 %>%
        select(2)
      
      other_columns <- input_file1 %>%
        select(6:ncol(input_file1))
      
      # Step 4: Filter out columns with non-numeric values in the other columns
      filtered_columns <- other_columns %>%
        select(where(~ all(!is.na(as.numeric(.)))))
      
      # Step 5: Combine the first column with the filtered columns
      selected_data <- bind_cols(first_column, filtered_columns)
      
      # Step 6: Remove the eighth semicolon and everything after it
      trimmed_semicolons <- selected_data %>%
        rowwise() %>%
        mutate(Taxonomy = {
          parts <- str_split(Taxonomy, ";")[[1]]
          if (length(parts) > 8) {
            str_c(parts[1:8], collapse = ";")
          } else {
            Taxonomy
          }
        }) %>%
        ungroup()
      
      # Step 7: Replace slashes in OTU names with underscores
      trimmed_slashes <- trimmed_semicolons %>%
        mutate(across(1, ~ str_replace_all(., "/", "_")))
      
      # Step 8: Replace semicolons with slashes
      replaced_semicolons <- trimmed_slashes %>%
        mutate(across(1, ~ str_replace_all(., ";", "/")))
      
      # Step 9: Remove the first backslash in the first column and everything before it
      removed_first_backslash <- replaced_semicolons %>%
        mutate(across(1, ~ str_replace(., "^.*?/", "")))
      
      # Step 10: Remove taxonomic prefixes
      removed_prefixes <- removed_first_backslash %>%
        mutate(across(1, ~ str_replace_all(., "(k__|p__|c__|o__|f__|g__|s__)", "")))
      
      # Step 11: Remove underscores
      removed_underscores <- removed_prefixes %>%
        mutate(across(1, ~ str_replace_all(., "_", "")))
      
      # Step 12: Trim trailing slashes in the Taxonomy column
      trimmed_taxonomy <- removed_underscores %>%
        mutate(Taxonomy = ifelse(str_detect(Taxonomy, "/$"), str_sub(Taxonomy, 1, -2), Taxonomy))
      
      # Step 13: Rename the Taxonomy column to O_T_U
      renamed_columns <- trimmed_taxonomy %>%
        rename(O_T_U = Taxonomy)
      
      # Drop OTUs where every sample has < 3 counts (no robust signal in any sample)
      filtered_data <- renamed_columns %>%
        filter(if_any(-O_T_U, ~ . >= 3))
      
      # Step 16: Remove rows with "Chloroplast" or "Mitochondria" in the O_T_U column
      final_data1 <- filtered_data %>%
        filter(!str_detect(O_T_U, "Chloroplast|Mitochondria"))
      
      # Create metadata file
      sample_names <- colnames(final_data1)[-1]
      metadata <- data.frame(
        Sample = sample_names,
        Treatment = case_when(
          str_detect(sample_names, "BjSa|BSM")  ~ "BjSa",
          str_detect(sample_names, "PAST|Past") ~ "Pasteurized",
          str_detect(sample_names, "RR|Rep")    ~ "Replant",
          str_detect(sample_names, "Std|STAND") ~ "Standard",
          str_detect(sample_names, "Mul|MUL")   ~ "Mulch",
          str_detect(sample_names, "Org|ORG")   ~ "Organic",
          str_detect(sample_names, "Carb|CARB") ~ "Carbon",
          str_detect(sample_names, "Plus|PLUS") ~ "Plus",
          str_detect(sample_names, "ASD|asd")   ~ "ASD",
          str_detect(sample_names, "NTC|ntc")   ~ "NTC",
          str_detect(sample_names, "FUM|fum")   ~ "Fumigated",
          str_detect(sample_names, "Comp|comp") ~ "Compost",
          TRUE                                  ~ "Enter Manually"
        ),
        Rep = vapply(
          str_extract_all(sample_names, "\\d+"),
          function(x) if (length(x) > 0) tail(x, 1) else NA_character_,
          character(1)
        ),
        stringsAsFactors = FALSE
      )

      # Update reactive value with metadata instead of immediate file write
      metadata_val(metadata)

      # Update the reactive value with the processed data
      processed_data1(final_data1)

      incProgress(0.5, detail = "OTU reformatting complete.")
    })

    # Taxonomy correction (if selected)
    tax_summary <- NULL
    if (input$tax_method != "none") {
      df <- processed_data1()

      all_genera <- vapply(df$O_T_U, function(otu) {
        parts <- str_split(otu, "/")[[1]]
        if (length(parts) >= 6) trimws(parts[6]) else NA_character_
      }, character(1))

      unique_genera <- sort(unique(all_genera[!is.na(all_genera) & all_genera != ""]))
      n_genera      <- length(unique_genera)
      taxonomy_cache <- list()
      log_rows       <- list()

      if (n_genera > 0) {
        withProgress(message = "Correcting taxonomy...", value = 0, {

          taxid_to_genus <- list()
          kingdom_filter <- if (input$datatype == "16S") {
            " AND (txid2[Subtree] OR txid2157[Subtree])"
          } else {
            " AND txid4751[Subtree]"
          }

          for (i in seq_along(unique_genera)) {
            genus <- unique_genera[i]
            incProgress(
              amount = 0.5 / n_genera,
              detail = paste0("Searching: ", i, " of ", n_genera, ": ", genus)
            )
            result <- tryCatch(
              rentrez::entrez_search(
                db     = "taxonomy",
                term   = paste0('"', genus, '"[Scientific Name] AND "genus"[Rank]',
                                kingdom_filter),
                retmax = 5
              ),
              error = function(e) list(ids = character(0))
            )
            if (length(result$ids) > 0) taxid_to_genus[[result$ids[1]]] <- genus
          }

          all_taxids <- names(taxid_to_genus)

          if (length(all_taxids) > 0) {
            fetch_batches <- split(all_taxids, ceiling(seq_along(all_taxids) / 200))
            n_fb <- length(fetch_batches)

            for (i in seq_along(fetch_batches)) {
              incProgress(
                amount = 0.3 / n_fb,
                detail = paste0("Fetching taxonomy: batch ", i, " of ", n_fb)
              )
              xml_text <- tryCatch(
                rentrez::entrez_fetch(db = "taxonomy", id = fetch_batches[[i]],
                                      rettype = "xml"),
                error = function(e) NULL
              )
              if (is.null(xml_text)) next
              doc <- tryCatch(xml2::read_xml(xml_text), error = function(e) NULL)
              if (is.null(doc)) next

              for (taxon in xml2::xml_find_all(doc, "/TaxaSet/Taxon")) {
                taxid <- xml2::xml_text(xml2::xml_find_first(taxon, "TaxId"))
                if (!taxid %in% names(taxid_to_genus)) next
                genus_name <- taxid_to_genus[[taxid]]

                ln_nodes <- xml2::xml_find_all(taxon, "LineageEx/Taxon")
                if (length(ln_nodes) == 0) next
                taxonomy_cache[[genus_name]] <- data.frame(
                  rank = sapply(ln_nodes, function(n)
                    xml2::xml_text(xml2::xml_find_first(n, "Rank"))),
                  name = sapply(ln_nodes, function(n)
                    xml2::xml_text(xml2::xml_find_first(n, "ScientificName"))),
                  stringsAsFactors = FALSE
                )
              }
            }
          } else {
            incProgress(0.3, detail = "No TaxIDs found in NCBI search.")
          }

          incProgress(0.2, detail = "Applying corrections to OTU table...")

          get_rank_val <- function(tax_df, rank_name) {
            val <- tax_df$name[tax_df$rank == rank_name]
            if (length(val) > 0 && !is.na(val[1])) tolower(trimws(val[1])) else NA_character_
          }

          rank_map       <- c("superkingdom", "phylum", "class", "order", "family", "genus")
          corrected_otus <- df$O_T_U

          for (k in seq_along(corrected_otus)) {
            otu     <- corrected_otus[k]
            parts   <- str_split(otu, "/")[[1]]
            n_parts <- length(parts)
            if (n_parts < 6) next

            original_genus <- trimws(parts[6])
            if (!original_genus %in% names(taxonomy_cache)) next

            tax      <- taxonomy_cache[[original_genus]]
            new_vals <- vapply(rank_map, function(r) get_rank_val(tax, r), character(1))

            if (input$datatype == "16S") {
              if (!is.na(new_vals["superkingdom"]) && !new_vals["superkingdom"] %in% c("bacteria", "archaea")) next
            } else {
              if (!any(tolower(tax$name) == "fungi")) next
            }

            orig_order  <- if (n_parts >= 4) parts[4] else NA_character_
            orig_family <- if (n_parts >= 5) parts[5] else NA_character_

            changed <- FALSE
            j_start <- if (input$datatype == "ITS") 2L else 1L
            for (j in j_start:6) {
              if (!is.na(new_vals[j]) && new_vals[j] != parts[j]) {
                parts[j] <- new_vals[j]
                changed  <- TRUE
              }
            }

            # Incertae sedis: SILVA/UNITE repeat the order name in the family field
            # as a placeholder when family is unresolved. If NCBI has no family rank
            # and the original family matched the original order, mark it explicitly.
            if (n_parts >= 5 && !is.na(orig_order) && !is.na(orig_family) &&
                orig_family == orig_order && is.na(new_vals["family"])) {
              parts[5] <- "incertae sedis"
              changed  <- TRUE
            }

            corrected_genus <- new_vals[6]
            if (n_parts >= 7 && !is.na(corrected_genus) &&
                corrected_genus != original_genus &&
                str_starts(parts[7], original_genus)) {
              parts[7] <- paste0(corrected_genus,
                                 str_sub(parts[7], nchar(original_genus) + 1))
              changed <- TRUE
            }

            if (changed) {
              new_otu <- paste(parts, collapse = "/")
              corrected_otus[k] <- new_otu
              log_rows[[length(log_rows) + 1]] <- data.frame(
                Original  = otu,
                Corrected = new_otu,
                stringsAsFactors = FALSE
              )
            }
          }

          df$O_T_U <- corrected_otus
          processed_data1(df)
        })

        db_label      <- "NCBI Entrez"
        n_resolved    <- length(taxonomy_cache)
        n_not_found   <- n_genera - n_resolved
        n_corrections <- length(log_rows)

        correction_log_val(list(
          timestamp     = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          datatype      = input$datatype,
          database      = db_label,
          n_genera      = n_genera,
          n_resolved    = n_resolved,
          n_not_found   = n_not_found,
          n_corrections = n_corrections,
          corrections   = if (n_corrections > 0) do.call(rbind, log_rows) else NULL
        ))

        tax_summary <- paste0(
          db_label, " correction: ", n_corrections, " OTU entries updated. ",
          n_resolved, "/", n_genera, " genera resolved."
        )
      } else {
        tax_summary <- "No genera found for taxonomy correction."
      }
    } else {
      correction_log_val(NULL)
    }

    # Display processed data
    output$outputTable <- renderTable({
      req(processed_data1())
      head(processed_data1(), 100)
    })

    progress_msg <- paste0(
      "OTU reformatting complete.",
      if (!is.null(tax_summary)) paste0(" ", tax_summary) else ""
    )
    output$progressText <- renderText({
      progress_msg
    })
  })
  
  # Allow user to download processed file
  output$downloadProcessed1 <- downloadHandler(
    filename = function() {
      req(datatype_val())
      paste0(format(Sys.time(), "%m%d%Y_%H%M%S"), "_", datatype_val(), "_preprocessed_for_explicet.txt")
    },
    content = function(file) {
      req(processed_data1())
      write.table(processed_data1(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )
  
  # Download handler for the metadata file
  output$downloadMetadata <- downloadHandler(
    filename = function() {
      req(datatype_val())
      paste0(format(Sys.time(), "%m%d%Y_%H%M%S"), "_", datatype_val(), "_metadata.txt")
    },
    content = function(file) {
      req(metadata_val())
      write.table(metadata_val(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )
  
  # Metadata verification UI
  output$metadata_verification_ui <- renderUI({
    if (is.null(metadata_val())) {
      return(div(style = "padding: 20px; text-align: center;",
                 p("Please upload and process OTU data first to generate metadata for verification.")))
    }

    tagList(
      div(
        style = "margin-bottom: 20px; padding: 15px; background-color: #fcfcfc; border: 1px solid #ddd; border-radius: 5px;",
        h4("Bulk Update Tools"),
        fluidRow(
          column(4,
                 textInput("bulkTreatment", "Treatment Name", ""),
                 actionButton("applyTreatment", "Set Treatment for Selected", class = "btn-primary")
          ),
          column(4,
                 textInput("bulkRep", "Rep Value", ""),
                 actionButton("applyRep", "Set Rep for Selected", class = "btn-primary"),
                 br(), br(),
                 actionButton("fillReps", "Fill Rep Sequence (1, 2, 3...)", class = "btn-info")
          ),
          column(4,
                 br(), br(),
                 helpText("Select rows in the table below to apply bulk updates."),
                 br(),
                 helpText("Set Rep for Selected applies one value to all selected rows (e.g. group all replicate-1 samples). Fill Rep Sequence numbers the selected rows 1, 2, 3...")
          )
        )
      ),
      div(
        style = "margin-bottom: 20px; padding: 15px; background-color: #fcfcfc; border: 1px solid #ddd; border-radius: 5px;",
        h4("Sample Name Tools"),
        helpText("Clean up sample names (e.g. strip a sequencing-ID tag). You can also edit any Sample cell directly in the table below. All edits rename the matching column in the Processed OTU File and update the preview automatically."),
        fluidRow(
          column(4, textInput("findSample", "Find", "")),
          column(4, textInput("replaceSample", "Replace with", "")),
          column(4,
                 checkboxInput("sampleRegex", "Use regular expression", FALSE),
                 actionButton("replaceSelected", "Replace in Selected", class = "btn-primary"),
                 br(), br(),
                 actionButton("replaceAll", "Replace in All Rows", class = "btn-warning")
          )
        ),
        uiOutput("sample_rename_status"),
        div(style = "max-height: 220px; overflow-y: auto; margin-top: 8px;",
            tableOutput("sample_rename_preview"))
      ),
      div(
        style = "margin-bottom: 10px; padding: 10px; background-color: #eaf3fb; border-left: 4px solid #2c7fb8; border-radius: 3px;",
        strong("Tip: "), "Double-click a cell to edit individual entries."
      ),
      DTOutput("metadataTable")
    )
  })

  # Editable Metadata Table
  output$metadataTable <- renderDT({
    req(metadata_val())
    datatable(
      metadata_val(),
      editable = TRUE,
      selection = 'multiple',
      options = list(pageLength = 25),
      rownames = FALSE
    )
  })

  metadata_proxy <- dataTableProxy("metadataTable")

  # Rename a sample everywhere: metadata Sample column + matching OTU table column.
  # Returns TRUE on success, FALSE (with notification) if the change is invalid.
  rename_samples <- function(old_names, new_names) {
    df <- metadata_val()

    # Build the full proposed Sample vector so we can validate the whole set.
    proposed <- df$Sample
    idx <- match(old_names, proposed)
    proposed[idx] <- new_names

    if (any(is.na(proposed) | trimws(proposed) == "")) {
      showNotification("Sample names cannot be empty. No changes made.", type = "error")
      DT::replaceData(metadata_proxy, metadata_val(), rownames = FALSE, resetPaging = FALSE)
      return(FALSE)
    }
    if (any(duplicated(proposed))) {
      showNotification("That would create duplicate sample names. No changes made.", type = "error")
      DT::replaceData(metadata_proxy, metadata_val(), rownames = FALSE, resetPaging = FALSE)
      return(FALSE)
    }

    otu <- processed_data1()
    if (!is.null(otu)) {
      cn <- colnames(otu)
      for (i in seq_along(old_names)) {
        if (!is.na(old_names[i]) && old_names[i] != new_names[i]) {
          cn[cn == old_names[i]] <- new_names[i]
        }
      }
      colnames(otu) <- cn
      processed_data1(otu)
    }

    df$Sample <- proposed
    metadata_val(df)
    TRUE
  }

  # Handle cell edits
  observeEvent(input$metadataTable_cell_edit, {
    info    <- input$metadataTable_cell_edit
    df      <- metadata_val()
    col_idx <- info$col + 1

    if (names(df)[col_idx] == "Sample") {
      old_name <- df[info$row, col_idx]
      # rename_samples validates and reverts the display on rejection
      if (rename_samples(old_name, info$value)) {
        showNotification("Sample name updated (OTU table column renamed).", type = "message")
      }
      return()
    }

    df[info$row, col_idx] <- info$value
    metadata_val(df)
  })

  # Find & Replace across sample names
  do_sample_replace <- function(rows) {
    if (is.null(input$findSample) || input$findSample == "") {
      showNotification("Enter a 'Find' value.", type = "warning")
      return()
    }
    df  <- metadata_val()
    old <- df$Sample[rows]
    pattern <- if (isTRUE(input$sampleRegex)) input$findSample else stringr::fixed(input$findSample)
    new <- str_replace_all(old, pattern, input$replaceSample)

    if (all(old == new)) {
      showNotification("No sample names matched the 'Find' value.", type = "message")
      return()
    }
    if (rename_samples(old, new)) {
      showNotification(paste0("Renamed ", sum(old != new), " sample name(s)."), type = "message")
    }
  }

  observeEvent(input$replaceSelected, {
    req(metadata_val())
    rows <- input$metadataTable_rows_selected
    if (is.null(rows) || length(rows) == 0) {
      showNotification("Select one or more rows first.", type = "warning")
      return()
    }
    do_sample_replace(rows)
  })

  observeEvent(input$replaceAll, {
    req(metadata_val())
    do_sample_replace(seq_len(nrow(metadata_val())))
  })

  # Live preview of what Find & Replace would do to sample names.
  # Cheap: operates only on the Sample name vector, not the OTU table.
  sample_rename_calc <- reactive({
    md <- metadata_val()
    if (is.null(md)) return(NULL)
    find <- input$findSample
    if (is.null(find) || find == "") return(list(status = "empty"))
    repl <- input$replaceSample
    if (is.null(repl)) repl <- ""
    pattern <- if (isTRUE(input$sampleRegex)) find else stringr::fixed(find)
    new <- tryCatch(str_replace_all(md$Sample, pattern, repl),
                    error = function(e) NULL)
    if (is.null(new)) return(list(status = "error"))
    list(status = "ok", old = md$Sample, new = new,
         changed = which(md$Sample != new))
  })

  output$sample_rename_status <- renderUI({
    calc <- sample_rename_calc()
    if (is.null(calc) || identical(calc$status, "empty")) return(NULL)
    if (identical(calc$status, "error")) {
      return(div(style = "color:#a00; font-size:12px; margin-top:6px;",
                 "Invalid regular expression — check the pattern."))
    }
    n <- length(calc$changed)
    if (n == 0) {
      return(div(style = "color:#666; font-size:12px; margin-top:6px;",
                 "No sample names match — nothing would change."))
    }
    div(style = "color:#333; font-size:12px; margin-top:6px;",
        sprintf("%d name(s) would change%s:", n,
                if (n > 3) " (showing first 3)" else ""))
  })

  output$sample_rename_preview <- renderTable({
    calc <- sample_rename_calc()
    req(calc, identical(calc$status, "ok"), length(calc$changed) > 0)
    show_n <- head(calc$changed, 3)
    data.frame(Current = calc$old[show_n],
               Becomes = calc$new[show_n],
               check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  # Bulk Treatment Update
  observeEvent(input$applyTreatment, {
    req(metadata_val(), input$metadataTable_rows_selected)
    df <- metadata_val()
    selected_rows <- input$metadataTable_rows_selected
    df$Treatment[selected_rows] <- input$bulkTreatment
    metadata_val(df)
    showNotification("Treatment updated for selected samples", type = "message")
  })

  # Bulk Rep Update (same value for all selected rows)
  observeEvent(input$applyRep, {
    req(metadata_val(), input$metadataTable_rows_selected)
    df <- metadata_val()
    selected_rows <- input$metadataTable_rows_selected
    df$Rep[selected_rows] <- input$bulkRep
    metadata_val(df)
    showNotification("Rep value updated for selected samples", type = "message")
  })

  # Fill Rep Sequence
  observeEvent(input$fillReps, {
    req(metadata_val(), input$metadataTable_rows_selected)
    df <- metadata_val()
    selected_rows <- input$metadataTable_rows_selected
    df$Rep[selected_rows] <- seq_along(selected_rows)
    metadata_val(df)
    showNotification("Rep sequence applied to selected samples", type = "message")
  })

  output$log_download_ui <- renderUI({
    req(correction_log_val())
    downloadButton("downloadLog", "Download Correction Log", class = "btn-info")
  })

  output$downloadLog <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%m%d%Y_%H%M%S"), "_", datatype_val(), "_taxonomy_correction_log.txt")
    },
    content = function(file) {
      req(correction_log_val())
      log <- correction_log_val()
      header <- c(
        "DataPrepper Taxonomy Correction Log",
        paste(rep("=", 40), collapse = ""),
        paste0("Run Date/Time:       ", log$timestamp),
        paste0("Data Type:           ", log$datatype),
        paste0("Database:            ", log$database),
        paste0("Genera queried:      ", log$n_genera),
        paste0("Genera resolved:     ", log$n_resolved),
        paste0("Genera not found:    ", log$n_not_found, " (left unchanged)"),
        paste0("OTU entries updated: ", log$n_corrections),
        "",
        "NOTE: Species epithets (e.g., gender agreement, emended spellings) were NOT",
        "evaluated by this tool. Verify species names independently using LPSN.",
        "",
        "Corrected Entries (tab-delimited, Original vs Corrected taxonomy string):",
        paste(rep("-", 40), collapse = "")
      )
      writeLines(header, file)
      if (!is.null(log$corrections) && nrow(log$corrections) > 0) {
        write.table(log$corrections, file, sep = "\t", row.names = FALSE,
                    col.names = TRUE, quote = FALSE, append = TRUE)
      } else {
        write("\nNo taxonomy corrections were applied.", file, append = TRUE)
      }
    }
  )

  observeEvent(input$process2, {
    req(input$file2)
    
    # Immediately update the UI to indicate processing has started.
    output$progressText <- renderText({
      "Processing..."
    })
    
    # Retrieve the input file's base name and build the output file name.
    base_filename <- tools::file_path_sans_ext(input$file2$name)
    output_filename2 <- paste0(base_filename, "_processed.txt")
    output_path2 <- file.path(session_dir, output_filename2)
    
    withProgress(message = "Processing file...", value = 0, {
      incProgress(0.2, detail = "Reading file...")
      
      # Read the input file. Adjust 'skip' if necessary.
      input_file2 <- tryCatch({
        read.delim(input$file2$datapath, header = TRUE, sep = "\t",
                   colClasses = "character", skip = 1)
      }, error = function(e) {
        showNotification(paste("Error reading file:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(input_file2)) return()

      ##### STEP 1: Filter out unwanted rows #####
      filtered_data <- tryCatch({
        input_file2 %>%
          filter(!is.na(OTU_Name),
                 OTU_Name != "",
                 !OTU_Name %in% c(";All_Libs", "root"))
      }, error = function(e) {
        showNotification(paste("Error filtering data:", e$message), type = "error")
        return(NULL)
      })
      
      ##### STEP 2: Remove the unwanted columns (OTU_Name and Total) #####
      cleaned_data <- tryCatch({
        filtered_data %>%
          select(-any_of(c("OTU_Name", "Total")))
      }, error = function(e) {
        showNotification(paste("Error cleaning data:", e$message), type = "error")
        return(NULL)
      })
      
      incProgress(0.3, detail = "Cleaning data...")
      
      ##### STEP 3: Clean cells: Remove % signs and round numbers #####
      processed_temp <- tryCatch({
        cleaned_data %>%
          mutate(across(everything(), function(x) {
            # Remove percent signs
            no_percent <- str_replace_all(x, "%", "")
            # Try converting to numeric; if successful, round to two decimal places.
            num_val <- suppressWarnings(as.numeric(no_percent))
            ifelse(!is.na(num_val),
                   format(round(num_val, 2), nsmall = 2),
                   no_percent)
          }))
      }, error = function(e) {
        showNotification(paste("Error processing data:", e$message), type = "error")
        return(NULL)
      })
      
      incProgress(0.2, detail = "Formatting data...")
      
      ##### STEP 4: Transpose the data and rename columns #####
      transposed <- tryCatch({
        # Transpose the data so that sample names (original column headers) become row names.
        t_data <- as.data.frame(t(processed_temp), stringsAsFactors = FALSE)
        # Convert the row names into a proper column called "Sample" so they are not dropped
        t_data <- tibble::rownames_to_column(t_data, var = "Sample")
        # Rename the remaining columns (which correspond to OTUs) as otu1, otu2, etc.
        n_otu <- ncol(t_data) - 1
        colnames(t_data)[-1] <- paste0("otu", seq_len(n_otu))
        t_data
      }, error = function(e) {
        showNotification(paste("Error transposing data:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(transposed)) return()

      incProgress(0.3, detail = "Saving file...")
      
      ##### STEP 5: Save the processed data #####
      # Write out the file. row.names = FALSE because we've converted sample names to a column.
      write.table(transposed, file = output_path2, sep = "\t", 
                  row.names = FALSE, col.names = TRUE, quote = FALSE)
      
      processed_data2(transposed)
      
      incProgress(0.1, detail = "Processing complete.")
    })
    
    # After processing, the output progress text is updated.
    output$progressText <- renderText({
      "File processing complete."
    })
    
    #Display processed data in the Shiny app
    output$outputTable <- renderTable({
      req(processed_data2())
      head(processed_data2(), 100)
    })
  })
  
  output$downloadProcessed2 <- downloadHandler(
    filename = function() {
      base_filename <- tools::file_path_sans_ext(input$file2$name)
      paste0(base_filename, "_processed.txt")
    },
    content = function(file) {
      req(processed_data2())
      write.table(processed_data2(), file, sep = "\t",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)