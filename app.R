# =============================================================================
# DepMap Cell Death Explorer
#
# Project context
# ---------------
# This app explores DepMap CRISPR gene dependency data to identify cancer cell
# lines that depend on cell death regulators. The focus is on separating death
# suppressors (BCL2L1, GPX4, MCL1), whose loss makes cells die, from death
# executioners (CASP3, BAX, MLKL), which drive cell death when active.
#
# The ccRCC kidney use case is the worked example. It uses a WT vs altered
# split as a positive control, validated against TP53 (Cohen's d = 1.83) and
# PTEN (d = 0.63). "Altered" means a mutation call is present OR log2 copy
# number is at or below DepMap's loss threshold (0.731).
#
# Data: all four tables are loaded from Bioconductor through ExperimentHub
# using the depmap package. See data/README.md.
#
# Gene effect convention: CRISPR scores are DepMap Chronos gene effects.
# Below -0.5 = essential (dependency). Above 0 = growth benefit.
# Cohen's d is computed as (altered - WT) / pooled SD, so a positive d means
# WT lines are more dependent (lower scores) than altered lines.
# =============================================================================

library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(DT)
library(shinycssloaders)
library(effsize)

ALL_CANCER <- "__all__"
DEFAULT_CANCER <- ALL_CANCER
DEFAULT_GENE <- "BCL2L1"
ESSENTIAL_CUT <- -0.5
# DepMap's discretized copy-number calls (log2(CN ratio + 1) scale, diploid ~ 1):
# heterozygous loss is (0.521, 0.731], deep deletion is <= 0.521. Any loss call
# is therefore <= 0.7311832. Source: DepMap forum, "Defining deep deletions and
# amplifications" and "Classifying copy number alterations".
CN_LOSS_CUT <- 0.7311832
COL_ESSENTIAL <- "#e05d5d"
COL_NEUTRAL <- "#9aa0a6"
COL_GROWTH <- "#6f9fd8"
COL_BG <- "#15171c"
COL_FG <- "#e4e2dd"

pretty_label <- function(x) tools::toTitleCase(gsub("_", " ", x))

format_p <- function(p) formatC(p, format = "e", digits = 1)

theme_academic <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = COL_BG, colour = NA),
      panel.background = element_rect(fill = COL_BG, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#2a2e36"),
      text = element_text(colour = COL_FG),
      axis.text = element_text(colour = "#b8b4aa"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "#b8b4aa"),
      legend.background = element_rect(fill = COL_BG, colour = NA),
      legend.key = element_rect(fill = COL_BG, colour = NA)
    )
}

# Two-group comparison of WT vs altered lines. Returns NULL when either group
# has fewer than two lines, since neither test is meaningful then.
compare_groups <- function(df) {
  wt <- df$dependency[df$group == "WT"]
  alt <- df$dependency[df$group == "Altered"]
  if (length(wt) < 2 || length(alt) < 2) return(NULL)
  list(
    n_wt = length(wt),
    n_alt = length(alt),
    p = suppressWarnings(wilcox.test(alt, wt)$p.value),
    d = effsize::cohen.d(alt, wt)$estimate
  )
}

interpret_result <- function(cmp, gene, df) {
  p_txt <- format_p(cmp$p)
  if (cmp$p >= 0.05) {
    sprintf(
      "No significant difference in %s dependency between WT and altered lines (p = %s, d = %.2f).",
      gene, p_txt, cmp$d
    )
  } else if (cmp$d > 0) {
    tail <- if (any(df$cn_loss)) ", consistent with copy number loss relieving the dependency" else ""
    sprintf(
      "WT lines are significantly more dependent on %s than altered lines (p = %s, d = %.2f)%s.",
      gene, p_txt, cmp$d, tail
    )
  } else {
    sprintf(
      "Altered lines are significantly more dependent on %s than WT lines (p = %s, d = %.2f).",
      gene, p_txt, cmp$d
    )
  }
}

group_stats <- function(df, label) {
  x <- df$dependency
  if (length(x) == 0) {
    return(tibble(Group = label, n = 0L, Median = NA_real_, Mean = NA_real_, Min = NA_real_, Max = NA_real_))
  }
  tibble(
    Group = label,
    n = length(x),
    Median = round(median(x), 3),
    Mean = round(mean(x), 3),
    Min = round(min(x), 3),
    Max = round(max(x), 3)
  )
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

app_theme <- bs_theme(
  version = 5,
  bg = COL_BG,
  fg = COL_FG,
  primary = "#c9a86a",
  secondary = "#5b6470",
  base_font = bslib::font_collection("Georgia", "serif"),
  heading_font = bslib::font_collection("Georgia", "serif")
)

ui <- page_sidebar(
  title = "DepMap Cell Death Explorer",
  theme = app_theme,
  sidebar = sidebar(
    width = 300,
    selectInput(
      "cancer", "Cancer type",
      choices = c("Loading..." = ""),
      selected = ""
    ),
    selectizeInput(
      "gene", "Gene",
      choices = NULL,
      selected = DEFAULT_GENE,
      options = list(placeholder = "Search gene symbol")
    ),
    checkboxInput("overlay", "Split by WT vs Altered", value = TRUE),
    actionButton("reset", "Reset to defaults", class = "btn-outline-secondary w-100"),
    tags$hr(),
    tags$small(
      paste0("Altered = mutation call present OR log2 copy number <= ", CN_LOSS_CUT, "."),
      "Scores: below −0.5 essential, above 0 growth benefit."
    )
  ),
  navset_card_tab(
    id = "tabs",

    nav_panel(
      "Distribution",
      div(
        class = "callout mb-3 p-3",
        style = "border-left: 3px solid #c9a86a; background: #1d2026;",
        tags$strong("Pan-cancer context — BCL2L1: "),
        "a death suppressor. Across 1,086 cell lines the median gene effect is −0.75 and 73% of lines are below −0.5. ",
        "Among WT vs altered lines (mutation or copy loss), altered lines are less dependent (Wilcoxon p = 3.5e-02, d = 0.45, n altered = 23)."
      ),
      downloadButton("dl_png", "Download PNG", class = "btn-sm mb-2"),
      withSpinner(plotOutput("dist_plot", height = "520px"), color = "#c9a86a")
    ),

    nav_panel(
      "Cell Line Table",
      downloadButton("dl_csv", "Download CSV", class = "btn-sm mb-2"),
      withSpinner(DTOutput("cell_table"), color = "#c9a86a")
    ),

    nav_panel(
      "Summary Statistics",
      withSpinner(tableOutput("stats_table"), color = "#c9a86a"),
      uiOutput("stats_text")
    ),

    nav_panel(
      "Gene Summary",
      withSpinner(uiOutput("gene_summary"), color = "#c9a86a")
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {

  # Load the four depmap tables once per session. The reactive caches its value
  # for the life of the session, so changing inputs never triggers a reload.
  datasets <- reactive({
    withProgress(message = "Loading DepMap data from ExperimentHub", value = 0, {
      incProgress(0.1, detail = "CRISPR gene effects")
      crispr <- depmap::depmap_crispr() %>%
        select(depmap_id, gene_name, dependency)

      incProgress(0.25, detail = "Cell line metadata")
      meta <- depmap::depmap_metadata() %>%
        select(depmap_id, cell_line_name, primary_disease, lineage) %>%
        distinct(depmap_id, .keep_all = TRUE)

      incProgress(0.25, detail = "Copy number (large table)")
      cn <- depmap::depmap_copyNumber() %>%
        select(depmap_id, gene_name, log_copy_number)

      incProgress(0.25, detail = "Mutation calls")
      mut <- depmap::depmap_mutationCalls() %>%
        select(depmap_id, gene_name)

      list(crispr = crispr, meta = meta, cn = cn, mut = mut)
    })
  })

  # Fill the dropdown and gene search once the data is available.
  observe({
    d <- datasets()

    lineages <- sort(unique(stats::na.omit(d$meta$lineage)))
    updateSelectInput(
      session, "cancer",
      choices = c(
        setNames(ALL_CANCER, "All Cancer Types"),
        setNames(lineages, pretty_label(lineages))
      ),
      selected = DEFAULT_CANCER
    )

    genes <- sort(unique(d$crispr$gene_name))
    updateSelectizeInput(session, "gene", choices = genes, selected = DEFAULT_GENE, server = TRUE)
  })

  # One row per cell line with CRISPR data for the selected gene, joined with
  # metadata, copy number and mutation status. Depends only on the gene, so
  # changing the cancer type does not recompute it.
  gene_frame <- reactive({
    d <- datasets()
    gene <- req(input$gene)

    cn <- d$cn %>%
      filter(gene_name == !!gene) %>%
      select(depmap_id, log_copy_number) %>%
      distinct(depmap_id, .keep_all = TRUE)

    mut_ids <- d$mut %>%
      filter(gene_name == !!gene) %>%
      pull(depmap_id) %>%
      unique()

    d$crispr %>%
      filter(gene_name == !!gene) %>%
      left_join(d$meta, by = "depmap_id") %>%
      left_join(cn, by = "depmap_id") %>%
      mutate(
        mutated = depmap_id %in% mut_ids,
        cn_loss = !is.na(log_copy_number) & log_copy_number <= CN_LOSS_CUT,
        altered = mutated | cn_loss,
        group = if_else(altered, "Altered", "WT"),
        status = case_when(
          mutated & cn_loss ~ "Mutated + CN loss",
          mutated ~ "Mutated",
          cn_loss ~ "CN loss",
          TRUE ~ "WT"
        )
      )
  })

  # Gene data restricted to the selected cancer type.
  cohort <- reactive({
    df <- gene_frame()
    cancer <- req(input$cancer)
    if (cancer != ALL_CANCER) {
      df <- filter(df, lineage == cancer)
    }
    df
  })

  cohort_label <- reactive({
    if (input$cancer == ALL_CANCER) "All Cancer Types" else pretty_label(input$cancer)
  })

  # -- Tab 1: distribution -----------------------------------------------------

  dist_plot <- reactive({
    df <- cohort()
    gene <- input$gene
    validate(need(nrow(df) > 0, "No cell lines with CRISPR data for this gene and cancer type."))

    df <- df %>%
      mutate(band = case_when(
        dependency < ESSENTIAL_CUT ~ "Essential (< −0.5)",
        dependency > 0 ~ "Growth benefit (> 0)",
        TRUE ~ "Neutral"
      ))

    band_colours <- setNames(
      c(COL_ESSENTIAL, COL_NEUTRAL, COL_GROWTH),
      c("Essential (< −0.5)", "Neutral", "Growth benefit (> 0)")
    )

    if (input$overlay) {
      df$group <- factor(df$group, levels = c("WT", "Altered"))
      cmp <- compare_groups(df)
      subtitle <- if (is.null(cmp)) {
        "Not enough WT or altered lines for a statistical comparison."
      } else {
        sprintf(
          "Wilcoxon p = %s    Cohen's d = %.2f    (WT n = %d, altered n = %d)",
          format_p(cmp$p), cmp$d, cmp$n_wt, cmp$n_alt
        )
      }

      p <- ggplot(df, aes(x = group, y = dependency)) +
        geom_violin(fill = "#3a3f47", colour = NA, trim = TRUE) +
        geom_jitter(aes(colour = band, shape = group), width = 0.15, size = 2, alpha = 0.85) +
        labs(x = NULL, shape = "Group")
    } else {
      df$cancer_label <- cohort_label()
      subtitle <- NULL
      p <- ggplot(df, aes(x = cancer_label, y = dependency)) +
        geom_violin(fill = "#3a3f47", colour = NA, trim = TRUE) +
        geom_jitter(aes(colour = band), width = 0.15, size = 1.9, alpha = 0.85) +
        labs(x = NULL)
    }

    p +
      geom_hline(yintercept = c(ESSENTIAL_CUT, 0), linetype = "dashed", colour = "#5b6470") +
      scale_colour_manual(values = band_colours, name = "Gene effect") +
      labs(
        title = paste(gene, "CRISPR gene effect in", cohort_label()),
        subtitle = subtitle,
        y = "Gene effect score"
      ) +
      theme_academic()
  })

  output$dist_plot <- renderPlot({
    dist_plot()
  }, bg = COL_BG)

  output$dl_png <- downloadHandler(
    filename = function() paste0(input$gene, "_", input$cancer, "_distribution.png"),
    content = function(file) {
      ggsave(file, plot = dist_plot(), width = 9, height = 6, dpi = 200, bg = COL_BG)
    }
  )

  # -- Tab 2: cell line table --------------------------------------------------

  table_df <- reactive({
    cohort() %>%
      transmute(
        `Cell Line Name` = cell_line_name,
        `Cancer Type` = primary_disease,
        `CRISPR Gene Effect` = round(dependency, 3),
        `Mutation Status` = status,
        `Copy Number (log2)` = round(log_copy_number, 2),
        Lineage = lineage
      ) %>%
      arrange(`CRISPR Gene Effect`)
  })

  output$cell_table <- renderDT({
    datatable(
      table_df(),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE)
    ) %>%
      formatStyle(
        "CRISPR Gene Effect",
        color = styleInterval(c(ESSENTIAL_CUT, 0), c(COL_ESSENTIAL, COL_NEUTRAL, COL_GROWTH))
      )
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0(input$gene, "_", input$cancer, "_cell_lines.csv"),
    content = function(file) {
      write.csv(table_df(), file, row.names = FALSE)
    }
  )

  # -- Tab 3: summary statistics -----------------------------------------------

  output$stats_table <- renderTable({
    df <- cohort()
    validate(need(nrow(df) > 0, "No cell lines with CRISPR data for this gene and cancer type."))

    if (input$overlay) {
      bind_rows(
        group_stats(df, "All lines"),
        group_stats(filter(df, group == "WT"), "WT"),
        group_stats(filter(df, group == "Altered"), "Altered")
      )
    } else {
      group_stats(df, "All lines")
    }
  }, striped = TRUE, digits = 3)

  output$stats_text <- renderUI({
    df <- cohort()
    if (!input$overlay || nrow(df) == 0) return(NULL)

    cmp <- compare_groups(df)
    if (is.null(cmp)) {
      return(tags$p(class = "mt-3", "Not enough WT or altered lines for a statistical comparison."))
    }

    tags$div(
      class = "mt-3",
      tags$p(
        tags$strong("Wilcoxon p = "), format_p(cmp$p), tags$br(),
        tags$strong("Cohen's d = "), sprintf("%.2f", cmp$d), tags$br(),
        tags$strong("Lines: "), sprintf("WT n = %d, altered n = %d", cmp$n_wt, cmp$n_alt)
      ),
      tags$p(tags$em(interpret_result(cmp, input$gene, df)))
    )
  })

  # -- Tab 4: gene summary -----------------------------------------------------

  output$gene_summary <- renderUI({
    df <- gene_frame()
    validate(need(nrow(df) > 0, "No CRISPR data for this gene."))

    n <- nrow(df)
    pct_essential <- mean(df$dependency < ESSENTIAL_CUT) * 100
    pct_growth <- mean(df$dependency > 0) * 100

    top <- df %>%
      filter(!is.na(lineage)) %>%
      group_by(lineage) %>%
      summarise(n = n(), median = median(dependency), .groups = "drop") %>%
      filter(n >= 3) %>%
      arrange(median) %>%
      slice_head(n = 5)

    tags$div(
      class = "p-3",
      style = "background: #1d2026; border-radius: 6px;",
      tags$h5(paste(input$gene, "across all cell lines")),
      tags$ul(
        tags$li(tags$strong("Cell lines with data: "), n),
        tags$li(tags$strong("Essential (score < −0.5): "), sprintf("%.1f%%", pct_essential)),
        tags$li(tags$strong("Growth benefit (score > 0): "), sprintf("%.1f%%", pct_growth))
      ),
      tags$h6("Strongest median dependency by cancer type (top 5)"),
      tags$p(class = "small text-muted", "Lineages with at least 3 cell lines."),
      tags$ol(
        lapply(seq_len(nrow(top)), function(i) {
          tags$li(sprintf(
            "%s: median %.2f (n = %d)",
            pretty_label(top$lineage[i]), top$median[i], top$n[i]
          ))
        })
      )
    )
  })

  # -- Reset -------------------------------------------------------------------

  observeEvent(input$reset, {
    updateSelectInput(session, "cancer", selected = DEFAULT_CANCER)
    updateSelectizeInput(session, "gene", selected = DEFAULT_GENE)
    updateCheckboxInput(session, "overlay", value = TRUE)
  })
}

shinyApp(ui, server)
