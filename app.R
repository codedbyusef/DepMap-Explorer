# =============================================================================
# DepMap Explorer
#
# Project context
# ---------------
# This app explores DepMap CRISPR gene dependency data to identify cancer cell
# lines that depend on death regulators. The focus is on separating death
# suppressors (BCL2L1, GPX4, MCL1), whose loss makes cells die, from death
# executioners (CASP3, BAX, MLKL), which drive apoptosis when active.
#
# The WT vs altered split is used as a positive control against TP53 and PTEN.
# The project spec gives Cohen's d = 1.83 (TP53) and 0.63 (PTEN). On the
# pan-cancer DepMap data with the definition below, the values are d = -1.76
# (TP53) and -0.40 (PTEN) on 24Q4. The spec values are NOT reproduced; the source of
# the spec values is unknown. "Altered" means a mutation call is present OR
# log2 copy number is at or below DepMap's loss threshold (0.731).
#
# Data: preprocess.R converts the DepMap 24Q4 CSV files to Parquet files in data/ once.
# The app reads them lazily with arrow and never touches the raw files
# at runtime. See data/README.md.
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
library(arrow)

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

DATA_DIR <- "data"

# Cell death gene panel for the Dependency Plot tab. Roles follow the project
# figure: death suppressors, executioners, and growth suppressors / other.
# Colours are lightened from the figure so they read on the dark background.
ROLE_COLOURS <- c(
  "death suppressor" = "#4da3e0",
  "executioner" = "#b8b4aa",
  "growth suppressor / other" = "#f0a030",
  "added gene" = "#8f98a3"
)
# Role groups offered as presets. All are selected by default.
PRESET_ROLES <- names(ROLE_COLOURS)[1:3]
DEFAULT_ROLES <- PRESET_ROLES
GENE_PANEL <- local({
  roles <- PRESET_ROLES
  data.frame(
    gene_name = c(
      "MCL1", "BCL2", "SLC7A11", "BCL2L1", "GPX4",
      "CASP7", "APAF1", "MLKL", "CASP8", "BID", "RIPK3", "CASP1", "RIPK1", "BAX",
      "BAK1", "GSDMD", "AIFM1", "CASP3", "CASP9", "PARP1",
      "VHL", "TP53", "PTEN"
    ),
    role = roles[c(rep(1, 5), rep(2, 15), rep(3, 3))],
    stringsAsFactors = FALSE
  )
})

pretty_label <- function(x) tools::toTitleCase(gsub("_", " ", x))

format_p <- function(p) formatC(p, format = "e", digits = 1)

# Open the preprocessed Parquet files. The big tables (crispr, copy number,
# mutations) stay on disk and are filtered through arrow per gene; the small
# ones are read into memory once when the app starts and shared by all sessions.
open_data <- function(dir = DATA_DIR) {
  path <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) {
      stop("Missing ", p, ". Run `Rscript get_data.R` to download the data files, or `Rscript preprocess.R` to build them.", call. = FALSE)
    }
    p
  }
  crispr <- open_dataset(path("crispr.parquet"))
  list(
    crispr = crispr,
    cn = open_dataset(path("copy_number.parquet")),
    mut = open_dataset(path("mutations.parquet")),
    meta = read_parquet(path("meta.parquet")),
    effects = read_parquet(path("lineage_effects.parquet")),
    go = read_parquet(path("go_enrichment.parquet")),
    genes = sort(unique(read_parquet(path("lineage_effects.parquet"), col_select = "gene_name")$gene_name))
  )
}

DATA <- open_data()

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
  title = "DepMap Explorer",
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
      options = list(
        placeholder = "Type to search gene symbol",
        openOnFocus = FALSE,
        maxOptions = 50
      )
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
      conditionalPanel(
        condition = "input.gene && input.gene !== ''",
        downloadButton("dl_png", "Download PNG", class = "btn-sm mb-2")
      ),
      withSpinner(plotOutput("dist_plot", height = "520px"), color = "#c9a86a")
    ),

    nav_panel(
      "Cell Line Table",
      conditionalPanel(
        condition = "input.gene && input.gene !== ''",
        downloadButton("dl_csv", "Download CSV", class = "btn-sm mb-2")
      ),
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
    ),

    nav_panel(
      "Dependency Plot",
      checkboxGroupInput(
        "dep_roles", "Gene groups",
        choices = PRESET_ROLES, selected = DEFAULT_ROLES, inline = TRUE
      ),
      selectizeInput(
        "dep_genes", "Genes on the plot",
        choices = NULL, multiple = TRUE, width = "100%",
        options = list(placeholder = "Type to add any gene", maxOptions = 50)
      ),
      tags$p(
        class = "small text-muted",
        "Ticking a group selects its genes and keeps any extra genes you added.",
        "Significance is corrected across the genes shown."
      ),
      downloadButton("dl_dep_png", "Download PNG", class = "btn-sm mb-2"),
      withSpinner(plotOutput("dep_plot", height = "auto"), color = "#c9a86a")
    ),

    nav_panel(
      "Cell Death GO",
      radioButtons(
        "go_backend", NULL,
        choices = c("Static (ggplot2)" = "ggplot2", "Interactive (plotly)" = "plotly"),
        inline = TRUE
      ),
      conditionalPanel(
        condition = "input.go_backend === 'ggplot2'",
        withSpinner(plotOutput("go_plot_static", height = "640px"), color = "#c9a86a")
      ),
      conditionalPanel(
        condition = "input.go_backend === 'plotly'",
        withSpinner(plotly::plotlyOutput("go_plot_interactive", height = "640px"), color = "#c9a86a")
      ),
      tags$p(
        class = "small text-muted mt-2",
        "GO Biological Process terms under cell death (GO:0008219) enriched among genes this lineage",
        "depends on more than other lineages (Cohen's d <= -0.5, q < 0.05). Colour is -log10 adjusted p."
      )
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {

  # Fill the dropdown and gene search once the data is available.
  observe({
    d <- DATA

    lineages <- sort(unique(stats::na.omit(d$meta$lineage)))
    updateSelectInput(
      session, "cancer",
      choices = c(
        setNames(ALL_CANCER, "All Cancer Types"),
        setNames(lineages, pretty_label(lineages))
      ),
      selected = DEFAULT_CANCER
    )

    genes <- d$genes
    updateSelectizeInput(session, "gene", choices = genes, selected = DEFAULT_GENE, server = FALSE)
  })

  # One row per cell line with CRISPR data for the selected gene, joined with
  # metadata, copy number and mutation status. Depends only on the gene, so
  # changing the cancer type does not recompute it.
  gene_frame <- reactive({
    d <- DATA
    gene <- req(input$gene)

    # Each filter runs inside arrow, so only this gene's rows are read from disk.
    cn <- d$cn %>%
      filter(gene_name == !!gene) %>%
      select(depmap_id, log_copy_number) %>%
      collect() %>%
      distinct(depmap_id, .keep_all = TRUE)

    mut_ids <- d$mut %>%
      filter(gene_name == !!gene) %>%
      select(depmap_id) %>%
      collect() %>%
      pull(depmap_id) %>%
      unique()

    d$crispr %>%
      filter(gene_name == !!gene) %>%
      select(depmap_id, gene_name, dependency) %>%
      collect() %>%
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
      )
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

  # -- Tab 5: dependency plot (cell death panel, lineage vs all others) ---------

  # Gene picker. Fill the choices once; ticking role groups rewrites the
  # selection as (extra genes the user added) + (genes of the ticked groups).
  observe({
    updateSelectizeInput(
      session, "dep_genes", choices = DATA$genes, server = FALSE,
      selected = GENE_PANEL$gene_name[GENE_PANEL$role %in% DEFAULT_ROLES]
    )
  })

  observeEvent(input$dep_roles, ignoreNULL = FALSE, ignoreInit = TRUE, {
    extra <- setdiff(isolate(input$dep_genes), GENE_PANEL$gene_name)
    updateSelectizeInput(
      session, "dep_genes",
      selected = c(GENE_PANEL$gene_name[GENE_PANEL$role %in% input$dep_roles], extra)
    )
  })

  # Lineage vs all other lineages for each selected gene. Cohen's d is lineage minus
  # others (negative = lineage more dependent). Significance is a Wilcoxon test
  # with Benjamini-Hochberg correction across the genes shown, the same family
  # of tests the plot shows. Scores are read from the crispr Parquet file for the
  # selected genes only.
  panel_effects <- reactive({
    cancer <- req(input$cancer)
    validate(need(
      cancer != ALL_CANCER,
      "Choose a cancer type in the sidebar. This plot compares one lineage against all other lineages."
    ))
    genes <- req(input$dep_genes)
    pc <- DATA$crispr %>%
      filter(gene_name %in% !!genes) %>%
      select(depmap_id, gene_name, dependency) %>%
      collect() %>%
      inner_join(select(DATA$meta, depmap_id, lineage), by = "depmap_id") %>%
      filter(!is.na(lineage))
    res <- lapply(split(pc, pc$gene_name), function(g) {
      a <- g$dependency[g$lineage == cancer]
      b <- g$dependency[g$lineage != cancer]
      if (length(a) < 5 || length(b) < 2) return(NULL)
      data.frame(
        gene_name = g$gene_name[1], n_lineage = length(a), n_other = length(b),
        cohen_d = effsize::cohen.d(a, b)$estimate,
        p = suppressWarnings(wilcox.test(a, b)$p.value)
      )
    })
    df <- bind_rows(res)
    validate(need(nrow(df) > 0, "Not enough cell lines in this lineage (at least 5 are needed)."))
    df$q <- p.adjust(df$p, method = "BH")
    df <- left_join(df, GENE_PANEL, by = "gene_name")
    df$role[is.na(df$role)] <- "added gene"
    df
  })

  dep_plot <- reactive({
    df <- panel_effects()
    df <- df[order(df$cohen_d), ]
    df$gene_name <- factor(df$gene_name, levels = df$gene_name)
    df$sig <- df$q < 0.05
    df$role <- droplevels(factor(df$role, levels = names(ROLE_COLOURS)))

    ggplot(df, aes(x = cohen_d, y = gene_name)) +
      geom_vline(xintercept = 0, colour = "#5b6470") +
      geom_segment(aes(x = 0, xend = cohen_d, yend = gene_name, colour = role), alpha = 0.55, linewidth = 0.7) +
      geom_point(data = df[!df$sig, ], aes(colour = role), size = 2.6) +
      geom_point(data = df[df$sig, ], aes(fill = role), shape = 21, colour = COL_FG, size = 4, stroke = 1) +
      scale_colour_manual(values = ROLE_COLOURS, name = NULL) +
      scale_fill_manual(values = ROLE_COLOURS, guide = "none") +
      scale_x_continuous(breaks = seq(-2, 3, by = 0.5)) +
      labs(
        title = paste("Dependency by gene:", cohort_label()),
        subtitle = "Negative d: lineage more dependent than the others.\nOutlined points: q < 0.05 (Wilcoxon, BH across the genes shown)",
        x = sprintf(
          "Cohen's d, %s (n=%d) vs all other lineages (n=%s)",
          tolower(cohort_label()), max(df$n_lineage), format(max(df$n_other), big.mark = ",")
        ),
        y = NULL
      ) +
      theme_academic() +
      theme(
        panel.grid.major.x = element_line(colour = "#2a2e36"),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(
          colour = unname(ROLE_COLOURS[as.character(df$role)]),
          face = ifelse(df$sig, "bold", "plain")
        ),
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0)
      )
  })

  # Plot height grows with the number of genes.
  dep_height <- function(n, per = 26, base = 150) max(380, base + per * n)

  output$dep_plot <- renderPlot({
    dep_plot()
  }, bg = COL_BG, height = function() {
    n <- tryCatch(nrow(panel_effects()), error = function(e) 10)
    dep_height(n)
  })

  output$dl_dep_png <- downloadHandler(
    filename = function() paste0(input$cancer, "_gene_dependency.png"),
    content = function(file) {
      ggsave(file, plot = dep_plot(), width = 9, height = dep_height(nrow(panel_effects())) / 100,
             dpi = 200, bg = COL_BG)
    }
  )

  # -- Tab 6: cell death GO sunburst (GOfan) -------------------------------------

  go_input <- reactive({
    cancer <- req(input$cancer)
    validate(need(
      cancer != ALL_CANCER,
      "Choose a cancer type in the sidebar. The enrichment is computed per lineage."
    ))
    df <- DATA$go %>%
      filter(lineage == cancer) %>%
      select(ID, p.adjust, qvalue, count) %>%
      as.data.frame()
    validate(need(nrow(df) > 0, "No significantly enriched cell death terms for this lineage."))
    df
  })

  go_sunburst <- function(backend) {
    df <- go_input()
    org <- org.Hs.eg.db::org.Hs.eg.db
    g <- GOfan::getGraph(df, org = org, onto = "BP")
    GOfan::sunburstGO(
      df, org = org, g = g, fill = "qvalue", onto = "BP",
      filterNodesByEdgeNumber = 0, plotBy = backend
    )
  }

  go_static <- function() {
    go_sunburst("ggplot2") +
      scale_fill_gradient(low = "#d4e6f5", high = "#3f86bd", name = "-log10\nadj. p") +
      labs(title = paste("Cell death GO terms:", cohort_label())) +
      theme(
        plot.background = element_rect(fill = COL_BG, colour = NA),
        panel.background = element_rect(fill = COL_BG, colour = NA),
        legend.background = element_rect(fill = COL_BG, colour = NA),
        text = element_text(colour = COL_FG),
        plot.title = element_text(face = "bold")
      )
  }

  output$go_plot_static <- renderPlot({
    req(input$go_backend == "ggplot2")
    go_static()
  }, bg = COL_BG)

  output$go_plot_interactive <- plotly::renderPlotly({
    req(input$go_backend == "plotly")
    go_sunburst("plotly")
  })

  # -- Reset -------------------------------------------------------------------

  observeEvent(input$reset, {
    updateSelectInput(session, "cancer", selected = DEFAULT_CANCER)
    updateSelectizeInput(session, "gene", selected = DEFAULT_GENE)
    updateCheckboxInput(session, "overlay", value = TRUE)
  })
}

shinyApp(ui, server)
