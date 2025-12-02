library(shiny)
library(rhandsontable)
library(yaml)
library(zip)
library(shinythemes)

# --- 0. CUSTOM CSS ---
genopred_css <- "
  /* --- GENOPRED THEME --- */
  body {
    background-color: #2c3e50;
    color: #ecf0f1;
    font-family: 'Source Sans Pro', Calibri, Candara, Arial, sans-serif;
  }
  
  a { color: #5dade2; } /* Brighter link color */
  h1, h2, h3, h5, h6 { color: #ecf0f1; }
  h4 { color: #ffffff; font-weight: 500; font-size: 16px; } /* Pure white for h4 */
  h1 { font-weight: 700; font-size: 26px; }
  
  p, li { 
    font-weight: 300 !important; 
    color: #f8f9fa; /* Brighter text for readability */
  }

  /* Bootstrap Text Utilities Overrides for Dark Theme */
  .text-success {
    color: #2ecc71 !important; /* Brighter green for 'Active' status */
  }
  .text-muted {
    color: #bdc3c7 !important; /* Much lighter grey for 'Empty' status */
  }
  .text-info {
    color: #5dade2 !important; /* Brighter blue */
  }

  /* Navbar styling (if we used navbarPage, but applied to title for consistency) */
  .container-fluid > .row > div > h2 {
    color: #ecf0f1;
    font-weight: 700;
    padding-bottom: 10px;
    border-bottom: 1px solid #121212;
  }

  /* Sidebar and Wells */
  .well {
    background-color: #34495e;
    border: 1px solid #2c3e50;
    color: #ecf0f1;
    box-shadow: none;
    border-radius: 5px;
  }

  /* Inputs */
  .form-control {
    background-color: #ecf0f1;
    color: #2c3e50;
    border: 1px solid #bdc3c7;
    border-radius: 5px;
  }
  .selectize-input, .selectize-dropdown {
    background-color: #ecf0f1;
    color: #2c3e50;
    border-radius: 5px;
  }
  /* Selected options in multi-select inputs (Custom Style) */
  .selectize-control.multi .selectize-input > div {
    cursor: pointer;
    margin: 0 3px 3px 0;
    padding: 1px 5px;
    background: #c0cfdd !important;
    color: #333 !important;
    border: 0 solid rgba(0,0,0,0);
  }
  
  /* Buttons */
  .btn {
    white-space: normal; /* Ensure text wraps on small screens */
    height: auto;        /* Allow button to grow in height */
    border-radius: 5px;
  }
  .btn-default {
    background-color: #007bff;
    color: white;
    border: none;
  }
  .btn-default:hover {
    background-color: #0056b3;
    color: white;
  }
  .btn-success {
    background-color: #2780e3;
    border-color: #2780e3;
  }
  .btn-info {
    background-color: #007bff; /* Matching your inline_button */
    border-color: #007bff;
  }

  /* Note Box (Replaces Alerts) */
  .note-box {
    border: 1px solid #2780e3;
    padding: 10px;
    margin: 20px 0;
    background-color: #2a495f !important;
    color: #9bcaef;
    border-radius: 5px;
  }
  .note-box strong { color: #4793d5; }
  
  /* Warning specifically - Orange version of note-box */
  .alert-warning {
    border: 1px solid #e67e22;
    background-color: #4e3629 !important; /* Muted dark orange background */
    color: #f5b041; /* Light orange text */
    padding: 10px;
    margin: 20px 0;
    border-radius: 5px;
  }
  .alert-warning strong { color: #e67e22; } /* Strong text matches border */  
  
  /* Help Block Text */
  .help-block {
    display: block;
    margin-top: 5px;
    margin-bottom: 10px;
    color: #cbcbcb;
  }
  
  /* Handsontable Overrides for Dark Mode */
  .handsontable {
    color: #000000; /* Black text for maximum contrast on white cells */
    overflow: hidden;
  }
  .handsontable th {
    background-color: #1d2935;
    color: #ecf0f1;
  }
  /* Fix for unreadable header when column is selected */
  .handsontable th.ht__highlight {
    background-color: #34495e;
    color: #24516b;
  }

  /* Tabs */
  .nav-tabs>li>a {
    color: #ecf0f1;
    white-space: normal; /* Allow text to wrap */
    height: auto;
    border-radius: 5px;
  }
  .nav-tabs>li>a:hover, .nav-tabs>li>a:focus {
    background-color: #1f2c39;
    color: #ecf0f1;
  }
  .nav-tabs>li.active>a, .nav-tabs>li.active>a:focus, .nav-tabs>li.active>a:hover {
    color: #2c3e50;
    background-color: #ecf0f1;
  }
  
  /* --- Styling for code blocks (User Request) --- */
  pre, code {
    border-radius: 5px;
    font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
  }
  
  pre {
    padding: 10px;
    overflow-x: auto;
    background-color: #34495e; /* Dark background to match theme */
    color: #ecf0f1;            /* Light text */
    border: 1px solid #4e5a5b; /* Subtle border */
  }
  
  /* Styling for inline code */
  p code, li code {
    padding: 2px 4px;
    border-radius: 5px;
    border: 1px solid #7f8c8d;
    background-color: #34495e; /* Dark background */
    color: #82ccdd;            /* Orange-red to make variables pop */
  }
"

# --- 1. DEFINE COLUMN STRUCTURES ---
# Changed 'filename' to 'path' as requested
cols_gwas <- c("name", "path", "population", "n", "label", "sampling", "prevalence", "mean", "sd")
cols_target <- c("name", "path", "type", "indiv_report", "unrel") 
cols_score <- c("name", "path", "label") 
cols_groups <- c("name", "gwas", "label")

# --- 2. CONFIGURATION DICTIONARIES ---

pgs_methods_single <- c(
  "P+T Clumping (ptclump)" = "ptclump", 
  "DBSLMM" = "dbslmm", 
  "PRS-CS" = "prscs", 
  "SBayesR" = "sbayesr", 
  "SBayesRC" = "sbayesrc",
  "Lassosum" = "lassosum", 
  "Lassosum2" = "lassosum2", 
  "LDpred2" = "ldpred2", 
  "MegaPRS" = "megaprs",
  "QuickPRS" = "quickprs"
)

pgs_methods_multi <- c(
  "X-Wing" = "xwing",
  "PRS-CSx" = "prscsx"
)

pgs_methods_all <- c(pgs_methods_single, pgs_methods_multi)

# --- 3. HELPER FUNCTIONS ---
make_template <- function(cols, rows=5) {
  data_list <- list()
  for(col in cols) {
    if(col %in% c("n", "sampling", "prevalence", "mean", "sd")) {
      data_list[[col]] <- as.numeric(rep(NA, rows))
    } else {
      data_list[[col]] <- rep("", rows)
    }
  }
  df <- data.frame(data_list, stringsAsFactors = FALSE)
  df <- df[, cols, drop=FALSE]
  return(df)
}

# --- 4. UI LAYOUT ---
ui <- fluidPage(
  tags$head(tags$style(HTML(genopred_css))), # Inject Custom CSS
  br(),
  theme = shinytheme("cosmo"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Project Status"),
      uiOutput("status_panel"),
      
      hr(),
      h4("Finalize"),
      p("Once tabs are filled:"),
      
      # --- Conditional Download Button ---
      # 1. Active Button (Shown when outdir is NOT empty)
      conditionalPanel(
        condition = "input.outdir.trim() != ''",
        downloadButton("download_bundle", "Download Config Bundle (.zip)", class = "btn-success btn-block")
      ),
      # 2. Disabled Button (Shown when outdir IS empty)
      conditionalPanel(
        condition = "input.outdir.trim() == ''",
        tags$button("Download Config Bundle (.zip)", id = "download_disabled", class = "btn btn-success btn-block disabled", type = "button"),
        div(class = "text-danger", style = "font-size: 0.8em; margin-top: 5px;", icon("exclamation-circle"), " Output Directory (tab 5) is required.")
      )
      
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        
        # --- 0. INSTRUCTIONS (New) ---
        tabPanel("0. Instructions", icon = icon("info-circle"),
                 h3("Welcome to the GenoPred Configuration Wizard"),
                 p("This tool helps you generate the required configuration files to run the GenoPred pipeline."),
                 hr(),
                 h4("Workflow Steps"),
                 tags$ol(
                   tags$li(strong("Define Inputs:"), " Fill in the GWAS List (Tab 1) and Target List (Tab 2)."),
                   tags$li(strong("Optional Inputs:"), " Add External Scores (Tab 3) or define GWAS Groups (Tab 4) for meta-analysis."),
                   tags$li(strong("Select Methods:"), " Choose your PGS methods in Basic Parameters (Tab 5)."),
                   tags$li(strong("Configure Resources:"), " Set paths and advanced options in Advanced Parameters (Tab 6)."),
                   tags$li(strong("Download:"), " Click 'Download Config Bundle' in the sidebar.")
                 ),
                 hr(),
                 h4("Important: File Paths"),
                 p("All file paths provided in this wizard (e.g., for GWAS summary stats, Target data, Output directory) must be either:"),
                 tags$ul(
                   tags$li(strong("Absolute paths:"), " e.g., ", code("/home/user/data/study1/sumstats.gz")),
                   tags$li(strong("Relative paths:"), " Relative to the ", code("GenoPred/pipeline/"), " folder.")
                 ),
                 hr(),
                 h4("Running the Pipeline"),
                 p("1. Unzip the bundle into your project directory:"),
                 pre("unzip genopred_config_YYYYMMDD.zip"),
                 p("2. Run the pipeline from the ", code("GenoPred/pipeline"), " folder, pointing to your config file:"),
                 pre("snakemake --profile slurm --use-conda --configfile=/path/to/your/config.yaml output_all"),
                 hr(),
                 p("For detailed documentation, please visit the ", 
                   a(href="https://opain.github.io/GenoPred/pipeline_readme.html", target="_blank", "GenoPred Website"), ".")
        ),
        
        # --- 1. GWAS LIST ---
        tabPanel("1. GWAS List", icon = icon("dna"),
                 br(),
                 p(class="lead", "Define the GWAS summary statistics to be used by the pipeline."),
                 radioButtons("gwas_mode", "Source:",
                              choices = c("Create new table" = "create", "Use existing file" = "existing", "None (NA)" = "none"),
                              inline = TRUE),
                 hr(),
                 conditionalPanel("input.gwas_mode == 'create'",
                                  rHandsontableOutput("hot_gwas"),
                                  br(),
                                  div(class = "note-box", 
                                      strong("Note:"), " The prevalence and sampling values are used to estimate the SNP-based heritability on the liability scale. Furthermore, prevalence and sampling, or mean and sd values are used to interpret the polygenic scores on the absolute scale."),
                                  
                                  # Column Legend / Help
                                  tags$details(
                                    tags$summary(class = "btn btn-info btn-sm", "Column Guide & File Format"),
                                    br(), br(),
                                    tags$ul(
                                      tags$li(strong("name:"), "ID for the GWAS sumstats. Cannot contain spaces (' ') or hyphens ('-')."),
                                      tags$li(strong("path:"), "File path to the GWAS summary statistics (uncompressed or gzipped)."),
                                      tags$li(strong("population:"), "Reference population (AFR, AMR, EAS, EUR, CSA, MID). If mixed, choose majority."),
                                      tags$li(strong("n:"), "Total sample size. Required if not in sumstats (else NA)."),
                                      tags$li(strong("sampling:"), "Proportion of cases in GWAS sample (binary traits, else NA)."),
                                      tags$li(strong("prevalence:"), "Population prevalence of phenotype (binary traits, else NA)."),
                                      tags$li(strong("mean:"), "Phenotype mean in general population (continuous traits, else NA)."),
                                      tags$li(strong("sd:"), "Phenotype sd in general population (continuous traits, else NA)."),
                                      tags$li(strong("label:"), "Human readable name.")
                                    ),
                                    hr(),
                                    strong("Required GWAS Sumstat Columns (Header):"),
                                    p("The pipeline accepts various header formats using a dictionary. Ensure columns are interpreted correctly by checking sumstat QC logs."),
                                    tags$ul(
                                      tags$li("Must contain: ", code("RSID"), " OR ", code("Chromosome"), " & ", code("Position")),
                                      tags$li("Must contain Effect Size: ", code("BETA"), ", ", code("OR"), ", ", code("log(OR)"), ", or ", code("Z-score")),
                                      tags$li("Must contain: ", code("P-value"), " OR ", code("Standard Error")),
                                      tags$li("Recommended: ", code("N"), " (per variant), ", code("EAF"), " (Effect Allele Freq), ", code("INFO"))
                                    )
                                  )
                 ),
                 conditionalPanel("input.gwas_mode == 'existing'",
                                  textInput("gwas_ext_path", "Path to existing gwas_list file:")
                 )
        ),
        
        # --- 2. TARGET LIST ---
        tabPanel("2. Target List", icon = icon("users"),
                 br(),
                 p(class="lead", "Define the target genotype datasets for polygenic scoring."),
                 radioButtons("target_mode", "Source:",
                              choices = c("Create new table" = "create", "Use existing file" = "existing", "None (NA)" = "none"),
                              inline = TRUE),
                 hr(),
                 conditionalPanel("input.target_mode == 'create'",
                                  rHandsontableOutput("hot_target"),
                                  br(),
                                  div(class = "note-box", 
                                      strong("Note:"), " If the prefix of your target genetic data files does not meet the requirements of GenoPred (see 'path' below), you can create symlinks (like a shortcut) to the original genetic data, and then specify these symlinks in the target_list."),
                                  uiOutput("indiv_report_warning"),
                                  
                                  tags$details(
                                    tags$summary(class = "btn btn-info btn-sm", "Column Guide"),
                                    br(), br(),
                                    tags$ul(
                                      tags$li(strong("name:"), "ID for the target dataset. Cannot contain spaces (' ') or hyphens ('-')."),
                                      tags$li(strong("path:"), "Path to the target genotype data. For type '23andMe', provide full file path. For type 'plink1', 'plink2', 'bgen', and 'vcf', provide the filename PREFIX (e.g. data.chr1-22)."),
                                      tags$li(strong("type:"), "Format of the target genotype dataset.", 
                                              tags$ul(
                                                tags$li("23andMe: Formatted data for an individual."),
                                                tags$li("plink1: Preimputed PLINK1 binary (.bed/.bim/.fam)."),
                                                tags$li("plink2: Preimputed PLINK2 binary (.pgen/.pvar/.psam)."),
                                                tags$li("bgen: Preimputed Oxford format (.bgen/.sample)."),
                                                tags$li("vcf: Preimputed gzipped VCF format (.vcf.gz).")
                                              )
                                      ),
                                      tags$li(strong("indiv_report:"), "Logical indicating whether reports for each individual should be generated. Use with caution if target data contains many individuals."),
                                      tags$li(strong("unrel:"), "Optional path to list of unrelated individuals.")
                                    )
                                  )
                 ),
                 conditionalPanel("input.target_mode == 'existing'",
                                  textInput("target_ext_path", "Path to existing target_list file:")
                 )
        ),
        
        # --- 3. SCORE LIST ---
        tabPanel("3. Score List", icon = icon("list"),
                 br(),
                 p(class="lead", "Define external scoring files (e.g., from PGS Catalog)."),
                 radioButtons("score_mode", "Source:",
                              choices = c("Create new table" = "create", "Use existing file" = "existing", "None (NA)" = "none"),
                              inline = TRUE),
                 hr(),
                 conditionalPanel("input.score_mode == 'create'",
                                  rHandsontableOutput("hot_score"),
                                  br(),
                                  div(class = "note-box", 
                                      strong("Note:"), " Externally derived PGS score files may have a poor variant overlap with the default GenoPred reference data (HapMap3). Score files with <75% of variants present in the reference are excluded from downstream target scoring."),
                                  
                                  tags$details(
                                    tags$summary(class = "btn btn-info btn-sm", "Column Guide & File Format"),
                                    br(), br(),
                                    tags$ul(
                                      tags$li(strong("name:"), "ID for the score (or PGS Catalog ID e.g., PGS000001)."),
                                      tags$li(strong("path:"), "Full path to score file. Leave blank for PGS Catalog auto-download (if name is a PGS ID)."),
                                      tags$li(strong("label:"), "Display name.")
                                    ),
                                    hr(),
                                    strong("Required Score File Columns (Header):"),
                                    p("GenoPred allows only one column of effect sizes per score file. It requires RSIDs or Chromosome/Position."),
                                    tags$ul(
                                      tags$li(code("rsID"), " or ", code("hm_rsID"), " - RSID"),
                                      tags$li(code("chr_name"), " or ", code("hm_chr"), " - Chromosome number"),
                                      tags$li(code("chr_position"), " or ", code("hm_pos"), " - Basepair position"),
                                      tags$li(code("effect_allele"), " - Allele corresponding to effect_weight"),
                                      tags$li(code("other_allele"), " - The other allele"),
                                      tags$li(code("effect_weight"), " - The effect size")
                                    )
                                  )
                 ),
                 conditionalPanel("input.score_mode == 'existing'",
                                  textInput("score_ext_path", "Path to existing score_list file:")
                 )
        ),
        
        # --- 4. GWAS GROUPS ---
        tabPanel("4. GWAS Groups", icon = icon("object-group"),
                 br(),
                 p(class="lead", "Define groups of GWAS for multi-source methods or meta-analysis."),
                 radioButtons("groups_mode", "Source:",
                              choices = c("Create new table" = "create", "Use existing file" = "existing", "None (NA)" = "none"),
                              inline = TRUE),
                 hr(),
                 conditionalPanel("input.groups_mode == 'create'",
                                  rHandsontableOutput("hot_groups"),
                                  br(),
                                  div(class = "note-box", 
                                      strong("Note:"), " Combine multiple GWAS into a single Group ID. 'gwas' should be a comma-separated list of names defined in Tab 1."),
                                  uiOutput("groups_warning"), # Warning for duplicate/invalid names
                                  
                                  tags$details(
                                    tags$summary(class = "btn btn-info btn-sm", "Column Guide"),
                                    br(), br(),
                                    tags$ul(
                                      tags$li(strong("name:"), "Unique ID for this group (distinct from single GWAS names)."),
                                      tags$li(strong("gwas:"), "Comma-separated list of GWAS names (e.g. 'height_ukb,height_bbj')."),
                                      tags$li(strong("label:"), "Display name for the group.")
                                    )
                                  )
                 ),
                 conditionalPanel("input.groups_mode == 'existing'",
                                  textInput("groups_ext_path", "Path to existing gwas_groups file:")
                 )
        ),
        
        # --- 5. BASIC PARAMETERS ---
        tabPanel("5. Basic Parameters", icon = icon("sliders-h"),
                 br(),
                 fluidRow(
                   column(6,
                          wellPanel(
                            h4("1. Output Settings"),
                            textInput("outdir", "Output Directory (required)", placeholder = "e.g. /home/user/project/outputs"),
                            helpText("The directory where outputs of the pipeline will be saved."),
                            hr(),
                            textInput("config_dir", "Configuration Directory", placeholder = "e.g. /home/user/project/inputs"),
                            helpText("The directory where you will unzip these config files (gwas_list.txt etc.)."),
                            helpText("If this is not specified, these configuration files must be in the same folder as the Snakefile (GenoPred/pipeline).")
                          )
                   ),
                   column(6,
                          wellPanel(
                            h4("2. Testing Mode"),
                            radioButtons("testing_mode", NULL,
                                         choices = c("Full Run (Production)" = "NA", 
                                                     "Test Run (Chromosome 22 only)" = "chr22"),
                                         inline = TRUE)
                          ),
                          wellPanel(
                            h4("3. Single-Source PGS Methods"),
                            helpText("Standard methods using one GWAS."),
                            selectInput("pgs_methods_basic", NULL,
                                        choices = pgs_methods_single,
                                        selected = c("sbayesrc"),
                                        multiple = TRUE)
                          )
                   )
                 )
        ),
        
        # --- 6. ADVANCED PARAMETERS ---
        tabPanel("6. Advanced Parameters", icon = icon("cogs"),
                 br(),
                 fluidRow(
                   column(6,
                          wellPanel(
                            h4("1. Resources & Reference"),
                            textInput("resdir", "Resources Directory", value = "NA"),
                            helpText("Set to 'NA' to use default (GenoPred/pipeline/resources/)."),
                            hr(),
                            textInput("refdir", "Alternative Reference Path", placeholder = "/path/to/reference_plink_prefix")
                          ),
                          wellPanel(
                            h4("2. Ancestry Settings"),
                            
                            selectInput("ancestry_adjustment", "Ancestry Adjustment Approach", 
                                        choices = c("Continuous Correction (Recommended)" = "continuous", 
                                                    "Discrete Correction" = "discrete"),
                                        multiple = TRUE,
                                        selected = "continuous"),
                            
                            # Only show threshold if 'discrete' is selected
                            conditionalPanel(
                              condition = "input.ancestry_adjustment && input.ancestry_adjustment.indexOf('discrete') > -1",
                              numericInput("ancestry_threshold", "Ancestry Probability Threshold", value = 0.95, min = 0, max = 1, step = 0.05),
                              helpText("Threshold applies only to Discrete Correction.")
                            )
                          )
                   ),
                   column(6,
                          # Conditional Multi-Source Methods
                          conditionalPanel(
                            condition = "input.groups_mode != 'none'",
                            wellPanel(
                              h4("3. Multi-Source Methods"),
                              div(class = "note-box", strong("Enabled:"), " GWAS Groups are active."),
                              
                              strong("Multi-Source PGS Methods (Jointly Optimised)"),
                              helpText("Methods that use multiple GWAS (e.g., PRS-CSx, X-Wing)."),
                              div(class = "alert alert-warning", style="font-size: 0.9em;",
                                  icon("clock"), "Warning: 'Currently implemented' jointly optimised methods are computationally intensive and slow. We recommend using independently optimised methods (Single-Source) combined with LEOPARD+QuickPRS below."),
                              selectInput("pgs_methods_advanced", NULL,
                                          choices = pgs_methods_multi,
                                          multiple = TRUE),
                              hr(),
                              
                              strong("LEOPARD + QuickPRS Combination"),
                              helpText("Combine PGS from single source methods for GWAS groups using weights derived using LEOPARD + QuickPRS."),
                              helpText(em("Only methods selected in 'Single-Source' (Tab 5) are available here.")),
                              selectInput("leopard_methods", NULL,
                                          choices = NULL, # Populated by server
                                          multiple = TRUE)
                            )
                          ),
                          
                          conditionalPanel(
                            condition = "input.groups_mode == 'none'",
                            div(class = "alert alert-warning", 
                                "Multi-Source PGS methods and Leopard are disabled because 'GWAS Groups' is set to None.")
                          ),
                          
                          wellPanel(
                            h4("4. Additional Parameters"),
                            p("For a full list of parameters, refer to the ", 
                              a(href="https://opain.github.io/GenoPred/pipeline_readme.html#Additional_parameters", target="_blank", "GenoPred Documentation"), "."),
                            helpText("Add any other keys here (YAML format)."),
                            textAreaInput("custom_yaml", NULL, rows = 3, 
                                          placeholder = "# memory_limit: 30000")
                          )
                   )
                 )
        )
      )
    )
  )
)

# --- 5. SERVER LOGIC ---
server <- function(input, output, session) {
  
  values <- reactiveValues()
  values$gwas <- make_template(cols_gwas)
  values$target <- make_template(cols_target)
  values$score <- make_template(cols_score)
  values$groups <- make_template(cols_groups)
  
  # --- RENDER TABLES ---
  
  output$hot_gwas <- renderRHandsontable({
    rhandsontable(values$gwas, rowHeaders = NULL) %>%
      hot_table(minSpareRows = 1) %>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>%
      hot_col("name", placeholder = "e.g. height_ukb") %>%
      hot_col("path", placeholder = "/path/to/sumstats.gz") %>%
      hot_col("population", type = "dropdown", source = c("EUR", "AFR", "AMR", "EAS", "CSA", "MID"), placeholder = "EUR") %>%
      hot_col("n", type = "numeric", placeholder = "10000") %>%
      hot_col("label", placeholder = "\"Height (UKB)\"") %>%
      hot_col("sampling", format = "0.00", placeholder = "0.5") %>%
      hot_col("prevalence", format = "0.00", placeholder = "0.1") %>%
      hot_col("mean", placeholder = "170") %>%
      hot_col("sd", placeholder = "10") %>%
      hot_cols(manualColumnResize = TRUE)
  })
  
  output$hot_groups <- renderRHandsontable({
    rhandsontable(values$groups, rowHeaders = NULL, stretchH = "all") %>%
      hot_table(minSpareRows = 1) %>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>%
      hot_col("name", placeholder = "e.g. height_meta") %>%
      hot_col("gwas", placeholder = "height_ukb,height_bbj") %>%
      hot_col("label", placeholder = "\"Height (Meta-analysis)\"") %>%
      hot_cols(manualColumnResize = TRUE)
  })
  
  output$hot_target <- renderRHandsontable({
    rhandsontable(values$target, rowHeaders = NULL, stretchH = "all") %>%
      hot_table(minSpareRows = 1) %>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>%
      hot_col("name", placeholder = "e.g. target_cohort") %>%
      hot_col("path", placeholder = "/path/to/plink_prefix") %>%
      hot_col("type", type = "dropdown", source = c("plink1", "plink2", "vcf", "bgen", "23andMe"), placeholder = "plink2") %>%
      hot_col("indiv_report", type = "dropdown", source = c("TRUE", "FALSE"), placeholder = "FALSE") %>%
      hot_col("unrel", placeholder = "/path/to/unrelated_ids.txt") %>%
      hot_cols(manualColumnResize = TRUE)
  })
  
  output$hot_score <- renderRHandsontable({
    rhandsontable(values$score, rowHeaders = NULL, stretchH = "all") %>%
      hot_table(minSpareRows = 1) %>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>%
      hot_col("name", placeholder = "PGS000001") %>%
      hot_col("path", placeholder = "/path/to/score.txt") %>%
      hot_col("label", placeholder = "\"Height (PGS Catalog)\"") %>%
      hot_cols(manualColumnResize = TRUE)
  })
  
  # --- DYNAMIC LEOPARD CHOICES ---
  observe({
    req(input$pgs_methods_basic)
    selected_methods <- pgs_methods_single[pgs_methods_single %in% input$pgs_methods_basic]
    updateSelectInput(session, "leopard_methods",
                      choices = selected_methods,
                      selected = input$leopard_methods)
  })
  
  # --- INDIVIDUAL REPORT WARNING ---
  output$indiv_report_warning <- renderUI({
    req(input$hot_target)
    target_df <- hot_to_r(input$hot_target)
    
    # Check if ANY row has indiv_report set to TRUE
    # Use na.rm = TRUE to handle the empty/NA value in the spare row
    if (any(target_df$indiv_report == "TRUE", na.rm = TRUE)) {
      div(class = "alert alert-warning", icon("exclamation-triangle"), 
          "Warning: Setting 'indiv_report' to TRUE generates reports for EVERY participant (slow for large N).")
    } else {
      NULL
    }
  })
  
  # --- GROUPS VALIDATION ---
  output$groups_warning <- renderUI({
    req(input$hot_gwas, input$hot_groups)
    gwas_df <- hot_to_r(input$hot_gwas)
    groups_df <- hot_to_r(input$hot_groups)
    
    gwas_df <- gwas_df[gwas_df$name != "", ]
    groups_df <- groups_df[groups_df$name != "", ]
    
    msg <- list()
    
    # 1. Check for Duplicate Group Names vs GWAS Names
    conflict <- intersect(gwas_df$name, groups_df$name)
    if(length(conflict) > 0) {
      msg[[length(msg)+1]] <- paste0("Error: Group Name(s) '", paste(conflict, collapse=", "), "' overlap with names in GWAS List.")
    }
    
    # 2. Check each group
    if(nrow(groups_df) > 0) {
      for(i in 1:nrow(groups_df)) {
        raw_list <- unlist(strsplit(groups_df$gwas[i], ","))
        gwas_in_group <- trimws(raw_list)
        
        # 2a. Check for duplicates in the comma-separated list itself
        if(anyDuplicated(gwas_in_group)) {
          msg[[length(msg)+1]] <- paste0("Error in Group '", groups_df$name[i], "': Contains duplicate GWAS names.")
        }
        
        # 2b. Check if GWAS exists in list
        missing <- gwas_in_group[!gwas_in_group %in% gwas_df$name]
        if(length(missing) > 0) {
          msg[[length(msg)+1]] <- paste0("Error in Group '", groups_df$name[i], "': GWAS ID(s) '", paste(missing, collapse=", "), "' not found in GWAS List.")
        }
        
        # 2c. Check Populations (must be distinct)
        if(nrow(gwas_df) > 0) {
          pops <- gwas_df$population[match(gwas_in_group, gwas_df$name)]
          pops <- pops[!is.na(pops)] # remove missing matches
          if(anyDuplicated(pops)) {
            msg[[length(msg)+1]] <- paste0("Error in Group '", groups_df$name[i], "': Contains multiple GWAS from the same population (", paste(unique(pops[duplicated(pops)]), collapse=","), "). Multi-source groups typically require distinct populations.")
          }
        }
      }
    }
    
    if(length(msg) > 0) {
      div(class = "alert alert-danger", icon("exclamation-circle"), HTML(paste(msg, collapse="<br>")))
    } else {
      NULL
    }
  })
  
  # --- STATUS PANEL ---
  has_data <- function(hot_input) {
    if(is.null(hot_input)) return(FALSE)
    df <- hot_to_r(hot_input)
    any(df[[1]] != "" & !is.na(df[[1]]))
  }
  
  # Status helper to handle file path inputs
  has_path <- function(path_input) {
    !is.null(path_input) && path_input != ""
  }
  
  output$status_panel <- renderUI({
    
    # Check GWAS
    if(input$gwas_mode == 'create') {
      gwas_stat <- if(has_data(input$hot_gwas)) tags$li(class="text-success", icon("check"), "GWAS List: Active") else tags$li(class="text-muted", "GWAS List: Empty")
    } else if (input$gwas_mode == 'existing') {
      gwas_stat <- if(has_path(input$gwas_ext_path)) tags$li(class="text-success", icon("check"), "GWAS List: File Linked") else tags$li(class="text-muted", "GWAS List: Empty")
    } else {
      gwas_stat <- tags$li(class="text-muted", "GWAS List: None")
    }
    
    # Check Target
    if(input$target_mode == 'create') {
      target_stat <- if(has_data(input$hot_target)) tags$li(class="text-success", icon("check"), "Target List: Active") else tags$li(class="text-muted", "Target List: Empty")
    } else if (input$target_mode == 'existing') {
      target_stat <- if(has_path(input$target_ext_path)) tags$li(class="text-success", icon("check"), "Target List: File Linked") else tags$li(class="text-muted", "Target List: Empty")
    } else {
      target_stat <- tags$li(class="text-muted", "Target List: None")
    }
    
    # Check Score
    if(input$score_mode == 'create') {
      score_stat <- if(has_data(input$hot_score)) tags$li(class="text-success", icon("check"), "Score List: Active") else tags$li(class="text-muted", "Score List: Empty")
    } else if (input$score_mode == 'existing') {
      score_stat <- if(has_path(input$score_ext_path)) tags$li(class="text-success", icon("check"), "Score List: File Linked") else tags$li(class="text-muted", "Score List: Empty")
    } else {
      score_stat <- tags$li(class="text-muted", "Score List: None")
    }
    
    # Check Groups
    if(input$groups_mode == 'create') {
      groups_stat <- if(has_data(input$hot_groups)) tags$li(class="text-success", icon("check"), "GWAS Groups: Active") else tags$li(class="text-muted", "GWAS Groups: Empty")
    } else if (input$groups_mode == 'existing') {
      groups_stat <- if(has_path(input$groups_ext_path)) tags$li(class="text-success", icon("check"), "GWAS Groups: File Linked") else tags$li(class="text-muted", "GWAS Groups: Empty")
    } else {
      groups_stat <- tags$li(class="text-muted", "GWAS Groups: None")
    }
    
    tags$ul(class="list-unstyled", gwas_stat, target_stat, score_stat, groups_stat)
  })
  
  # --- PROCESS DATA FUNCTION ---
  process_data <- function(hot_input, type="standard") {
    if(is.null(hot_input)) return(NULL)
    df <- hot_to_r(hot_input)
    key_col <- "name" 
    
    df <- df[!is.na(df[[key_col]]) & df[[key_col]] != "", ]
    if(nrow(df) == 0) return(NULL)
    
    if("label" %in% colnames(df)) {
      clean_labels <- gsub('^"|"$', '', as.character(df$label))
      df$label <- paste0('"', clean_labels, '"')
    }
    
    if(type == "groups") return(df)
    
    # Path Logic for Score (Auto-download if empty path)
    if(type == "score" && "path" %in% colnames(df)) {
      df$path <- ifelse(df$path == "" | is.na(df$path), NA, df$path)
    }
    
    return(df)
  }
  
  # --- DOWNLOAD HANDLER ---
  output$download_bundle <- downloadHandler(
    filename = function() { paste0("genopred_config_", format(Sys.Date(), "%Y%m%d"), ".zip") },
    content = function(file) {
      
      if (trimws(input$outdir) == "") {
        showNotification("Error: Output Directory must be specified!", type = "error")
        stop("Output Directory is mandatory.")
      }
      
      tmpdir <- tempdir()
      setwd(tmpdir)
      files <- c()
      
      combined_methods <- input$pgs_methods_basic
      if(input$groups_mode != "none" && !is.null(input$pgs_methods_advanced)) {
        combined_methods <- c(combined_methods, input$pgs_methods_advanced)
      }
      
      config <- list(
        outdir = input$outdir,
        config_file = paste0(input$config_dir,'/config.yaml'),
        resdir = if(input$resdir == "") "NA" else input$resdir, 
        pgs_methods = combined_methods,
        testing = if(input$testing_mode == "NA") "NA" else input$testing_mode,
        ancestry_threshold = input$ancestry_threshold,
        ancestry_adjustment = input$ancestry_adjustment
      )
      
      if(input$refdir != "") {
        config$refdir <- input$refdir
      }
      
      if (!is.null(input$leopard_methods) && input$groups_mode != "none") {
        config$leopard_methods <- input$leopard_methods
      }
      
      # Process Tables
      if(input$gwas_mode == "create") {
        df <- process_data(input$hot_gwas)
        if(!is.null(df)) {
          write.table(df, "gwas_list.txt", sep=" ", quote=FALSE, row.names=FALSE, na="NA")
          files <- c(files, "gwas_list.txt")
          config$gwas_list <- "gwas_list.txt"
        } else { config$gwas_list <- "NA" }
      } else if(input$gwas_mode == "existing") {
        config$gwas_list <- input$gwas_ext_path
      } else { config$gwas_list <- "NA" }
      
      if(input$groups_mode == "create") {
        df <- process_data(input$hot_groups, type="groups")
        if(!is.null(df)) {
          write.table(df, "gwas_groups.txt", sep=" ", quote=FALSE, row.names=FALSE, na="NA")
          files <- c(files, "gwas_groups.txt")
          config$gwas_groups <- "gwas_groups.txt"
        } else { config$gwas_groups <- "NA" }
      } else if(input$groups_mode == "existing") {
        config$gwas_groups <- input$groups_ext_path
      } else { config$gwas_groups <- "NA" }
      
      if(input$target_mode == "create") {
        df <- process_data(input$hot_target)
        if(!is.null(df)) {
          write.table(df, "target_list.txt", sep=" ", quote=FALSE, row.names=FALSE, na="NA")
          files <- c(files, "target_list.txt")
          config$target_list <- "target_list.txt"
        } else { config$target_list <- "NA" }
      } else if(input$target_mode == "existing") {
        config$target_list <- input$target_ext_path
      } else { config$target_list <- "NA" }
      
      if(input$score_mode == "create") {
        df <- process_data(input$hot_score, type="score")
        if(!is.null(df)) {
          write.table(df, "score_list.txt", sep=" ", quote=FALSE, row.names=FALSE, na="NA")
          files <- c(files, "score_list.txt")
          config$score_list <- "score_list.txt"
        } else { config$score_list <- "NA" }
      } else if(input$score_mode == "existing") {
        config$score_list <- input$score_ext_path
      } else { config$score_list <- "NA" }
      
      # Write Config
      format_list_string <- function(items) {
        if(is.null(items) || length(items) == 0) return("[]")
        quoted_items <- paste0("'", items, "'")
        paste0("[", paste(quoted_items, collapse = ", "), "]")
      }
      
      # Helper to prepend config directory if file was generated by app
      apply_config_dir <- function(filename) {
        if (filename == "NA") return("NA")
        
        # Only modify the specific files generated by this app
        generated_files <- c("gwas_list.txt", "gwas_groups.txt", "target_list.txt", "score_list.txt")
        
        if (filename %in% generated_files) {
          # Get dir from input
          dir_path <- input$config_dir
          
          # If user provided a path
          if (!is.null(dir_path) && dir_path != "") {
            # Ensure trailing slash
            if (substr(dir_path, nchar(dir_path), nchar(dir_path)) != "/") {
              dir_path <- paste0(dir_path, "/")
            }
            return(paste0(dir_path, filename))
          }
        }
        
        # Return original (either no config_dir set, or it's an existing absolute path provided by user)
        return(filename)
      }
      
      # --- Build YAML Content Conditionally ---
      # Start with mandatory parameters
      yaml_lines <- c(
        paste0("outdir: ", config$outdir)
      )
      
      yaml_lines <- c(
        c(yaml_lines, paste0("config_file: ", config$config_file))
      )
      
      # Optional: Resources Directory (Default NA)
      if(config$resdir != "NA" && config$resdir != "") {
        yaml_lines <- c(yaml_lines, paste0("resdir: ", config$resdir))
      }
      
      # PGS Methods (Always needed if pipeline runs)
      yaml_lines <- c(yaml_lines, paste0("pgs_methods: ", format_list_string(config$pgs_methods)))
      
      # Optional: Leopard Methods
      if(!is.null(config$leopard_methods)) {
        yaml_lines <- c(yaml_lines, paste0("leopard_methods: ", format_list_string(config$leopard_methods)))
      }
      
      # Optional: Testing Mode (Default NA)
      if(config$testing != "NA") {
        yaml_lines <- c(yaml_lines, paste0("testing: ", config$testing))
      }
      
      # Optional: File Lists (Default NA)
      g_list <- apply_config_dir(config$gwas_list)
      if(g_list != "NA") yaml_lines <- c(yaml_lines, paste0("gwas_list: ", g_list))
      
      g_groups <- apply_config_dir(config$gwas_groups)
      if(g_groups != "NA") yaml_lines <- c(yaml_lines, paste0("gwas_groups: ", g_groups))
      
      t_list <- apply_config_dir(config$target_list)
      if(t_list != "NA") yaml_lines <- c(yaml_lines, paste0("target_list: ", t_list))
      
      s_list <- apply_config_dir(config$score_list)
      if(s_list != "NA") yaml_lines <- c(yaml_lines, paste0("score_list: ", s_list))
      
      # Optional: Ancestry Threshold (Default 0.95)
      if(config$ancestry_threshold != 0.95) {
        yaml_lines <- c(yaml_lines, paste0("ancestry_threshold: ", config$ancestry_threshold))
      }
      
      # Optional: Ancestry Adjustment (Default ['continuous'])
      # Check if it differs from default single value "continuous"
      is_default_adj <- length(config$ancestry_adjustment) == 1 && config$ancestry_adjustment[1] == "continuous"
      if(!is_default_adj) {
        yaml_lines <- c(yaml_lines, paste0("ancestry_adjustment: ", format_list_string(config$ancestry_adjustment)))
      }
      
      # Optional: Reference Data
      if(!is.null(config$refdir)) {
        yaml_lines <- c(yaml_lines, paste0("refdir: ", config$refdir))
      }
      
      writeLines(yaml_lines, "config.yaml")
      
      if(trimws(input$custom_yaml) != "") {
        cat("\n# --- Custom Parameters ---\n", file="config.yaml", append=TRUE)
        cat(input$custom_yaml, file="config.yaml", append=TRUE)
      }
      files <- c(files, "config.yaml")
      
      zip(zipfile = file, files = files)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)