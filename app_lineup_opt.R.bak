# ============================================================================
#  app_lineup_opt.R - OFFENSE: best batting order for our 9 vs a given pitcher.
#  Standalone app; shares data/{hitter_object,pitcher_object,run_values}.rds.
#  run: shiny::runApp("app_lineup_opt.R")
# ============================================================================
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(tidyverse); library(reactable); library(plotly); library(gtools)
})
source("continuous_common.R")
source("matchup_kernel_continuous.R")
source("lineup_sim_continuous.R")
source("lineup_optimizer.R")
SHAPE_OK <- tryCatch({ source("build_batter_response.R"); source("matchup_shape.R")
  source("shape_pa.R"); TRUE }, error=function(e){ message("shape model unavailable: ", e$message); FALSE })
km <- load_cont_kernel("data")
CURRENT_SEASON <- "2026"

PITCHERS <- km$po$buckets %>% filter(season==CURRENT_SEASON) %>%
  distinct(PitcherId, Pitcher, PitcherTeam) %>%
  filter(!is.na(PitcherId), !is.na(Pitcher), Pitcher!="") %>%
  mutate(label=paste0(Pitcher," (",PitcherTeam,")")) %>% arrange(Pitcher)
pid_of <- function(l) PITCHERS$PitcherId[match(l, PITCHERS$label)]
TEAMS  <- sort(unique(km$ho$meta$BatterTeam[km$ho$meta$season==CURRENT_SEASON]))

navy <- bs_theme(version=5, bg="#0F172A", fg="#F1F5F9", primary="#2563EB",
  base_font=font_google("Inter"), heading_font=font_google("Inter"))
rt <- reactableTheme(backgroundColor="#1E293B", color="#F1F5F9", borderColor="#334155",
  stripedColor="#273448", highlightColor="#2D3F58",
  headerStyle=list(backgroundColor="#0F172A", color="#94A3B8"))

ui <- page_navbar(title="NECBL Lineup Optimizer", theme=navy,
  nav_panel("Best Batting Order", icon=icon("bolt"),
    layout_sidebar(sidebar=sidebar(width=360,
      selectInput("team","Team", choices=TEAMS,
                  selected=if("North Shore Navigators"%in%TEAMS)"North Shore Navigators" else TEAMS[1]),
      uiOutput("hitters_ui"),
      selectInput("opp","Opposing pitcher", choices=PITCHERS$label),
      hr(),
      selectizeInput("battle","Position battle (optional): competing hitters",
                     choices=NULL, multiple=TRUE),
      numericInput("battle_max","How many of them play", value=1, min=0, step=1),
      p("Enter 10+ hitters to let it pick the best 9. Flag guys competing for a spot and cap how many play.",
        class="text-muted small"),
      hr(),
      uiOutput("lock_ui"),
      p("Lock a hitter to a slot; the rest optimize around him.", class="text-muted small"),
      hr(),
      sliderInput("innings","Innings", 1, 9, 9),
      checkboxInput("shape", "Shape-aware matchups (pitch shape + location)", value = FALSE),
      sliderInput("sims","MC sims (for clicked lineup's chart)", 500, 3000, 1000, step=500),
      actionButton("go","Optimize", class="btn-primary w-100")),
    tagList(
      layout_columns(col_widths=c(7,5),
        card(card_header("Top orders \u2014 click one for the full breakdown"),
             reactableOutput("tbl")),
        card(card_header("Hitter pool"), reactableOutput("pool_tbl"))),
      uiOutput("detail_card")))))

server <- function(input, output, session){
  output$hitters_ui <- renderUI({
    req(input$team)
    pool <- km$ho$meta %>% filter(season==CURRENT_SEASON, BatterTeam==input$team) %>%
      arrange(desc(pa)) %>% distinct(Batter) %>% pull(Batter)
    selectizeInput("batters","Our hitters (add 9+; top-9 auto-filled, editable)",
      choices=pool, selected=head(pool,9), multiple=TRUE)
  })
  # keep the battle dropdown limited to the hitters currently in the pool
  observeEvent(input$batters, {
    updateSelectizeInput(session, "battle", choices=input$batters,
                         selected=intersect(input$battle, input$batters), server=FALSE)
  })
  output$pool_tbl <- renderReactable({
    req(input$batters)
    d <- resolve_lineup_c(km, input$batters, CURRENT_SEASON) %>%
      mutate(N=row_number()) %>% select(N, Batter, Side=BatterSide)
    reactable(d, theme=rt, defaultPageSize=12)
  })

  output$lock_ui <- renderUI({
    req(input$batters)
    tagList(
      selectizeInput("lock1_b","Lock a hitter", choices=c("", input$batters), selected=""),
      selectInput("lock1_s","to slot", choices=1:9, selected=1),
      selectizeInput("lock2_b","Lock a hitter", choices=c("", input$batters), selected=""),
      selectInput("lock2_s","to slot", choices=1:9, selected=2)
    )
  })

  res <- eventReactive(input$go, {
    req(length(input$batters)>=9)
    options(necbl.shape_mode = isTRUE(input$shape) && SHAPE_OK)
    ctx <- list(pid=pid_of(input$opp), season=CURRENT_SEASON, innings=input$innings, sims=input$sims)
    battle <- if (length(input$battle)>=1) input$battle else NULL
    bmax   <- if (!is.null(battle)) as.integer(input$battle_max) else NULL
    locked <- list()
    if (!is.null(input$lock1_b) && nzchar(input$lock1_b)) locked[[as.character(input$lock1_s)]] <- input$lock1_b
    if (!is.null(input$lock2_b) && nzchar(input$lock2_b)) locked[[as.character(input$lock2_s)]] <- input$lock2_b
    orders <- optimize_lineup(km, input$batters, ctx$pid, ctx$season,
                              innings=ctx$innings, n_sims=ctx$sims, method="local",
                              battle=battle, battle_max=bmax,
                              locked=if(length(locked)) locked else NULL)
    xrv_of <- function(ord){ cache <- build_pa_cache(km, ord, ctx$pid, ctx$season)
      order_xrv_c(cache, ord, ctx$innings) }
    tbl <- imap_dfr(orders, function(o,i) tibble(Rank=i,
      Order=paste(o$order, collapse=" \u2192 "),
      `Mean runs`=round(o$mean_runs,2),
      xRV=round(xrv_of(o$order),2),
      P10=round(o$p10,1), P90=round(o$p90,1), `5+ runs %`=round(100*o$big,1)))
    list(orders=orders, tbl=tbl, ctx=ctx)
  })
  output$tbl <- renderReactable({
    o <- res()
    reactable(o$tbl, theme=rt, striped=TRUE, highlight=TRUE,
              selection="single", onClick="select", defaultPageSize=10,
              columns=list(
                Rank=colDef(maxWidth=55),
                Order=colDef(name="Lineup (best first)", minWidth=300),
                `Mean runs`=colDef(maxWidth=95),
                xRV=colDef(maxWidth=80)))
  })

  sel <- reactive({ getReactableState("tbl","selected") })
  detail <- reactive({ o <- res(); i <- sel(); if (is.null(o) || is.null(i)) return(NULL); o$orders[[i]] })
  output$detail_card <- renderUI({
    if (is.null(detail())) return(NULL)
    card(card_header("Selected lineup"),
         uiOutput("detail_head"), reactableOutput("detail_tbl"),
         plotlyOutput("detail_dist", height="300px"))
  })
  output$detail_head <- renderUI({
    d <- detail(); if (is.null(d)) return(NULL)
    tags$div(
      tags$table(class="table table-dark table-sm", style="margin-bottom:8px;",
        tags$tr(tags$td("Mean runs"), tags$td(strong(round(d$mean_runs,2)))),
        tags$tr(tags$td("P10 / P90"), tags$td(paste0(round(d$p10,1)," / ",round(d$p90,1)))),
        tags$tr(tags$td("5+ runs %"), tags$td(paste0(round(100*d$big,1),"%")))),
      tags$div(style="font-size:15px; line-height:1.7;",
        lapply(seq_along(d$order), function(k)
          tags$div(tags$b(paste0(k, ". ")), d$order[k]))))
  })
  output$detail_tbl <- renderReactable({
    d <- detail(); if (is.null(d)) return(NULL)
    tb <- d$rows %>% filter(tto==1) %>% transmute(Slot=slot, Batter, xRV=round(xrv,3),
      OB=map_dbl(line,~sum(.x[c("BB","HBP","single","double","triple","home_run")])),
      HR=map_dbl(line,~.x["home_run"])) %>%
      mutate(`OB%`=round(100*OB,1), `HR%`=round(100*HR,1)) %>% select(Slot,Batter,xRV,`OB%`,`HR%`)
    reactable(tb, theme=rt, defaultPageSize=9)
  })
  output$detail_dist <- renderPlotly({
    d <- detail(); if (is.null(d)) return(NULL)
    plot_ly(x=d$runs, type="histogram", marker=list(color="#16A34A")) %>%
      layout(paper_bgcolor="#1E293B", plot_bgcolor="#1E293B", font=list(color="#F1F5F9"),
             xaxis=list(title="runs scored", color="#94A3B8"),
             yaxis=list(title="", color="#94A3B8"), bargap=0.05, margin=list(t=6,b=30))
  })
}
shinyApp(ui, server)
