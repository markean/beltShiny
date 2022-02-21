library(shiny)
library(melt)
library(ggplot2)
beta <- c("1 / m", "m / n", "1")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Near Bayes"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("theta", "True mean",
                  value = 0, min = ,-3, max = 3, step = 0.05),
      sliderInput("sigma", "True sd", value = 1, min = 0.5, max = 3),
      sliderInput("theta0", "Prior mean",
                  value = -1, min = ,-3, max = 3, step = 0.05),
      sliderInput("sigma0", "Prior sd", value = 1, min = 0.5, max = 3),
      sliderInput("grid", "Grid range",
                  value = c(-1, 1), min = -4, max = 4, step = 0.05),
      numericInput("seed", "Seed",
                   value = 1, min = 1, max = .Machine$integer.max),
      numericInput("n", "n: sample size (10 <= n <= 1000)", value = 100,
                   min = 10, max = 1000),
      numericInput("m", "m: pseudo sample size (2 <= m <= 1000)",
                   value = 2, min = 10, max = 1000),
      selectInput("beta", "beta (fractional weight of one pseudo observation)",
                  beta),
      sliderInput("plotlim", "Plot x-axis range",
                  value = c(-1, 1), min = -4, max = 4, round = F, step = 0.05),
      checkboxGroupInput("wel", "Weighted EL",
                         choices = c("wel1" = "wel1", "wel2" = "wel2",
                                     "wel3" = "wel3", "wel4" = "wel4"),
                         selected = c("wel1"))),
    mainPanel(textOutput("summary"),
              verbatimTextOutput("code"),
              plotOutput("plot1"),
              plotOutput("plot2"))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(warn = -1)
  # parameters 1
  theta <- reactive(input$theta)
  sigma <- reactive(input$sigma)
  theta0 <- reactive(input$theta0)
  sigma0 <- reactive(input$sigma0)
  # parameters 2
  n <- reactive(input$n)
  m <- reactive(input$m)
  grid <- reactive(input$grid)
  plotlim <- reactive(input$plotlim)
  # parameters 3
  seed <- reactive(input$seed)
  beta <- reactive(input$beta)
  wel <- reactive(input$wel)

  # summary statistics
  output$summary <- renderText("Summary of sample data")
  output$code <- renderPrint({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    summary(x)
    }
  )
  # posterior density plot (discrete, normalized)
  output$plot1 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # generate grid
    theta_grid <- seq(grid()[1], grid()[2], length.out = 1000)
    # posterior mean
    posterior_mean <-
      (theta0() / sigma0()^2 + sum(x) / sigma()^2) /
      (1 / sigma0()^2 + n() / sigma()^2)
    # posterior sd
    posterior_sd <- (1 / sigma0()^2 + n() / sigma()^2)^(-1 / 2)

    # weights
    if (beta() == "1 / m") {
      beta <- 1 / m()
    } else if (beta() == "m / n") {
      beta <- m() / n()
    } else {
      beta <- 1
    }
    w <- c(rep(1, n()), rep(beta, m()))

    true_log_pd <- function(par) {
      dnorm(par, mean = posterior_mean, sd = posterior_sd, log = T)
    }
    EL_logLR <- function(par) {
      z <- el_mean(par, x, rep(1, n()),
                   list(threshold = 1e+10, maxit = 100, abstol = 1e-06))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }

    pseudo_sample <- function(par) {
      pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
      while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
        pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
      }
      pseudo_x
    }
    # logLR
    logLR <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                       abstol = 1e-06))$optim$log.prob[1:n()])
    }
    logWLR <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                       abstol = 1e-06))$optim$log.wprob[1:n()])
    }
    # logLR with augmented data
    logLR_aug <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                   abstol = 1e-06))$optim$logLR
    }
    logWLR_aug <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                   abstol = 1e-06))$optim$logWLR
    }
    # true posterior density (discrete, normalized)
    true_log_pd <- sapply(theta_grid, true_log_pd)
    true_pd <- exp(true_log_pd) / sum(exp(true_log_pd))
    # EL posterior density (discrete, normalized)
    EL_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
      sapply(theta_grid, EL_logLR)
    EL_pd <- exp(EL_log_pd) / sum(exp(EL_log_pd), na.rm = T)
    df_true <- data.frame(x = theta_grid, y = true_pd, type = "true")
    df_EL <- data.frame(x = theta_grid, y = EL_pd, type = "EL")

    # posterior density with real data (discrete, normalized)
    if ("wel1" %in% wel()) {
      set.seed(seed())
      WEL1_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logLR)
      WEL1_pd <- exp(WEL1_log_pd) / sum(exp(WEL1_log_pd))
      df_WEL1 <- data.frame(x = theta_grid, y = WEL1_pd, type = "WEL1")
    } else {
      df_WEL1 <- NULL
    }
    if ("wel2" %in% wel()) {
      set.seed(seed())
      WEL2_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logWLR)
      WEL2_pd <- exp(WEL2_log_pd) / sum(exp(WEL2_log_pd))
      df_WEL2 <- data.frame(x = theta_grid, y = WEL2_pd, type = "WEL2")
    } else {
      df_WEL2 <- NULL
    }
    # posterior density with augmented data (discrete, normalized)
    if ("wel3" %in% wel()) {
      set.seed(seed())
      WEL3_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logLR_aug)
      WEL3_pd <- exp(WEL3_log_pd) / sum(exp(WEL3_log_pd))
      df_WEL3 <- data.frame(x = theta_grid, y = WEL3_pd, type = "WEL3")
    } else {
      df_WEL3 <- NULL
    }
    if ("wel4" %in% wel()) {
      set.seed(seed())
      WEL4_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logWLR_aug)
      WEL4_pd <- exp(WEL4_log_pd) / sum(exp(WEL4_log_pd))
      df_WEL4 <- data.frame(x = theta_grid, y = WEL4_pd, type = "WEL4")
    } else {
      df_WEL4 <- NULL
    }

    ggplot(rbind(df_true, df_EL, df_WEL1, df_WEL2, df_WEL3, df_WEL4),
           aes(x, y, color = type)) +
      geom_path(alpha = 0.5, na.rm = TRUE) +
      labs(x = "theta", y = "posterior density (normalized)",
           title = "Posterior Densities") +
      geom_vline(xintercept =  c(min(x), max(x)), linetype = "dashed",
                 alpha = 0.5) +
      xlim(plotlim()[1], plotlim()[2])
  })

  # LR plot (discrete, normalized)
  output$plot2 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # generate grid
    theta_grid <- seq(grid()[1], grid()[2], length.out = 1000)
    # weights
    if (beta() == "1 / m") {
      beta <- 1 / m()
    } else if (beta() == "m / n") {
      beta <- m() / n()
    } else {
      beta <- 1
    }
    w <- c(rep(1, n()), rep(beta, m()))

    EL_logLR <- function(par) {
      z <- el_mean(par, x, rep(1, n()),
                   list(threshold = 1e+10, maxit = 100, abstol = 1e-06))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }
    pseudo_sample <- function(par) {
      pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
      while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
        pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
      }
      pseudo_x
    }
    # logLR
    logLR <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                       abstol = 1e-06))$optim$log.prob[1:n()])
    }
    logWLR <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                       abstol = 1e-06))$optim$log.wprob[1:n()])
    }
    # logLR with augmented data
    logLR_aug <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                   abstol = 1e-06))$optim$logLR
    }
    logWLR_aug <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                   abstol = 1e-06))$optim$logWLR
    }

    # EL density (discrete, normalized)
    EL_log_d <- sapply(theta_grid, EL_logLR)
    EL_d <- exp(EL_log_d) / sum(exp(EL_log_d), na.rm = T)
    df_EL <- data.frame(x = theta_grid, y = EL_d, type = "EL")
    # WEL1 density
    if ("wel1" %in% wel()) {
      set.seed(seed())
      WEL1_log_d <- sapply(theta_grid, logLR)
      WEL1_d <- exp(WEL1_log_d) / sum(exp(WEL1_log_d))
      df_WEL1 <- data.frame(x = theta_grid, y = WEL1_d, type = "WEL1")
    } else {
      df_WEL1 <- NULL
    }
    # WEL2 density
    if ("wel2" %in% wel()) {
      set.seed(seed())
      WEL2_log_d <- sapply(theta_grid, logWLR)
      WEL2_d <- exp(WEL2_log_d) / sum(exp(WEL2_log_d))
      df_WEL2 <- data.frame(x = theta_grid, y = WEL2_d, type = "WEL2")
    } else {
      df_WEL2 <- NULL
    }
    # WEL3 density
    if ("wel3" %in% wel()) {
      set.seed(seed())
      WEL3_log_d <- sapply(theta_grid, logLR_aug)
      WEL3_d <- exp(WEL3_log_d) / sum(exp(WEL3_log_d))
      df_WEL3 <- data.frame(x = theta_grid, y = WEL3_d, type = "WEL3")
    } else {
      df_WEL3 <- NULL
    }
    # WEL3 density
    if ("wel4" %in% wel()) {
      set.seed(seed())
      WEL4_log_d <- sapply(theta_grid, logWLR_aug)
      WEL4_d <- exp(WEL4_log_d) / sum(exp(WEL4_log_d))
      df_WEL4 <- data.frame(x = theta_grid, y = WEL4_d, type = "WEL4")
    } else {
      df_WEL4 <- NULL
    }

    # density plot (discrete, normalized)
    ggplot(rbind(df_EL, df_WEL1, df_WEL2, df_WEL3, df_WEL4),
           aes(x, y, color = type)) +
      geom_path(alpha = 0.5, na.rm = TRUE) +
      labs(x = "theta", y = "density (normalized)",
           title = "Densities") +
      geom_vline(xintercept =  c(min(x), max(x)), linetype = "dashed",
                 alpha = 0.5) +
      xlim(plotlim()[1], plotlim()[2])
  })
}

# Run the application
shinyApp(ui = ui, server = server)
