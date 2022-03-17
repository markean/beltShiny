library(shiny)
library(tidyverse)
# library(nloptr)
library(melt)
library(ggplot2)
library(ggsci)
ui <- fluidPage(
  titlePanel("Bayesian Weighted EL with Fractional Pseudo Observations (Normal-Normal Model)"),
  tabsetPanel(
    tabPanel("LR and Posterior Density", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 sliderInput("theta", "sampling mean",
                             value = 0, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma", "sampling sd",
                             value = 1, min = 0.5, max = 3),
                 sliderInput("theta0", "prior mean",
                             value = -1, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma0", "prior sd",
                             value = 1, min = 0.5, max = 3),
                 sliderInput("grid", "grid range",
                             value = c(-1, 3), min = -5, max = 5, step = 0.1),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 100,
                             min = 10, max = 200, step = 10),
                 sliderInput("m", "m: pseudo sample size",
                             value = 2, min = 2, max = 100, step = 1),
                 sliderInput("strength", "c: strength of pseudo sample",
                             value = 1, min = 1, max = 5000),
                 checkboxGroupInput("type", "weighted EL to use",
                                    choices = c("WEL (w/o pseudo obs. & discard weights)" = "wel",
                                                "WETEL (w/o pseudo obs. & discard weights)" = "et"),
                                    selected = c("et"),
                                    inline = FALSE)),

               mainPanel(
                 textOutput("summary1"),
                 verbatimTextOutput("code1"),
                 plotOutput("plot1"),
                 plotOutput("plot2")
               )
             )
    )
  )
)

server <- function(input, output) {
  options(warn = -1)
  # parameters for densities
  theta <- reactive(input$theta)
  sigma <- reactive(input$sigma)
  theta0 <- reactive(input$theta0)
  sigma0 <- reactive(input$sigma0)
  n <- reactive(input$n)
  m <- reactive(input$m)
  grid <- reactive(input$grid)
  seed <- reactive(input$seed)
  strength <- reactive(input$strength)
  type <- reactive(input$type)

  # summary statistics
  output$summary1 <- renderText("Summary of sample data")
  output$code1 <- renderPrint({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    summary(x)
  }
  )

  # LR plot
  output$plot1 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # generate grid
    theta_grid <- seq(grid()[1], grid()[2], length.out = 500)
    # weights
    # beta <- 1 / m()
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w

    EL_logLR <- function(par) {
      z <- el_mean(par, x, control = list(abstol = 1e-04,
                                          threshold = 1e+10))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }
    logWLR <- function(par) {
      x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
      pp <- el_mean(par, c(x, x_pseudo), w, control =
                      list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
      p1 <- pp[seq_len(n())]
      p2 <- pp[-seq_len(n())]
      sum(p1) + n() * (log(n() + strength()) - log(n() + m()))
    }
    opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-08)
    logET <- function(par) {
      x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
      g <- c(x, x_pseudo) - par

      # eval_f <- function(l) mean(w * exp(l * g))
      # eval_grad_f <- function(l) mean(w * exp(l * g) * g)
      # lambda <- nloptr(x0 = 0, eval_f = eval_f, eval_grad_f = eval_grad_f,
      #                  opts = opts)$solution

      lambda_finder <- function(l) mean(w * exp(l * g) * g)
      lambda <- uniroot(lambda_finder, extendInt = "yes",
                        lower = -1e+10, upper = 1e+10)$root

      unnormalized_prob <- w * exp(lambda * g)
      log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
      sum(log_prob[seq_len(n())]) + n() * log(n() + strength())
    }

    # EL logLR
    EL_log_d <- sapply(theta_grid, EL_logLR)
    df_EL <- data.frame(x = theta_grid, y = EL_log_d, type = "EL")
    dt <- df_EL
    # WEL logLR
    if ("wel" %in% type()) {
      # set.seed(seed())
      WEL_logLR <- sapply(theta_grid, logWLR)
      df_WEL <- data.frame(x = theta_grid, y = WEL_logLR, type = "WEL")
      dt <- rbind(dt, df_WEL)
    }
    # ET logLR
    if ("et" %in% type()) {
      # set.seed(seed())
      ET_logLR <- sapply(theta_grid, logET)
      df_ET <- data.frame(x = theta_grid, y = ET_logLR, type = "WETEL")
      dt <- rbind(dt, df_ET)
    }

    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("EL", "WEL", "WETEL"))
    # logLR plot
    ggplot(dt,
           aes(x, y, color = type, linetype = type)) +
      geom_path(alpha = 0.5, na.rm = TRUE) +
      labs(x = expression(theta), y = "logLR",
           title = "logLR") +
      geom_vline(xintercept =  c(min(x), mean(x), max(x)), linetype = "dashed",
                 alpha = 0.5) +
      xlim(grid()[1], grid()[2]) +
      theme(axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, size = 1),
            panel.grid = element_blank(),
            legend.text = element_text(size = 10, color = "black"),
            legend.background = element_rect(fill = alpha("white", 0)),
            legend.key = element_rect(fill = alpha("white", 1)),
            legend.title = element_blank()) +
      scale_color_npg()
  })
  # posterior density plot (discrete, normalized)
  output$plot2 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # generate grid
    theta_grid <- seq(grid()[1], grid()[2], length.out = 500)
    # posterior mean
    posterior_mean <-
      (theta0() / sigma0()^2 + sum(x) / sigma()^2) /
      (1 / sigma0()^2 + n() / sigma()^2)
    # posterior sd
    posterior_sd <- (1 / sigma0()^2 + n() / sigma()^2)^(-1 / 2)
    # weights
    # beta <- 1 / m()
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w


    true_log_pd <- function(par) {
      dnorm(par, mean = posterior_mean, sd = posterior_sd, log = T)
    }
    EL_logLR <- function(par) {
      z <- el_mean(par, x, control = list(abstol = 1e-04,
                                          threshold = 1e+10))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }
    logWLR <- function(par) {
      x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
      pp <- el_mean(par, c(x, x_pseudo), w, control =
                      list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
      p1 <- pp[seq_len(n())]
      p2 <- pp[-seq_len(n())]
      sum(p1) + n() * (log(n() + strength()) - log(n() + m()))
    }
    opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-08)
    logET <- function(par) {
      x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
      g <- c(x, x_pseudo) - par

      # eval_f <- function(l) mean(w * exp(l * g))
      # eval_grad_f <- function(l) mean(w * exp(l * g) * g)
      # lambda <- nloptr(x0 = 0, eval_f = eval_f, eval_grad_f = eval_grad_f,
      #                  opts = opts)$solution

      lambda_finder <- function(l) mean(w * exp(l * g) * g)
      lambda <- uniroot(lambda_finder, extendInt = "yes",
                        lower = -1e+10, upper = 1e+10)$root

      unnormalized_prob <- w * exp(lambda * g)
      log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
      sum(log_prob[seq_len(n())]) + n() * log(n() + strength())
    }

    # true posterior density (discrete, normalized)
    true_log_pd <- sapply(theta_grid, true_log_pd)
    true_pd <- exp(true_log_pd) / sum(exp(true_log_pd))
    # EL posterior density (discrete, normalized)
    EL_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
      sapply(theta_grid, EL_logLR)
    EL_pd <- exp(EL_log_pd) / sum(exp(EL_log_pd), na.rm = T)
    df_true <- data.frame(x = theta_grid, y = true_pd, type = "True")
    df_EL <- data.frame(x = theta_grid, y = EL_pd, type = "EL")
    dt <- rbind(df_true, df_EL)
    # posterior density with real data (discrete, normalized)
    if ("wel" %in% type()) {
      set.seed(seed())
      WEL_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logWLR)
      WEL_pd <- exp(WEL_log_pd) / sum(exp(WEL_log_pd))
      df_WEL1 <- data.frame(x = theta_grid, y = WEL_pd, type = "WEL")
      dt <- rbind(dt, df_WEL1)
    }
    if ("et" %in% type()) {
      set.seed(seed())
      ET_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logET)
      ET_pd <- exp(ET_log_pd) / sum(exp(ET_log_pd))
      df_ET <- data.frame(x = theta_grid, y = ET_pd, type = "WETEL")
      dt <- rbind(dt, df_ET)
    }

    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("EL", "WEL", "WETEL", "True"))
    ggplot(dt, aes(x, y, color = type, linetype = type)) +
      geom_path(alpha = 0.5, na.rm = TRUE) +
      labs(x = expression(theta), y = "value",
           title = "Normalized Posterior Densities") +
      geom_vline(xintercept =  c(min(x), mean(x), max(x)), linetype = "dashed",
                 alpha = 0.5) +
      xlim(grid()[1], grid()[2]) +
      theme(axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, size = 1),
            panel.grid = element_blank(),
            legend.text = element_text(size = 10, color = "black"),
            legend.background = element_rect(fill = alpha("white", 0)),
            legend.key = element_rect(fill = alpha("white", 1)),
            legend.title = element_blank()) +
      scale_color_npg()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
