library(tidyverse)
library(melt)
library(VGAM)
library(ggplot2)
library(ggsci)

ui <- fluidPage(
  titlePanel("Bayesian Weighted EL with Fractional Pseudo Observations (Normal-Normal Model)"),
  tabsetPanel(
    tabPanel("LR and Posterior Density", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 sliderInput("shape", "alpha (skew normal)",
                             value = 10, min = 0, max = 20, step = 0.1),
                 sliderInput("theta0", "prior mean",
                             value = 1, min = ,-1, max = 1, step = 0.05),
                 sliderInput("sigma0", "prior sd",
                             value = 1, min = 0.1, max = 2),
                 sliderInput("grid", "grid range",
                             value = c(0, 2.3), min = -1, max = 5, step = 0.1),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 20,
                             min = 10, max = 100, step = 1),
                 sliderInput("m", "m: pseudo sample size",
                             value = 10, min = 2, max = 100, step = 1),
                 checkboxGroupInput("wel", "weighted EL to use",
                                    choices = c("WEL1 (qnorm)" = "wel1",
                                                "WEL2 (qunif)" = "wel2"),
                                    selected = c("wel1"),
                                    inline = FALSE)),
               mainPanel(
                 textOutput("summary1"),
                 verbatimTextOutput("code1"),
                 plotOutput("plot_data"),
                 plotOutput("plot1"),
                 plotOutput("plot2")
               )
             )
    ),
    # tabPanel("Metropolis Sampling", fluid = TRUE,
    #          sidebarLayout(
    #            sidebarPanel(
    #              sliderInput("B", "B: iteration number", value = 1000,
    #                          min = 1000, max = 5000, step = 100),
    #              sliderInput("sigma_p", "proposal sd",
    #                          value = 1 / 2, min = 0.1, max = 2),
    #              sliderInput("shape_mcmc", "sampling shape (skew normal)",
    #                          value = 1, min = 0, max = 10, step = 0.1),
    #              sliderInput("theta_mcmc", "model mean",
    #                          value = 0, min = ,-3, max = 3, step = 0.05),
    #              sliderInput("sigma_mcmc", "model sd",
    #                          value = 1, min = 0.5, max = 3),
    #              sliderInput("theta0_mcmc", "prior mean",
    #                          value = -1, min = ,-3, max = 3, step = 0.05),
    #              sliderInput("sigma0_mcmc", "prior sd",
    #                          value = 1, min = 0.5, max = 3),
    #              numericInput("seed_mcmc", "seed",
    #                           value = 1, min = 1, max = .Machine$integer.max),
    #              sliderInput("n_mcmc", "n: sample size", value = 100,
    #                          min = 10, max = 500, step = 10),
    #              sliderInput("m_mcmc", "m: pseudo sample size",
    #                          value = 2, min = 2, max = 100, step = 1),
    #              selectInput("wel_mcmc", "weighted EL to use",
    #                          choices = c("WEL1 (w/o pseudo obs. & discard weights)" = "wel1",
    #                                      "WEL2 (w/ pseudo obs. & retain weights)" = "wel2"))),
    #            mainPanel(
    #              plotOutput("plot3"),
    #              # plotOutput("plot4")
    #            )
    #          )
    # )
  )
)


server <- function(input, output) {
  options(warn = -1)
  # parameters for densities
  shape <- reactive(input$shape)
  theta0 <- reactive(input$theta0)
  sigma0 <- reactive(input$sigma0)
  n <- reactive(input$n)
  m <- reactive(input$m)
  grid <- reactive(input$grid)
  seed <- reactive(input$seed)
  wel <- reactive(input$wel)

  # summary statistics
  output$summary1 <- renderText("Summary of sample data")
  output$code1 <- renderPrint({
    set.seed(seed())
    x <- rskewnorm(n(), shape = shape())
    summary(x)
  }
  )
  # density plot
  output$plot_data <- renderPlot({
    set.seed(seed())
    x <- rskewnorm(n(), shape = shape())
    ggplot(data.frame(x), aes(x)) + geom_density()
  })

  # LR plot
  output$plot1 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rskewnorm(n(), shape = shape())
    # generate grid
    grid <- seq(grid()[1], grid()[2], length.out = 1000)
    # weights
    beta <- 1 / m()
    w <- c(rep(1, n()), rep(beta, m()))

    EL_logLR <- function(par) {
      z <- el_mean(par, x, control = list(abstol = 1e-05,
                                          threshold = 1e+10))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }
    AEL_logLR <- function(par) {
      g <- x - par
      g_pseudo <- -log(n()) / 2 * mean(g) * 0.5
      el_eval(c(g, g_pseudo),
              control = list(abstol = 1e-05, threshold = 1e+10))$optim$logLR
    }
    logWLR1 <- function(par) {
      x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = sd(x))
      pp <- el_mean(par, c(x, x_pseudo), w, control =
                      list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
      p1 <- pp[seq_len(n())]
      p2 <- pp[-seq_len(n())]
      sum(p1) + n() * (log(n() + 1) - log(n() + m()))
    }
    logWLR2 <- function(par) {
      x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
      pp <- el_mean(par, c(x, x_pseudo), w, control =
                      list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
      p1 <- pp[seq_len(n())]
      p2 <- pp[-seq_len(n())]
      sum(p1) + n() * (log(n() + 1) - log(n() + m()))
    }

    # EL logLR
    df_EL <- data.frame(x = grid, y = sapply(grid, EL_logLR), type = "EL")
    df_AEL <- data.frame(x = grid, y = sapply(grid, AEL_logLR), type = "AEL")
    dt <- rbind(df_EL, df_AEL)
    # WEL1 density
    if ("wel1" %in% wel()) {
      df_WEL1 <- data.frame(x = grid, y = sapply(grid, logWLR1), type = "WEL1")
      dt <- rbind(dt, df_WEL1)
    }
    # WEL2 density
    if ("wel2" %in% wel()) {
      df_WEL2 <- data.frame(x = grid, y = sapply(grid, logWLR2), type = "WEL2")
      dt <- rbind(dt, df_WEL2)
    }

    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("EL", "AEL", "WEL1", "WEL2"))
    # density plot (discrete, normalized)
    ggplot(dt,
           aes(x, y, color = type)) +
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
    x <- rskewnorm(n(), shape = shape())
    # generate grid
    grid <- seq(grid()[1], grid()[2], length.out = 1000)

    var <- 1 - 2 / pi * (shape()^2 / (1 + shape()^2))
    # posterior mean
    posterior_mean <-
      (theta0() / sigma0()^2 + sum(x) / var) /
      (1 / sigma0()^2 + n() / var)
    # posterior sd
    posterior_sd <- (1 / sigma0()^2 + n() / var)^(-1 / 2)

    # weights
    beta <- 1 / m()
    w <- c(rep(1, n()), rep(beta, m()))

    normal_log_pd <- function(par) {
      dnorm(par, mean = posterior_mean, sd = posterior_sd, log = T)
    }
    # normal model posterior density (discrete, normalized)
    normal_model_log_pd <- sapply(grid, normal_log_pd)
    normal_model_pd <- exp(normal_model_log_pd) / sum(exp(normal_model_log_pd))
    df_true <- data.frame(x = grid, y = normal_model_pd, type = "Normal model")
    # EL posterior density (discrete, normalized)
    EL_logLR <- function(par) {
      z <- el_mean(par, x, control = list(abstol = 1e-05,
                                          threshold = 1e+10))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }
    EL_log_pd <- dnorm(grid, theta0(), sigma0(), log = T) +
      sapply(grid, EL_logLR)
    EL_pd <- exp(EL_log_pd) / sum(exp(EL_log_pd), na.rm = T)
    df_EL <- data.frame(x = grid, y = EL_pd, type = "EL")

    logWLR1 <- function(par) {
      x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = sd(x))
      pp <- el_mean(par, c(x, x_pseudo), w, control =
                      list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
      p1 <- pp[seq_len(n())]
      p2 <- pp[-seq_len(n())]
      sum(p1) + n() * (log(n() + 1) - log(n() + m()))
    }
    logWLR2 <- function(par) {
      x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
      pp <- el_mean(par, c(x, x_pseudo), w, control =
                      list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
      p1 <- pp[seq_len(n())]
      p2 <- pp[-seq_len(n())]
      sum(p1) + n() * (log(n() + 1) - log(n() + m()))
    }

    # EL logLR
    # df_AEL <- data.frame(x = grid, y = sapply(grid, AEL_logLR), type = "AEL")
    # adjust factor levels for plot
    dt <- rbind(df_true, df_EL)
    # posterior density with real data (discrete, normalized)
    if ("wel1" %in% wel()) {
      set.seed(seed())
      WEL1_log_pd <- dnorm(grid, theta0(), sigma0(), log = T) +
        sapply(grid, logWLR1)
      WEL1_pd <- exp(WEL1_log_pd) / sum(exp(WEL1_log_pd))
      df_WEL1 <- data.frame(x = grid, y = WEL1_pd, type = "WEL1")
      dt <- rbind(dt, df_WEL1)
    }
    if ("wel2" %in% wel()) {
      set.seed(seed())
      WEL2_log_pd <- dnorm(grid, theta0(), sigma0(), log = T) +
        sapply(grid, logWLR2)
      WEL2_pd <- exp(WEL2_log_pd) / sum(exp(WEL2_log_pd))
      df_WEL2 <- data.frame(x = grid, y = WEL2_pd, type = "WEL2")
      dt <- rbind(dt, df_WEL2)
    }
    # adjust factor levels for plot
    dt$type <- factor(dt$type,
                      levels = c("Normal model", "EL", "WEL1", "WEL2", "True"))
    ggplot(dt, aes(x, y, color = type)) +
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

  # # parameters for MCMC
  # B <- reactive(input$B)
  # sigma_p <- reactive(input$sigma_p)
  # shape_mcmc <- reactive(input$shape_mcmc)
  # theta_mcmc <- reactive(input$theta_mcmc)
  # sigma_mcmc <- reactive(input$sigma_mcmc)
  # theta0_mcmc <- reactive(input$theta0_mcmc)
  # sigma0_mcmc <- reactive(input$sigma0_mcmc)
  #
  # seed_mcmc <- reactive(input$seed_mcmc)
  # beta_mcmc <- reactive(input$beta_mcmc)
  # n_mcmc <- reactive(input$n_mcmc)
  # m_mcmc <- reactive(input$m_mcmc)
  # wel_mcmc <- reactive(input$wel_mcmc)
  # # MCMC samples
  # output$plot3 <- renderPlot({
  #   # seed for sampling
  #   set.seed(seed_mcmc())
  #   # simulate data
  #   x <- rskewnorm(n_mcmc(), location = 0, scale = 1, shape = shape_mcmc())
  #   # weights
  #   beta <- 1 / m_mcmc()
  #   w <- c(rep(1, n_mcmc()), rep(beta, m_mcmc()))
  #
  #   # theta_sample_true <- vector("numeric", B())
  #   # theta_sample_true[1] <- theta0_mcmc()
  #   theta_sample_model <- vector("numeric", B())
  #   theta_sample_model[1] <- theta0_mcmc()
  #   theta_sample_el <- vector("numeric", B())
  #   theta_sample_el[1] <- theta0_mcmc()
  #   theta_sample_wel <- vector("numeric", B())
  #   theta_sample_wel[1] <- theta0_mcmc()
  #
  #   # set.seed(seed_mcmc())
  #   # for (i in 2:B()) {
  #   #   # sample proposal value
  #   #   theta_proposal <- rnorm(1, theta_sample_true[i - 1], sigma_p())
  #   #   # compute log ratio posterior densities
  #   #   logr <-
  #   #     (dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc(), log = TRUE) -
  #   #        dnorm(theta_sample_true[i - 1],
  #   #              theta0_mcmc(), sigma0_mcmc(), log = TRUE)) +
  #   #     (dskewnorm(theta_proposal, shape = shape_mcmc(), log = TRUE) -
  #   #        dskewnorm(theta_sample_true[i - 1], shape = shape_mcmc(), log = TRUE))
  #   #   # sample uniform random variable
  #   #   u <- runif(1)
  #   #   # accept or reject
  #   #   if (log(u) < logr) {
  #   #     theta_sample_true[i] <- theta_proposal
  #   #   } else {
  #   #     theta_sample_true[i] <- theta_sample_true[i - 1]
  #   #   }
  #   # }
  #
  #   set.seed(seed_mcmc())
  #   for (i in 2:B()) {
  #     # sample proposal value
  #     theta_proposal <- rnorm(1, theta_sample_model[i - 1], sigma_p())
  #     # compute log ratio posterior densities
  #     logr <-
  #       (dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc(), log = TRUE) -
  #          dnorm(theta_sample_model[i - 1],
  #                theta0_mcmc(), sigma0_mcmc(), log = TRUE)) +
  #       (sum(dnorm(x, theta_proposal, sigma_mcmc(), log = TRUE)) -
  #          sum(dnorm(x, theta_sample_model[i - 1], sigma_mcmc(), log = TRUE)))
  #     # sample uniform random variable
  #     u <- runif(1)
  #     # accept or reject
  #     if (log(u) < logr) {
  #       theta_sample_model[i] <- theta_proposal
  #     } else {
  #       theta_sample_model[i] <- theta_sample_model[i - 1]
  #     }
  #   }
  #
  #
  #   logWLR1 <- function(par) {
  #     x_pseudo <- par + qnorm(1:m_mcmc() / (m_mcmc() + 1), mean = 0, sd = sd(x))
  #     pp <- el_mean(par, c(x, x_pseudo), w, control =
  #                     list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
  #     p1 <- pp[seq_len(n_mcmc())]
  #     p2 <- pp[-seq_len(n_mcmc())]
  #     sum(p1) + n_mcmc() * (log(n_mcmc() + 1) - log(n_mcmc() + m_mcmc()))
  #   }
  #   logWLR2 <- function(par) {
  #     x_pseudo <- par + qunif(1:m_mcmc() / (m_mcmc() + 1), min = -IQR(x) / 2,
  #                             max = IQR(x) / 2)
  #     pp <- el_mean(par, c(x, x_pseudo), w, control =
  #                     list(abstol = 1e-05, threshold = 1e+10))$optim$log.prob
  #     p1 <- pp[seq_len(n_mcmc())]
  #     p2 <- pp[-seq_len(n_mcmc())]
  #     sum(p1) + n_mcmc() * (log(n_mcmc() + 1) - log(n_mcmc() + m_mcmc()))
  #   }
  #   set.seed(seed_mcmc())
  #   for (i in 2:B()) {
  #     # sample proposal value
  #     theta_proposal <- rnorm(1, theta_sample_wel[i - 1], sigma_p())
  #     # compute log ratio posterior densities
  #     if (wel_mcmc() == "wel1") {
  #       LR1 <- logWLR1(theta_proposal)
  #       LR2 <- logWLR1(theta_sample_wel[i - 1])
  #     } else {
  #       LR1 <- logWLR2(theta_proposal)
  #       LR2 <- logWLR2(theta_sample_wel[i - 1])
  #     }
  #     logr <-
  #       (dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc(), log = TRUE) -
  #          dnorm(theta_sample_model[i - 1],
  #                theta0_mcmc(), sigma0_mcmc(), log = TRUE)) +
  #       (LR1 - LR2)
  #     # sample uniform random variable
  #     u <- runif(1)
  #     # accept or reject
  #     if (log(u) < logr) {
  #       theta_sample_wel[i] <- theta_proposal
  #     } else {
  #       theta_sample_wel[i] <- theta_sample_wel[i - 1]
  #     }
  #   }
  #
  #
  #
  #
  #
  #
  #   df_model <- data.frame(x = seq_len(B()), y = theta_sample_model,
  #                          type = "model")
  #   df_wel <- data.frame(x = seq_len(B()), y = theta_sample_wel,
  #                        type = toupper(wel_mcmc()))
  #
  #   ggplot(rbind(df_model, df_wel), aes(x, y, color = type)) +
  #     geom_path() +
  #     theme(axis.text = element_text(size = 12, color = "black"),
  #           axis.title = element_text(size = 12),
  #           panel.background = element_blank(),
  #           panel.border = element_rect(fill = NA, size = 1),
  #           panel.grid = element_blank(),
  #           legend.text = element_text(size = 10, color = "black"),
  #           legend.background = element_rect(fill = alpha("white", 0)),
  #           legend.key = element_rect(fill = alpha("white", 1)),
  #           legend.title = element_blank()) +
  #     labs(x = "Iteration", y = expression(theta)) +
  #     scale_color_npg()
  # })
  #
  # # MCMC samples
  # output$plot4 <- renderPlot({
  #   # seed for sampling
  #   set.seed(seed_mcmc())
  #   # simulate data
  #   x <- rnorm(n_mcmc(), mean = theta_mcmc(), sd = sigma_mcmc())
  #   # weights
  #   if (beta_mcmc() == "1 / m") {
  #     beta <- 1 / m_mcmc()
  #   } else if (beta_mcmc() == "m / n") {
  #     beta <- m_mcmc() / n_mcmc()
  #   } else {
  #     beta <- 1
  #   }
  #   w <- c(rep(1, n_mcmc()), rep(beta, m_mcmc()))
  #
  #   pseudo_sample <- function(par) {
  #     par + qnorm(1:m_mcmc() / (m_mcmc() + 1))
  #   }
  #   logWLR <- function(par) {
  #     x_aug <- c(x, pseudo_sample(par))
  #     sum(el_mean(par, x_aug, w, list(threshold = 1e+10, maxit = 100,
  #                                     abstol = 1e-06))$optim$log.prob[1:n_mcmc()]) +
  #       n_mcmc() * (log(n_mcmc()) - log(n_mcmc() + m_mcmc()))
  #   }
  #   logWLR_aug <- function(par) {
  #     x_aug <- c(x, pseudo_sample(par))
  #     el_mean(par, x_aug, w, list(threshold = 1e+10, maxit = 100,
  #                                 abstol = 1e-06))$optim$logWLR
  #   }
  #
  #   theta_sample_el <- vector("numeric", B())
  #   theta_sample_el[1] <- theta0_mcmc()
  #   theta_sample_wel <- vector("numeric", B())
  #   theta_sample_wel[1] <- theta0_mcmc()
  #
  #   for (i in 2:B()) {
  #     # sample proposal value
  #     theta_proposal <- rnorm(1, theta_sample_el[i - 1], sigma_p())
  #
  #     # compute log ratio posterior densities
  #     LR1 <- el_mean(theta_proposal, x)$optim$logLR
  #     LR2 <- el_mean(theta_sample_el[i - 1], x)$optim$logLR
  #     logr <-
  #       (log(dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc())) -
  #          log(dnorm(theta_sample_el[i - 1], theta0_mcmc(), sigma0_mcmc()))) +
  #       (LR1 - LR2)
  #     # sample uniform random variable
  #     u <- runif(1)
  #
  #     # accept or reject
  #     if (log(u) < logr) {
  #       theta_sample_el[i] <- theta_proposal
  #     } else {
  #       theta_sample_el[i] <- theta_sample_el[i - 1]
  #     }
  #   }
  #
  #   for (i in 2:B()) {
  #     # sample proposal value
  #     theta_proposal <- rnorm(1, theta_sample_wel[i - 1], sigma_p())
  #     if (wel_mcmc() == "wel1") {
  #       LR1 <- logWLR(theta_proposal)
  #       LR2 <- logWLR(theta_sample_wel[i - 1])
  #     } else {
  #       LR1 <- logWLR_aug(theta_proposal)
  #       LR2 <- logWLR_aug(theta_sample_wel[i - 1])
  #     }
  #     logr <-
  #       (log(dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc())) -
  #          log(dnorm(theta_sample_wel[i - 1], theta0_mcmc(), sigma0_mcmc()))) +
  #       (LR1 - LR2)
  #     # sample uniform random variable
  #     u <- runif(1)
  #
  #     # accept or reject
  #     if (log(u) < logr) {
  #       theta_sample_wel[i] <- theta_proposal
  #     } else {
  #       theta_sample_wel[i] <- theta_sample_wel[i - 1]
  #     }
  #   }
  #
  #   df_el <- data.frame(x = seq_len(B()), y = theta_sample_el, type = "EL")
  #   df_wel <- data.frame(x = seq_len(B()), y = theta_sample_wel,
  #                        type = toupper(wel_mcmc()))
  #
  #   # posterior mean
  #   posterior_mean <-
  #     (theta0_mcmc() / sigma0_mcmc()^2 + sum(x) / sigma_mcmc()^2) /
  #     (1 / sigma0_mcmc()^2 + n_mcmc() / sigma_mcmc()^2)
  #   # posterior sd
  #   posterior_sd <- (1 / sigma0_mcmc()^2 + n_mcmc() / sigma_mcmc()^2)^(-1 / 2)
  #
  #   ggplot(rbind(df_el, df_wel), aes(y, fill = type, color = type)) +
  #     geom_histogram(aes(y = ..density..), position = "identity",
  #                    bins = 50, alpha = 0.2) +
  #     theme(axis.text = element_text(size = 12, color = "black"),
  #           axis.title = element_text(size = 12),
  #           panel.background = element_blank(),
  #           panel.border = element_rect(fill = NA, size = 1),
  #           panel.grid = element_blank(),
  #           legend.text = element_text(size = 10, color = "black"),
  #           legend.background = element_rect(fill = alpha("white", 0)),
  #           legend.key = element_rect(fill = alpha("white", 1)),
  #           legend.title = element_blank()) +
  #     labs(x = expression(theta)) +
  #     stat_function(fun = dnorm,
  #                   args = list(mean = posterior_mean,
  #                               sd = posterior_sd)) +
  #     scale_color_npg()
  # })
}

# Run the application
shinyApp(ui = ui, server = server)
