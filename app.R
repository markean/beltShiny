library(shiny)
library(melt)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Near Bayes (Normal - Normal)"),
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
                             value = c(-1, 1), min = -4, max = 4, step = 0.05),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 100,
                             min = 10, max = 500, step = 10),
                 sliderInput("m", "m: pseudo sample size",
                             value = 2, min = 2, max = 100, step = 1),
                 selectInput("beta", "beta (fractional weight of one pseudo observation)",
                             choices = c("1 / m", "m / n", "1")),
                 sliderInput("plotlim", "Plot x-axis range",
                             value = c(-1, 1), min = -4, max = 4,
                             round = F, step = 0.05),
                 selectInput("method", "pseudo sample method:",
                             choices = c("qnorm" = "qnorm", "rnorm" = "rnorm")),
                 checkboxGroupInput("wel", "Weighted EL",
                                    choices = c("wel1" = "wel1", "wel2" = "wel2",
                                                "wel3" = "wel3", "wel4" = "wel4"),
                                    selected = c("wel1"))),

               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", verbatimTextOutput("code"), plotOutput("plot1"),
                            plotOutput("plot2"))
                 )
                 # textOutput("summary"),
                 # verbatimTextOutput("code"),
                 # plotOutput("plot1"),
                 # plotOutput("plot2")
               )
             )
             ),
    tabPanel("Metropolis Sampling", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 sliderInput("B", "B: iteration number", value = 1000,
                             min = 1000, max = 5000, step = 100),
                 sliderInput("sigma_p", "proposal sd",
                             value = 1 / 2, min = 0.1, max = 2),
                 sliderInput("theta", "sampling mean",
                             value = 0, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma", "sampling sd",
                             value = 1, min = 0.5, max = 3),
                 sliderInput("theta0", "prior mean",
                             value = -1, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma0", "Prior sd",
                             value = 1, min = 0.5, max = 3),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 100,
                             min = 10, max = 500, step = 10),
                 sliderInput("m", "m: pseudo sample size",
                             value = 2, min = 2, max = 100, step = 1),
                 selectInput("beta", "beta (fractional weight of one pseudo observation)",
                             choices = c("1 / m", "m / n", "1")),
                 selectInput("method", "pseudo sample method:",
                             choices = c("qnorm" = "qnorm", "rnorm" = "rnorm")),
                 checkboxGroupInput("wel", "Weighted EL",
                                    choices = c("wel1" = "wel1", "wel2" = "wel2",
                                                "wel3" = "wel3", "wel4" = "wel4"),
                                    selected = c("wel1"))),

               mainPanel(
                 plotOutput("plot3"),
                 plotOutput("plot4")
               )
             )

             )
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
  method <- reactive(input$method)
  wel <- reactive(input$wel)

  # parameters 4
  B <- reactive(input$B)
  sigma_p <- reactive(input$sigma_p)

  # summary statistics
  output$summary <- renderText("Summary of sample data")
  output$code <- renderPrint({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    summary(x)
    })

  # LR plot
  output$plot1 <- renderPlot({
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
      z <- el_mean(par, x, control = list(threshold = 1e+10))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }
    pseudo_sample <- function(par) {
      if (method() == "qnorm") {
        pseudo_x <- par + qnorm(1:m() / (m() + 1))
      } else {
        pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
          pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        }
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
    df_EL <- data.frame(x = theta_grid, y = EL_log_d, type = "EL")
    # WEL1 density
    if ("wel1" %in% wel()) {
      set.seed(seed())
      WEL1_log_d <- sapply(theta_grid, logLR)
      WEL1_d <- exp(WEL1_log_d) / sum(exp(WEL1_log_d))
      df_WEL1 <- data.frame(x = theta_grid, y = WEL1_log_d, type = "WEL1")
    } else {
      df_WEL1 <- NULL
    }
    # WEL2 density
    if ("wel2" %in% wel()) {
      set.seed(seed())
      WEL2_log_d <- sapply(theta_grid, logWLR)
      WEL2_d <- exp(WEL2_log_d) / sum(exp(WEL2_log_d))
      df_WEL2 <- data.frame(x = theta_grid, y = WEL2_log_d, type = "WEL2")
    } else {
      df_WEL2 <- NULL
    }
    # WEL3 density
    if ("wel3" %in% wel()) {
      set.seed(seed())
      WEL3_log_d <- sapply(theta_grid, logLR_aug)
      WEL3_d <- exp(WEL3_log_d) / sum(exp(WEL3_log_d))
      df_WEL3 <- data.frame(x = theta_grid, y = WEL3_log_d, type = "WEL3")
    } else {
      df_WEL3 <- NULL
    }
    # WEL3 density
    if ("wel4" %in% wel()) {
      set.seed(seed())
      WEL4_log_d <- sapply(theta_grid, logWLR_aug)
      WEL4_d <- exp(WEL4_log_d) / sum(exp(WEL4_log_d))
      df_WEL4 <- data.frame(x = theta_grid, y = WEL4_log_d, type = "WEL4")
    } else {
      df_WEL4 <- NULL
    }

    # density plot (discrete, normalized)
    ggplot(rbind(df_EL, df_WEL1, df_WEL2, df_WEL3, df_WEL4),
           aes(x, y, color = type)) +
      geom_path(alpha = 0.5, na.rm = TRUE) +
      labs(x = expression(theta), y = "logLR",
           title = "logLR") +
      geom_vline(xintercept =  c(min(x), mean(x), max(x)), linetype = "dashed",
                 alpha = 0.5) +
      xlim(plotlim()[1], plotlim()[2])
  })
  # posterior density plot (discrete, normalized)
  output$plot2 <- renderPlot({
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
      z <- el_mean(par, x, control = list(threshold = 1e+10))$optim
      if (!isTRUE(z$convergence)) {
        return(NA)
      }
      z$logLR
    }

    pseudo_sample <- function(par) {
      if (method() == "qnorm") {
        pseudo_x <- par + qnorm(1:m() / (m() + 1))
      } else {
        pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
          pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        }
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
  # MCMC samples
  output$plot3 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # weights
    if (beta() == "1 / m") {
      beta <- 1 / m()
    } else if (beta() == "m / n") {
      beta <- m() / n()
    } else {
      beta <- 1
    }
    w <- c(rep(1, n()), rep(beta, m()))

    pseudo_sample <- function(par) {
      if (method() == "qnorm") {
        pseudo_x <- par + qnorm(1:m() / (m() + 1))
      } else {
        pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
          pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        }
      }
      pseudo_x
    }

    theta_sample_el <- vector("numeric", B())
    theta_sample_el[1] <- mean(x)
    theta_sample_wel <- vector("numeric", B())
    theta_sample_wel[1] <- mean(x)

    for (i in 2:B()) {
      # sample proposal value
      theta_proposal <- rnorm(1, theta_sample_el[i - 1], sigma_p())

      # compute log ratio posterior densities
      LR1 <- el_mean(theta_proposal, x)$optim$logLR
      LR2 <- el_mean(theta_sample_el[i - 1], x)$optim$logLR
      logr <-
        (log(dnorm(theta_proposal, theta0(), sigma0())) -
           log(dnorm(theta_sample_el[i - 1], theta0(), sigma0()))) +
        (LR1 - LR2)
      # sample uniform random variable
      u <- runif(1)

      # accept or reject
      if (log(u) < logr) {
        theta_sample_el[i] <- theta_proposal
      } else {
        theta_sample_el[i] <- theta_sample_el[i - 1]
      }
    }

    for (i in 2:B()) {
      # sample proposal value
      theta_proposal <- rnorm(1, theta_sample_wel[i - 1], sigma_p())

      # pseudo sample
      x_aug <- c(x, pseudo_sample(theta_proposal))
      # compute log ratio posterior densities
      LR1 <- sum(el_mean(theta_proposal, x_aug, w)$optim$log.prob[1:n()])
      LR2 <- sum(el_mean(theta_sample_wel[i - 1], x_aug, w)$optim$log.prob[1:n()])
      logr <-
        (log(dnorm(theta_proposal, theta0(), sigma0())) -
           log(dnorm(theta_sample_wel[i - 1], theta0(), sigma0()))) +
        (LR1 - LR2)
      # sample uniform random variable
      u <- runif(1)

      # accept or reject
      if (log(u) < logr) {
        theta_sample_wel[i] <- theta_proposal
      } else {
        theta_sample_wel[i] <- theta_sample_wel[i - 1]
      }
    }

    df_el <- data.frame(ss = 1:length(theta_sample_el), x = theta_sample_el, type = "EL")
    df_wel <- data.frame(ss = 1:length(theta_sample_wel), x = theta_sample_wel, type = "WEL")

    ggplot(rbind(df_el, df_wel), aes(x = ss, y = x, color = type)) +
      geom_path() +
      theme(axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, size = 1),
            panel.grid = element_blank(),
            legend.text = element_text(size = 10, color = "black"),
            legend.background = element_rect(fill = alpha("white", 0)),
            legend.key = element_rect(fill = alpha("white", 1)),
            legend.title = element_blank()) +
      labs(x = "Iteration", y = expression(theta)) +
      scale_color_npg()
  })
  # MCMC samples
  output$plot4 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(input$n, mean = input$theta, sd = input$sigma)
    # weights
    if (beta() == "1 / m") {
      beta <- 1 / m()
    } else if (beta() == "m / n") {
      beta <- m() / n()
    } else {
      beta <- 1
    }
    w <- c(rep(1, n()), rep(beta, m()))

    pseudo_sample <- function(par) {
      if (method() == "qnorm") {
        pseudo_x <- par + qnorm(1:m() / (m() + 1))
      } else {
        pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
          pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        }
      }
      pseudo_x
    }

    theta_sample_el <- vector("numeric", B())
    theta_sample_el[1] <- mean(x)
    theta_sample_wel <- vector("numeric", B())
    theta_sample_wel[1] <- mean(x)

    aa <- function(i) {
      # sample proposal value
      theta_proposal <- rnorm(1, theta_sample_el[i - 1], input$sigma_p)

      # compute log ratio posterior densities
      LR1 <- el_mean(theta_proposal, x)$optim$logLR
      LR2 <- el_mean(theta_sample_el[i - 1], x)$optim$logLR
      logr <-
        (log(dnorm(theta_proposal, input$theta0, input$sigma0)) -
           log(dnorm(theta_sample_el[i - 1], input$theta0), input$sigma0)) +
        (LR1 - LR2)
      # sample uniform random variable
      u <- runif(1)

      # accept or reject
      if (log(u) < logr) {
        theta_sample_el[i] <- theta_proposal
      } else {
        theta_sample_el[i] <- theta_sample_el[i - 1]
      }
    }
    sapply(2:input$B, aa)

    # for (i in 2:input$B) {
    #   # sample proposal value
    #   theta_proposal <- rnorm(1, theta_sample_el[i - 1], input$sigma_p)
    #
    #   # compute log ratio posterior densities
    #   LR1 <- el_mean(theta_proposal, x)$optim$logLR
    #   LR2 <- el_mean(theta_sample_el[i - 1], x)$optim$logLR
    #   logr <-
    #     (log(dnorm(theta_proposal, input$theta0, input$sigma0)) -
    #        log(dnorm(theta_sample_el[i - 1], input$theta0), input$sigma0)) +
    #     (LR1 - LR2)
    #   # sample uniform random variable
    #   u <- runif(1)
    #
    #   # accept or reject
    #   if (log(u) < logr) {
    #     theta_sample_el[i] <- theta_proposal
    #   } else {
    #     theta_sample_el[i] <- theta_sample_el[i - 1]
    #   }
    # }
    #
    bb <- function(i) {
      # sample proposal value
      theta_proposal <- rnorm(1, theta_sample_wel[i - 1], input$sigma_p)

      # pseudo sample
      x_aug <- c(x, pseudo_sample(theta_proposal))
      # compute log ratio posterior densities
      LR1 <- sum(el_mean(theta_proposal, x_aug, w)$optim$log.prob[1:input$n])
      LR2 <- sum(el_mean(theta_sample_wel[i - 1], x_aug, w)$optim$log.prob[1:input$n])
      logr <-
        (log(dnorm(theta_proposal, input$theta0, input$sigma0)) -
           log(dnorm(theta_sample_wel[i - 1], input$theta0, input$sigma0))) +
        (LR1 - LR2)
      # sample uniform random variable
      u <- runif(1)

      # accept or reject
      if (log(u) < logr) {
        theta_sample_wel[i] <- theta_proposal
      } else {
        theta_sample_wel[i] <- theta_sample_wel[i - 1]
      }

    }
    sapply(2:B(), bb)

    # for (i in 2:input$B) {
    #   # sample proposal value
    #   theta_proposal <- rnorm(1, theta_sample_wel[i - 1], input$sigma_p)
    #
    #   # pseudo sample
    #   x_aug <- c(x, pseudo_sample(theta_proposal))
    #   # compute log ratio posterior densities
    #   LR1 <- sum(el_mean(theta_proposal, x_aug, w)$optim$log.prob[1:input$n])
    #   LR2 <- sum(el_mean(theta_sample_wel[i - 1], x_aug, w)$optim$log.prob[1:input$n])
    #   logr <-
    #     (log(dnorm(theta_proposal, input$theta0, input$sigma0)) -
    #        log(dnorm(theta_sample_wel[i - 1], input$theta0, input$sigma0))) +
    #     (LR1 - LR2)
    #   # sample uniform random variable
    #   u <- runif(1)
    #
    #   # accept or reject
    #   if (log(u) < logr) {
    #     theta_sample_wel[i] <- theta_proposal
    #   } else {
    #     theta_sample_wel[i] <- theta_sample_wel[i - 1]
    #   }
    # }

    df_el <- data.frame(ss = 1:length(theta_sample_el), x = theta_sample_el, type = "EL")
    df_wel <- data.frame(ss = 1:length(theta_sample_wel), x = theta_sample_wel, type = "WEL")

    ggplot(rbind(df_el, df_wel), aes(x, fill = type, color = type)) +
      geom_histogram(aes(y = ..density..), position = "identity",
                     bins = 50, alpha = 0.2) +
      theme(axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, size = 1),
            panel.grid = element_blank(),
            legend.text = element_text(size = 10, color = "black"),
            legend.background = element_rect(fill = alpha("white", 0)),
            legend.key = element_rect(fill = alpha("white", 1)),
            legend.title = element_blank()) +
      labs(x = expression(theta)) +
      scale_color_npg()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
