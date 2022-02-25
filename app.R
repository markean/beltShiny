library(shiny)
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
                             value = c(-1, 1), min = -4, max = 4, step = 0.05),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 100,
                             min = 10, max = 500, step = 10),
                 sliderInput("m", "m: pseudo sample size",
                             value = 2, min = 2, max = 100, step = 1),
                 selectInput("beta", "beta (fractional weight for each pseudo observation)",
                             choices = c("1 / m", "m / n", "1")),
                 #             round = F, step = 0.05),
                 # selectInput("method", "pseudo sample method:",
                 #             choices = c("qnorm" = "qnorm", "rnorm" = "rnorm")),
                 checkboxGroupInput("wel", "weighted EL to use",
                                    choices = c("WEL1 (w/o pseudo obs. & discard weights)" = "wel1",
                                                "WEL2 (w/ pseudo obs. & retain weights)" = "wel2"),
                                    selected = c("wel1"),
                                    inline = FALSE)),

               mainPanel(
                 # tabsetPanel(
                 #   tabPanel("Summary", verbatimTextOutput("code"), plotOutput("plot1"),
                 #            plotOutput("plot2"))
                 # )
                 textOutput("summary1"),
                 verbatimTextOutput("code1"),
                 plotOutput("plot1"),
                 plotOutput("plot2")
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
                 sliderInput("theta_mcmc", "sampling mean",
                             value = 0, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma_mcmc", "sampling sd",
                             value = 1, min = 0.5, max = 3),
                 sliderInput("theta0_mcmc", "prior mean",
                             value = -1, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma0_mcmc", "prior sd",
                             value = 1, min = 0.5, max = 3),
                 numericInput("seed_mcmc", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n_mcmc", "n: sample size", value = 100,
                             min = 10, max = 500, step = 10),
                 sliderInput("m_mcmc", "m: pseudo sample size",
                             value = 2, min = 2, max = 100, step = 1),
                 selectInput("beta_mcmc", "beta (fractional weight for each pseudo observation)",
                             choices = c("1 / m", "m / n", "1")),
                 selectInput("wel_mcmc", "weighted EL to use",
                             choices = c("WEL1 (w/o pseudo obs. & discard weights)" = "wel1",
                                         "WEL2 (w/ pseudo obs. & retain weights)" = "wel2"))),
               mainPanel(
                 textOutput("summary2"),
                 verbatimTextOutput("code2"),
                 plotOutput("plot3"),
                 plotOutput("plot4")
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
  # plotlim <- reactive(input$plotlim)
  seed <- reactive(input$seed)
  beta <- reactive(input$beta)
  # method <- reactive(input$method)
  wel <- reactive(input$wel)

  # parameters for MCMC
  B <- reactive(input$B)
  sigma_p <- reactive(input$sigma_p)
  theta_mcmc <- reactive(input$theta_mcmc)
  sigma_mcmc <- reactive(input$sigma_mcmc)
  theta0_mcmc <- reactive(input$theta0_mcmc)
  sigma0_mcmc <- reactive(input$sigma0_mcmc)
  seed_mcmc <- reactive(input$seed_mcmc)
  beta_mcmc <- reactive(input$beta_mcmc)
  n_mcmc <- reactive(input$n_mcmc)
  m_mcmc <- reactive(input$m_mcmc)
  wel_mcmc <- reactive(input$wel_mcmc)

  # summary statistics
  output$summary1 <- renderText("Summary of sample data")
  output$code1 <- renderPrint({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    summary(x)
    }
  )
  output$summary2 <- renderText("Summary of sample data")
  output$code2 <- renderPrint({
    set.seed(seed_mcmc())
    x <- rnorm(n_mcmc(), mean = theta_mcmc(), sd = sigma_mcmc())
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
        # pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        # while (par <= range(pseudo_x)[1] || par >= range(pseudo_x)[2]) {
        #   pseudo_x <- rnorm(m(), mean = par, sd = sigma0())
        # }
      par + qnorm(1:m() / (m() + 1))
    }
    # logLR
    # logLR <- function(par) {
    #   pseudo_x <- pseudo_sample(par)
    #   data_f <- c(x, pseudo_x)
    #   sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
    #                                    abstol = 1e-06))$optim$log.prob[1:n()])
    # }
    logWLR <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      # sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
      #                                  abstol = 1e-06))$optim$log.wprob[1:n()])
      sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                       abstol = 1e-06))$optim$log.prob[1:n()]) +
        n() * (log(n()) - log(n() + m()))
    }
    # logLR with augmented data
    # logLR_aug <- function(par) {
    #   pseudo_x <- pseudo_sample(par)
    #   data_f <- c(x, pseudo_x)
    #   el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
    #                                abstol = 1e-06))$optim$logLR
    # }
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
    dt <- df_EL
    # WEL1 density
    if ("wel1" %in% wel()) {
      set.seed(seed())
      WEL1_log_d <- sapply(theta_grid, logWLR)
      WEL1_d <- exp(WEL1_log_d) / sum(exp(WEL1_log_d))
      df_WEL1 <- data.frame(x = theta_grid, y = WEL1_log_d, type = "WEL1")
      dt <- rbind(dt, df_WEL1)
    }
    # WEL2 density
    if ("wel2" %in% wel()) {
      set.seed(seed())
      WEL2_log_d <- sapply(theta_grid, logWLR_aug)
      WEL2_d <- exp(WEL2_log_d) / sum(exp(WEL2_log_d))
      df_WEL2 <- data.frame(x = theta_grid, y = WEL2_log_d, type = "WEL2")
      dt <- rbind(dt, df_WEL2)
    }

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
            legend.title = element_blank())
      # scale_color_npg()
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
      par + qnorm(1:m() / (m() + 1))
    }
    # logLR
    # logLR <- function(par) {
    #   pseudo_x <- pseudo_sample(par)
    #   data_f <- c(x, pseudo_x)
    #   sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
    #                                    abstol = 1e-06))$optim$log.prob[1:n()])
    # }
    logWLR <- function(par) {
      pseudo_x <- pseudo_sample(par)
      data_f <- c(x, pseudo_x)
      # sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
      #                                  abstol = 1e-06))$optim$log.wprob[1:n()])
      sum(el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
                                       abstol = 1e-06))$optim$log.prob[1:n()]) +
        n() * (log(n()) - log(n() + m()))
    }
    # logLR with augmented data
    # logLR_aug <- function(par) {
    #   pseudo_x <- pseudo_sample(par)
    #   data_f <- c(x, pseudo_x)
    #   el_mean(par, data_f, w, list(threshold = 1e+10, maxit = 100,
    #                                abstol = 1e-06))$optim$logLR
    # }
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
    df_true <- data.frame(x = theta_grid, y = true_pd, type = "True")
    df_EL <- data.frame(x = theta_grid, y = EL_pd, type = "EL")
    dt <- rbind(df_true, df_EL)
    # posterior density with real data (discrete, normalized)
    if ("wel1" %in% wel()) {
      set.seed(seed())
      WEL1_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logWLR)
      WEL1_pd <- exp(WEL1_log_pd) / sum(exp(WEL1_log_pd))
      df_WEL1 <- data.frame(x = theta_grid, y = WEL1_pd, type = "WEL1")
      dt <- rbind(dt, df_WEL1)
    }
    if ("wel2" %in% wel()) {
      set.seed(seed())
      WEL2_log_pd <- dnorm(theta_grid, theta0(), sigma0(), log = T) +
        sapply(theta_grid, logWLR_aug)
      WEL2_pd <- exp(WEL2_log_pd) / sum(exp(WEL2_log_pd))
      df_WEL2 <- data.frame(x = theta_grid, y = WEL2_pd, type = "WEL2")
      dt <- rbind(dt, df_WEL2)
    }
    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("EL", "WEL1", "WEL2", "True"))

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
  # MCMC samples
  output$plot3 <- renderPlot({
    # seed for sampling
    set.seed(seed_mcmc())
    # simulate data
    x <- rnorm(n_mcmc(), mean = theta_mcmc(), sd = sigma_mcmc())
    # weights
    if (beta_mcmc() == "1 / m") {
      beta <- 1 / m_mcmc()
    } else if (beta_mcmc() == "m / n") {
      beta <- m_mcmc() / n_mcmc()
    } else {
      beta <- 1
    }
    w <- c(rep(1, n_mcmc()), rep(beta, m_mcmc()))

    pseudo_sample <- function(par) {
      par + qnorm(1:m_mcmc() / (m_mcmc() + 1))
    }
    logWLR <- function(par) {
      x_aug <- c(x, pseudo_sample(par))
      sum(el_mean(par, x_aug, w, list(threshold = 1e+10, maxit = 100,
                                      abstol = 1e-06))$optim$log.prob[1:n_mcmc()]) +
        n_mcmc() * (log(n_mcmc()) - log(n_mcmc() + m_mcmc()))
    }
    logWLR_aug <- function(par) {
      x_aug <- c(x, pseudo_sample(par))
      el_mean(par, x_aug, w, list(threshold = 1e+10, maxit = 100,
                                  abstol = 1e-06))$optim$logWLR
    }

    theta_sample_el <- vector("numeric", B())
    theta_sample_el[1] <- theta0_mcmc()
    theta_sample_wel <- vector("numeric", B())
    theta_sample_wel[1] <- theta0_mcmc()

    for (i in 2:B()) {
      # sample proposal value
      theta_proposal <- rnorm(1, theta_sample_el[i - 1], sigma_p())

      # compute log ratio posterior densities
      LR1 <- el_mean(theta_proposal, x)$optim$logLR
      LR2 <- el_mean(theta_sample_el[i - 1], x)$optim$logLR
      logr <-
        (log(dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc())) -
           log(dnorm(theta_sample_el[i - 1], theta0_mcmc(), sigma0_mcmc()))) +
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
      # x_aug <- c(x, pseudo_sample(theta_proposal))
      # compute log ratio posterior densities
      if (wel_mcmc() == "wel1") {
        LR1 <- logWLR(theta_proposal)
        LR2 <- logWLR(theta_sample_wel[i - 1])
      } else {
        LR1 <- logWLR_aug(theta_proposal)
        LR2 <- logWLR_aug(theta_sample_wel[i - 1])
      }
      # LR1 <- sum(el_mean(theta_proposal, x_aug, w)$optim$log.prob[1:n_mcmc()])
      # LR2 <- sum(el_mean(theta_sample_wel[i - 1], x_aug, w)$optim$log.prob[1:n_mcmc()])
      logr <-
        (log(dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc())) -
           log(dnorm(theta_sample_wel[i - 1], theta0_mcmc(), sigma0_mcmc()))) +
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

    df_el <- data.frame(x = seq_len(B()), y = theta_sample_el,
                        type = "EL")
    df_wel <- data.frame(x = seq_len(B()), y = theta_sample_wel,
                         type = toupper(wel_mcmc()))

    ggplot(rbind(df_el, df_wel), aes(x, y, color = type)) +
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
    set.seed(seed_mcmc())
    # simulate data
    x <- rnorm(n_mcmc(), mean = theta_mcmc(), sd = sigma_mcmc())
    # weights
    if (beta_mcmc() == "1 / m") {
      beta <- 1 / m_mcmc()
    } else if (beta_mcmc() == "m / n") {
      beta <- m_mcmc() / n_mcmc()
    } else {
      beta <- 1
    }
    w <- c(rep(1, n_mcmc()), rep(beta, m_mcmc()))

    pseudo_sample <- function(par) {
      par + qnorm(1:m_mcmc() / (m_mcmc() + 1))
    }
    logWLR <- function(par) {
      x_aug <- c(x, pseudo_sample(par))
      sum(el_mean(par, x_aug, w, list(threshold = 1e+10, maxit = 100,
                                      abstol = 1e-06))$optim$log.prob[1:n_mcmc()]) +
        n_mcmc() * (log(n_mcmc()) - log(n_mcmc() + m_mcmc()))
    }
    logWLR_aug <- function(par) {
      x_aug <- c(x, pseudo_sample(par))
      el_mean(par, x_aug, w, list(threshold = 1e+10, maxit = 100,
                                  abstol = 1e-06))$optim$logWLR
    }

    theta_sample_el <- vector("numeric", B())
    theta_sample_el[1] <- theta0_mcmc()
    theta_sample_wel <- vector("numeric", B())
    theta_sample_wel[1] <- theta0_mcmc()

    for (i in 2:B()) {
      # sample proposal value
      theta_proposal <- rnorm(1, theta_sample_el[i - 1], sigma_p())

      # compute log ratio posterior densities
      LR1 <- el_mean(theta_proposal, x)$optim$logLR
      LR2 <- el_mean(theta_sample_el[i - 1], x)$optim$logLR
      logr <-
        (log(dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc())) -
           log(dnorm(theta_sample_el[i - 1], theta0_mcmc(), sigma0_mcmc()))) +
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
      if (wel_mcmc() == "wel1") {
        LR1 <- logWLR(theta_proposal)
        LR2 <- logWLR(theta_sample_wel[i - 1])
      } else {
        LR1 <- logWLR_aug(theta_proposal)
        LR2 <- logWLR_aug(theta_sample_wel[i - 1])
      }
      logr <-
        (log(dnorm(theta_proposal, theta0_mcmc(), sigma0_mcmc())) -
           log(dnorm(theta_sample_wel[i - 1], theta0_mcmc(), sigma0_mcmc()))) +
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

    df_el <- data.frame(x = seq_len(B()), y = theta_sample_el, type = "EL")
    df_wel <- data.frame(x = seq_len(B()), y = theta_sample_wel,
                         type = toupper(wel_mcmc()))

    # posterior mean
    posterior_mean <-
      (theta0_mcmc() / sigma0_mcmc()^2 + sum(x) / sigma_mcmc()^2) /
      (1 / sigma0_mcmc()^2 + n_mcmc() / sigma_mcmc()^2)
    # posterior sd
    posterior_sd <- (1 / sigma0_mcmc()^2 + n_mcmc() / sigma_mcmc()^2)^(-1 / 2)

    ggplot(rbind(df_el, df_wel), aes(y, fill = type, color = type)) +
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
      stat_function(fun = dnorm,
                    args = list(mean = posterior_mean,
                                sd = posterior_sd)) +
      scale_color_npg()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
