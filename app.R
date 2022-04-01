library(shiny)
library(tidyverse)
library(plotly)
# library(nloptr)
library(melt)
library(ggplot2)
library(ggsci)
ui <- fluidPage(
  titlePanel("040121"),
  tabsetPanel(
    tabPanel("proportion", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 sliderInput("par", "theta",
                             value = 10, min = ,-10, max = 10, step = 0.1),
                 sliderInput("theta", "sampling mean",
                             value = -3, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma", "sampling sd",
                             value = 1, min = 0.5, max = 3),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 100,
                             min = 10, max = 200, step = 10),
                 sliderInput("m", "m: pseudo sample size",
                             value = 2, min = 2, max = 200, step = 1),
                 sliderInput("strength", "c: strength of pseudo sample",
                             value = 1, min = 1, max = 1000)),

               mainPanel(
                 textOutput("summary1"),
                 verbatimTextOutput("code1"),
                 textOutput("summary2"),
                 tableOutput("code2"),
                 textOutput("summary3"),
                 tableOutput("code3"),
                 textOutput("summary4"),
                 tableOutput("code4"),
                 textOutput("summary5"),
                 tableOutput("code5"),
               )
             )
    ),
    tabPanel("Logistic Regression", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 sliderInput("b0", "beta0",
                             value = 1, min = ,-3, max = 3, step = 0.1),
                 sliderInput("b1", "beta1",
                             value = 1, min = -3, max = 3, step = 0.1),
                 sliderInput("mean_lm", "x mean",
                             value = 0, min = ,-3, max = 3, step = 0.1),
                 sliderInput("sd_lm", "x sd",
                             value = 1, min = 0.5, max = 3, step = 0.1),
                 sliderInput("grid_lm", "grid length",
                             value = 10, min = 1, max = 20, step = 0.1),
                 sliderInput("logLR_threshold_lm", "logLR threshold",
                             value = -50, min = -500, max = -1, step = 1),
                 numericInput("seed_lm", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n_lm", "n: sample size", value = 100,
                             min = 10, max = 200, step = 1)),
               mainPanel(
                 plotlyOutput("plot1"),
                 plotOutput("plot2")
               )
             )
    )
  )
)

server <- function(input, output) {
  options(warn = -1)
  # parameters for densities
  par <- reactive(input$par)
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
  output$summary2 <- renderText("WETEL (real observations)")
  output$code2 <- renderTable({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # weights
    # beta <- 1 / m()
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w

    # pseudo obervations
    # x_pseudo <- par() + qnorm(seq_len(m()) / (m() + 1), sd = sd(x))
    x_pseudo <- par() + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
    g <- c(x, x_pseudo) - par()
    lambda_finder <- function(l) mean(w * exp(l * g) * g)
    lambda <- uniroot(lambda_finder, extendInt = "yes",
                      lower = -1e+10, upper = 1e+10)$root

    unnormalized_prob <- w * exp(lambda * g)
    log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
    prob <- exp(log_prob)
    p_real <- prob[seq_len(n())]

    percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
    }
    aa <- data.frame("sum.pi" = sum(p_real),
                     "max.pi" = max(p_real),
                     "max.pi proportion" = percent(max(p_real) / sum(p_real)),
                     "logLR" = sum(log_prob[seq_len(n())]) + n() * log(n() + strength()))
    aa
  }
  )
  output$summary3 <- renderText("WETEL (pseudo observations)")
  output$code3 <- renderTable({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # weights
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w

    # pseudo obervations
    # x_pseudo <- par() + qnorm(seq_len(m()) / (m() + 1), sd = sd(x))
    x_pseudo <- par() + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
    g <- c(x, x_pseudo) - par()
    lambda_finder <- function(l) mean(w * exp(l * g) * g)
    lambda <- uniroot(lambda_finder, extendInt = "yes",
                      lower = -1e+10, upper = 1e+10)$root

    unnormalized_prob <- w * exp(lambda * g)
    log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
    prob <- exp(log_prob)
    p_real <- prob[seq_len(n())]
    p_pseudo <- prob[-seq_len(n())]
    percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
    }

    aa <- data.frame("sum.pi" = sum(p_pseudo),
                     "max.pi" = max(p_pseudo),
                     "max.pi.proportion" = percent(max(p_pseudo) /
                                                     sum(p_pseudo)))
    aa
  }
  )
  output$summary4 <- renderText("WEL (real observations)")
  output$code4 <- renderTable({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # weights
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w

    # pseudo obervations
    # x_pseudo <- par() + qnorm(seq_len(m()) / (m() + 1), sd = sd(x))
    x_pseudo <- par() + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)

    g <- c(x, x_pseudo) - par()
    pp <- el_eval(g, w, control = list(tol = 1e-05, th = 1e+10))
    prob <- exp(pp$log.prob)
    p_real <- prob[seq_len(n())]

    percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
    }

    aa <- data.frame("sum.pi" = sum(p_real),
                     "max.pi" = max(p_real),
                     "max.pi.proportion" = percent(max(p_real) /
                                                     sum(p_real)),
                     "logLR" = sum(pp$log.prob[seq_len(n())]) + n() * log(n() + strength()))
    aa
  }
  )
  output$summary5 <- renderText("WEL (pseudo observations)")
  output$code5 <- renderTable({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # weights
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w

    # pseudo obervations
    # x_pseudo <- par() + qnorm(seq_len(m()) / (m() + 1), sd = sd(x))
    x_pseudo <- par() + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
                              max = IQR(x) / 2)
    g <- c(x, x_pseudo) - par()
    pp <- el_eval(g, w, control = list(tol = 1e-05, th = 1e+10))
    prob <- exp(pp$log.prob)

    p_real <- prob[seq_len(n())]
    p_pseudo <- prob[-seq_len(n())]
    percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
    }
    aa <- data.frame("sum.pi" = sum(p_pseudo),
                     "max.pi" = max(p_pseudo),
                     "max.pi.proportion" = percent(max(p_pseudo) /
                                                     sum(p_pseudo)))
    aa
  }
  )



  # parameters for logistic regression
  b0 <- reactive(input$b0)
  b1 <- reactive(input$b1)
  mean_lm <- reactive(input$mean_lm)
  sd_lm <- reactive(input$sd_lm)
  n_lm <- reactive(input$n_lm)
  m_lm <- reactive(input$m_lm)
  grid_lm <- reactive(input$grid_lm)
  logLR_threshold_lm <- reactive(input$logLR_threshold_lm)
  seed_lm <- reactive(input$seed_lm)

  output$plot1 <- renderPlotly({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n_lm(), mean = mean_lm(), sd = sd_lm())

    l <- b0() + b1() * x
    mu <- 1 / (1 + exp(-l))
    ff <- function(x) {
      sample(c(1, 0), 1, replace = F, prob = c(x, 1 - x))
    }
    y <- sapply(mu, ff)
    df <- data.frame(y, x)

    # extract estimate
    fit <- el_glm(y ~., df, control = list(th = 1e+3, maxit = 300),
                  family = binomial)
    est <- fit$coefficients


    lhs <- matrix(c(1, 0, 0, 1), nrow = 2)
    grid <- tibble(
      x = seq(est[1] - grid_lm() / 2, est[1] + grid_lm() / 2, length.out = 70),
      y = seq(est[2] - grid_lm() / 2, est[2] + grid_lm() / 2, length.out = 70)
    )

    # EL
    z <- matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
    for (i in seq_len(nrow(z))) {
      for (j in seq_len(nrow(z))) {
        z[i, j] <- lht(fit, lhs = lhs, rhs = c(grid$x[i], grid$y[j]),
                       control = list(th = 1e+3, maxit = 300))$optim$logLR
      }
    }
    plot_ly(x = ~grid$x, y = ~grid$y, z = ~t(z)) %>%
      layout(title = "EL",
             scene = list(
               xaxis = list(title = "b0"),
               yaxis = list(title = "b1"),
               zaxis = list(title = "logLR", range = c(logLR_threshold_lm(), 0))
             )) %>%
      add_surface()


  })
  output$plot2 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n_lm(), mean = mean_lm(), sd = sd_lm())

    l <- b0() + b1() * x
    mu <- 1 / (1 + exp(-l))
    ff <- function(x) {
      sample(c(1, 0), 1, replace = F, prob = c(x, 1 - x))
    }
    y <- sapply(mu, ff)
    df <- data.frame(y, x)

    # extract estimate
    fit <- el_glm(y ~., df, control = list(th = 1e+3, maxit = 300),
                  family = binomial)
    est <- fit$coefficients


    lhs <- matrix(c(1, 0, 0, 1), nrow = 2)



    grid2 <- expand_grid(
      x = seq(est[1] - grid_lm() / 2, est[1] + grid_lm() / 2, length.out = 70),
      y = seq(est[2] - grid_lm() / 2, est[2] + grid_lm() / 2, length.out = 70))
    grid2$value <- NA

    for (i in seq_len(nrow(grid2))) {
      grid2$value[i] <- lht(fit, lhs = lhs, rhs = as.numeric(grid2[i, -3]),
                            control = list(th = 1e+3, maxit = 300))$optim$logLR
    }
    ggplot(grid2, aes(x, y, z = value)) +
      geom_contour_filled(bins = 20) +
      geom_point(data = data.frame(x = est[1], y = est[2], value = NA),
                 aes(x, y)) +
      labs(x = "b0", y = "b1")
  })

}

# Run the application
shinyApp(ui = ui, server = server)
