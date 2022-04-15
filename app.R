library(shiny)
library(tidyverse)
library(melt)
library(ggplot2)
library(ggsci)
ui <- fluidPage(
  titlePanel("041522"),
  tabsetPanel(
    tabPanel("LR", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 sliderInput("theta", "sampling mean",
                             value = 0, min = ,-3, max = 3, step = 0.05),
                 sliderInput("sigma", "sampling sd", value = 1, min = 0.5, max = 3),
                 sliderInput("grid", "grid range", value = c(-1, 3), min = -5, max = 5,
                             step = 0.1),
                 numericInput("seed", "seed",
                              value = 1, min = 1, max = .Machine$integer.max),
                 sliderInput("n", "n: sample size", value = 50,
                             min = 10, max = 200, step = 10),
                 sliderInput("m", "m: pseudo sample size",
                             value = 2, min = 2, max = 100, step = 1),
                 sliderInput("strength", "s1: strength of pseudo sample",
                             value = 1, min = 1, max = 5000),
                 sliderInput("stc", "s2: strength of continuous mixing",
                             value = 1, min = 0.1, max = 100),
                 checkboxGroupInput("type", "methods",
                                    choices = c("true" = "true",
                                                "EL" = "el",
                                                "ETEL" = "etel",
                                                "WEL (w/o pseudo obs.)" = "wel",
                                                "WETEL (w/o pseudo obs. & discard weight)" = "wetel",
                                                "TEST (N(theta, 1) mixing added to ETEL)" = "test"),
                                    selected = c("true", "wetel", "test"),
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
  n <- reactive(input$n)
  m <- reactive(input$m)
  grid <- reactive(input$grid)
  seed <- reactive(input$seed)
  strength <- reactive(input$strength)
  stc <- reactive(input$stc)
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


    dt <- data.frame()
    # true logLR
    if ("true" %in% type()) {
      true <- function(par) {
        sum(dnorm(x, mean = par, sd = sigma(), log = TRUE)) -
          sum(dnorm(x, mean = mean(x), sd = sigma(), log = TRUE))
      }
      df_true <- data.frame(x = theta_grid, y = sapply(theta_grid, true),
                            type = "true")
      dt <- rbind(dt, df_true)
    }

    # EL logLR
    if ("el" %in% type()) {
      EL <- function(par) {
        el_mean(par, x, control = melt_control(th = 500))$optim$logLR
      }
      df_EL <- data.frame(x = theta_grid, y = sapply(theta_grid, EL), type = "EL")
      dt <- rbind(dt, df_EL)
    }
    # ELEL logLR
    if ("etel" %in% type()) {
      ETEL <- function(par) {
        g <- x - par
        # eval_f <- function(l) mean(exp(l * g))
        # eval_grad_f <- function(l) mean(exp(l * g) * g)
        # lambda <- nloptr(x0 = 0, eval_f = eval_f, eval_grad_f = eval_grad_f,
        #                  opts = opts)$solution
        lambda_finder <- function(l) mean(exp(l * g) * g)
        lambda <- uniroot(lambda_finder, extendInt = "yes",
                          lower = -1e+10, upper = 1e+10)$root

        unnormalized_prob <- exp(lambda * g)
        log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
        sum(log_prob) + n() * log(n())
      }
      df_ETEL <- data.frame(x = theta_grid, y = sapply(theta_grid, ETEL),
                            type = "ETEL")
      dt <- rbind(dt, df_ETEL)
    }
    # WEL logLR
    if ("wel" %in% type()) {
      WEL <- function(par) {
        # x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
        #                         max = IQR(x) / 2)
        x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = sd(x))
        pp <- el_mean(par, c(x, x_pseudo), w, control = melt_control(th = 500))$log.prob
        p1 <- pp[seq_len(n())]
        p2 <- pp[-seq_len(n())]
        sum(p1) + n() * (log(n() + strength()))
      }
      df_WEL <- data.frame(x = theta_grid, y = sapply(theta_grid, WEL),
                           type = "WEL")
      dt <- rbind(dt, df_WEL)
    }
    # WETEL logLR
    if ("wetel" %in% type()) {
      WETEL <- function(par) {
        # x_pseudo <- par + qunif(1:m() / (m() + 1), min = -IQR(x) / 2,
        #                         max = IQR(x) / 2)
        x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = sd(x))
        g <- c(x, x_pseudo) - par

        eval_f <- function(l) mean(w * exp(l * g))
        eval_grad_f <- function(l) mean(w * exp(l * g) * g)
        # lambda <- nloptr(x0 = 0, eval_f = eval_f, eval_grad_f = eval_grad_f,
        #                  opts = opts)$solution
        lambda_finder <- function(l) mean(w * exp(l * g) * g)
        lambda <- uniroot(lambda_finder, extendInt = "yes",
                          lower = -1e+10, upper = 1e+10)$root

        unnormalized_prob <- w * exp(lambda * g)
        log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
        sum(log_prob[seq_len(n())]) + n() * log(n() + strength())
      }
      df_WETEL <- data.frame(x = theta_grid, y = sapply(theta_grid, WETEL),
                             type = "WETEL")
      dt <- rbind(dt, df_WETEL)
    }

    # test logLR
    if ("test" %in% type()) {
      test <- function(par) {
        g <- x - par
        test_l <- function(l) sum(exp(l * g) * g) +
          stc() * l * exp(l^2 / 2 - l * par)
        l_test <- uniroot(test_l, extendInt = "yes", lower = -1e+10,
                          upper = 1e+10)$root
        # const <- sum(exp(l_test * g)) +
        #   strength() * exp(l_test^2 / 2 - l_test * par)
        const <- sum(exp(l_test * g))
        logp <- l_test * g - log(const)
        sum(logp) + n() * log(n())
      }
      df_test <- data.frame(x = theta_grid, y = sapply(theta_grid, test),
                            type = "TEST")
      dt <- rbind(dt, df_test)
    }
    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("true", "EL", "ETEL", "WEL", "WETEL",
                                          "TEST"))
    # logLR plot
    ggplot(dt,
           aes(x, y, color = type, linetype = type)) +
      geom_path(alpha = 0.5, na.rm = TRUE) +
      labs(x = expression(theta), y = "logLR") +
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


    # test logLR
    if ("test" %in% type()) {
      test <- function(par) {
        g <- x - par
        test_l <- function(l) sum(exp(l * g) * g) +
          stc() * l * exp(l^2 / 2 - l * par)
        l_test <- uniroot(test_l, extendInt = "yes", lower = -1e+10,
                          upper = 1e+10)$root
        l_test
      }
      dt <-  data.frame(x = theta_grid, y = sapply(theta_grid, test))
      ggplot(dt, aes(x, y)) +
        geom_path(alpha = 0.5, na.rm = TRUE) +
        labs(x = expression(theta), y = expression(lambda),
             title = "Monotonicity of lambda (exponential tilting)") +
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
    }

    # adjust factor levels for plot
    # dt$type <- factor(dt$type, levels = c("test"))
    #

  })
}

# Run the application
shinyApp(ui = ui, server = server)
