library(shiny)
library(tidyverse)
library(melt)
library(ggplot2)
library(ggsci)
ui <- fluidPage(
  titlePanel("042222"),
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
                 selectInput(
                   "level", "level", c("95%", "99%", "99.9%"),
                   multiple = FALSE
                 ),
                 checkboxGroupInput("type", "methods",
                  choices = c("true" = "true",
                              "EL" = "el",
                              "ETEL" = "etel",
                              "WEL (w/o pseudo obs.)" = "wel",
                              "WETEL (w/o pseudo obs. & discard weight)" = "wetel",
                              "TEST (N(theta, 1) added to ETEL)" = "test",
                              "TEST2 (N(0, 1) added to ETEL)" = "test2"),
                  selected = c("true", "test", "test2"),
                                    inline = FALSE)),
               mainPanel(
                 textOutput("summary1"),
                 verbatimTextOutput("code1"),
                 plotOutput("plot1"),
                 textOutput("summary2"),
                 verbatimTextOutput("code2"),
                 plotOutput("plot2"),
                 plotOutput("plot3")
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
  level <- reactive(input$level)
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
        el_mean(par, x, control = control_el(th = 500))$optim$logLR
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
        x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = 1)
        pp <- el_mean(par, c(x, x_pseudo), w, control = control_el(th = 500))$log.prob
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
        x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = 1)
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
        const <- sum(exp(l_test * g)) +
          stc() * exp(l_test^2 / 2 - l_test * par)
        logp <- l_test * g - log(const)
        sum(logp) + n() * log(n())

        # g <- x - par
        # test_l <- function(l) sum(exp(l * g) * g) +
        #   stc() * l * exp(l^2 / 2 - l * par)
        # l_test <- uniroot(test_l, extendInt = "yes", lower = -1e+10,
        #                   upper = 1e+10)$root
        # const <- sum(exp(l_test * g)) +
        #   stc() * exp(l_test^2 / 2 - l_test * par)
        # logp <- l_test * g - log(const)
        # sum(logp) + n() * log(n() + stc()) +
        #   log(n() + stc()) * (l_test^2 / 2 - l_test * par)
      }
      df_test <- data.frame(x = theta_grid, y = sapply(theta_grid, test),
                            type = "TEST")
      dt <- rbind(dt, df_test)
    }
    if ("test2" %in% type()) {
      test <- function(par) {
        g <- x - par
        test_l <- function(l) sum(exp(l * g) * g) +
          stc() * l * exp(l^2 / 2)
        l_test <- uniroot(test_l, extendInt = "yes", lower = -1e+10,
                          upper = 1e+10)$root
        const <- sum(exp(l_test * g)) + stc() * exp(l_test^2 / 2)
        logp <- l_test * g - log(const)
        sum(logp) + n() * log(n())
      }
      df_test <- data.frame(x = theta_grid, y = sapply(theta_grid, test),
                            type = "TEST2")
      dt <- rbind(dt, df_test)
    }
    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("true", "EL", "ETEL", "WEL", "WETEL",
                                          "TEST", "TEST2"))
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

  output$summary2 <- renderText("EL Confidence Interval")
  output$code2 <- renderPrint({
    set.seed(seed())
    x <- rnorm(n(), mean = theta(), sd = sigma())
    if (level() == "95%") {
      l <- 0.95
    } else if (level() == "99%") {
      l <- 0.99
    } else {
      l <- 0.999
    }
    fit <- el_mean(0, x)
    names(fit$coefficients) <- "theta"
    confint(fit, level = l)
  }
  )
  output$plot2 <- renderPlot({
    # seed for sampling
    set.seed(seed())
    # simulate data
    x <- rnorm(n(), mean = theta(), sd = sigma())
    # generate grid
    theta_grid <- seq(mean(x) - 1 * sd(x),
                      mean(x) + 1 * sd(x), length.out = 500)
    # weights
    # beta <- 1 / m()
    beta <- strength() / m()
    w <- c(rep(1, n()), rep(beta, m()))
    w <- (length(w) / sum(w)) * w

    dt <- data.frame()
    # true logLR
    if ("true" %in% type()) {
      true <- function(par) {
        -2 * (sum(dnorm(x, mean = par, sd = sigma(), log = TRUE)) -
          sum(dnorm(x, mean = mean(x), sd = sigma(), log = TRUE)))
      }
      df_true <- data.frame(x = theta_grid, y = sapply(theta_grid, true),
                            type = "true")
      dt <- rbind(dt, df_true)
    }

    # EL logLR
    if ("el" %in% type()) {
      EL <- function(par) {
        -2 * el_mean(par, x, control = control_el(th = 500))$optim$logLR
      }
      df_EL <- data.frame(x = theta_grid, y = sapply(theta_grid, EL), type = "EL")
      dt <- rbind(dt, df_EL)
    }
    # ELEL logLR
    if ("etel" %in% type()) {
      ETEL <- function(par) {
        g <- x - par
        lambda_finder <- function(l) mean(exp(l * g) * g)
        lambda <- uniroot(lambda_finder, extendInt = "yes",
                          lower = -1e+10, upper = 1e+10)$root
        unnormalized_prob <- exp(lambda * g)
        log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
        - 2*(sum(log_prob) + n() * log(n()))
      }
      df_ETEL <- data.frame(x = theta_grid, y = sapply(theta_grid, ETEL),
                            type = "ETEL")
      dt <- rbind(dt, df_ETEL)
    }
    # WEL logLR
    if ("wel" %in% type()) {
      WEL <- function(par) {
        x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = 1)
        pp <- el_mean(par, c(x, x_pseudo), w,
                      control = control_el(th = 500))$log.prob
        p1 <- pp[seq_len(n())]
        p2 <- pp[-seq_len(n())]
        - 2*(sum(p1) + n() * (log(n() + strength())))
      }
      df_WEL <- data.frame(x = theta_grid, y = sapply(theta_grid, WEL),
                           type = "WEL")
      dt <- rbind(dt, df_WEL)
    }
    # WETEL logLR
    if ("wetel" %in% type()) {
      WETEL <- function(par) {
        x_pseudo <- par + qnorm(1:m() / (m() + 1), mean = 0, sd = 1)
        g <- c(x, x_pseudo) - par
        eval_f <- function(l) mean(w * exp(l * g))
        eval_grad_f <- function(l) mean(w * exp(l * g) * g)
        lambda_finder <- function(l) mean(w * exp(l * g) * g)
        lambda <- uniroot(lambda_finder, extendInt = "yes",
                          lower = -1e+10, upper = 1e+10)$root
        unnormalized_prob <- w * exp(lambda * g)
        log_prob <- log(unnormalized_prob) - log(sum(unnormalized_prob))
        -2 * (sum(log_prob[seq_len(n())]) + n() * log(n() + strength()))
      }
      df_WETEL <- data.frame(x = theta_grid,
                             y = sapply(theta_grid, WETEL),
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
        const <- sum(exp(l_test * g)) +
          stc() * exp(l_test^2 / 2 - l_test * par)
        logp <- l_test * g - log(const)
        sum(logp) + n() * log(n())
        - 2 * (sum(logp) + n() * log(n()))

        # g <- x - par
        # test_l <- function(l) sum(exp(l * g) * g) +
        #   stc() * l * exp(l^2 / 2 - l * par)
        # l_test <- uniroot(test_l, extendInt = "yes", lower = -1e+10,
        #                   upper = 1e+10)$root
        # # const <- sum(exp(l_test * g)) +
        # #   strength() * exp(l_test^2 / 2 - l_test * par)
        # const <- sum(exp(l_test * g))
        # logp <- l_test * g - log(const)
        # - 2 * (sum(logp) + n() * log(n()))
      }
      df_test <- data.frame(x = theta_grid, y = sapply(theta_grid, test),
                            type = "TEST")
      dt <- rbind(dt, df_test)
    }

    if ("test2" %in% type()) {
      test <- function(par) {
        g <- x - par
        test_l <- function(l) sum(exp(l * g) * g) +
          stc() * l * exp(l^2 / 2)
        l_test <- uniroot(test_l, extendInt = "yes", lower = -1e+10,
                          upper = 1e+10)$root
        const <- sum(exp(l_test * g)) + stc() * exp(l_test^2 / 2)
        logp <- l_test * g - log(const)
        # sum(logp) + n() * log(n())
        - 2 * (sum(logp) + n() * log(n()))
      }
      df_test <- data.frame(x = theta_grid, y = sapply(theta_grid, test),
                            type = "TEST2")
      dt <- rbind(dt, df_test)
    }

    if (level() == "95%") {
      cutoff <- qchisq(0.95, 1)
    } else if (level() == "99%") {
      cutoff <- qchisq(0.99, 1)
    } else {
      cutoff <- qchisq(0.999, 1)
    }

    # adjust factor levels for plot
    dt$type <- factor(dt$type, levels = c("true", "EL", "ETEL", "WEL", "WETEL",
                                          "TEST", "TEST2"))
    dt <- dt |> group_by(type) |> filter(y < cutoff) |>
      slice(which.min(x), which.max(x)) |> select(type, x) |>
      mutate(y = max(x)) |> slice(1)

    ggplot(dt, aes(x = x, y = type, color = type)) +
      geom_errorbar(aes(xmin = x, xmax = y), width = 0.2) +
      labs(x = expression(theta), y = NULL, title = "Confidence Intervals") +
      geom_vline(xintercept = mean(x), linetype = "dashed", alpha = 0.5) +
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
  output$plot3 <- renderPlot({
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
