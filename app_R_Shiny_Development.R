# app.R
library(shiny)
library(ggplot2)

# ---------- Theory helpers ----------

# Exact density of S(t) for Exp(lambda) arrivals and Exp(beta) jumps
# Compound Poisson with exponential jumps
# f_S(s; t, lambda, beta) = exp(-lambda t - beta s) * sqrt(lambda t beta / s) * I1(2 * sqrt(lambda t beta s)), s > 0
# Atom at 0: P(S(t) = 0) = exp(-lambda t)

f_S_density <- function(s, t, lambda, beta) {
  dens <- numeric(length(s))
  idx  <- which(s > 0)
  if (length(idx) > 0) {
    x <- s[idx]
    z <- 2 * sqrt(lambda * t * beta * x)
    term <- exp(-lambda * t - beta * x) *
      sqrt(lambda * t * beta / x) *
      besselI(z, nu = 1)
    dens[idx] <- term
  }
  dens
}

p0_mass <- function(t, lambda) exp(-lambda * t)

# Theoretical mean and variance for compound Poisson with Exp(beta) jumps
# E[S(t)]  = lambda t / beta
# Var[S(t)] = 2 lambda t / beta^2
S_mean <- function(t, lambda, beta) (lambda * t) / beta
S_var  <- function(t, lambda, beta) (2 * lambda * t) / (beta^2)

# Normal approximation at large t: N(mean, var)
f_normal_approx <- function(s, t, lambda, beta) {
  mu <- S_mean(t, lambda, beta)
  v  <- S_var(t, lambda, beta)
  dnorm(s, mean = mu, sd = sqrt(v))
}

# ---------- Simulation helpers ----------

# Simulate one path of S(t) on [0, T] over a given grid of times
simulate_path <- function(T, lambda, beta, time_grid) {
  # Number of arrivals by time T
  nmax <- rpois(1, lambda * T)
  if (nmax == 0) return(rep(0, length(time_grid)))
  
  # Arrival epochs via cumulative sum of exponentials
  interarrivals <- rexp(nmax, rate = lambda)
  arrivals <- cumsum(interarrivals)
  arrivals <- arrivals[arrivals <= T]
  k <- length(arrivals)
  if (k == 0) return(rep(0, length(time_grid)))
  
  jumps   <- rexp(k, rate = beta)
  cumsums <- cumsum(jumps)
  
  # Piecewise-constant right-continuous process
  S  <- numeric(length(time_grid))
  ai <- 1
  for (i in seq_along(time_grid)) {
    t <- time_grid[i]
    while (ai <= k && arrivals[ai] <= t) ai <- ai + 1
    if (ai == 1) {
      S[i] <- 0
    } else {
      S[i] <- cumsums[ai - 1]
    }
  }
  S
}

# Sample S(t*) by compounding: draw N ~ Pois(lambda t*), then sum N Exp(beta)
sample_S_at_t <- function(tstar, lambda, beta, n_sims) {
  N     <- rpois(n_sims, lambda * tstar)
  Svals <- numeric(n_sims)
  positive_idx <- which(N > 0)
  if (length(positive_idx) > 0) {
    for (j in positive_idx) {
      Svals[j] <- sum(rexp(N[j], rate = beta))
    }
  }
  Svals
}

# ---------- UI ----------

ui <- fluidPage(
  titlePanel("Compound Poisson S(t): sensitivity to exponential interarrivals (lambda) and jumps (beta)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "Interarrival rate lambda:", min = 0.1, max = 5,
                  value = 1, step = 0.1),
      sliderInput("beta",   "Jump rate beta:",          min = 0.1, max = 5,
                  value = 1, step = 0.1),
      sliderInput("T",      "Time horizon T:",          min = 10,  max = 1000,
                  value = 200, step = 10),
      sliderInput("nPaths", "Number of simulated paths:", min = 1, max = 200,
                  value = 50, step = 1),
      sliderInput("nGrid",  "Points in time grid:",     min = 100, max = 3000,
                  value = 1000, step = 100),
      hr(),
      sliderInput("tstar",  "Histogram time t*:",       min = 1, max = 1000,
                  value = 100, step = 1),
      sliderInput("nSims",  "Histogram samples:",       min = 1000, max = 50000,
                  value = 10000, step = 1000),
      checkboxInput("showExact",  "Overlay exact density", TRUE),
      checkboxInput("showNormal", "Overlay normal approximation", TRUE),
      checkboxInput("logY",       "Log y-scale for histogram density", FALSE),
      hr(),
      checkboxInput("showMoments", "Show theoretical mean/variance curves", TRUE),
      helpText("Tip: Increase lambda for more frequent jumps; decrease beta for larger jumps.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Trajectories",
                 plotOutput("trajPlot", height = "420px"),
                 verbatimTextOutput("trajStats")),
        tabPanel("Histogram at t*",
                 plotOutput("histPlot", height = "420px"),
                 verbatimTextOutput("histStats")),
        tabPanel("Moments over time",
                 plotOutput("momentsPlot", height = "420px"))
      )
    )
  )
)

# ---------- Server ----------

server <- function(input, output, session) {
  
  time_grid <- reactive(seq(0, input$T, length.out = input$nGrid))
  
  # Simulate paths
  paths <- reactive({
    tg     <- time_grid()
    lambda <- input$lambda
    beta   <- input$beta
    nP     <- input$nPaths
    
    Smat <- matrix(0, nrow = length(tg), ncol = nP)
    for (j in 1:nP) {
      Smat[, j] <- simulate_path(T = input$T, lambda = lambda, beta = beta,
                                 time_grid = tg)
    }
    list(t = tg, S = Smat)
  })
  
  output$trajPlot <- renderPlot({
    pr <- paths()
    df <- data.frame(
      t    = rep(pr$t, times = ncol(pr$S)),
      S    = as.vector(pr$S),
      path = factor(rep(seq_len(ncol(pr$S)), each = length(pr$t)))
    )
    
    ggplot(df, aes(x = t, y = S, group = path, color = path)) +
      geom_line(linewidth = 0.7, alpha = 0.7, show.legend = FALSE) +
      labs(
        title = "Simulated trajectories of S(t)",
        x     = "Time t",
        y     = "S(t)"
      ) +
      theme_minimal()
  })
  
  output$trajStats <- renderText({
    lambda <- input$lambda
    beta   <- input$beta
    T      <- input$T
    
    muT  <- S_mean(T, lambda, beta)
    varT <- S_var(T, lambda, beta)
    p0   <- p0_mass(T, lambda)
    
    paste0(
      "Theoretical at T = ", T, ":\n",
      "  Mean E[S(T)] = ", round(muT, 4), "\n",
      "  Var Var[S(T)] = ", round(varT, 4), "\n",
      "  Atom at zero P{S(T)=0} = ", format(p0, digits = 4)
    )
  })
  
  # Histogram at t*
  samples <- reactive({
    sample_S_at_t(input$tstar, input$lambda, input$beta, input$nSims)
  })
  
  output$histPlot <- renderPlot({
    svals  <- samples()
    tstar  <- input$tstar
    lambda <- input$lambda
    beta   <- input$beta
    
    df <- data.frame(S = svals)
    
    p <- ggplot(df, aes(x = S)) +
      geom_histogram(aes(y = ..density..),
                     bins = 50,
                     fill = "#4C78A8",
                     color = "white",
                     alpha = 0.9) +
      labs(
        title = paste0("Histogram of S(t*) at t* = ", tstar, " (density scale)"),
        x     = "S(t*)",
        y     = ifelse(input$logY, "log-density", "density")
      ) +
      theme_minimal()
    
    # Overlay grids
    sgrid_max <- max(
      max(svals),
      S_mean(tstar, lambda, beta) + 6 * sqrt(S_var(tstar, lambda, beta))
    )
    sgrid <- seq(1e-6, sgrid_max, length.out = 1000)
    
    if (input$showExact) {
      dens_exact <- f_S_density(sgrid, tstar, lambda, beta)
      p <- p + geom_line(aes(x = sgrid, y = dens_exact),
                         linewidth = 1.1,
                         color = "#F58518")
    }
    
    if (input$showNormal) {
      dens_norm <- f_normal_approx(sgrid, tstar, lambda, beta)
      p <- p + geom_line(aes(x = sgrid, y = dens_norm),
                         linewidth = 1.1,
                         linetype = "dashed",
                         color = "#54A24B")
    }
    
    if (input$logY) {
      p <- p + scale_y_log10()
    }
    
    p
  })
  
  output$histStats <- renderText({
    tstar  <- input$tstar
    lambda <- input$lambda
    beta   <- input$beta
    
    mu <- S_mean(tstar, lambda, beta)
    v  <- S_var(tstar, lambda, beta)
    p0 <- p0_mass(tstar, lambda)
    
    paste0(
      "Theoretical at t* = ", tstar, ":\n",
      "  Mean E[S(t*)] = ", round(mu, 4), "\n",
      "  Var Var[S(t*)] = ", round(v, 4), "\n",
      "  Atom at zero P{S(t*)=0} = ", format(p0, digits = 4)
    )
  })
  
  # Moments over time
  output$momentsPlot <- renderPlot({
    # If checkbox is off, don't draw anything
    if (!input$showMoments) return(NULL)
    
    lambda <- input$lambda
    beta   <- input$beta
    tgrid  <- seq(0, input$T, length.out = 300)
    
    df <- data.frame(
      t     = rep(tgrid, 2),
      value = c(S_mean(tgrid, lambda, beta),
                S_var(tgrid,  lambda, beta)),
      kind  = factor(rep(c("Mean", "Variance"),
                         each = length(tgrid)))
    )
    
    ggplot(df, aes(x = t, y = value, color = kind)) +
      geom_line(linewidth = 1.1) +
      labs(
        title = "Theoretical mean and variance of S(t)",
        x     = "Time t",
        y     = "Value"
      ) +
      scale_color_manual(values = c("Mean" = "#E45756",
                                    "Variance" = "#72B7B2")) +
      theme_minimal()
  })
}

shinyApp(ui, server)