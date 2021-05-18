source('play/VaccToImmunity.R')

discretify <- function(dfun, tsum = 1, size = 20) {
  interm <- (0:(size - 1)) %>% sapply(function(x) integrate(dfun, x + 0.5, x + 1.5)$value)
  interm + (tsum - sum(interm)) / length(interm)
}

convolve.ts <- function(ts.hist, dist.discrete) {
  sum(tail(ts.hist, length(dist.discrete)) * dist.discrete)
}

convolve.sim <- function(ts.hist, dist.discrete, horizon) {
  ts.interm <- ts.hist
  for (n in 1:horizon) {
    ts.interm <- c(ts.interm, convolve.ts(ts.interm, dist.discrete))
  }
  ts.interm
}

combine.dists <- function(distAB, distBC) {
  # Key intuition: split AB flow into each day. On day N, distAB[N] of A
  # will follow distBC shifted by N+1.
  # 6th: AB[1] * BC[5] + AB[2] * BC[4] + AB[3] * BC[3] + ...
  # because weibull was fit to discrete distribution (and lim x->0 = Inf),
  # we use the the int_-Inf^Inf's discrete analog where all inputs are natural numbers
  distComb <- function(n) {
    sum(sapply(1:(n-1), function(x)
      distAB[x] * distBC[n - x]), na.rm = TRUE)
  }
  c(0, sapply(2:(length(distAB) + length(distBC)), distComb))
}

uncombine.dists <- function(distAB, distAC) {
  # Outputs BC
  # Differential example: BC[1] is probability that distAC - distAB = 1
  stopifnot(length(distAB) == length(distAC))
  dist.size <- length(distAB)
  BC.fun <- function (x) {
    # if x > 0, then start at 1, else 1 - x
    AB.params <- max(1 - x, 1):min(dist.size - x, dist.size)
    sum(sapply(AB.params, function (ab.param)
      distAC[ab.param] * distAB[ab.param - x]))
  }
  sapply((-dist.size+1):(dist.size-1), BC.fun)
}

optim_delay <- function(state0, leader, lagger) {
  sdlag <- function(x) {
    df <- indicators %>% filter(state == state0) %>% arrange(date)
    vec <- df[[leader]] / lead(df[[lagger]], round(x))
    quants <- quantile(vec[!is.na(vec)], c(.15, .85))
    quants[2] - quants[1]
  }
  
  try(plot(sapply(0:50, sdlag)))
  
  res <- optimize(sdlag, c(0, 50))
  return(res)
}

# https://www.medrxiv.org/content/10.1101/2020.09.04.20188516v1.full.pdf
initial_params = list(
  # Generation time (2020.09.04.20188516v1, supplementary materials)
  # Multiplied by Transmission Rate Rt
  SE = c(discretify(function(x)
    dweibull(x, shape = 3.2862, scale = 6.1244), size = 15), rep(0, 5)),
  # https://bmjopen.bmj.com/content/bmjopen/10/8/e039652.full.pdf
  EI = discretify(function(x)
    dlnorm(x, meanlog = 1.63, sdlog = 0.5)),
  # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0246772
  # https://github.com/JungsikNoh/COVID19_Estimated-Size-of-Infectious-Population/blob/main/output/state_summary/2021-04-15/regns2.csv
  
  # Interval from symptom onset to case confirmation:
  # Extrapolated from Belgium---may need state-specific weighted mean based on
  # age distribution
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589278/
  IQ = discretify(function(x)
    dweibull(x, shape = 0.75, scale = 5.657)), # (previously shape = 0.9)
  
  # Fixed delay of 9-12 days derived from SE
  # Multiply by proportion recovering independently (~99% of not hospitalized)
  # https://arxiv.org/pdf/2006.01283.pdf (Compendium)
  ER = discretify(function(x)
    dweibull(x - 9, shape = 1.5, scale = 3)),
  
  # https://www.thelancet.com/cms/10.1016/S1473-3099(20)30287-5/attachment/fde0cd38-f985-499f-9595-6578e367930c/mmc2.pdf
  IH = discretify(function(x)
    dlnorm(x, meanlog = 1.23, sdlog = 0.79)),
  
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01726-3#Sec8
  # https://annalsofintensivecare.springeropen.com/articles/10.1186/s13613-020-00749-6#Sec12
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589278/
  #   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589278/figure/ijerph-17-07560-f0A3/?report=objectonly
  HR = discretify(function(x)
    dlnorm(x, meanlog = 2.3, sdlog = 0.77), size = 40),
  HD = discretify(function(x)
    dweibull(x, shape = 1/0.842, scale = exp(2.4)), size = 40)
)

initial_params$EQ <- combine.dists(initial_params$EI, initial_params$IQ)
initial_params$EH <- combine.dists(initial_params$EI, initial_params$IH)
# Small limitation: Assuming hospitalization has
# no effect on distribution of duration between exposure and death, only proportion
# initial_params$QD satisfies EQ * QD = EH * HD
initial_params$QD <- discretify(function(x) dlnorm(x - 1.358, 1.662, 1.113), size=40)

search.deconvolve <- function(AB, AC, BC.pdf, BC.pdf.params, errfn = function(a, b) mean(abs(a - b))) {
  BC.length <- length(AC) - length(AB)
  param.error <- function(params) {
    BC.disc <- sapply(1:BC.length, function(x) do.call(BC.pdf, as.list(c(x, params))))
    AC.candidate <- combine.dists(AB, BC.disc)
    # MAE
    err <- errfn(AC.candidate, AC)
    # plot(BC.disc)
    err
  }
  optim(BC.pdf.params, param.error)
}

# parameters
# compendium value is 2 to 4%
prop.exposed.hosp <- 0.01
# https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2773971
# https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2769387/
prop.hosp.die <- 0.055
# https://www.cdc.gov/nchs/covid19/mortality-overview.htm
# 65.3% of all deaths hospitalized
prop.die.hosp <- 0.653
prop.exposed.die <- prop.exposed.hosp * prop.hosp.die / prop.die.hosp
prop.nonhosp.die <- prop.exposed.die * (1 - prop.die.hosp)
prop.exposed.detected <- 0.4

# state-specific IFR from age distribution
# https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
prop.exposed.die.states <- data.ages %>% ungroup %>%
  mutate(total = `(-1,17]` + `(17,49]` + `(49,64]` + `(64,200]`,
         prop.1 = `(-1,17]` / total,
         prop.2 = `(17,49]` / total,
         prop.3 = `(49,64]` / total,
         prop.4 = `(64,200]` / total,
         # ifr = (prop.1 * 14 / 1000000 +
         #          prop.2 * 350 / 1000000 +
         #          prop.3 * 2500 / 1000000 +
         #          prop.4 * 40000 / 1000000))
         ifr = (prop.1 * 20 / 1000000 +
           prop.2 * 500 / 1000000 +
           prop.3 * 6000 / 1000000 +
           prop.4 * 90000 / 1000000))

state.props <- function(abbr) {
  # TODO apply age distribution
  props <- initial_params
  props$ifr <- filter(prop.exposed.die.states, state == abbr)$ifr
  # Hospitalized death rate
  props$hdr <- prop.hosp.die
  # Nonhospitalized death rate
  props$ehr <- prop.exposed.hosp
  props$ndr <- (props$ifr - prop.exposed.hosp * prop.hosp.die) / 1
  props$total <- filter(prop.exposed.die.states, state == abbr)$total
  props
}

library(stringr)
run.simulation <- function(params, queues, compartments, Rt, Vt, start.date) {
  Q.SE <- queues$Q.SE
  Q.EQ <- queues$Q.EQ
  Q.ER <- queues$Q.ER
  Q.EH <- queues$Q.EH
  Q.HR <- queues$Q.HR
  Q.HD <- queues$Q.HD
  Q.QD <- queues$Q.QD
  S <- compartments$S
  E <- compartments$E
  Q <- compartments$Q
  H <- compartments$H
  D <- compartments$D
  R <- compartments$R
  Compartment.History <- tibble(
    # S=numeric(), E=numeric(), Q=numeric(), R=numeric(), H=numeric(), D=numeric(), date=as.Date(numeric())
    S=S,
    E=E,
    Q=Q,
    R=R,
    H=H,
    D=D,
    date=start.date-days(1)
  )
  the.date <- start.date
  
  for (i in 1:length(Rt)) {
    Ri <- Rt[i]
    V <- Vt[i]
    IPCi <- params$ipc[i]
    # Step: Generate deltas
    # Special case for QR.SE: QR.SE is Q.SE - Q.EQ - Q.EH - Q.ER
    dSE <- convolve.ts(Q.SE * 0.75, params$SE) * Ri * max(0, S) / params$total
    dER <- convolve.ts(Q.SE, params$ER) * (1 - prop.exposed.hosp) * (1 - params$ndr)
    dEQ <- convolve.ts(Q.SE, params$EQ) * (1 - prop.exposed.hosp) / IPCi
    dEH <- convolve.ts(Q.SE, params$EH) * prop.exposed.hosp
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7920817/
    # https://www.medrxiv.org/content/10.1101/2021.04.21.21255473v1.full-text
    dHR <- convolve.ts(Q.EH, params$HR) * (1 - prop.hosp.die)
    dHD <- convolve.ts(Q.EH, params$HD) * prop.hosp.die
    dQD <- convolve.ts(Q.EQ, params$QD) * params$ndr * IPCi
    dSR <- if (is.na(V)) 0 else V
    # Step: Apply deltas
    Q.SE <- append(Q.SE, dSE)
    #QR.SE <- append(QR.SE, dSE - dEQ - dEH) # rethink this
    Q.ER <- append(Q.ER, dER)
    Q.EQ <- append(Q.EQ, dEQ)
    Q.EH <- append(Q.EH, dEH)
    Q.HR <- append(Q.HR, dHR)
    Q.HD <- append(Q.HD, dHD)
    Q.QD <- append(Q.QD, dQD)
    S <- max(0, S - dSE - dSR)
    E <- E - dER - dEH + dSE
    Q <- Q - dQD + dEQ
    R <- R + dER + dHR + dSR
    H <- H - dHD - dHR + dEH
    D <- D + dHD + dQD
    the.date <- the.date + days(1)
    #print(str_glue("dSE={dSE} drSE={dSE - dEQ - dEH} dER={dER} dEH={dEH} dQD={dQD} dEQ={dEQ} dHR={dHR} dHD={dHD}"))
    Compartment.History <- Compartment.History %>%
      add_row(S=S, E=E, Q=Q, R=R, H=H, D=D, date=the.date)
  }
  list(history=Compartment.History, queues=list(
    Q.SE=Q.SE,
    Q.EQ=Q.EQ,
    Q.ER=Q.ER,
    Q.EH=Q.EH,
    Q.HR=Q.HR,
    Q.HD=Q.HD,
    Q.QD=Q.QD
  ))
}

library(lubridate)
library(memoise)

# We prime before vaccinations???
prime.simulation <- function(state_a, indicators, start.date) {
  start.date <- as.Date(start.date)
  end.date <- start.date + days(40)
  indicators <- indicators %>% group_by(state) %>% arrange(date) %>%
    mutate(cases.smooth = rollmean(cases, k=7, na.pad=TRUE),
           hosp.smooth = rollmean(hosp, k=7, na.pad=TRUE),
           deaths.smooth = rollmean(deaths, k=7, na.pad=TRUE))
  ind.all <- indicators %>% filter(state == state_a, date >= end.date)
  ind <- indicators %>% filter(state == state_a, date >= start.date, date <= end.date)
  # Find start H/D, end H/D, then fit shape of curve to get dEH - dHD - dHR
  # and dHD + dQD
  # dH + dD = dEH + dHR + dQD
  # dHD / dQD = 0.635 / (1 - 0.635)
  dH <- c(0, diff(ind$hosp.smooth))
  dD <- c(0, diff(ind$deaths.smooth))
  dHD <- prop.die.hosp * dD
  dQD <- (1 - prop.die.hosp) * dD
  dHR <- dH * (1 - prop.hosp.die) / 10 # Rough estimate as sum(params$HR[1:10]) ~ 0.54
  dEH <- dH + dHD + dHR
  # Median of 18 days ED; SE leads k*ED by 18 days
  # Flow & conservation breakdown: dED / prop.exposed.die = lead(dSE, 18)
  ind.later <- filter(indicators,
                      state == state_a,
                      date >= start.date + days(18),
                      date <= end.date + days(18)) %>% arrange(date)
  # dHL <- c(0, diff(ind.later$hosp.smooth))
  dEDL <- c(0, diff(ind.later$deaths))
  # dEHL <- dHL + prop.die.hosp * dDL + ind.later$hosp.smooth * (1 - prop.hosp.die) / 20
  dSE <- dEDL / state.props(state_a)$ifr * 3 # was 4
  dER <- dEH / prop.exposed.hosp * (1 - prop.exposed.hosp)
  dEQ <- c(0, diff(ind$cases.smooth))
  drSE <- dSE - dEQ - dEH
  props <- state.props(state_a)
  last.day <- ind[nrow(ind),]
  # infections per case
  props$ipc <- ind.all$deaths / ind.all$cases / props$ifr
  # props$rSE <- pmax(props$SE - pmax(initial_params$EQ[1:20] / props$ipc, initial_params$ER), 0)
  
  list(
    queues = list(
      Q.SE = dSE,
      Q.EQ = dEQ,
      Q.ER = dER,
      Q.EH = dEH,
      Q.HD = dHD,
      Q.HR = dHR,
      Q.QD = dQD),
    compartments = list(
      S = props$total - last.day$deaths / props$ifr - sum(dSE, na.rm = TRUE),
      E = sum(dSE, na.rm=TRUE),
      Q = last.day$cases,
      R = last.day$deaths / props$ifr - last.day$deaths,
      H = last.day$hosp,
      D = last.day$deaths),
    start.date = end.date,
    params = props
  )
}

fitting.end <- as.Date('2021-03-30')
horizon.size <- 45 # max 45

# run.simulation(pd$params, pd$Qi.SE, pd$QRi.SE, pd$Qi.EQ, pd$Qi.ER, pd$Qi.EH, pd$Qi.HR, pd$Qi.HD, pd$Qi.QD, pd$S, pd$E, pd$Q, pd$R, pd$H, pd$D, rep(1, 50), rep(0, 50), pd$start.date) -> res
prepare.indicators.can <- function(fitting.start) {
  fitting.start <- as.Date(fitting.start)
  immunization.start <- min(data.vac.owid$date) - days(28)
  prevac.vector <- rep(0, max(0, (immunization.start - fitting.start) / ddays(1)))
  state.ci.column <- function(state_abb) {
    state_name <- if (state_abb == 'DC') 'District of Columbia' else state.name[match(state_abb, state.abb)]
    print(state_name)
    if (state_name == 'New York') state_name <- 'New York State' # OWID oddities
    efc <- state.eff.curve(state_name, date.min = immunization.start, date.max = fitting.end + days(horizon.size))
    c(prevac.vector, efc$ci)
  }
  indicators.can %>%
    filter(date >= fitting.start, date <= fitting.end + days(horizon.size)) %>%
    filter((state %in% state.abb) | state == 'DC') %>%
    mutate(ci = state.ci.column(first(state)))
}

# library(optimParallel)
fit.simulation <- function(prime, indicators, sta, start.date, Rt) {
  step.size <- 1
  # cluster <- makeCluster(8)
  
  # Expected in signals: cases.smooth, hosp.smooth, deaths.smooth, ci
  fit.horizon <- function(prime, signals, Rpast, date) {
    Rt.point.slope <- function(slope_a) {
      intercept <- slope_a[1]
      slope <- slope_a[2]
      output <- rep(intercept, horizon.size)
      for (i in 1:length(output)) {
        output[i] <- output[i] + (i-1) * slope
      }
      #(0:(horizon.size-1)) * slope + Rpast
      # for (order in 1:(length(slope)-1)) {
      #   for (i in 1:length(output)) {
      #     output[i] <- output[i] + (i-1)^order * slope[order+1] / factorial(order)
      #   }
      # }
      output
    }
    
    Rt.slope.error <- function(slope, show = FALSE) {
      Rt <- Rt.point.slope(slope)
      res <- run.simulation(prime$params, prime$queues, prime$compartments,
                     Rt, signals$ci, date)
      smae <- function(a, b) {
        diff.mae <- mean(abs(diff(tail(a, horizon.size)) - diff(tail(b, horizon.size))))
        abs.mae <- mean(abs(tail(a, horizon.size) - tail(b, horizon.size)))
        diff.mae + abs.mae*3
      }
      #plot(ggplot(tibble(a=1:(horizon.size-1),
     #                    b=tail(res$history$Q, horizon.size) %>% diff,
      #                   q=tail(signals$cases.smooth, horizon.size) %>% diff)) +
      #  geom_line(aes(a, b)) + geom_line(aes(a,q), color='red'))
      
      # if (readline('continue?') == 'no') stop()
      cases.error <- smae(res$history$Q, signals$cases.smooth) * 0.004
      hosp.error <- smae(res$history$H, signals$hosp.smooth) * 0.13
      death.error <- smae(res$history$D, signals$deaths.smooth)
      # deaths roughly 0.02 of cases. hosps roughly 0.1 of cases.
      if (show)
        print(str_glue("hosp.error={hosp.error} death.error={death.error} cases.error={cases.error}"))
      composite.error <- hosp.error + death.error + cases.error
      if (is.nan(composite.error) || is.infinite(composite.error) || is.na(composite.error))
        1e10 else composite.error
    }
    #sapply((20:30), Rt.slope.error)
    optres <- optim(par=c(min(max(0.1, Rpast), 12), 0), Rt.slope.error)
                    # lower=c(0.1, 0.2), upper=c(10, 5), method = 'L-BFGS-B')
    sim <- run.simulation(prime$params, prime$queues, prime$compartments,
                   Rt.point.slope(optres$par), c(0, diff(signals$ci)), date)
    
    pstyle <- theme_gray(base_family = 'Linux Libertine O', base_size = 11)
    pcases <- ggplot(signals) + geom_line(aes(date, cases.smooth), color='red') +
      geom_line(data = sim$history, aes(date, Q)) +
      ylab('Total') +
      ggtitle('Cases') + pstyle
    phosp <- ggplot(signals) + geom_line(aes(date, hosp), color='red') +
      geom_line(data = sim$history, aes(date, H)) +
      ylab('Total') +
      ggtitle('Hospitalizations') + pstyle
    pdeaths <- ggplot(signals) + geom_line(aes(date, deaths.smooth), color='red') +
      geom_line(data = sim$history, aes(date, D)) +
      ylab('Total') +
      ggtitle('Deaths') + pstyle
    pdcases <- ggplot(signals) + geom_line(aes(date, c(NA, diff(cases.smooth))), color='red') +
      geom_line(data = sim$history, aes(date, diff(Q))) +
      ylab('Daily Increase') +
      ggtitle('Cases') + pstyle
    pdhosp <- ggplot(signals) + geom_line(aes(date, c(NA, diff(hosp))), color='red') +
      geom_line(data = sim$history, aes(date, diff(H))) +
      ylab('Daily Increase') +
      ggtitle('Hospitalizations') + pstyle
    pddeaths <- ggplot(signals) + geom_line(aes(date, c(NA, diff(deaths.smooth))), color='red') +
      geom_line(data = sim$history, aes(date, diff(D))) +
      ylab('Daily Increase') +
      ggtitle('Deaths') + pstyle
    grid.arrange(pcases, phosp, pdeaths, nrow=3)
    
    print(optres$par)
    Rt.slope.error(optres$par, show = TRUE)
    list(sim=sim, par=optres$par, Rt=Rt.point.slope(optres$par))
  }
  
  compartment.history <- tibble()
  Rt.history <- c()
  
  while (start.date <= fitting.end) {
    end.date <- min(max(indicators$date), start.date + days(horizon.size))
    fitting <- fit.horizon(prime, filter(indicators, state == sta,
                              date >= start.date, date <= end.date), Rt,
                           start.date)
    compartment.history <- rbind(compartment.history, fitting$sim$history %>% slice(2:(step.size+1)))
    Rt.history <- append(Rt.history, fitting$Rt[1:step.size])
    Rt <- fitting$Rt[step.size+1]
    prime$start.date <- start.date + days(step.size)
    prime$compartments <- as.list(slice(fitting$sim$history, step.size+1))
    prime$queues <- fitting$sim$queues
    prime$queues <- lapply(prime$queues, function(d) {
      length(d) <- length(d) - (horizon.size - step.size)
      d
    })
    prime$params$ipc <- tail(prime$params$ipc, length(prime$params$ipc) - step.size)
    print(str_glue("S={prime$compartments$S} Processed horizon period from {start.date} to {end.date} Rt={Rt}"))
    start.date <- start.date + days(step.size)
  }
  
  # stopCluster(cluster)
  compartment.history$Rt <- Rt.history
  list(history=compartment.history, Rt=Rt, prime=prime)
}

try.model <- function(state_abb, start.date, Rt=1.0) {
  ind.prep <- prepare.indicators.can(start.date - days(30))
  primed <- prime.simulation(state_abb, ind.prep, start.date)
  finfit <- fit.simulation(primed, ind.prep, state_abb, primed$start.date, Rt)
  diagplots <- ggplot(inner_join(finfit$history, filter(ind.prep, state == state_abb), by='date'), aes(x=date)) +
    geom_line(aes(y=cases.smooth), color='red') + geom_line(aes(y=Q))
  plot(diagplots)
  list(ind.prep, primed, finfit, state_abb)
}

library(gridExtra)
plot.diagnosis <- function(model.results) {
  global.theme <- theme_bw(base_family = 'Linux Libertine O') + theme(plot.background = element_rect('#f1f2f4'))
  joint.table <- model.results[[3]]$history %>%
    inner_join(filter(model.results[[1]], state == model.results[[4]]))
    
  plot.cases <- ggplot(joint.table) +
    geom_line(aes(date, Q)) +
    geom_line(aes(date, cases.smooth), color = 'red') +
    ggtitle('Fitted Cases') + global.theme
  
  plot.hosp <- ggplot(joint.table) +
    geom_line(aes(date, H)) +
    geom_line(aes(date, hosp), color = 'red') +
    ggtitle('Fitted Hospitalizations') + global.theme
  
  plot.deaths <- ggplot(joint.table) +
    geom_line(aes(date, D)) +
    geom_line(aes(date, deaths.smooth), color = 'red') +
    ggtitle('Fitted Deaths') + global.theme
  
  plot.Rt <- ggplot(joint.table) + geom_line(aes(date, Rt)) +
    ggtitle('Modeled Transmission Rate') + global.theme

  grid.arrange(plot.cases, plot.hosp, plot.deaths, plot.Rt, nrow=2, ncol=2)
}

continue.sim <- function(fit.output, ci, Rt) {
  run.simulation(fit.output[[2]]$params, fit.output[[2]]$queues, fit.output[[2]]$compartments,
                        Rt, ci, max(fit.output[[3]]$history$date) + days())
}

plot.sim <- function(sim.output) {
  plot.Q <- ggplot(sim.output$history, aes(date, Q)) + geom_line() + theme_bw()
  plot.S <- ggplot(sim.output$history, aes(date, S)) + geom_line() + theme_bw()
  plot.D <- ggplot(sim.output$history, aes(date, D)) + geom_line() + theme_bw()
  plot.E <- ggplot(sim.output$history, aes(date, E)) + geom_line() + theme_bw()
  plot.R <- ggplot(sim.output$history, aes(date, R)) + geom_line() + theme_bw()
  plot.H <- ggplot(sim.output$history, aes(date, H)) + geom_line() + theme_bw()
  
  
  grid.arrange(plot.Q, plot.S, plot.D, plot.E, plot.R, plot.H, nrow=2, ncol=3)
}

plot.sims <- function(sims, states) {
  output <- tibble()
  for (i in 1:length(states)) {
    sims[[i]]$history$state <- states[i]
    output <- rbind(output, sims[[i]]$history)
  }
  plot.Q <- ggplot(output, aes(date, Q, color = state)) + geom_line() + theme_bw()
  plot.S <- ggplot(output, aes(date, S, color = state)) + geom_line() + theme_bw()
  plot.D <- ggplot(output, aes(date, D, color = state)) + geom_line() + theme_bw()
  plot.E <- ggplot(output, aes(date, E, color = state)) + geom_line() + theme_bw()
  plot.R <- ggplot(output, aes(date, R, color = state)) + geom_line() + theme_bw()
  plot.H <- ggplot(output, aes(date, H, color = state)) + geom_line() + theme_bw()
  
  
  grid.arrange(plot.Q, plot.S, plot.D, plot.E, plot.R, plot.H, nrow=2, ncol=3)
}

# len.1 <- 300; plot.sims(list(continue.sim(try.ma.final, rep(1000, len.1), rep(1, len.1)), continue.sim(try.ma.final, rep(0, len.1), rep(1, len.1)), continue.sim(try.ma.final, rep(10000, len.1), rep(1, len.1))), c('1000', '0', '10000'))