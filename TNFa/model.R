## Library dependencies and plot theme --------------------------------------------

library(dMod) # Version 30 Oct 2017

runscript <- FALSE

## Model Definition ------------------------------------------------------

# Read in model csv
reactionlist <- read.csv("topology_full_publication.csv")
myodelist <- as.eqnlist(reactionlist,
                        volumes = c(nIkb = "Vnuc", pnNfk = "Vnuc", nNfk = "Vnuc", nNfkIkb = "Vnuc"))
myode <- as.eqnvec(myodelist)

# Read reduced model csv
reactionlist_red <- read.csv("topology_reduced_publication.csv")
myodelist_red <- as.eqnlist(reactionlist_red,
                        volumes = c(nIkb = "Vnuc", pnNfk = "Vnuc", nNfk = "Vnuc", nNfkIkb = "Vnuc"))
myode_red <- as.eqnvec(myodelist_red)

# Conserved quantities
conserved <- c(NfkIkb = "1*NfkIkb+1*NfkpIkb+1*pNfkIkb+1*pNfkpIkb+1*pNfk+1*Nfk + Vnuc * (1*pnNfk+1*nNfk+1*nNfkIkb)", 
               Ikk = "Ikk + pIkk + ppIkk + iIkk")

# Set list of forcings
forcings <- NULL


# Define new observables based on ODE states
observables0 <- c(
  # Fluorescence microscopy data
  Nuclei = "(s_Nuclei*(pnNfk + nNfk + nNfkIkb) + off_Nuclei)*exp(-tau_Nuclei*time)",
  Cytoplasm = "(s_Cytoplasm*(NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb + pNfk + Nfk) + off_Cytoplasm)*exp(-tau_Cytoplasm*time)",
  A20 = "(s_A20 * A20 + off_A20)*exp(-tau_A20*time)",
  # Immunoblot data
  IkBa = "s_IkBa*(NfkIkb + pNfkIkb + NfkpIkb + pIkb + Ikb) + off_IkBa",
  pp65="s_pp65*(pNfkpIkb + pNfkIkb + pNfk) + off_pp65",
  pIkBa="s_pIkBa*(NfkpIkb + pNfkpIkb + pIkb) + off_pIkBa",
  pIKK = "s_pIKK*(pIkk + ppIkk) + off_pIKK",
  # Immunoblot data from IP experiments
  IkBalpha_IP = "s_IkBa_IP*(NfkIkb + pNfkIkb + NfkpIkb + pIkb + Ikb) + off_IkBa_IP",
  p65_IP = "s_p65_IP*(NfkIkb + pNfkIkb + NfkpIkb + pNfkpIkb + pNfk + Nfk) + off_p65_IP",
  pIkBalpha_tp65_IP = "s_pIkBalpha_tp65_IP*(NfkpIkb + pNfkpIkb) + off_pIkBalpha_tp65_IP",
  pp65_IP = "s_pp65_IP*(pNfkpIkb + pNfkIkb + pNfk) + off_pp65_IP",
  tIkBalpha_pp65_IP = "s_tIkBalpha_pp65_IP*(pNfkpIkb + pNfkIkb) + off_tIkBalpha_pp65_IP",
  tIkBalpha_tp65_IP = "s_tIkBalpha_tp65_IP*(NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb) + off_tIkBalpha_tp65_IP",
  ratio_Nfk_Ikb = "(NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb + pNfk + Nfk)/(NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb + pIkb + Ikb)"
)

# Log-transformed observables
observables <- paste0("log(", observables0, ")"); names(observables) <- names(observables0)


errors <- as.eqnvec(
  paste0("sigma_", names(observables)),
  names = names(observables)
)[c("Nuclei", "Cytoplasm", "A20")]

# Symmetries are broken by TNF, total_NFKB, total_IKK and prod_Ikb
# TNF is an input and unaffected by DCF.
# The second and third are not affected by DCF but prod_Ikb might be
# When seeing an impact on prod_mIkb_by_nNfk or the steady state level
# of mIkb, this can be due to a "hidden" impact on prod_Ikb.

g <- Y(observables, myode, compile = TRUE, modelname = "obsfn", attach.input = FALSE)
err <- Y(errors, c(observables, myode), states = names(observables), compile = TRUE, modelname = "errfn", attach.input = FALSE)

# Same for reduced model
g_red <- Y(observables, myode_red, compile = TRUE, modelname = "obsfn_red", attach.input = FALSE)
err_red <- Y(errors, c(observables, myode_red), states = names(observables), compile = TRUE, modelname = "errfn_red", attach.input = FALSE)


# Generate the model C files, compile them and return a list with func and extended.
fixed <- NULL
model0 <- odemodel(myode, forcings = forcings, compile = TRUE, jacobian = "inz.lsodes", modelname = "model", gridpoints = 361)
model0_red <- odemodel(myode_red, forcings = forcings, compile = TRUE, jacobian = "inz.lsodes", modelname = "model_red", gridpoints = 361)

# Prediction function
x <- Xs(model0, optionsSens = list(method = "lsodes", atol = 1e-9, rtol = 1e-9))
x_red <- Xs(model0_red, optionsSens = list(method = "lsodes", atol = 1e-9, rtol = 1e-9))

# Control the steady state by two concentrations
controlSS <- c("total_NFKB", "total_IKK")



## Prepare data ------------------------------------------------------------

datasheet <- readRDS("data.rds")

# Create datalist 
data <- as.datalist(datasheet, split.by = c("treatment", "compound", "cell", "scale"))

# Set up conditions
condition.grid <- attr(data, "condition.grid")
conditions <- rownames(condition.grid)


## Parameter specification ------------------------------------------------------


## Specific part

conditions <- c("no_no_HepG2_1", "TNF_no_HepG2_1", "no_DCF_HepG2_1", "TNF_DCF_HepG2_1")

setP <- function(conditions, doModelReduction = FALSE, penalized = FALSE, prediction = FALSE, PHH = FALSE) {

  ## General part 
  
  # Get inner parameters (symbols occurring in the equations except for forcings and time)
  if (doModelReduction) 
    innerpars <<- getSymbols(c(names(myode_red), myode_red, observables, controlSS, errors), exclude = c(forcings, "time"))
  else
    innerpars <<- getSymbols(c(names(myode), myode, observables, controlSS, errors), exclude = c(forcings, "time"))
  
  # Initialize transformation by identity
  trafo <- repar("x~x", x = innerpars)
  if (doModelReduction) {
    # Explicit steady states
    trafo <- repar("x~y", x = myodelist_red$states, y = "0", trafo)
    trafo[names(conserved)] <- controlSS
    trafo["nIkb"] <- "nIkb"
  } else {
    # Fix initials of SS-determined states
    trafo <- repar("x~y", x = setdiff(myodelist$states, controlSS), y = "1", trafo)
  }
  # Symmetries
  trafo <- repar("x~y", x = c("s_A20", "Vnuc", "total_NFKB", "total_IKK"), y = "1", trafo)
  # Log transform
  trafo <- repar("x~exp(x)", x = innerpars, trafo)
  
  
  Reduce("+", lapply(conditions, function(C) with(condition.grid[C,], {
    
    obspars <- getSymbols(observables, exclude = c(myodelist$states, "time"))
    errpars <- getSymbols(errors)
    modelpars <- getParameters(model0)
    if (doModelReduction) modelpars <- getParameters(model0_red)
    
    if (cell == "HepG2") {
      # Symmetries
      # trafo <- repar("x~y", x = c("total_NFKB", "total_IKK"), y = "0", trafo)
    }
    
    # Different error and observation parameters for donor experiments
    specific <- c(obspars, errpars)
    if (scale != "1")
      trafo <- repar("x ~ x + Delta_experiment_cell_scale_x", trafo, x = specific, cell = cell, scale = scale)
    
    # Treatment-specific parameters
    if (treatment == "no")
      trafo["TNF"] <- "0"
    if (treatment == "TNF")
      trafo["TNF"] <- "dose_TNF"
    
    # Different model parameters per compound
    if (doModelReduction) {
      subf <- myodelist_red
      specific <- getSymbols(as.eqnvec(subf), exclude = subf$states)
    } else {
      subf <- myodelist
      specific <- getSymbols(as.eqnvec(subf), exclude = subf$states)
    }
    
    if (penalized) {
      specific <- c("act_Ikk_by_TNF", "trigger_iIkk", "act_pIkk", "act_Ikb_by_Ikk", "prod_mIkb_by_nNfk", "degrad_mIkb", "shuttle_RnaA20")
    }
    if (compound != "no")
      trafo <- repar("x ~ x + Delta_compound_x", trafo, x = specific, compound = compound)
    
    if (prediction) {
      
      trafo <- repar("Delta_compound_act_Ikb_by_Ikk ~ 0", trafo, compound = compound)
      trafo <- repar("Delta_compound_act_Ikk_by_TNF ~ -effect_compound_Ikk", trafo, compound = compound)
      trafo <- repar("Delta_compound_act_pIkk ~ effect_compound_Ikk", trafo, compound = compound)
      trafo <- repar("Delta_compound_degrad_mIkb ~ -effect_compound_Ikb", trafo, compound = compound)
      trafo <- repar("Delta_compound_prod_mIkb_by_nNfk ~ -effect_compound_Ikb", trafo, compound = compound)
      trafo <- repar("Delta_compound_shuttle_RnaA20 ~ 0", trafo, compound = compound)
      trafo <- repar("Delta_compound_trigger_iIkk ~ 0", trafo, compound = compound)
      
    }
    
    # Use symmetries to reduce model
    specific <- c("uptake", "prod_Ikb", "build_RnaA20", "split_NfkIkb")
    trafo <- repar("x~y", x = specific, y = "0", trafo)
    
    # Model reduction
    if (doModelReduction) {
      trafo <- replaceSymbols("deact_TNFR", "log(1e-3)", trafo)
      trafo <- replaceSymbols("int_Nfk", "log(1e-2)", trafo) # couples eta_int* parameters
      trafo <- replaceSymbols("ext_nNfkIkb", "log(1000)", trafo)
      trafo["tau_A20"] <- 0 # does not couple to any other parameter
      trafo["off_IkBa_IP"] <- "0"
      if (!PHH) trafo["off_IkBa"] <- "0"
      
      trafo["nIkb"] <- 0
      trafo <- replaceSymbols("form_complex_nuc", "log(1000)", trafo)
      trafo <- replaceSymbols("deact_pnNfk", "log(1000)", trafo)
      
    }
    
    if (PHH) {
      
      timerates <- getSymbols(myodelist_red$rates, myodelist$states)
      timerates <- intersect(timerates, names(trafo))
      trafo[timerates] <- paste0("exp(timescale)*(", trafo[timerates], ")")
      
    }
    
    
    P(trafo, condition = C, attach.input = TRUE)
    
    
  })))
  
}
p <- setP(conditions)


## Steady state transformation
myode.SS <- as.eqnvec(subset(myodelist, !grepl("uptake*TNF", Rate, fixed = TRUE)))
myode.SS[names(conserved)] <- paste0(controlSS, " - (", conserved, ")")
pSS <- P(myode.SS, method = "implicit", modelname = "pSS", compile = TRUE)



## Prediction and objective Functions -------------------------------------------------------


# Model times, data times and all times
times <- 0:360
timesD <- sort(unique(lbind(data)$time))

# Initalize parameters 
initializeParameters <- function(p, reduced = FALSE) {
  
  fixedpars <<- c("dose_TNF")
  outerpars <<- setdiff(getParameters(p), fixedpars)
  
  obspars <- getSymbols(observables, exclude = c(myodelist$states, "time"))
  deltapars <<- outerpars[grep("^Delta", outerpars)]
  
  etapars <<- outerpars[grep("^eta", outerpars)]
  modelpars <<- setdiff(outerpars, c(deltapars, etapars))
  offpars <<- obspars[grepl("^off", obspars)]
  
  
  mu.model <<- structure(rep(-1, length(modelpars)), names = modelpars)
  mu.model[offpars] <<- -2
  sigma.model <<- structure(rep(4, length(modelpars)), names = modelpars)
  sigma.model[offpars] <<- 2
  
  mu.delta <<- structure(rep(0, length(deltapars)), names = deltapars)
  sigma.delta <<- paste("sigma", deltapars, sep = "_")
  
  mu.eta <<- structure(rep(0, length(etapars)), names = etapars)
  sigma.eta <<- structure(rep(4, length(etapars)), names = etapars)
  
  mu.sigma <<- structure(rep(0, length(sigma.eta) + length(sigma.delta)), names = c(sigma.delta, sigma.eta))
  
  prior.model <<- constraintL2(mu.model, sigma = sigma.model)
  prior.delta <<- constraintL2(mu.delta,  sigma = 4) 
  prior.eta <<- constraintL2(mu.eta, sigma = 4)
  
  pouter <<- c(mu.model, mu.delta, mu.eta)
  fixed <<- structure(c(1), names = fixedpars)
  
}



## Skript part (run step by step) -------------------------------------------

if (runscript) {

  ## Part 1: Setup for TNF/control conditions in the full model
  conditions <- c("no_no_HepG2_1", "TNF_no_HepG2_1") 
  p <- setP(conditions)
  initializeParameters(p)
  
  obj <- normL2(data[conditions], g*x*pSS*p, err) + prior.model + prior.eta
  
  bestfit <- readRDS("bestfit_step1.rds")
  plot((g*x*pSS*p)(times, bestfit, fixed = fixed), data[conditions]) + scale_x_sqrt()
  
  
  ## Part 2: TNF/control in the reduced model
  p <- setP(conditions, doModelReduction = TRUE)
  initializeParameters(p)
  
  p_ratios <- P(
    trafo = eqnvec(
      r1 = "RnaA20/A20" # All mRNA vs protein ratios are required to be small
    )
  )
  constr_ratios <- constraintL2(mu = c(r1 = 0), sigma = 0.1, attr.name = "data")
  
  obj <- normL2(data[conditions], g_red*x_red*p, err_red) + prior.model + prior.eta
  
  bestfit <- readRDS("bestfit_step2.rds")
  plot((g_red*x_red*p)(times, bestfit, fixed = fixed), data[conditions]) + scale_x_sqrt()
  
  
  ## Part 3: Reduced model with diclofenac parameters 
  conditions <- c("no_no_HepG2_1", "TNF_no_HepG2_1", "no_DCF_HepG2_1", "TNF_DCF_HepG2_1") 
  p <- setP(conditions, doModelReduction = TRUE)
  initializeParameters(p)
  obj <- normL2(data[conditions], g_red*x_red*p, err_red) + prior.model + prior.eta + prior.delta
  
  bestfit <- readRDS("bestfit_step3.rds")
  plot((g_red*x_red*p)(times, bestfit, fixed = fixed), data[conditions]) + scale_x_sqrt()
  
  
  # Step 4: Reduced model with reduced number of diclofenac parameters
  
  conditions <- c("no_no_HepG2_1", "TNF_no_HepG2_1", "no_DCF_HepG2_1", "TNF_DCF_HepG2_1") 
  p <- setP(conditions, doModelReduction = TRUE, penalized = TRUE)
  initializeParameters(p)
  obj <- normL2(data[conditions], g_red*x_red*p, err_red) + prior.model + prior.eta + prior.delta
  

  bestfit <- readRDS("bestfit_step4.rds")
  plot((g_red*x_red*p)(times, bestfit, fixed = fixed), data[conditions]) + scale_x_sqrt()
  
  
  ## Step 5: Compound fitting
  compounds <- c("DCF", "AMD", "XIM", "APAP", "FIAU")
  
  conditions <- c("TNF_no_HepG2_1", "TNF_DCF_HepG2_1",
                  "TNF_no_HepG2_3", "TNF_AMD_HepG2_3", "TNF_XIM_HepG2_3", "TNF_APAP_HepG2_3", "TNF_FIAU_HepG2_3")
  p <- setP(conditions, doModelReduction = TRUE, penalized = TRUE, prediction = TRUE)
  initializeParameters(p)
  

  
  bestfit <- readRDS("bestfit_step4.rds") 
  bestfit <- bestfit[setdiff(names(bestfit), deltapars), drop = TRUE]
  common <- intersect(names(pouter), names(bestfit))
  different <- setdiff(names(pouter), names(bestfit))
  pouter[common] <- bestfit[common]

 
  prior.delta <- constraintL2(mu = pouter[different], sigma = 4)
  
  
  subdata <- data[conditions]
  subdata$TNF_DCF_HepG2_1 <- subset(subdata$TNF_DCF_HepG2_1, name %in% c("pIKK", "pIkBa", "pp65", "IkBa"))
  obj <- normL2(subdata, g_red*x_red*p, err_red) + prior.delta
  myfit <- trust(obj, pouter[different], rinit = .1, rmax = 1, fixed = c(fixed, pouter[common]))
  plot((g_red*x_red*p)(times, myfit$argument, fixed = c(fixed, pouter[common])), data, 
       name %in% c("A20", "Nuclei", "pIkBa", "pIKK", "pp65") & condition %in% conditions[-(1)]) + 
    facet_wrap(~name, scales = "free", ncol = 1)
  plot((g_red*x_red*p)(times, myfit$argument, fixed = c(fixed, pouter[common])), data, 
       name %in% c("A20", "Nuclei", "pIkBa", "pIKK", "pp65") & condition %in% conditions[-(1)] & time <= 60) + 
    facet_wrap(~name, scales = "free", ncol = 1)
  
  
  ## Step 6: Fit primary human hepatocytes
  conditions <- c("TNF_no_PHH_2", "TNF_DCF_PHH_2")
  p <- setP(conditions, doModelReduction = TRUE, penalized = TRUE, PHH = TRUE, prediction = FALSE)
  initializeParameters(p)
  bestfit <- readRDS("bestfit_step4.rds")
  
  pouter[intersect(names(bestfit), names(pouter))] <- bestfit[intersect(names(bestfit), names(pouter))]
  different <- setdiff(names(pouter), names(bestfit))
  different <- union(different, names(pouter)[grepl("DCF", names(pouter))])
  common <- setdiff(names(pouter), different)
  
  subdata <- subset(data[conditions], !name %in% c("Cytoplasm", "Nuclei"))
  
  prior.delta <- constraintL2(mu = pouter[different], sigma = 4)
  obj <- normL2(subdata, g_red*x_red*p, err_red) + prior.delta
  myfit <- trust(obj, pouter[different], rinit = .1, rmax = 1, fixed = c(fixed, pouter[common]))
  
  plot((g_red*x_red*p)(times, myfit$argument, fixed = c(fixed, pouter[common])), data[conditions]) + scale_x_sqrt()
  
  
}
