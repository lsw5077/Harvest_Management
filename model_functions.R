# Model run functions for fish eco-evolutionary model

# Load some required packages

require(tidyverse)
require(truncnorm)

# Make an inverse function of the in function, which tells
# us if a value is NOT in a vector. 

`%!in%` = Negate(`%in%`)

# Model run function:

# Individual growth rate function:
# Commonly called 'von Bertalanffy's'
# Although the common name is a misnomer b/c
# really this version is from Beverton and Holt


growVB <- function(lAsym, # Asymptotic length
                   k, # growth rate
                   t, # age
                   t0) { # theoretical age at which length = 0
  
  # rounding to the mm to prevent some discontinuities later
  
  lt = round(lAsym*(1-exp(-k*(t-t0))))
  
  return(lt)
  
}

# Find age, inverse of individual growth function
# This is statistically fraught in the real world, but it works here.

findAge <-  function(t0,
                     lt,
                     lAsym,
                     k) {
  
  t = -1*(log(1-lt/lAsym))/k + t0
  
  return(t)
  
  
}



# function to find mass from length
# Standard mass eq, using base10 standard weight model
# I hate that this is base10 too. Unfortunately, this model is for fish biologists, who use base10 for some reason.

find_mass <- function(massConstant,
                      massMultiplier,
                      lt){
  
  logmass = massConstant + (massMultiplier*log10(lt))
  mass = 10^logmass
  
  return(mass)
}

# Find constant for calculating k

findKonstant <- function(lAsymMean,
                  kStart) {
  
  c = log(kStart)/log(lAsymMean)
  
  return(c)
  
}


# Start the population

makePop <- function(cohortSize, # size of initial cohort
                    ageMin, # minimum age of starting fish
                    ageMax, # maximum age of starting fish
                    lAsymMean, # mean starting asymptotic length
                    lAsymSD, # sd of starting asymptotic length
                    lReproMean, # mean Length at first reproduction 
                    t0, # hypothetical time when length = 0
                    kStart, # growth rate
                    massMultiplier, # multiplier for standard mass eq
                    massConstant, # added constant for standard mass eq
                    harvestFrame, # regulation df
                    mortRate, # linear mortality rate for big/old fish
                    mortMultiplier, # parameters for age-specific mort function
                    mortExponent,
                    ageCut, # The age cutoff where fish transition from varying to constant mortality
                    vLength # length at which fish become vulnerable to catch
                    ) { 
  
  matRat = lReproMean/lAsymMean # calculate the length at maturity/total length
  
  Konstant = findKonstant(lAsymMean = lAsymMean,
                          kStart = kStart) # calculate the constant for the asymptotic size-growth rate function
  
  #startPop <- data.frame(age = seq(from = ageMin, to = ageMax)) %>% # Make dataframe of starting age range
   #           mutate(cohortSurvive = round(cohortSize*((1-mortRate)^age))) # predict survival
  
  # Make a dataframe of fish at age composition
  
  
  startPop <- data.frame(age = seq(ageMin,ageMax,1),
                         cohortSurvive = NA) %>%
    rowwise() %>%
    mutate(mortRate = ifelse(age < ageCut, 
                             mortMultiplier*exp(mortExponent*age),      
                             mortRate),
           cohortSurvive = as.numeric(cohortSurvive)) %>%
    as.data.frame()
  
  startPop$cohortSurvive[1] <- cohortSize
  
  
  for (w in 2:nrow(startPop)) {
    
    startPop$cohortSurvive[w] <- floor(startPop$cohortSurvive[w-1]*
                                         (1-startPop$mortRate[w-1]))
    
  }
  
  # Expand population dataframe based on age, mortality, and predicted size
  
  
  startPop <- startPop[rep(seq(nrow(startPop)), as.numeric(startPop$cohortSurvive)),] %>% # expand by cohort number
    mutate(ID = paste("0", row_number(), sep = "_")) %>% rowwise() %>% # assign ID
    mutate(lAsymGeno = rnorm(1, lAsymMean, lAsymSD), # Draw an asymptotic length
           lAsymPheno = lAsymGeno, # For the starting population, their genotypic and phenotypic sizes are the same
           k = exp(Konstant*log(lAsymPheno)), # Calculate k growth constant
           lt = growVB(t = age,
                        k = k,
                        t0 = t0,
                        lAsym = lAsymPheno), # Calculate length
           lRepro = matRat*lAsymPheno, # Calculate reproductive length
           tRepro = findAge(t0 = t0,
                            lAsym = lAsymPheno,
                            k = k,
                            lt = lRepro), # calculate age at reproduction
           stage = ifelse (lt < lRepro, "Juvenile", "Adult"), # assign stage
           mass = find_mass(massConstant = massConstant, 
                            massMultiplier = massMultiplier,
                            lt = lt), # find mass given length
           vulnerable = ifelse(lt >= vLength, 1, 0) # assign vulnerability to catch based on length
    ) %>% 
    select(ID, age, lt, lRepro, lAsymGeno, lAsymPheno, k, tRepro, mass, stage, vulnerable) %>%
    rowwise() %>%
    mutate(legal = harvestFrame$legal[harvestFrame$binMin <= round(lt) &
                                        harvestFrame$binMax > round(lt)], 
           vulnerable = ifelse(lt >= vLength, 1,0))  # assign legality of harvest and vulnerability to catch
  
  return(startPop) # return initiated population
  
}


# growth

growPop <- function(pop, # population passed from initation function
                    harvestFrame, # regulation dataframe
                    t0, # parameters from the individual growth model
                    massConstant, # standard mass params
                    massMultiplier,
                    vLength) { # length at which vulnerability starts
  
  
  growPop <- pop %>%
    mutate(age = age + 1, # advance age
           lt = growVB(t = age,
                        k = k,
                        t0 = t0,
                        lAsym = lAsymPheno), #  Find length
           stage = ifelse (lt < lRepro, "Juvenile", "Adult"), # update stage if necessary
           mass = find_mass(massConstant = massConstant, 
                            massMultiplier = massMultiplier,
                            lt = lt)) %>% # update mass per newlength
    rowwise() %>%
    mutate(legal = harvestFrame$legal[harvestFrame$binMin <= round(lt) &
                                        harvestFrame$binMax > round(lt)],
           vulnerable = ifelse(lt >= vLength, 1,0)) %>% # update legality and vulnerability
    select(ID, age, lt, lRepro, lAsymGeno, lAsymPheno, k, tRepro, mass, stage, vulnerable, legal)
  
  
  return(growPop)
  
}


# Reproduction, w/ summer survival

reproduce <- function(pop, # population passed from previous step
                      rickerA,
                      rickerB,
                      g,
                      eggMass, # mass of each egg in grams
                      recruitZ, # expected survival from egg to first fall
                      vLength, # Minimum vulnerable length
                      h2lAsym, # heritability of asymptotic length
                      lReproMean,
                      t0, # theoretical time when length = 0
                      kStart, # Starting Brody growth constant
                      lAsymMean, # starting asymptotic length mean
                      lAsymSD, # starting asymptotic length sd
                      timeStep, # timestep within simulation
                      sim, # simulation number
                      massConstant, # mass constant for standard mass eq
                      massMultiplier, # mass multiplier for standard mass eq
                      harvestFrame, # harvest dataframe
                      lowerLengthLimit, # species minimum length
                      upperLengthLimit) { # species maximum length
  
   # calculate how many offspring each adult fish should generate based on their size
  # Correcting by a factor of 0.5 b/c we only want females
  # Final final revision. Can we change logistic part so that it just does recruitment?
  
  
  matRat = lReproMean/lAsymMean
  
  Konstant = findKonstant(lAsymMean = lAsymMean,
                          kStart = kStart)
  
  

  repro_pop <- pop %>%
    filter(stage == "Adult") %>%
     mutate( # calculate expected number of offspring
    nOffspring = floor((0.5* # sex ratio
                        g* # gonadal-somatic index
                        mass* # fish's mass
                        recruitZ # expected survival from egg to first fall
                        )/eggMass) # number of fall recruits (mass allocation*mass*first summer survival)
    ) %>% # divide total repro investment by egg mass to get total recruits
    filter(nOffspring > 0) # Keep only individuals who produce at least one offspring.
  
   rickerRecruits <- round(nrow(pop[pop$stage == "Adult",])*exp(rickerA + rickerB*nrow(pop[pop$stage == "Adult",])))
  
  
  # This dummy row keeps the function from throwing an error in 
  # iterations where 0 offspring are produced (b/c you can't rep by 0!)
  # We drop it in the next step.
  
  dummy <-  data.frame(ID = NA,
                       age = NA,
                       lt = NA,
                       lRepro = NA,
                       lAsymGeno = NA,
                       lAsymPheno = NA,
                       k = NA,
                       tRepro = NA,
                       mass = NA,
                       stage = NA,
                       vulnerable = NA,
                       legal = NA,
                       nOffspring = 1)
  
  repro_pop <- rbind(repro_pop, dummy)
  
  repro_pop <- repro_pop[rep(seq(nrow(repro_pop)), as.numeric(repro_pop$nOffspring)),] %>% # replicate by # of offspring
    as.data.frame() %>% # drop any rows w/ NA length, i.e. the dummy row
    mutate(ID = paste(sim, timeStep, row_number(), sep = "_"), # Make a new ID
           age = 0) %>% rowwise() %>% # Assign an age, then for each new individual:
    mutate(lAsymGeno = # Calculate asymptotic length "genotype." Really a quantitative trait.
    ifelse(is.na(lt), NA, rtruncnorm(n = 1, # Because this is designed for people who think more about fish than evolutionary
                                     a = lowerLengthLimit, # biology, it uses a truncated normal distribution so you can set
                                     b = upperLengthLimit, # species upper and lower possible length limits.
                                     mean = lAsymGeno, # Draw a value from a truncated normal whose mean = the parent's mean
                                     sd = sd(pop$lAsymGeno))), # and whose sd = the population sd
           lAsymPheno = lAsymGeno*h2lAsym + # Now, calculate the phenotype by summing the product of heritability and "genotype" 
                        (1-h2lAsym)*(rnorm(n = 1, mean = lAsymMean, sd = lAsymSD)), # And the product of 1-heritability and the phenotypic modifier, here the initial condition
           # calculate all resulting life history params based on phenotype: 
           k = exp(Konstant*log(lAsymPheno)), # k growth constant
           g = g,
           lRepro = matRat*lAsymPheno, # Length at reproduction
           tRepro = findAge(t0 = t0,
                            lt = lRepro,
                            lAsym = lAsymPheno,
                            k = k), # Age at reproduction
           lt = growVB(t = 0,
                       k = k,
                       t0 = t0,
                       lAsym = lAsymPheno), # find length at age 0
           stage = "Juvenile", # Assign stage = juvenile
           mass = find_mass(massConstant = massConstant,
                            massMultiplier = massMultiplier,
                            lt = lt)) %>% filter(!is.na(lt)) %>%  # calculate mass of new age0 fish
          filter(!is.na(lt)) %>% # drop the dummy row
          mutate(vulnerable = ifelse(lt >= vLength, 1,0), # assign vulnerability based on length
                        legal = harvestFrame$legal[harvestFrame$binMin <= round(lt) &
                                              harvestFrame$binMax > round(lt)] # assign legality based on regulations
                        ) %>% 
    sample_n(rickerRecruits) %>%
    select(ID, age, lt, lRepro, lAsymGeno, lAsymPheno, k, tRepro, mass, stage, vulnerable, legal)
  
  

  pop <- rbind(pop, repro_pop) 
  
  repro_out <- list(pop = pop,
                    nRecruits = nrow(pop[pop$age == 0,]))
  
  return(repro_out)
  

}

# Harvest


harvest <- function(pop, # population being harvested
                    probCatch, # probability of being caught given vulnerability
                    probHarv, # probability of being harvested given legal capture
                    prefLength, # preferred/ "trophy" fish length
                    discardMortRate) { # probability of dying if released given capture
  
  # Harvested 
  
  harvPop <- pop %>%
    filter(vulnerable == 1, # filter pop to legal and vulnerable fish
           legal == 1) %>%
    rowwise() %>% # for each fish:
    mutate(caught = rbinom(n = 1, 
                           size = 1,
                           prob = probCatch), # draw a catch outcome from a binomial dist
           harvested = ifelse(caught == 0, 0,
                              rbinom(n = 1,
                                     size = 1,
                                     prob = probHarv)), # draw a harvest outcome for those caught
           discardMort = ifelse(caught == 1 &  
                                  harvested == 0,
                                rbinom(n = 1,
                                       size = 1,
                                       prob = discardMortRate), 0)) %>% # draw a discard mortality outcome for those released
    filter(harvested == 1 | discardMort == 1)
  
  # Discard Mort from catch of vulnerable illegal fish
  
  illegalDiscards <- pop %>% 
    filter(vulnerable == 1,
           legal == 0) %>% # filter to fish vulnerable to catch but not legal to harvest
    rowwise() %>%
    mutate(caught = rbinom(n = 1,
                           size = 1,
                           prob = probCatch), # draw catch outcome given catch prob
           discardMort = ifelse(caught == 1,
                                rbinom(n = 1,
                                       size = 1,
                                       prob = discardMortRate), 0)) %>% # draw discard mortality outcome given discardMort rate
    filter(discardMort == 1) # keep only those killed by angling discard
  
  
  pop <- pop %>% filter(ID %!in% harvPop$ID,
                        ID %!in% illegalDiscards$ID) # filter the population to only those who survived catch and harvest
  
  
  massHarv <- sum(harvPop$mass[harvPop$harvested == 1]) # sum up total biomass harvested
  nDiscard <- nrow(illegalDiscards) # count discard morts
  massDiscard <- sum(illegalDiscards$mass,
                     harvPop$mass[harvPop$discardMort == 1]) # sum up discard mass
  meanHarvLength <- mean(harvPop$lt) # calculate mean length of fish harvested
  
  nHarv <- nrow(harvPop %>% filter(harvested == 1)) # count up harvested fish 
  nPref <- nrow(harvPop %>% filter(harvested == 1, lt >= prefLength)) # count up number of preferred-size fish harvested
  
  
  nHarv <- 
  
  harvOut <- list(pop = pop, 
                  nHarv = ifelse(is.na(nHarv), 0, nHarv),
                  massHarv = ifelse(is.na(massHarv), 0, massHarv),
                  nDiscard = ifelse(is.na(nDiscard), 0, nDiscard),
                  massDiscard = ifelse(is.na(massDiscard), 0, massDiscard),
                  nPref = ifelse(is.na(nPref), 0, nPref),
                  meanHarvLength = meanHarvLength)
  
  return(harvOut) # output results as list.
  
  
}

# natural death

naturalDeath <- function(pop, # population passed
                         mortMultiplier, # parameters for age-specific mort function
                         mortExponent,
                         mortRate, # linear mortality rate for big/old fish
                         ageCut, # The age cutoff where fish transition from varying to constant mortality
                         compMort) { # whether mortality is compensatory

  
  # Fish undergo natural mortality
  survive_pop <- pop %>%
    rowwise() %>%
    mutate (probDeath = ifelse(age < ageCut,
                               mortMultiplier*exp(mortExponent*age),
                               mortRate),
            dead = rbinom(n = 1, size = 1,
                          prob = probDeath)) %>%
    filter(dead == 0) %>% # Draw a death value from a binomial distribution and exclude dead fish
    select(ID, age, lt, lRepro, lAsymGeno, lAsymPheno, k, tRepro, mass, stage, vulnerable, legal)
  

  return(survive_pop)
  
}
                                                                           

# Run model function to bundle everything up and make nice outputs

modelSpringSpawn <- function(nSim,
                             nSteps,
                             compMort,
                             massConstant,
                             massMultiplier,
                             cohortSize,
                             ageMin,
                             ageMax,
                             lAsymMean,
                             lAsymSD,
                             vLength,
                             recruitZ,
                             eggMass,
                             g,
                             rickerA,
                             rickerB,
                             h2lAsym,
                             mortRate,
                             mortMultiplier,
                             mortExponent,
                             ageCut,
                             t0,
                             kStart,
                             lReproMean,
                             discardMortRate,
                             probCatch,
                             probHarv,
                             harvestFrame,
                             prefLength,
                             lowerLengthLimit, 
                             upperLengthLimit,
                             keepPops,
                             popPath = NULL){
  
  # Make output dataframe
  
  OUT <- expand.grid(seq(1:nSim),
                     seq(1:nSteps)) %>%
    rename("sim" = Var1, "t" = Var2) %>%
    mutate(N = NA,
           nJuvenile = NA,
           nAdult = NA,
           meanlRepro = NA,
           meantRepro = NA,
           meanLength = NA, 
           meanK = NA,
           sdK = NA,
           sdlRepro = NA,
           sdtRepro = NA,
           nPref = NA,
           nHarv = NA,
           meanHarvLength = NA,
           massHarv = NA,
           nDiscard = NA,
           massDiscard = NA,
           recruits = NA)
  
  
  # For each simulation, 
  
  for (i in 1:nSim) {
    
    # Make population
    
    pop <- makePop (cohortSize = cohortSize, # size of initial cohort
                    ageMin = ageMin, # minimum age of starting fish
                    ageMax = ageMax, # maximum age of starting fish
                    lAsymMean = lAsymMean, # mean starting asymptotic length
                    lAsymSD = lAsymSD, # sd of starting asymptotic length
                    lReproMean = lReproMean, # ratio of total length when fish matures 
                    t0 = t0, # hypothetical time when length = 0
                    kStart = kStart, # growth rate
                    massMultiplier = massMultiplier, # multiplier for standard mass eq
                    massConstant = massConstant, # added constant for standard mass eq
                    harvestFrame = harvestFrame, # regulation df
                    mortRate = mortRate, # mortality rate 
                    mortMultiplier = mortMultiplier,
                    mortExponent = mortExponent,
                    ageCut = ageCut,
                    vLength = vLength) # vulnerability length
    
    
    # Fill in timestep 1 population values
    
    OUT$N[OUT$sim == i & OUT$t == 1] <- nrow(pop)
    OUT$nJuvenile[OUT$sim == i & OUT$t == 1] <- nrow(pop[pop$stage == 'Juvenile',])
    OUT$nAdult[OUT$sim == i & OUT$t == 1] <- nrow(pop[pop$stage == 'Adult',])
    OUT$recruits[OUT$sim == i & OUT$t == 1] = nrow(pop[pop$age == 1,])
    
    OUT$meanlRepro[OUT$sim == i & OUT$t == 1] <- mean(pop$lRepro)
    OUT$meanlAsymGeno[OUT$sim == i & OUT$t == 1] <- mean(pop$lAsymGeno)
    OUT$meanlAsymPheno[OUT$sim == i & OUT$t == 1] <- mean(pop$lAsymPheno)
  
    OUT$meanK[OUT$sim == i & OUT$t == 1] <- mean(pop$k)
    OUT$sdK[OUT$sim == i & OUT$t == 1] <- sd(pop$k)
    OUT$sdlRepro[OUT$sim == i & OUT$t == 1] <- sd(pop$lRepro)
    OUT$sdtRepro[OUT$sim == i & OUT$t == 1] <- sd(pop$tRepro)
    OUT$sdlAsymGeno[OUT$sim == i & OUT$t == 1] <- sd(pop$lAsymGeno)
    OUT$sdlAsymPheno[OUT$sim == i & OUT$t == 1] <- sd(pop$lAsymPheno)
    OUT$meantRepro[OUT$sim == i & OUT$t == 1] = mean(pop$tRepro)
    OUT$meanLength[OUT$sim == i & OUT$t == 1] = mean(pop$lt)
    
    for (j in 2:nSteps) {
      
      print(paste("Sim number", i, "of", nSim, 
                  ",", "Timestep", j, "of", nSteps,
                  "N = ", nrow(pop)))
      
      # Fish Grow  

      pop <- growPop(pop = pop,
                     harvestFrame = harvestFrame,
                     t0 = t0,
                     massConstant = massConstant,
                     massMultiplier = massMultiplier,
                     vLength = vLength)
     
      Nt0 = nrow(pop)
      
      
      repro_out <- reproduce(pop = pop,
                             rickerA = rickerA,
                             rickerB = rickerB,
                             eggMass = eggMass,
                             recruitZ = recruitZ,
                             g = g,
                             vLength = vLength,
                             h2lAsym = h2lAsym,
                             lAsymMean = lAsymMean,
                             lAsymSD = lAsymSD,
                             t0 = t0,
                             lReproMean = lReproMean,
                             kStart = kStart,
                             timeStep = j,
                             sim = i, 
                             massConstant = massConstant,
                             massMultiplier = massMultiplier,
                             harvestFrame = harvestFrame,
                             lowerLengthLimit = lowerLengthLimit,
                             upperLengthLimit = upperLengthLimit)
      
      nRecruits <- repro_out$nRecruits
      pop <- repro_out$pop
      
      # fish are harvested
      
      harv <- harvest(pop,
                      probHarv = probHarv,
                      probCatch = probCatch,
                      discardMortRate = discardMortRate,
                      prefLength = prefLength)
      
      pop <- harv$pop
      nHarv <- harv$nHarv
      
      # fish die naturally
      
      pop <- naturalDeath(pop = pop,
                          mortMultiplier = mortMultiplier,
                          mortExponent = mortExponent,
                          ageCut = ageCut,
                          mortRate = mortRate,
                          compMort = compMort)
      
    # write pop to file so we can inspect if model crashes      

    if (keepPops == TRUE) {write.csv(pop, paste0(popPath, i,"-", j,".csv"))}
      
      
      OUT$N[OUT$sim == i & OUT$t == j] <- nrow(pop)
      OUT$nJuvenile[OUT$sim == i & OUT$t == j] <- nrow(pop[pop$stage == 'Juvenile',])
      OUT$nAdult[OUT$sim == i & OUT$t == j] <- nrow(pop[pop$stage == 'Adult',])
      OUT$recruits[OUT$sim == i & OUT$t == j] = nrow(pop[pop$age == 0,])
      OUT$nHarv[OUT$sim == i & OUT$t == j] = harv$nHarv
      
      OUT$meanlRepro[OUT$sim == i & OUT$t == j] <- mean(pop$lRepro)
      OUT$sdlRepro[OUT$sim == i & OUT$t == j] <- sd(pop$lRepro)
      OUT$sdtRepro[OUT$sim == i & OUT$t == j] <- sd(pop$tRepro)
      OUT$sdlAsymGeno[OUT$sim == i & OUT$t == j] <- sd(pop$lAsymGeno)
      OUT$sdlAsymPheno[OUT$sim == i & OUT$t == j] <- sd(pop$lAsymPheno)
      
      OUT$meanlAsymGeno[OUT$sim == i & OUT$t == j] <- mean(pop$lAsymGeno)
      OUT$meanlAsymPheno[OUT$sim == i & OUT$t == j] <- mean(pop$lAsymPheno)
      
      OUT$meantRepro[OUT$sim == i & OUT$t == j] = mean(pop$tRepro)
      OUT$meanLength[OUT$sim == i & OUT$t == j] = mean(pop$lt)
      OUT$meanHarvLength[OUT$sim == i & OUT$t == j] = harv$meanHarvLength
      
      
      OUT$meanK[OUT$sim == i & OUT$t == j] <- mean(pop$k)
      OUT$sdK[OUT$sim == i & OUT$t == j] <- sd(pop$k)
      OUT$massHarv[OUT$sim == i & OUT$t == j] = harv$massHarv
      OUT$nDiscard[OUT$sim == i & OUT$t == j] = harv$nDiscard
      OUT$massDiscard[OUT$sim == i & OUT$t == j] = harv$massDiscard
      OUT$nPref[OUT$sim == i & OUT$t == j] = harv$nPref
      
      
      if(nrow(pop) == 0) {
        break
      }
      
      
    }
    
  }  
  
  return(OUT)
  
}



