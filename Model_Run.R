
library(tidyverse)
library(egg)

# Minimum harvest length

source("~functions_09152020.R")


floorFrame <- data.frame(binMin = c(0, 381),
                         binMax = c(380.999, 2000),
                         legal = c(0, 1))

# Protected slot regulation

protectedSlot <- data.frame(binMin = c(0, 381, 508),
                            binMax = c(380.999, 507.999, 2000),
                            legal = c(1, 0, 1))

# Harvest slot regulation

harvestSlot <- data.frame(binMin = c(0, 381, 508),
                          binMax = c(380.999, 507.999, 2000),
                          legal = c(0, 1, 0))



nSim = 3
nSteps = 100
massConstant = -5.453
massMultiplier = 3.180
cohortSize = 1000
ageMin = 0
ageMax = 12

h2lAsym = 0.3
lAsymMean = 700
lAsymSD = 100

vLength = 300
recruitZ = 0.0001
eggMass = 0.0001
g = 0.24

mortRate = 0.3
t0 = -1.23
kStart = 0.13

discardMortRate = 0.1
probHarv = 1
prefLength = 508
compMort = FALSE
lReproMean = 420

upperLengthLimit = 1200
lowerLengthLimit = 50

mortMultiplier = 0.7
mortExponent = -0.114
ageCut = 0

rickerA = 3.392
rickerB = -0.001176

keepPops = TRUE


floorHi <- modelSpringSpawn(nSim = 3,
                           nSteps = 100,
                           compMort = compMort,
                           massConstant = massConstant,
                           massMultiplier = massMultiplier,
                           cohortSize = cohortSize,
                           ageMin = ageMin,
                           ageMax = ageMax,
                           lAsymMean = lAsymMean,
                           lAsymSD = lAsymSD,
                           vLength = vLength,
                           recruitZ = recruitZ,
                           eggMass = eggMass,
                           g = g,
                           rickerA = rickerA,
                           rickerB = rickerB,
                           h2lAsym = 0.3,
                           mortRate = mortRate,
                           mortMultiplier = mortMultiplier,
                           mortExponent = mortExponent,
                           ageCut = ageCut,
                           t0 = t0,
                           kStart = kStart,
                           lReproMean = lReproMean,
                           discardMortRate = discardMortRate,
                           probCatch = 0.6,
                           probHarv = probHarv,
                           harvestFrame = floorFrame,
                           prefLength = prefLength,
                           upperLengthLimit = upperLengthLimit,
                           lowerLengthLimit = lowerLengthLimit, 
                           keepPops = TRUE,
                           popPath = "floorHi/")

write_rds(floorHi, path = "floorHi.rds")

floorHiNoEvo <- modelSpringSpawn(nSim = 3,
                                nSteps = 100,
                                compMort = compMort,
                                massConstant = massConstant,
                                massMultiplier = massMultiplier,
                                cohortSize = cohortSize,
                                ageMin = ageMin,
                                ageMax = ageMax,
                                lAsymMean = lAsymMean,
                                lAsymSD = lAsymSD,
                                vLength = vLength,
                                recruitZ = recruitZ,
                                eggMass = eggMass,
                                g = g,
                                rickerA = rickerA,
                                rickerB = rickerB,
                                h2lAsym = 0,
                                mortRate = mortRate,
                                mortMultiplier = mortMultiplier,
                                mortExponent = mortExponent,
                                ageCut = ageCut,
                                t0 = t0,
                                kStart = kStart,
                                lReproMean = lReproMean,
                                discardMortRate = discardMortRate,
                                probCatch = 0.6,
                                probHarv = probHarv,
                                harvestFrame = floorFrame,
                                prefLength = prefLength,
                                upperLengthLimit = upperLengthLimit,
                                lowerLengthLimit = lowerLengthLimit,
                                keepPops = TRUE,
                                popPath = "floorHiNoEvo/")

#write_rds(floorHiNoEvo, path = "floorHiNoEvo.rds")


floorLo <- modelSpringSpawn(nSim = 3,
                           nSteps = 100,
                           compMort = compMort,
                           massConstant = massConstant,
                           massMultiplier = massMultiplier,
                           cohortSize = cohortSize,
                           ageMin = ageMin,
                           ageMax = ageMax,
                           lAsymMean = lAsymMean,
                           lAsymSD = lAsymSD,
                           vLength = vLength,
                           recruitZ = recruitZ,
                           eggMass = eggMass,
                           g = g,
                           h2lAsym = 0.3,
                           rickerA = rickerA,
                           rickerB = rickerB,
                           mortRate = mortRate,
                           mortMultiplier = mortMultiplier,
                           mortExponent = mortExponent,
                           ageCut = ageCut,
                           t0 = t0,
                           kStart = kStart,
                           lReproMean = lReproMean,
                           discardMortRate = discardMortRate,
                           probCatch = 0.1,
                           probHarv = probHarv,
                           harvestFrame = floorFrame,
                           prefLength = prefLength,
                           upperLengthLimit = upperLengthLimit,
                           lowerLengthLimit = lowerLengthLimit,
                           keepPops = TRUE,
                           popPath = "floorLo/")

#write_rds(floorLo, path = "floorLo.rds")


floorLoNoEvo <- modelSpringSpawn(nSim = 3,
                           nSteps = 100,
                           compMort = compMort,
                           massConstant = massConstant,
                           massMultiplier = massMultiplier,
                           cohortSize = cohortSize,
                           ageMin = ageMin,
                           ageMax = ageMax,
                           lAsymMean = lAsymMean,
                           lAsymSD = lAsymSD,
                           vLength = vLength,
                           recruitZ = recruitZ,
                           eggMass = eggMass,
                           g = g,
                           h2lAsym = 0,
                           rickerA = rickerA,
                           rickerB = rickerB,
                           mortRate = mortRate,
                           mortMultiplier = mortMultiplier,
                           mortExponent = mortExponent,
                           ageCut = ageCut,
                           t0 = t0,
                           kStart = kStart,
                           lReproMean = lReproMean,
                           discardMortRate = discardMortRate,
                           probCatch = 0.1,
                           probHarv = probHarv,
                           harvestFrame = floorFrame,
                           prefLength = prefLength,
                           upperLengthLimit = upperLengthLimit,
                           lowerLengthLimit = lowerLengthLimit,
                           keepPops = TRUE,
                           popPath = "floorLoNoEvo/")

#write_rds(floorLoNoEvo2, path = "floorLoNoEvo2.rds")


# Harvest slot



HiHslot <- modelSpringSpawn(nSim = 3,
                            nSteps = 100,
                            compMort = compMort,
                            massConstant = massConstant,
                            massMultiplier = massMultiplier,
                            cohortSize = cohortSize,
                            ageMin = ageMin,
                            ageMax = ageMax,
                            lAsymMean = lAsymMean,
                            lAsymSD = lAsymSD,
                            vLength = vLength,
                            recruitZ = recruitZ,
                            eggMass = eggMass,
                            rickerA = rickerA,
                            rickerB = rickerB,
                            h2lAsym = 0.3,
                            g = g,
                            mortRate = mortRate,
                            mortMultiplier = mortMultiplier,
                            mortExponent = mortExponent,
                            ageCut = ageCut,
                            t0 = t0,
                            kStart = kStart,
                            lReproMean = lReproMean,
                            discardMortRate = discardMortRate,
                            probCatch = 0.6,
                            probHarv = probHarv,
                            harvestFrame = harvestSlot,
                            prefLength = prefLength,
                            upperLengthLimit = upperLengthLimit,
                            lowerLengthLimit = lowerLengthLimit,
                            keepPops = TRUE,
                            popPath = "hSlotHi/")


#write_rds(HiHslot, path = "HiHSlot.rds")


LoHslot <- modelSpringSpawn(nSim = 3,
                            nSteps = 100,
                            compMort = compMort,
                            massConstant = massConstant,
                            massMultiplier = massMultiplier,
                            cohortSize = cohortSize,
                            ageMin = ageMin,
                            ageMax = ageMax,
                            lAsymMean = lAsymMean,
                            lAsymSD = lAsymSD,
                            vLength = vLength,
                            recruitZ = recruitZ,
                            eggMass = eggMass,
                            rickerA = rickerA,
                            rickerB = rickerB,
                            h2lAsym = 0.3,
                            g = g,
                            mortRate = mortRate,
                            mortMultiplier = mortMultiplier,
                            mortExponent = mortExponent,
                            ageCut = ageCut,
                            t0 = t0,
                            kStart = kStart,
                            lReproMean = lReproMean,
                            discardMortRate = discardMortRate,
                            probCatch = 0.1,
                            probHarv = probHarv,
                            harvestFrame = harvestSlot,
                            prefLength = prefLength,
                            upperLengthLimit = upperLengthLimit,
                            lowerLengthLimit = lowerLengthLimit,
                            keepPops = TRUE,
                            popPath = "hSlotLo/")

#write_rds(LoHslot, path = "LoHSlot.rds")


HiHslotNoEvo <- modelSpringSpawn(nSim = 3,
                                 nSteps = 100,
                                 compMort = compMort,
                                 massConstant = massConstant,
                                 massMultiplier = massMultiplier,
                                 cohortSize = cohortSize,
                                 ageMin = ageMin,
                                 ageMax = ageMax,
                                 lAsymMean = lAsymMean,
                                 lAsymSD = lAsymSD,
                                 vLength = vLength,
                                 recruitZ = recruitZ,
                                 eggMass = eggMass,
                                 rickerA = rickerA,
                                 rickerB = rickerB,
                                 h2lAsym = 0,
                                 g = g,
                                 mortRate = mortRate,
                                 mortMultiplier = mortMultiplier,
                                 mortExponent = mortExponent,
                                 ageCut = ageCut,
                                 t0 = t0,
                                 kStart = kStart,
                                 lReproMean = lReproMean,
                                 discardMortRate = discardMortRate,
                                 probCatch = 0.6,
                                 probHarv = probHarv,
                                 harvestFrame = harvestSlot,
                                 prefLength = prefLength,
                                 upperLengthLimit = upperLengthLimit,
                                 lowerLengthLimit = lowerLengthLimit,
                                 keepPops = TRUE,
                                 popPath = "hSlotHiNoEvo/")

#write_rds(HiHslotNoEvo, path = "HiHslotNoEvo.rds")


LoHslotNoEvo <- modelSpringSpawn(nSim = 3,
                                 nSteps = 100,
                                 compMort = compMort,
                                 massConstant = massConstant,
                                 massMultiplier = massMultiplier,
                                 cohortSize = cohortSize,
                                 ageMin = ageMin,
                                 ageMax = ageMax,
                                 lAsymMean = lAsymMean,
                                 lAsymSD = lAsymSD,
                                 h2lAsym = 0,
                                 g = g,
                                 vLength = vLength,
                                 recruitZ = recruitZ,
                                 eggMass = eggMass,
                                 rickerA = rickerA,
                                 rickerB = rickerB,
                                 mortRate = mortRate,
                                 mortMultiplier = mortMultiplier,
                                 mortExponent = mortExponent,
                                 ageCut = ageCut,
                                 t0 = t0,
                                 kStart = kStart,
                                 lReproMean = lReproMean,
                                 discardMortRate = discardMortRate,
                                 probCatch = 0.1,
                                 probHarv = probHarv,
                                 harvestFrame = harvestSlot,
                                 prefLength = prefLength,
                                 upperLengthLimit = upperLengthLimit,
                                 lowerLengthLimit = lowerLengthLimit,
                                 keepPops = TRUE,
                                 popPath = "hSlotLoNoEvo/")

write_rds(LoHslotNoEvo, path = "LoHSlotNoEvo.rds")

# Protected slot



HiPslot <- modelSpringSpawn(nSim = 3,
                            nSteps = 100,
                            compMort = compMort,
                            massConstant = massConstant,
                            massMultiplier = massMultiplier,
                            cohortSize = cohortSize,
                            ageMin = ageMin,
                            ageMax = ageMax,
                            lAsymMean = lAsymMean,
                            lAsymSD = lAsymSD,
                            vLength = vLength,
                            recruitZ = recruitZ,
                            eggMass = eggMass,
                            g = g,
                            rickerA = rickerA,
                            rickerB = rickerB,
                            h2lAsym = 0.3,
                            mortRate = mortRate,
                            mortMultiplier = mortMultiplier,
                            mortExponent = mortExponent,
                            ageCut = ageCut,
                            t0 = t0,
                            kStart = kStart,
                            lReproMean = lReproMean,
                            discardMortRate = discardMortRate,
                            probCatch = 0.6,
                            probHarv = probHarv,
                            harvestFrame = protectedSlot,
                            prefLength = prefLength,
                            upperLengthLimit = upperLengthLimit,
                            lowerLengthLimit = lowerLengthLimit,
                            keepPops = TRUE,
                            popPath = "pSlotHi/")

write_rds(HiPslot, path = "HiPSlot.rds")




HiPslotNoEvo <- modelSpringSpawn(nSim = 3,
                                 nSteps = 100,
                                 compMort = compMort,
                                 massConstant = massConstant,
                                 massMultiplier = massMultiplier,
                                 cohortSize = cohortSize,
                                 ageMin = ageMin,
                                 ageMax = ageMax,
                                 lAsymMean = lAsymMean,
                                 lAsymSD = lAsymSD,
                                 vLength = vLength,
                                 recruitZ = recruitZ,
                                 eggMass = eggMass,
                                 g = g,
                                 rickerA = rickerA,
                                 rickerB = rickerB,
                                 h2lAsym = 0,
                                 mortRate = mortRate,
                                 mortMultiplier = mortMultiplier,
                                 mortExponent = mortExponent,
                                 ageCut = ageCut,
                                 t0 = t0,
                                 kStart = kStart,
                                 lReproMean = lReproMean,
                                 discardMortRate = discardMortRate,
                                 probCatch = 0.6,
                                 probHarv = probHarv,
                                 harvestFrame = protectedSlot,
                                 prefLength = prefLength,
                                 upperLengthLimit = upperLengthLimit,
                                 lowerLengthLimit = lowerLengthLimit,
                                 keepPops = TRUE,
                                 popPath = "pSlotHiNoEvo/")

write_rds(HiPslotNoEvo, path = "HiPSlotNoEvo.rds")



LoPslot <- modelSpringSpawn(nSim = 3,
                            nSteps = 100,
                            compMort = compMort,
                            massConstant = massConstant,
                            massMultiplier = massMultiplier,
                            cohortSize = cohortSize,
                            ageMin = ageMin,
                            ageMax = ageMax,
                            lAsymMean = lAsymMean,
                            lAsymSD = lAsymSD,
                            vLength = vLength,
                            recruitZ = recruitZ,
                            eggMass = eggMass,
                            g = g,
                            rickerA = rickerA,
                            rickerB = rickerB,
                            h2lAsym = 0.3,
                            mortRate = mortRate,
                            mortMultiplier = mortMultiplier,
                            mortExponent = mortExponent,
                            ageCut = ageCut,
                            t0 = t0,
                            kStart = kStart,
                            lReproMean = lReproMean,
                            discardMortRate = discardMortRate,
                            probCatch = 0.1,
                            probHarv = probHarv,
                            harvestFrame = protectedSlot,
                            prefLength = prefLength,
                            upperLengthLimit = upperLengthLimit,
                            lowerLengthLimit = lowerLengthLimit,
                            keepPops = TRUE,
                            popPath = "pSlotLo/")

#write_rds(LoPslot, path = "LoPSlot.rds")



LoPslotNoEvo <- modelSpringSpawn(nSim = 3,
                                 nSteps = 100,
                                 compMort = compMort,
                                 massConstant = massConstant,
                                 massMultiplier = massMultiplier,
                                 cohortSize = cohortSize,
                                 ageMin = ageMin,
                                 ageMax = ageMax,
                                 lAsymMean = lAsymMean,
                                 lAsymSD = lAsymSD,
                                 vLength = vLength,
                                 recruitZ = recruitZ,
                                 eggMass = eggMass,
                                 g = g,
                                 rickerA = rickerA,
                                 rickerB = rickerB,
                                 h2lAsym = 0,
                                 mortRate = mortRate,
                                 mortMultiplier = mortMultiplier,
                                 mortExponent = mortExponent,
                                 ageCut = ageCut,
                                 t0 = t0,
                                 kStart = kStart,
                                 lReproMean = lReproMean,
                                 discardMortRate = discardMortRate,
                                 probCatch = 0.1,
                                 probHarv = probHarv,
                                 harvestFrame = protectedSlot,
                                 prefLength = prefLength,
                                 upperLengthLimit = upperLengthLimit,
                                 lowerLengthLimit = lowerLengthLimit,
                                 keepPops = TRUE,
                                 popPath = "pSlotLoNoEvo/")

write_rds(LoPslotNoEvo, path = "LoPSlotNoEvo.rds")
