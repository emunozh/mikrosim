#library(GREGWT)
setwd("/home/esteban/workspace/R/GREGWT/src/")
source("GREGWT.R")

setwd("/home/esteban/Documents/MAPS/Germany")

### Get the data from the micro census (2002)

mikro.raw = read.csv("./Data/Survey/mz02_cf.csv", sep=";")
# columns to keep for simulation:
# age, marital status, household size, weights
keep.simulation = c(
    "ef30",   # Age 
    "ef35",   # Marital status
    "ef521",  # Household size
    "ef750")  # Weights
mikro.simulation <- mikro.raw[names(mikro.raw) %in% keep.simulation]
# columns to keep for result:
# cold operating cost, warm operating cost
keep.result = c("ef464", "ef466")
mikro.result <- mikro.raw[names(mikro.raw) %in% keep.result]
mikro.result$ef464[mikro.result$ef464 == 9998]  <- NA
mikro.result$ef464[mikro.result$ef464 == 9999]  <- NA
mikro.result$ef464[mikro.result$ef464 == 8]  <- NA
mikro.result <- mikro.result$ef466 - mikro.result$ef464
# remove all observations with NaN values 
mikro.simulation <- mikro.simulation[complete.cases(mikro.result),]
mikro.result <- mikro.result[complete.cases(mikro.result)]

### Rearrange data to fit the Census 2011

# Age
age.01 <- vector(length=dim(mikro.simulation)[1])
age.02 <- vector(length=dim(mikro.simulation)[1])
age.03 <- vector(length=dim(mikro.simulation)[1])
age.04 <- vector(length=dim(mikro.simulation)[1])
age.05 <- vector(length=dim(mikro.simulation)[1])
age.06 <- vector(length=dim(mikro.simulation)[1])
age.07 <- vector(length=dim(mikro.simulation)[1])
age.08 <- vector(length=dim(mikro.simulation)[1])
age.09 <- vector(length=dim(mikro.simulation)[1])
age.10 <- vector(length=dim(mikro.simulation)[1])
age.11 <- vector(length=dim(mikro.simulation)[1])
age.01[mikro.simulation$ef30 < 3] = 1
age.02[mikro.simulation$ef30 < 6 & mikro.simulation$ef30 >= 3] = 1
age.03[mikro.simulation$ef30 < 15 & mikro.simulation$ef30 >= 6] = 1
age.04[mikro.simulation$ef30 < 18 & mikro.simulation$ef30 >= 15] = 1
age.05[mikro.simulation$ef30 < 25 & mikro.simulation$ef30 >= 18] = 1
age.06[mikro.simulation$ef30 < 30 & mikro.simulation$ef30 >= 25] = 1
age.07[mikro.simulation$ef30 < 40 & mikro.simulation$ef30 >= 30] = 1
age.08[mikro.simulation$ef30 < 50 & mikro.simulation$ef30 >= 40] = 1
age.09[mikro.simulation$ef30 < 65 & mikro.simulation$ef30 >= 50] = 1
age.10[mikro.simulation$ef30 < 75 & mikro.simulation$ef30 >= 65] = 1
age.11[mikro.simulation$ef30 >= 75] = 1

# Marital status
mst.01 <- vector(length=dim(mikro.simulation)[1])
mst.02 <- vector(length=dim(mikro.simulation)[1])
mst.03 <- vector(length=dim(mikro.simulation)[1])
mst.04 <- vector(length=dim(mikro.simulation)[1])
mst.01[mikro.simulation$ef35 == 1] = 1
mst.02[mikro.simulation$ef35 == 2] = 1
mst.03[mikro.simulation$ef35 == 3] = 1
mst.04[mikro.simulation$ef35 == 4] = 1

# Household size
hhs.01 <- vector(length=dim(mikro.simulation)[1])
hhs.02 <- vector(length=dim(mikro.simulation)[1])
hhs.03 <- vector(length=dim(mikro.simulation)[1])
hhs.04 <- vector(length=dim(mikro.simulation)[1])
hhs.05 <- vector(length=dim(mikro.simulation)[1])
hhs.06 <- vector(length=dim(mikro.simulation)[1])
hhs.01[mikro.simulation$ef521 == 1] = 1
hhs.02[mikro.simulation$ef521 == 2] = 1
hhs.03[mikro.simulation$ef521 == 3] = 1
hhs.04[mikro.simulation$ef521 == 4] = 1
hhs.05[mikro.simulation$ef521 == 5] = 1
hhs.06[mikro.simulation$ef521 >= 6] = 1

# put everything on a data frame
mikro.input = data.frame(
    age.01 = age.01,
    age.02 = age.02,
    age.03 = age.03,
    age.04 = age.04,
    age.05 = age.05,
    age.06 = age.06,
    age.07 = age.07,
    age.08 = age.08,
    age.09 = age.09,
    age.10 = age.10,
    age.11 = age.11,
    mst.01 = mst.01,
    mst.02 = mst.02,
    mst.03 = mst.03,
    mst.04 = mst.04,
    hhs.01 = hhs.01,
    hhs.02 = hhs.02,
    hhs.03 = hhs.03,
    hhs.04 = hhs.04,
    hhs.05 = hhs.05,
    hhs.06 = hhs.06)

# and the vector with the weights
dx <- mikro.simulation$ef750

### Get the data from the census (2011)
nan.strings = c('nan', '.')
gem.alt = read.csv("./Data/Gemeinden/ALTER_AF-all.csv",
    colClasses=c("character",rep("numeric",6)),
    na.strings = nan.strings)
gem.fam = read.csv("./Data/Gemeinden/FAMSTND_KURZ-all.csv",
    colClasses=c("character",rep("numeric",6)),
    na.strings = nan.strings)
gem.hhs = read.csv("./Data/Gemeinden/HHGROESS_KLASS-all.csv",
    colClasses=c("character",rep("numeric",7)),
    na.strings = nan.strings)

# Save the working space
save.image(file="Steps/01.RData")

# Select a single federal state (eg: 05 is the code for Nordrhein-Westfalen)
AGS.code = "05"
AGS.length = 2
gem.alt <- gem.alt[substr(gem.alt$X,1,AGS.length)==AGS.code, ]
gem.fam <- gem.fam[substr(gem.fam$X,1,AGS.length)==AGS.code, ]
gem.hhs <- gem.hhs[substr(gem.hhs$X,1,AGS.length)==AGS.code, ]

# Remove unwanted columns from the data
# age
drop <- c("Total")
gem.alt <- gem.alt[,!(names(gem.alt) %in% drop)]
# marital status
drop <- c("Total", "No.data")
gem.fam <- gem.fam[,!(names(gem.fam) %in% drop)]
# Household size
drop <- c("Total")
gem.hhs <- gem.hhs[,!(names(gem.hhs) %in% drop)]

# Merge all data into a big data frame
gem.input <- merge(gem.alt, gem.fam, by.x = "X", by.y = "X")
gem.input <- merge(gem.input, gem.hhs, by.x = "X", by.y = "X")

# define the number of areas to run
#areas.number = 4
areas.number = dim(gem.input)[1]

Result = data.frame(
    area=vector(length=areas.number),
    heat=vector(length=areas.number))

# Save the working space
save.image(file="Steps/02.RData")

# Prepare data for simulation
area.code <- gem.input[, 1]
Tx.s <- gem.input[, 2:dim(gem.input)[2]]
names(Tx.s) <- names(mikro.input)
Simulation.Data <- prepareData(mikro.input, Tx.s)
mikro.input.s <- Simulation.Data$X
Tx.s <- Simulation.Data$Tx

# loop through all areas 
for(i in seq(1, areas.number)){

    # Create a vector with the area totals
    Tx <- Tx.s[i,]

    # Store the area code
    acode <- area.code[i]

    # Get new weights with GREGWT
    Weights = GREGWT(mikro.input.s, dx, Tx, bounds=c(0,Inf))
    cat("estimated weights for:", acode, "\t", i, "/", areas.number, "\n")
    fw <- Weights$Final.Weights
    heat.expenditure <- sum(mikro.result * fw / sum(fw), na.rm=TRUE)
    Result[i,] <- c(acode, heat.expenditure)}

# Save the working space
save.image(file="Steps/03.RData")

# Load some data
# load("Steps/de.RData")

# Save the result
Result <- Result[Result$heat > 0, ]
write.csv(Result, file="SimulationResult.csv")

# Make some nice plots with the result
heat <- as.numeric(Result$heat)
jpeg(filename="HeatExpenditure.jpeg", width=600, height=600, quality=100)
plot(sort(heat),
     main="Heat expenditure in German municipalities",
     ylab="Monthly heat expediture in EUR",
     xlab="Sorted municipalities")
abline(h=mean(heat, na.rm=TRUE), col='red', lw=3)
dev.off()


jpeg(filename="HeatExpenditureHist.jpeg", width=600, height=600)
hist(heat, main="Histogram of heat expediture in German municipalities")
dev.off()
