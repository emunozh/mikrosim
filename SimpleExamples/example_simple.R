# load the functions
library('GREGWT')

## Data ######
# Variables:
#   1. age          (1=16-30 years and 0='otherwise')
#   2. sex          (1=female and 0=male)
#   3. employment   (1=unemployed and 0= ‘otherwise’)
#   4. income       (in real unit values 0, 1, 2, 3, 4 and 5)
#   5. location     (1= rural and 0= urban)
X <- data.frame( 
    age =        c(1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,0,0,1,0,1),
    sex =        c(1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,0,0,0,1,1,0),
    employment = c(0,1,1,1,0,1,0,1,0,0,1,1,1,1,1,0,1,0,1,1,0,1,0,1,0),
    income =     c(0,3,2,5,0,1,0,4,0,0,1,3,2,5,4,0,3,0,2,4,0,5,0,1,0),
    location =   c(0,1,1,0,1,0,1,0,1,1,0,1,1,1,0,0,1,0,1,0,1,1,0,0,1))

# Initial weights
dx <- c(4,5,6,5,3,4,6,4,5,3,5,4,3,6,4,5,6,3,6,4,5,3,5,4,3)

# True population totals
Tx <- data.frame(age=50,sex=45,employment=70,income=200,location=65)

# Get new weights with GREGWT
Weights.GREGWT = GREGWT(X, dx, Tx, bounds=c(0,Inf))

# display new weights
#print(Weights.New$TAD)

summary(Weights.GREGWT)

X["Initial.Weights"] <- dx
#View(dx)
X["Final.Weights"] <- Weights.GREGWT$Final.Weights
#View(Weights.New$Final.Weights)

#library("Hmisc")
#latex(X, file="simpleExample01.tex")
#latex(dx, file="simpleExample01dx.tex")
#latex(Tx, file="simpleExample01Tx.tex")
