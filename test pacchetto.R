rm(list = ls())
options(digits = 5)


# Get the current R file location
current_file_location <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Set the working directory to the current file location
setwd(current_file_location)

#installazione libreria dalla cartella
install.packages("gdpoisson_0.2.0.tar.gz", repos = NULL, type = "source")

library('gdpoisson')

?pgdpois
?dgdpois
#CDF
pgdpois(0:10,lambda=5, theta=3)

#PMF
dgdpois(0:10,lambda=5, theta=3)

dgdpois(1000:1010,lambda=5, theta=3)

#underdispersion
dgdpois(0:10,lambda=5.3, theta=0.2)

#poisson case
dgdpois(0:10,lambda=5, theta=1)
dpois(0:10,lambda=5)

#ideal underdispersion behavior for any mean as theta goes to 0
round(dgdpois(20:30, lambda=22.8, theta=0.0001),4)
round(dgdpois(1000:1010, lambda=1002.3, theta=0.0001),4)
round(dgdpois(1000:1010, lambda=1002.3, theta=0.00001),4)


#quantile
?qgdpois
qgdpois(p=0.99,lambda=5,theta=6)

#campione
?rgdpois
dgdpois(0:10,lambda=10, theta=0.4)

#sottodisperso
(sample1=rgdpois(100,10,0.4))
mean(sample1)
(sample3=rgdpois(100,1000,0.01))
mean(sample3)

#poisson
(sample2=rgdpois(100,10,1))
mean(sample2)
#sovradisperso
(sample3=rgdpois(100,10,8))
mean(sample3)




#glm from simulated data (simulated with distribution)
?glm.gdp


n <- 300
x1sim <- rnorm(n,14,3)
x2sim <- runif(n)
beta_true <- c(0.8, 0.3, -0.5)
theta_true <- 4

# Simulate response variable y
data <- data.frame(x1 = x1sim, x2 = x2sim)
X <- model.matrix(~ x1 + x2, data)
eta <- X %*% beta_true
lambda <- exp(eta)
lambda

y <- sapply(lambda, function(l) rgdpois(1, l , theta_true) )
data$y <- y
y

poisglm = glm(y ~ x1 + x2, data=data, family=poisson(link = "log"))
summary(poisglm)
sum((residuals(poisglm, type="pearson")^2))/ (n - 4)

gdpoismodel <- glm.gdp(y ~ x1 + x2, data = data, max_retries = 10)

summary(gdpoismodel)

#added variable
x3sim <- rnorm(n,14,3)
data$x3 <- x3sim
gdpoismodel2 <- glm.gdp(y ~ x1 + x2+x3, data = data, max_retries = 10)
summary(gdpoismodel2)


#null model
gdpoismodelnull <- glm.gdp(y ~ 1, data = data, max_retries = 10)

summary(gdpoismodelnull)


#other tests
gammamodel <- glm(y ~ x1 + x2, data = data, family =  Gamma(link = "log"))
summary(gammamodel)
normalmod <- glm(y ~ x1 + x2, data = data)
summary(normalmod)


COMPoisglm <- glm.cmp(y ~ x1 + x2, data = data)
summary(COMPoisglm)

?glm.cmp
library(boot)


poisglm = glm(y ~ x, data=cloth, family=poisson(link = "log"))
summary(poisglm)

qpoisglm = glm(y ~ x, data=cloth, family=quasipoisson(link = "log"))
summary(qpoisglm)
logLik(poisglm)

gdpoismodel <- glm.gdp(y ~ x, data=cloth, max_retries = 10)

summary(gdpoismodel)

sum((residuals(poisglm, type="pearson")^2))/ (n - 4)




#esempio con dati sovradispersi da un esame
load("dati_esame_m2.RData")

poisglm1 = glm(naffairs ~ .-kids-educ, data=Affairs, family=poisson(link = "log"))
summary(poisglm1)

sum((residuals(poisglm1, type="pearson")^2))/ (poisglm1$df.residual)

#la convergenza puo' richiedere alcuni secondi
gdglm1= glm.gd(naffairs ~ .-kids-educ, data=Affairs)
summary(gdglm1)

#simili livelli di significativita con quasi poisson
quasipoisglm= glm(naffairs ~ .-kids-educ, data=Affairs, family=quasipoisson(link = "log"))
summary(quasipoisglm)

#AIC
poisglm1$aic
gdglm1$aic


#compoisson

#install.packages("devtools")
#devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)
COMPoisglm <- glm.cmp(naffairs ~ . - kids - educ, data = Affairs)
summary(COMPoisglm)
?glm.cmp
?mpcmp

sum(pcomp(0:200, 20,5, lower.tail = FALSE))
dcomp(0:200, 20,0.3)

