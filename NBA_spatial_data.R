##The Data

nba <- nba[!is.na(nba$converted_x), ]
nba <- nba[!is.na(nba$converted_y), ]

options(digits=9)
dat <- cbind(nba$converted_x, nba$converted_y, nba$result, nba$shot_distance, nba$type)
dat <- as.data.frame(dat[-1,])
names(dat) <- c("x", "y", "result", "distance","type")

#target variable
vec <- dat$result
levels(vec) <- c("1", "0")
dat$result <- as.numeric(as.character(vec))

dat$x <- as.numeric(as.character(dat$x))
dat$y <- as.numeric(as.character(dat$y))

#colocated variables
dat$distance <- as.numeric(as.character(dat$distance)) 

vec2 <- dat$type
levels(vec2) <- c("1", "2", "3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                  "18","19","20","21","22")
dat$type <- as.numeric(as.character(vec2))

#getting rid of multiple data points
dat <- aggregate(dat$type, list(x = dat$x, y = dat$y, result=dat$result, distance=dat$distance), mean)
colnames(dat)[5] <- "type"

head(dat)


##Non-Spatial Analysis

####Result Variable --- Binary Variable
library(fitdistrplus)
fit.norm <- fitdist(dat$result, "norm")
plot(fit.norm)

par(mfrow=c(2,2))
hist(dat$result)
boxplot(dat$result)
plot(ecdf(dat$result))


####Distance Variable
fit.norm <- fitdist(dat$distance, "norm")
plot(fit.norm)

par(mfrow=c(2,2))
hist(dat$distance)
boxplot(dat$distance)
plot(ecdf(dat$distance))
junk <- dev.off()


#####Performing transformations to see if they give more normal distribution

#####Applying a square root transformation
dist12 <- (dat$distance)^(1/2)
fit.norm <- fitdist(dist12, "norm")
par(mfrow=c(1,1))
plot(fit.norm)

#####Transforming to the 1/3 power 
dist2 <- (dat$distance)^(1/3)
fit.norm <- fitdist(dist2, "norm")
plot(fit.norm)

#####Appyling a square root transformation
dist12 <- (dat$distance)^(2)
fit.norm <- fitdist(dist12, "norm")
plot(fit.norm)

dat$distance <- dat$distance^(1/2)


####Type Variable
fit.norm <- fitdist(dat$type, "norm")
plot(fit.norm)

par(mfrow=c(2,2))
hist(dat$type)
boxplot(dat$type)
plot(ecdf(dat$type))


#####Performing transformations to give more normal distribution 

#####Applying a log transformation 
type2 <- log(dat$type)
par(mfrow=c(2,2))
hist(type2)
boxplot(type2)
plot(ecdf(type2), xlim=c(0,1.5))


#####Appyling a square root transformation
type2 <- dat$type^(1/2)
fit.norm <- fitdist(type2, "norm")
par(mfrow=c(1,1))
plot(fit.norm)

dat$type <- log(dat$type)


##Variograms
library(gstat)
summary(dat)

###Create gstat object
g <- gstat(id="result", formula = result~1, locations=~x+y, data=dat)

plot(variogram(g, alpha=c(0,45,90,135)))

q <- variogram(g)
plot(q)
m <- vgm(psill= .22, "Gau", range=22, nugget=.20)

###Fit a model to the variogram
v.fit <- fit.variogram(q, m, fit.method = 6)
plot(q, v.fit)

#Append
g1 <- gstat(g, id="distance", formula = distance~1, locations=~x+y, data=dat)
g1 <- gstat(g1, id="type", formula = type~1, locations=~x+y,  data=dat)

#Plot the sample variograms
vm <- variogram(g1)
plot(vm)

#Fit a model
vm.fit <- fit.lmc(vm, g1, model=v.fit, correct.diagonal = 1.01, fit.lmc = TRUE)
#vm.fit
plot(vm, vm.fit)

vm.fit


#PREDICTIONS 

x.range <- as.integer(range(dat[,1]))
y.range <- as.integer(range(dat[,2])) 
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=0.05), 
                   y=seq(from=y.range[1], to=y.range[2], by=0.05))
#Ordinary Kriging 
ok <- krige(id="result", formula=result~1, locations = ~x+y, data= dat, newdata=grd, model=m)

#Universal Kriging 
uk <- krige(id="result", formula=result~x+y, locations = ~x+y, model=m, data=dat, newdata=grd)

#Indicator Kriging (same as Ordinary Kriging)
ik <- krige(id="result", formula=result~1, locations=~x+y, model=v.fit, data=dat, newdata=grd)


##Cross Validation
#Ordinary Kriging CV
cv_ok <- krige.cv(result~1, data=dat, locations =~x+y, model=v.fit, nfold=nrow(dat))
PRESS_ok <- sum(cv_ok$residual^2)

#Universal Kriging CV
cv_uk <- krige.cv(result~x+y, data=dat, locations = ~x+y, model=v.fit, nfold=nrow(dat))
PRESS_uk <- sum(cv_uk$residual^2)

#Co-Kriging CV
store <- capture.output(cv_ck <- gstat.cv(vm.fit))
PRESS_ck <- sum(cv_ck$residual^2)

#Indicator Kriging CV
cv_ik <- krige.cv(result~1, data=dat, locations = ~x+y, model=v.fit, nfold=nrow(dat))
PRESS_ik <- sum(cv_ik$residual^2)

####Comparing the errors for different kriging

#####Ordinary Kriging
PRESS_ok

#####Universal Kriging
PRESS_uk

#####Co-Kriging
PRESS_ck

#####Indicator Kriging
PRESS_ik


##Raster Map
#Collapse prediction values to matrix
qqq <- matrix(ik$result.pred,length(seq(from=x.range[1], to=x.range[2], by=0.05)),
              length(seq(from=y.range[1], to=y.range[2], by=0.05)))

image(seq(from=x.range[1], to=x.range[2], by=0.05), seq(from=y.range[1],to=y.range[2], by=0.05),
      qqq, xlab="West to East",ylab="South to North", main="Raster map of the predicted values") 
points(dat)

#Adding contours
contour(seq(from=x.range[1], to=x.range[2], by=0.05), 
        seq(from=y.range[1], to=y.range[2], by=0.05), qqq, 
        add=TRUE, col="black", levels=seq(0, 50, by=.2), labcex=1)


##H-Scatterplot
library(sp)

#Convert the data into spatial data
coordinates(dat) <- ~x+y 

hscat(result~1, dat, c(0,5,10,15,20,25,30))
hscat(distance~1, dat, c(0,5,10,15,20,25,30))
hscat(type~1, dat, c(0,5,10,15,20,25,30))
