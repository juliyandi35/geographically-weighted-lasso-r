#Statistika Deskriptif
library(readxl)
Data.IPM<-read_excel("IPM.xlsx")
library(dplyr)
IPM<- Data.IPM %>% select(-`Kabupaten/Kota`) #menghapus kolom
summary(IPM)

# Uji Multikolinearitas
library(car)
regresi<-lm(formula=Y~X1+X2+X3+X4,data=IPM)
vif(regresi)

#Uji Normalitas
library(nortest)
lillie.test(IPM$Y)

#Uji Heterogenitas (Breusch-Pagan)
library(lmtest)
bptest(regresi)

library(GWmodel)
library(spDataLarge)
library(spdep)
#Penskalaan Data ke Dalam Data Frame Spasial#
Data.m<-as.matrix(IPM)
Data.s<-scale(Data.m)
Coords<-as.matrix(cbind(IPM$longitude,IPM$latitude))
Data.spdf<-SpatialPointsDataFrame(Coords,as.data.frame(Data.s))
bw.e.fix<-bw.gwr(Y~ X1+X2+X3+X4, data=Data.spdf,
                 approach="CV",kernel="exponential",adaptive=FALSE,
                 p=2,longlat=FALSE)
#Matriks Pembobot model GWR#
library(fields)
dii<-rdist(Coords)
dii
W<-exp(-dii/bw.e.fix)
W
#Pemodelan GWR untuk pengujian Persamaan (11) #
GWRModel.e.fix<-gwr.basic(Y~ X1+X2+X3+X4,data=Data.spdf,
                          bw=bw.e.fix,kernel="exponential",adaptive=FALSE,p=2,
                          F123.test=T)
GWRModel.e.fix

#Pencarian nilai VIF lokal#
Diag.GWR.e.fix<-gwr.collin.diagno(Y~X1+X2+X3+X4,
                                  data=Data.spdf,
                                  bw=bw.e.fix,kernel="exponential",adaptive=FALSE, p=2, theta=0)
names(Diag.GWR.e.fix)
Diag.GWR.e.fix$VIF



gwl.est<-function (form, locs, data, kernel = "exp", cv.tol)
{
  lhs<- as.character(form)[2]
  rhs<- as.character(form)[3]
  rhs.v<- strsplit(rhs, " + ", fixed = TRUE)
  n.l<- length(rhs.v[[1]])
  db<- data
  y <- db[, lhs]
  N <- dim(db)[1]
  X <- rep(1, N)
  for (i in 1:n.l) X <- cbind(X, db[, rhs.v[[1]][i]])
  library(fields)
  S <- rdist(locs)
  rmspe<- NA
  library(lars)
  band.ub<- ceiling(max(S))
  band.lb <- min(S) + 0.01 * band.ub
  gwr.rmse<- function(y, yhat){
    N <- length(y)
    sqrt(sum((y - yhat)^2) / N)
  }
  if (missing(cv.tol)) {
    lm1 <- lm(form, data = db)
    lm.rmse<- gwr.rmse(y, lm1$fitted.values)
    cv.tol<- lm.rmse * 0.05
  }
  gwl.bw.cv <- function(lb, ub, eps, X, y, S, kernel){
    a <- lb
    b <- ub
    c <- (a+b)/2
    diff<- b - a
    N <- dim(X)[1]
    gwl.cv.err<- function(phi, X, y, S, N, kernel){
      LARGE <- 1000000
      gwr.Gauss<-function (dist2, bandwidth)
      {
        w <- exp((-0.5) * ((dist2^2)/(bandwidth^2)))
        w
      }
      gwr.Exp<-function (dist2, bandwidth)
      {
        w <- exp((-1) * ((abs(dist2))/(bandwidth)))
        w
      }
      if (kernel == "exp") W <- gwr.Exp(S, phi)
      if (kernel == "gauss") W <- gwr.Gauss(S, phi)
      W.sqrt<- sqrt(W)
      p <- dim(X)[2]
      sol.i<- array(0, N)
      err.i<- array(0, N)
      frac.i<- array(0, N)
      bool.i<- array(0, dim=c(N,p))

      for (i in 1:N){
        W.i<- diag(W.sqrt[i,])
        W.i<- W.i[-i,-i]
        X.W <- W.i %*% X[-i,]
        y.W<- W.i %*% y[-i]
        lars.obj <- lars(X.W,y.W)
        betas<- lars.obj$beta

        sol.n<- dim(betas)[1]
        err.min<- LARGE
        sol<- 0
        frac<- 1
        bool<- rep(1,p)
        for (k in 2:sol.n){
          yhat<- X[i,] %*% betas[k,]
          err<- (y[i] - yhat)^2
          f <- sum(abs(betas[k,])) / sum(abs(betas[sol.n,]))
          if (err <err.min){
            err.min<- err
            sol<- k
            frac<- f
            bool<- ifelse(betas[k,] == 0, 0, 1)
          }
        }
        sol.i[i] <- sol
        err.i[i] <- err.min
        frac.i[i] <- frac
        bool.i[i,] <- bool
      }
      RMSPE <- sqrt(sum(err.i) / N)
      params<- list(RMSPE, sol.i, frac.i, bool.i)
      names(params) <- c("err", "sol", "frac", "bool")
      params
    }
    while (diff >eps){
      a.c<- (a+c)/2
      c.b<- (c+b)/2
      CV <- gwl.cv.err(a.c, X, y, S, N, kernel)
      RMSE.a.c<- CV$err
      CV <- gwl.cv.err(c.b, X, y, S, N, kernel)
      RMSE.c.b<- CV$err

      if (RMSE.a.c<RMSE.c.b){
        b <- c.b
        RMSE.b<- RMSE.c.b
      }
      if (RMSE.a.c>RMSE.c.b){
        a <- a.c
        RMSE.a<- RMSE.a.c
      }
      c <- (a+b)/2
      diff<- abs(b - a)
    }
    CV <- gwl.cv.err(lb, X, y, S, N, kernel)
    RMSE.lb <- CV$err
    sol.lb <- CV$sol
    frac.lb <- CV$frac
    bool.lb <- CV$bool
    CV <- gwl.cv.err(ub, X, y, S, N, kernel)
    RMSE.ub<- CV$err
    sol.ub<- CV$sol
    frac.ub<- CV$frac
    bool.ub<- CV$bool
    CV <- gwl.cv.err(c, X, y, S, N, kernel)
    RMSE.c<- CV$err
    sol.c<- CV$sol
    frac.c<- CV$frac
    bool.c<- CV$bool
    # Check bounds
    if (RMSE.lb <RMSE.c){
      c <- lb
      sol.c<- sol.lb
      RMSE.c<- RMSE.lb
      frac.c<- frac.lb
      bool.c<- bool.lb
    }
    if (RMSE.ub<RMSE.c){c <- ub}
    params<- list(c, sol.c, RMSE.c, frac.c, bool.c)
    names(params) <- c("phi", "sol", "RMSPE", "frac", "bool")
    params
  }
  g.bw <- gwl.bw.cv(band.lb, band.ub, cv.tol, X, y, S, kernel)
  bw<- g.bw$phi
  sol<- g.bw$sol
  frac<- g.bw$frac
  bool<- g.bw$bool
  rmspe<- g.bw$RMSPE
  params<- list(frac,bw, rmspe)
  names(params) <- c("fraco","phi", "RMSPE")
  params
}

GWLModel.e<-gwl.est(Y~ X1+X2+X3+X4,locs=Coords,
                    data=Data.spdf,kernel="exp")
GWLModel.e

#Pembobot model GWL#
W<-exp(-dii/GWLModel.e$phi)
W
