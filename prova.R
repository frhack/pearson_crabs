library("dplyr")
library("ggpubr")
library(splitstackshape)
library(data.table)
library(ggprism)
library(MixtureInf)
library(mixtools)
library(GGally, quietly = TRUE)
library(tidyverse)
library(mclust)
x <- data(pearson)
#pearson['New Column'] = 123

y <- expandRows(pearson, "freq")

  my_data = setDT(expandRows(pearson, "freq"))[,][]

my_data$ratio = jitter(my_data$ratio + 0.0020, amount = 0.0020)

my_data$scaled = scale(my_data$ratio)

#my_data[1] <- data.frame(lapply(my_data[1], jitter( amount = 0.00005)))
#my_data <- ToothGrowth

#dplyr::sample_n(my_data, 10)

#barplot(table(cut(df1$Time, breaks = seq(0, 10, by = 2))))

shtest <-  shapiro.test(my_data$ratio)

vvv = shtest["p.value"]

pval = unlist(shtest["p.value"])
pmethod = unlist(shtest["method"])

shtest_tibble <- tibble::tribble(
  ~group1,~p.method,   ~p.format,
  "1",    pmethod,  format.pval(pval, digits=4 )
)


gghistogram(my_data$ratio,binwidth =0.0040, 
            add_density = TRUE,
            y= "..density..",
            add = "mean",
            color = "#00AFBB",
            fill = "#00AFBB", palette = c("#00AFBB", "#E7B800"),)


ggboxplot(my_data, y = "ratio",title = "Crabs", ylab = "Ratio", palette = "jco")


ggboxplot(my_data, y = "ratio",title = "Crabs", ylab = "Ratio", palette = "jco", width=0.5)+

  add_pvalue(
    shtest_tibble, y = 0.705,x = 1,
    label="{p.method}"
  )+
  add_pvalue(
    shtest_tibble,  y = 0.700,x = 1,
    label="p = {p.format}"
  )

library("ggpubr")
ggdensity(my_data$ratio, 
          main = "Density plot of tooth length",
          xlab = "Tooth length")


p = ggqqplot(my_data$ratio)

palette = c("#00AFBB", "#E7B800")

set_palette(p, palette)




x <- crabs.df$ratio
m <- c(mean(x), mean(x^2), mean(x^3), mean(x^4), mean(x^5), mean(x^6) )

PearsonMOM <- function(moms) {
  cums = toCumulants(moms)
  
  k1 = cums[1]
  
  k2 = cums[2]
  
  k3 = cums[3]
  
  k4 = cums[4]
  
  k5 = cums[5]
  
  k6 = cums[6]
  
  sols = 0
  
  if (k3 == 0 & k5 == 0) {
    if (k4 == 0) {
      print(
        cat(
          "Data resembles a single Gaussian with mean ",
          k1,
          "
and variance ",
          k2,
          ". No honest mixture.",
          "\n"
        )
      )
      
      sols = 1
      
    }
    else {
      if (k4 > 0) {
        print(cat("Data is consistent only with an equal means model."
                  , "\n"))
      }
      else{
        print(
          cat(
            "Data displays symmetry, different means alternative
still explored.",
            "\n"
          )
        )
      }
      e2 = k2 ^ 2 - k4 / 3 + k2 * k6 / (5 * k4)
      
      e1 = 2 * k2 + k6 / (5 * k4)
      
      d = e1 ^ 2 - 4 * e2
      
      if ((e1 > 0 & e2 > 0) & d > 0) {
        v1 = (e1 - sqrt(d)) / 2
        
        v2 = (e1 + sqrt(d)) / 2
        
        
        a = (v2 - k2) / (v2 - v1)
        
        print(cat(
          "(alpha,mu1,mu2,sigma1,sigma2)=
",
          c(a, k1, k1, sqrt(v1), sqrt(v2)),
          "\n"
        ))
        
        sols = 1
        
      }
      else {
        print(cat("Negative variance found, discarding equal means."
                  , "\n"))
      }
    }
  }
  rootsp = polyroot(
    c(
      -8 * k3 ^ 6,
      -32 * k3 ^ 4 * k4,
      -21 * k3 ^ 2 * k4 ^ 2 - 24 * k3 ^ 3 * k5,
      96 * k3 ^ 4 + 9 * k4 ^ 3 - 36 * k3 * k4 * k5,
      148 * k3 ^ 2 * k4 - 6 * k5 ^ 2,
      24 * k3 * k5 + 30 * k4 ^ 2,
      12 * k3 ^ 2,
      28 * k4,
      0,
      8
    )
  )
  
  rootsord = rootsp[sort.list(abs(Im(rootsp)))]
  
  cat("Pearson’s polynomial roots: ", rootsord, "\n")
  
  rootsreal = subset(rootsord, abs(Im(rootsord)) < 0.01 & Re(rootsord) < 0)
  
  cat(
    "Pearson’s polynomial appears to have ",
    length(rootsreal),
    "negative real roots.",
    "\n"
  )
  
  if (length(rootsreal) > 0) {
    p = Re(rootsreal)
    
    s = (2 * k3 ^ 3 + 6 * k3 * k4 * p + 3 * k5 * p ^ 2 - 8 * k3 * p ^ 3) / (p *
                                                                              (4 * k3 ^ 2 + 3 * k4 * p + 2 * p ^ 3))
    
    m1 = (s - sqrt(s ^ 2 - 4 * p)) / 2
    
    m2 = (s + sqrt(s ^ 2 - 4 * p)) / 2
    
    R1 = p + k2
    
    R2 = (k3 / p + s) / 3
    
    var1 = R1 - m1 * R2
    
    var2 = R1 - m2 * R2
    
    sigma1 = sqrt(ifelse(var1 >= 0, var1, NA))
    
    sigma2 = sqrt(ifelse(var2 >= 0, var2, NA))
    
    alpha = m2 / (m2 - m1)
    
    mu1 = m1 + k1
    
    mu2 = m2 + k1
    
    sixth = Inf
    
    
    print(sixth)
    for (i in 1:length(rootsreal)) {
      if (is.na(sigma1[i]) | is.na(sigma2[i])) {
        cat("Negative variance found, removing root.", "\n")
      }
      else{
        sixth[i] = abs(mixtmoms(c(
          alpha[i], mu1[i], mu2[i],
          sigma1[i], sigma2[i]
        ))[6] - moms[6])
        cat(
          "(i,alpha,mu1,mu2,sigma1,sigma2, sixth[i])=",
          c(i, alpha[i], mu1[i], mu2[i], sigma1[i], sigma2[i], sixth[i]),
          "\n"
        )
        
        
        
        
        sols = sols + 1
        
      }
    }
  }
  if (sols > 1) {
    j = which.min(sixth) 
    
    print(j)
    print(sixth)
    cat(
      "Of the ",
      sols,
      " statistically meaningful solutions, the closest to the sample’s sixth moment is",
      "(alpha,mu1,mu2,sigma1,sigma2)=
",
      c(alpha[j], mu1[j], mu2[j], sigma1[j], sigma2[j]),"\n");
    return(c(alpha[j], mu1[j], mu2[j], sigma1[j], sigma2[j]))
    
  } else{
    if (sols == 1) {
      print("Unique statistically meaningful solution found.")
    }
    else {
      print("No solutions, the data does not come from a mixture
of two Gaussians.")
    }
  }
}

toCumulants <- function(moms) {
  m1 = moms[1]
  m2 = moms[2]
  m3 = moms[3]
  m4 = moms[4]
  m5 = moms[5]
  m6 = moms[6]

  k1 = m1
  k2 = m2 - m1 ^ 2
  k3 = m3 - 3 * m1 * m2 + 2 * m1 ^ 3
  k4 = m4 - 4 * m1 * m3 - 3 * m2 ^ 2 + 12 * m1 ^ 2 * m2 - 6 * m1 ^ 4
  k5 = m5 - 5 * m1 * m4 - 10 * m2 * m3 + 20 * m1 ^ 2 * m3 + 
       30 * m1 * m2 ^ 2 - 60 * m1 ^ 3 * m2 + 24 * m1 ^ 5
  
  k6 = m6 - 6 * m1 * m5 - 15 * m2 * m4 + 30 * m1 ^ 2 * m4 - 10 * m3 ^ 2 + 
       120 * m1 * m2 * m3 - 120 * m1 ^ 3 * m3 + 30 * m2 ^ 3 - 270 * m1 ^ 2 * m2 ^ 2 + 360 * m1 ^ 4 * m2 - 120 * m1 ^ 6
  
  cums = c(k1, k2, k3, k4, k5, k6)
  
  return(cums)
}


mixtmoms <- function(params) {
  alpha = params[1]
  mu1 = params[2]
  mu2 = params[3]
  sigma1 = params[4]
  sigma2 = params[5]
  
  m1 = alpha * mu1 + (1 - alpha) * mu2
  
  m2 = alpha*(mu1^2 + sigma1^2) + (1 - alpha) * (mu2^2+ sigma2^2)
  
  m3 = alpha * (mu1 ^ 3 + 3 * mu1 * sigma1 ^ 2) + 
    (1 - alpha) * (mu2 ^ 3 + 3 * mu2 * sigma2 ^ 2)
  
  m4 = alpha * (mu1 ^ 4 + 6 * mu1 ^ 2 * sigma1 ^ 2 + 3 * sigma1 ^ 4) + 
    (1 - alpha) * (mu2 ^ 4 + 6 * mu2 ^ 2 * sigma2 ^ 2 + 3 * sigma2 ^ 4)
  
  m5 = alpha * (mu1^5 + 10*mu1^3 * sigma1^2 + 15 * mu1 * sigma1 ^ 4) + 
    (1 - alpha) * (mu2 ^ 5 + 10 * mu2 ^ 3 * sigma2 ^ 2 + 15 * mu2 * sigma2^4)
  
  m6 = alpha * (mu1 ^ 6 + 15 * mu1 ^ 4 * sigma1 ^ 2 + 45 * mu1 ^ 2 * sigma1 ^ 4 + 15 * sigma1 ^ 6) + 
    (1 - alpha) * (mu2 ^ 6 + 15 * mu2 ^ 4 * sigma2 ^ 2 + 45 * mu2 ^ 2 * sigma2 ^4 + 15 * sigma2 ^ 6)
  
  moms = c(m1, m2, m3, m4, m5, m6)
  
  return(moms)
}

m

rr = PearsonMOM(m)
rr

if (!require('pacman')) install.packages('pacman'); library('pacman')

lambda1  = rr[1]
mu1 = rr[2];
mu2 = rr[3];
sigma1 = rr[4];
sigma2 = rr[5];


ppp= ggplot(crabs.df, aes(x = ratio)) +
  #geom_line(stat="density")+  
  geom_line(stat="count")+
  #scale_y_continuous(labels=scales::)+
  #scale_y_continuous(labels = percent_format(accuracy = 1))+
  #scale_y_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24))+
  #geom_histogram(binwidth = 0.004) +
  #geom_histogram(aes(y = ..density..)) + geom_density()
  geom_density(aes(y = ..density.. * (1000 * 0.004)), col = 2)+
mapply(
  function(mean, sd, lambda,col, n, binwidth) {
    stat_function(
      fun = function(x) {
        (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda 
      },
      color = col, linetype = "dashed"
    )
  },
  mean = c(mu1,mu2), #mean
  sd = c(sigma1,sigma2), #standard deviation
  lambda = c(lambda1,1-lambda1), #amplitude
  col = c("red","blue"),
  n = length(crabs.df$ratio), #sample size
  binwidth = 0.004 #binwidth used for histogram
)



ppp1= ggplot(crabs.df, aes(x = ratio)) +
  #geom_line(aes(by = ratio),stat = "prop")+  
  #geom_line(stat="count", adjust = 0.5)+
  geom_bar(stat="count", binwidth = 0.004,aes(adjust = 0.5))+
  
  #scale_y_continuous(labels=scales::)+
  #scale_y_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24))+
  #geom_histogram(binwidth = 0.004) +
  #geom_histogram(aes(y = ..density..)) + geom_density()
  geom_density(aes(y = ..density.. , col = 2))+
  mapply(
    function(mean, sd, lambda,col, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd))  * lambda 
        },
        color = col, linetype = "dashed"
      )
    },
    mean = c(mu1,mu2), #mean
    sd = c(sigma1,sigma2), #standard deviation
    lambda = c(lambda1,1-lambda1), #amplitude
    col = c("red","blue"),
    n = length(crabs.df$ratio), #sample size
    binwidth = 0.004 #binwidth used for histogram
  )


ppp2 = ggplot(crabs.df, aes(x = ratio)) +
  #geom_line(aes(by = ratio),stat = "prop")+  
  #geom_line(stat="count", adjust = 0.5)+

  #scale_y_continuous(labels=scales::)+
  #scale_y_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24))+
  #geom_histogram(binwidth = 0.004) +
  #geom_histogram(aes(y = ..density..)) + geom_density()
  geom_density()+
  mapply(
    function(mean, sd, lambda,col, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd))  * lambda 
        },
        color = col, linetype = "dashed"
      )
    },
    mean = c(mu1,mu2), #mean
    sd = c(sigma1,sigma2), #standard deviation
    lambda = c(lambda1,1-lambda1), #amplitude
    col = c("red","blue"),
    n = length(crabs.df$ratio), #sample size
    binwidth = 0.004 #binwidth used for histogram
  )

cutp = seq(min(crabs.df$ratioj), max(crabs.df$ratioj), by=.001)
multdata = makemultdata(crabs.df$ratioj, cuts = cutp)

#kSelection = multmixmodel.sel(multdata,  comps = c(1,2,3,4,5,6,7))

#kSelection

mixplot = ggplot(crabs.df, aes(x = ratio)) +
  #geom_histogram(binwidth = 0.004) +
  mapply(
    function(mean, sd, lambda, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
        }
      )
    },
    mean = my_mix[["mu"]], #mean
    sd = my_mix[["sigma"]], #standard deviation
    lambda = my_mix[["lambda"]], #amplitude
    n = length(crabs.df$ratio), #sample size
    binwidth = 0.004 #binwidth used for histogram
  )



set.seed(12) #better make this reproducible
observations <- tibble(value = c(
  rnorm(n = 250, mean = 100, sd = 40), #the first normal distribution
  rnorm(n = 250, mean = 200, sd = 40),
  rnorm(n = 250, mean = 300, sd = 40)
)
)

#and a quick gander...
aaa = ggplot(observations, aes(x = value)) + 
  geom_histogram(binwidth = 0.05)

cut = seq(min(observations$value), max(observations$value), by=1)
multdata = makemultdata(observations$value, cuts = c(50,100,150,200,300))

kSelection = multmixmodel.sel(multdata,  comps = c(1,2,3,4,5,6))



x = crabs.df$ratioj
(hc1 <- hc(data = x, modelName = "V"))
## Call:
## hc(data = X, modelName = "VVV", use = "SVD") 
## 
## Model-Based Agglomerative Hierarchical Clustering 
## Model name        = VVV 
## Use               = SVD 
## Number of objects = 145

BIC1 <- mclustBIC(x, initialization = list(hcPairs = hc1)) # default 
plot(BIC1)


msEst <- mstep(modelName = "V", data = x,
               z = unmap(x))
names(msEst)
ret = em(modelName = msEst$modelName, data = x,
   parameters = msEst$parameters)





dens <- densityMclust(x,  modelNames = "V")
plot(dens$BIC)
print(summary(dens, parameters = TRUE))


br <- seq(min(x), max(x), length = 21)
plot(dens, what = "density", data = x, breaks = br)

print(dens$parameters)


irisBIC <- mclustBIC(crabs.df$ratioj)
plot(irisBIC)
# Clustering
mod1 =  me(modelName = "V",data = crabs.df$ratioj, z = unmap(crabs.df$ratioj))
summary(mod1)
#plot(mod1, what = c("classification"))


my_mix <- normalmixEM(observations$value, k = 3, epsilon = 1e-09) #telling it to find two gaussians in the observations
my_mix[["mu"]]

bbb = ggplot(observations, aes(x = value)) +
  geom_histogram(binwidth = 0.05) +
  mapply(
    function(mean, sd, lambda, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
        }
      )
    },
    mean = my_mix[["mu"]], #mean
    sd = my_mix[["sigma"]], #standard deviation
    lambda = my_mix[["lambda"]], #amplitude
    n = length(observations$value), #sample size
    binwidth = 0.05 #binwidth used for histogram
  )

#```{r}
densE = densityMclust(crabs.df$ratioj,  modelNames = "E", verbose = FALSE, plot = FALSE)
densV = densityMclust(crabs.df$ratioj,  modelNames = "V", verbose = FALSE, plot = FALSE)
#```
varianceV = sqrt(densV$parameters$variance$sigmasq)
meansV = densV$parameters$mean
probV = densV$parameters$pro
print(varianceV)
print(meansV)




varianceE = sqrt(densE$parameters$variance$sigmasq)
meansE = densE$parameters$mean
probE = densE$parameters$pro
print(varianceE)
print(meansE)

pppX= ggplot(crabs.df, aes(x = ratio)) +
  #geom_line(stat="density")+  
  #geom_line(stat="count")+
  #scale_y_continuous(labels=scales::)+
  #scale_y_continuous(labels = percent_format(accuracy = 1))+
  #scale_y_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24))+
  #geom_histogram(binwidth = 0.004) +
  #geom_histogram(aes(y = ..density..)) + geom_density()
  geom_density(aes(y = ..density.. * (1000 * 0.004)), col = 2)+
  mapply(
    function(mean, sd, lambda,col, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda 
        },
        color = col, linetype = "dashed"
      )
    },
    mean = c(mu1,mu2), #mean
    sd = c(sigma1,sigma2), #standard deviation
    lambda = c(lambda1,1-lambda1), #amplitude
    col = c("red","blue"),
    n = length(crabs.df$ratio), #sample size
    binwidth = 0.004 #binwidth used for histogram
  )+
  #geom_function(
  #  fun = function(x) {
  #    (lambda1*(dnorm(x, mean = mu1, sd = sigma1)) + (1-lambda1)* (dnorm(x, mean = mu2, sd = sigma2))) * 1000 * 0.004  
  #  },
  #  linetype = "dashed"
  #)+
geom_function(
  fun = function(x) {
    (0.25*(dnorm(x, mean = meansE[1], sd = varianceE)) + (0.75)* (dnorm(x, mean = mu2, sd = varianceE))) * 1000 * 0.004  
  },
  linetype = "dashed"
)

pppX

#ppp
#ppp1
#ppp2
#kSelection

#aaa
#bbb