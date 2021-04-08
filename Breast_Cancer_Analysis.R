# Import libraries
library(gRbase)
library(gRim)
library(ggplot2)
library(gridExtra)
library(stats)
library(igraph)
library(mgm)
library(arm)
        
# Read the dataset
setwd("C:/Users/Utente/Desktop/UniMi/Probabilistic Modeling/Project")
breast_cancer <- read.csv(file = 'breast_cancer.csv')
breast_cancer <- breast_cancer[c(10,1,2,3,8)]
        
breast_cancer$Classification[breast_cancer$Classification == 1] <- "Healthy"
breast_cancer$Classification[breast_cancer$Classification == 2] <- "Patients"
breast_cancer$Classification <- factor(breast_cancer$Classification)
        
# Density distribution of continuous variables with respect to Classification
d1 <- ggplot(breast_cancer, aes(x = Age)) + geom_density(aes(fill = Classification), 
     alpha= 0.5) +  scale_fill_manual(values = c("#228b22","#ff8c00")) + theme_bw()
d2 <- ggplot(breast_cancer, aes(x = BMI)) + geom_density(aes(fill = Classification),
     alpha= 0.5) +  scale_fill_manual(values = c("#228b22","#ff8c00")) + theme_bw()
d3 <- ggplot(breast_cancer, aes(x = Glucose)) + geom_density(aes(fill = Classification), 
     alpha= 0.5) +  scale_fill_manual(values = c("#228b22","#ff8c00")) + theme_bw()
d4 <- ggplot(breast_cancer, aes(x = Resistin)) + geom_density(aes(fill = Classification), 
      alpha= 0.5) + scale_fill_manual(values = c("#228b22","#ff8c00")) + theme_bw()

x11()
grid.arrange(d1,d2,d3,d4, nrow = 2, ncol=2) 
        
SS <- CGstats(breast_cancer)
SS
        
apply(SS$center,1,sd)/apply(SS$center,1,mean)
        
# Saturated Model 
ms <- mmod(~Classification*Glucose*Resistin, data = breast_cancer)
can.parms <- ms$fitinfo$parms
can.parms
        
apply(can.parms$h,1,sd) / apply(can.parms$h,1,mean)
        
# Partial correlation matrix
pc <- cov2pcor(solve(can.parms$K))
pc
        
# Score-based selection
sat <- mmod(~.^., data = breast_cancer)
forw <- mmod (~.^1 , data = breast_cancer)
        
m_aic <- stepwise(sat, data = breast_cancer)
m_bic <- stepwise(sat, k=log(nrow(breast_cancer)), details = 0)
m_ci <- stepwise(sat, "test",alpha=0.05, details = 0, headlong = TRUE)
f_bic <- stepwise (forw ,k=log ( nrow (breast_cancer)), direction ="forward", details = 0)
f_aic <- stepwise (forw, direction ="forward", detail = 0)
        
x11()
par(mfrow=c(3,2))
plot(as(m_aic, "igraph"), main = "AIC",
             vertex.shape = "circle", vertex.label.font = 50, vertex.color = "white",
             vertex.label.color = "black", vertex.size = 65, edge.color = "black")
        
plot(as(m_bic, "igraph"), main = "BIC",
             vertex.shape = "circle", vertex.label.font = 65, vertex.color = "white",
             vertex.label.color = "black", vertex.size = 65, edge.color = "black")
        
plot(as(m_ci, "igraph"), main = "Conditional Independence",
             vertex.shape = "circle", vertex.label.font = 65, vertex.color = "white",
             vertex.label.color = "black", vertex.size = 65, edge.color = "black")
        
plot(as(f_aic, "igraph"), main = "Forward AIC",
             vertex.shape = "circle", vertex.label.font = 60, vertex.color = "white",
             vertex.label.color =  "black", vertex.size = 50, edge.color = "black")
        
plot(as(f_bic, "igraph"), main = "Forward BIC",
             vertex.shape = "circle", vertex.label.font = 60, vertex.color = "white",
             vertex.label.color = "black", vertex.size = 50, edge.color = "black")
        
# Deviance comparison
pchisq(m_bic$fitinfo$dev, m_bic$fitinfo$dimension["df"], lower.tail = FALSE)
pchisq(m_aic$fitinfo$dev, m_aic$fitinfo$dimension["df"], lower.tail = FALSE)
pchisq(m_ci$fitinfo$dev, m_ci$fitinfo$dimension["df"], lower.tail=FALSE)
pchisq(f_bic$fitinfo$dev , f_bic$fitinfo$dimension["df"], lower.tail = FALSE)
pchisq(f_aic$fitinfo$dev , f_aic$fitinfo$dimension["df"], lower.tail = FALSE)
        
        
# Prediction
breast_pred <- read.csv(file = 'breast_cancer.csv')
breast_pred <- breast_pred[c(10,c(1,2,3,4,5,6,7,8,9))]
        
breast_pred$Classification[breast_pred$Classification == 1] <- "Healthy"
breast_pred$Classification[breast_pred$Classification == 2] <- "Patients"
breast_pred$Classification <- factor(breast_pred$Classification)
        
log_reg <- bayesglm(Classification~., data = breast_pred, family = "binomial")
summary(log_reg)
        
smp <- floor(0.75 * nrow(breast_pred))
set.seed(123)
train_ind <- sample(seq_len(nrow(breast_pred)), size = smp)
train <- breast_pred[train_ind, ]
test <- breast_pred[-train_ind, ]
        
train$Classification <- ifelse(train$Classification == "Healthy" ,0 ,1)
test$Classification <- ifelse(test$Classification == "Healthy" ,0 ,1)
        
fit_mgm<- mgm ( data =train ,
                        type = c("c", rep("g" ,9)) ,
                        level = c(2, rep (1 ,9)) ,
                        k=2,
                        lambdaSel = "CV",
                        lambdaFolds = 10,
                        ruleReg = "OR")
        
f <- (fit_mgm$pairwise$wadj)
rownames(f) <- colnames(f) <- names(breast_pred)
f[1 ,]
x11()
qgraph :: qgraph (fit_mgm$pairwise$wadj , layout ="spring", repulsion = 1.2 ,
                          edge.color =fit_mgm$pairwise$edgecolor ,
                          nodeNames = colnames(train), legend = TRUE )
        
        
        
p_mgm <- predict ( fit_mgm , data = test ,
                           errorCon = c("RMSE", "R2"),
                           errorCat = c("CC", "nCC"))
p_mgm$errors
        