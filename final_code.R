rm(list=ls(all=TRUE))

library(jsonlite)
require("ggplot2")
library("ggplot2")
require("sqldf")
require("readr")

{r}
df_bus_To = stream_in(file("https://www.scss.tcd.ie/~arwhite/Teaching/CS7DS3/Yelp%20data/Business_Toronto_Restaurant.json"))
df_rv_To = stream_in(file("https://www.scss.tcd.ie/~arwhite/Teaching/CS7DS3/Yelp%20data/Review_Toronto_Restaurant.json"))
names(df_bus_To)
names(df_rv_To)

# Look at categories
class(df_bus_To$categories)
head(df_bus_To$categories)

cat_total <- unlist(df_bus_To$categories)
class(cat_total)
length(cat_total)

## Repeated non-numeric values should have 'factor' class
cat_total <- factor(cat_total)
nlevels(cat_total)

## Look at all categories together 
cat_total <- unlist(df_bus_To$categories)

## Which categories are most popular?
cat_names_sort <- sort(table(cat_total), decreasing = TRUE)
head(cat_names_sort, n = 25)
tail(cat_names_sort, n = 25)

## Bin plot for cateorories except 'restaurant' label
cat_total <- as.data.frame(cat_total)
cat_noRes <- as.data.frame(cat_total[cat_total$cat_total != 'Restaurants',])
colnames(cat_noRes) <- 'label'
p1 <- ggplot(data= cat_noRes,aes(x = reorder(as.numeric(label), label, length))) + stat_count() + theme(axis.title.x=element_blank(),
                                                                                                        axis.text.x=element_blank())
p1

cat_names <- names(cat_names_sort)[2:21] ## 1 is Restaurants - we don't need this

cat_bus_ind_mat <- sapply(df_bus_To$categories, function(y) as.numeric(cat_names %in% y))
cat_bus_ind_mat <- t(cat_bus_ind_mat)
colnames(cat_bus_ind_mat) <- cat_names
df_To_tidy_cat <- cbind(df_bus_To, cat_bus_ind_mat)

# See all listed neighborhood names
unique(df_bus_To$neighborhood) 
# Select neighborhood
sel_nbhd <- df_bus_To$neighborhood == "Scarborough" | df_bus_To$neighborhood == "Etobicoke"
# Select opening status
isOpen <- df_bus_To$is_open == TRUE
# Select Indian Restaurant
num_res = nrow(df_bus_To)
isIndianRest = dim(num_res)
for(i in 1:num_res){
  isIndianRest[i] = "Restaurants" %in% df_bus_To$categories[[i]] & "Indian" %in% df_bus_To$categories[[i]]
}

# Create data subset
df_OpenIndRes <- subset(df_bus_To, isIndianRest & isOpen, select =  c("stars", "neighborhood","is_open","categories"))
df_ck <- subset(df_bus_To, sel_nbhd & isIndianRest & isOpen, select =  c("stars", "neighborhood","is_open","categories"))
dim(df_OpenIndRes)
dim(df_ck)

## Make a new indicator variable for neighborhood
df_ck$nbd_ind <- as.numeric(factor(df_ck$neighborhood))
sum(df_ck$neighborhood != "Scarborough")

p1 <- ggplot(df_ck) + geom_boxplot(aes(neighborhood, stars, fill = neighborhood)) + geom_jitter(aes(neighborhood, stars, shape = neighborhood))
p1

p2 <- ggplot(df_ck) + geom_violin(aes(neighborhood, stars, fill = neighborhood)) + geom_jitter(aes(neighborhood, stars, shape = neighborhood))
p2
ggsave("C:/Users/Neil-PC/Desktop/compare_2_nbh_1.jpg",last_plot(),device = "jpeg")

tapply(df_ck$stars, df_ck$neighborhood,mean)
tapply(df_ck$stars, df_ck$neighborhood,median)
tapply(df_ck$stars, df_ck$neighborhood,sd)
tapply(df_ck$stars, df_ck$neighborhood,function(x) {return(1/var(x))})
t.test(stars ~ nbd_ind, data=df_ck, var.equal = TRUE)

compare_2_gibbs <- function(y, ind, mu0 = 3, tau0 = 1/0.64, del0 = 0, gamma0 = 1, a0 = 2, b0 = 1, maxiter = 5000)
{
  y1 <- y[ind == 1]
  y2 <- y[ind == 2]
  
  n1 <- length(y1) 
  n2 <- length(y2)
  
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  
  for(s in 1 : maxiter) 
  {
    
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    
    ##update mu
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    
    ##update del
    gamman <-  gamma0 + tau*(n1 + n2)
    deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

require(MCMCpack)

fit <- compare_2_gibbs(df_ck$stars, df_ck$nbd_ind)

y1_sim <- rnorm(5000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
y2_sim <- rnorm(5000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))

ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))
mean(y2_sim>y1_sim)
ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) + geom_abline(slope = 1, intercept = 0)

plot(as.mcmc(fit))
acf(as.mcmc(fit))
raftery.diag(as.mcmc(fit))

apply(fit, 2, mean)
apply(fit, 2, sd)
mean(1/sqrt(fit[, 3])) 

cat_total <- unlist(df_bus_To$categories)
class(cat_total)
cat_total <- factor(cat_total)
nlevels(cat_total)

## Look at all categories together 
cat_total <- unlist(df_bus_To$categories)

## Which categories are most popular?
cat_names_sort <- sort(table(cat_total), decreasing = TRUE)
head(cat_names_sort, n = 25)

cat_names <- names(cat_names_sort)[2:11] ## 1 is Restaurants - we don't need this
cat_bus_ind_mat <- sapply(df_bus_To$categories, function(y) as.numeric(cat_names %in% y))
cat_bus_ind_mat <- t(cat_bus_ind_mat)
colnames(cat_bus_ind_mat) <- cat_names
df_To_tidy_cat <- cbind(df_bus_To, cat_bus_ind_mat)

install.packages("rapportools")
library(rapportools)
# Select opening status
isOpen <- df_bus_To$is_open == TRUE
# Select entries with neighborhood information
hasNbd <- !is.empty(df_bus_To$neighborhood)
# Create data subset for open restaurants in all neighborhoods
df_OpenRes <- subset(df_bus_To, isOpen & hasNbd, select =  c("stars", "neighborhood"))
dim(df_OpenRes)

## Make a new indicator variable for neighborhood
df_OpenRes$nbd_ind <- as.numeric(factor(df_OpenRes$neighborhood))

ggplot(df_OpenRes) + geom_boxplot(aes(x = reorder(nbd_ind, stars, median), stars, fill = reorder(nbd_ind, stars, median)), show.legend=FALSE)
ggplot(df_OpenRes, aes(x = reorder(nbd_ind, nbd_ind, length))) + stat_count()
ggplot(df_OpenRes, aes(stars)) + stat_bin()

ggplot(data.frame(size = tapply(df_OpenRes$stars, df_OpenRes$neighborhood, length), mean_score = tapply(df_OpenRes$stars, df_OpenRes$neighborhood, mean)), aes(size, mean_score)) + geom_point()
reorder(df_OpenRes$neighborhood,df_OpenRes$neighborhood, length)
tapply(df_OpenRes$neighborhood, df_OpenRes$neighborhood, length)


Gibbs sampler to model the diff b/w the mean stars of opening restaurants
compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
  ### weakly informative priors
  a0 <- 2 ; b0 <- 1 ## tau_w hyperparameters
  eta0 <- 2 ; t0 <- 1 ## tau_b hyperparameters
  mu0<- 3 ; gamma0 <- 1/0.64
  ###
  
  ### starting values
  m <- nlevels(factor(ind))
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <- 1/var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta-mu)^2)/2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}
fit2 <- compare_m_gibbs(df_OpenRes$stars, df_OpenRes$nbd_ind)

plot(as.mcmc(fit2))
raftery.diag(as.mcmc(fit2))

apply(fit2$params, 2, mean)
apply(fit2$params, 2, sd)
mean(1/sqrt(fit2$params[, 3]))
sd(1/sqrt(fit2$params[, 3]))


theta_hat=apply(fit2$theta,2,mean)

ggplot(data.frame(size = tapply(df_OpenRes$stars, df_OpenRes$neighborhood, length), theta_hat = theta_hat), aes(size, theta_hat)) + geom_point()

fit_n <- as.mcmc(fit2$params)
plot(fit_n)
acf(fit_n)
raftery.diag(fit_n)

rest_review_data <- stream_in(file("https://www.scss.tcd.ie/~arwhite/Teaching/CS7DS3/Yelp%20data/Review_Toronto_Restaurant.json"), flatten = TRUE)
rest_data <- stream_in(file("https://www.scss.tcd.ie/~arwhite/Teaching/CS7DS3/Yelp%20data/Business_Toronto_Restaurant.json"), flatten = TRUE)

head(rest_data)
dim(rest_data)
str(rest_data)


#Part 2
#preping data for model
imp_review_count <- sqldf('select business_id, sum(useful) as total_useful_reviews, count(useful) as total_reviews,
                          avg(useful) as average_total_Reviews
                          from rest_review_data
                          group by business_id')

rest_data = merge(rest_data, imp_review_count, by="business_id", all = T)


#removing identifier columns and dependent variable
pca_data = subset(rest_data, select = -c(business_id, name, address, is_open, categories))

#check available variables
colnames(pca_data)

#check variable class
str(pca_data)

library(dummies)
library(dplyr)

pca_data$review_count = as.numeric(pca_data$review_count)
pca_data$total_reviews = as.numeric(pca_data$total_reviews)
pca_data$total_useful_reviews = as.numeric(pca_data$total_useful_reviews)

logical_cols = colnames(pca_data[, sapply(pca_data, class) == "logical"])
pca_data[c(logical_cols)][is.na(pca_data[c(logical_cols)])] <- FALSE

#convert all logical to numeric
cols <- sapply(pca_data, is.logical)
pca_data[,cols] <- lapply(pca_data[,cols], as.numeric)
str(pca_data)
pca_data <- subset(pca_data, select = -c(hours.Monday, hours.Tuesday, hours.Wednesday, hours.Thursday,
                                        hours.Friday, hours.Saturday, hours.Sunday, attributes.AgesAllowed,
                                        attributes.Smoking, attributes.WiFi, attributes.RestaurantsAttire,
                                        attributes.NoiseLevel, attributes.Alcohol, neighborhood,
                                        city, state, postal_code))
str(pca_data)
pca_data$attributes.RestaurantsPriceRange2 = as.numeric(pca_data$attributes.RestaurantsPriceRange2)

library(zoo)
pca_data = na.aggregate(pca_data)
pca_data = Filter(var, pca_data)

#principal component analysis
prin_comp <- prcomp(pca_data, scale. = TRUE)
names(prin_comp)

#pLOTTING
biplot(prin_comp, scale = 1)

#Checking for dimensions
dim(prin_comp$x)
importance = as.array(prin_comp$center)
importance = sort(importance, decreasing = TRUE)
importance[1:25]

names_imp = names(importance[1:25])

is_open = rest_data$is_open
data_set = subset(pca_data, select = c(names_imp))
data_set = cbind(is_open, data_set)
data_set = data_set[complete.cases(data_set), ]
dat_new <- subset(data_set, select = -is_open)
head(dat_new)
dat_new <- as.data.frame(apply(dat_new, 2, scale))
dat_new <- cbind(is_open = data_set$is_open, dat_new)

#Using GLM
glm1 <- glm(is_open ~ ., data = dat_new, family = binomial())
head(predict(glm1)) ## values on log-odds scale
pred_glm <- plogis(predict(glm1)) ## this is logistic function, maps to [0,1]
boxplot(pred_glm  ~ dat_new$is_open)
table(pred_glm > 0.525, dat_new$is_open) ## 0.5 is an arbitrary threshold






library(jsonlite)
require("ggplot2")
library("ggplot2")
require("sqldf")
require("readr")
library("tidytext")
library("lattice")
library("udpipe")
library("ggrepel")
library("syuzhet")
library("reshape2")
library(dplyr)

reviews_from_data=sqldf("Select text from rest_review_data")

#assigning sentiment score
sentiment_score=get_sentiment(reviews_from_data$text)

#most positive and negative review
#most.positive_review=reviews_from_data$text[sentiment_score == max(sentiment_score)]
#most.negative_review=reviews_from_data$text[sentiment_score == min(sentiment_score)]


#categorizing wrt to positive, negative and neutral
#positive.reviews=reviews_from_data$text[sentiment_score>0]
#negative.reviews=reviews_from_data$text[sentiment_score<0]
#neutral.reviews=reviews_from_data$text[sentiment_score=0]

#categorizing reviews and adding the sentiment analysis to respective store
category_review_sentiments=ifelse(sentiment_score<0,"Negative",ifelse(sentiment_score>0,"Positive","Neutral"))
rest_review_data=cbind(rest_review_data,category_review_sentiments)

review_useful <- rest_review_data[c('business_id','useful','category_review_sentiments')]
review_useful$category_review_sentiments = ifelse(review_useful$category_review_sentiments == 'Negative', 0, 1)
review_useful$positive_useful <- review_useful$useful * review_useful$category_review_sentiments
review_useful$negative_useful <- ifelse(review_useful$category_review_sentiments == 0, review_useful$useful, 0)
summary(review_useful)

useful_pos_neg <- sqldf('select business_id, sum(positive_useful) as total_positive_useful, 
                          sum(negative_useful) as total_negative_useful, 
                          avg(positive_useful) as average_positive_useful,
                          avg(negative_useful) as average_negative_useful
                          from review_useful
                          group by business_id')

rest_data = merge(rest_data, useful_pos_neg, by="business_id", all = T)
 
ggplot(rest_data[(!is.na(rest_data$is_open) & !is.na(rest_data$total_positive_useful)),], aes(x = is_open, fill = total_positive_useful)) +
  geom_density(alpha=0.5, aes(fill=factor(total_positive_useful))) + labs(title=" ") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + theme_grey()

ggplot(rest_data[(!is.na(rest_data$is_open) & !is.na(rest_data$total_negative_useful)),], aes(x = is_open, fill = total_negative_useful)) +
  geom_density(alpha=0.5, aes(fill=factor(total_negative_useful))) + labs(title=" ") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + theme_grey()







#part 3
data_for_cluster = rest_data[c('latitude','longitude')]

library(mclust)
set.seed(100)
fit1 <- Mclust(data_for_cluster, G=1:9)
fit2 <- Mclust(data_for_cluster, G=10:19)
fit3 <- Mclust(data_for_cluster, G=20:29)

plot(fit1, what = "BIC")
plot(fit2, what = "BIC")
plot(fit3, what = "BIC")

fit1$BIC
#VVV,9    VVV,7    VVE,9 
#58460.61 57845.81 57778.13 

fit2$BIC
#VVV,16   VVV,19   VVV,18 
#63977.61 63515.55 63324.17 

fit3$BIC
#VVV,29   VVV,28   VVV,27 
#63645.02 63589.34 63538.96 


#plotting clusters for best one!
fit2 <- Mclust(data_for_cluster, G = 16, modelNames = "VVV")
plot(fit2, what = "classification")
plot(fit2, what = "uncertainty")

#data_closed = rest_data$text[sentiment_score=0]
clust_class = fit2$classification
clust_uncertainty = fit2$uncertainty
#clust_class
data_for_cluster$cluster_class = clust_class
data_for_cluster$cluster_uncertainty = clust_uncertainty

data_for_cluster$neighborhood = rest_data$neighborhood
data_for_cluster$restraurants = rest_data$business_id

Cluster_records = sqldf('select cluster_class as Cluster ,
                                count(distinct neighborhood) as Number_neighborhoods,
                                count(distinct restraurants) as Number_Restaurants
                                from data_for_cluster
                                group by cluster_class')
Cluster_records
rest_in_cluster = sqldf('select cluster_class as Cluster ,
                                count(distinct restraurants) as Number_Restaurants
                                from data_for_cluster
                                group by cluster_class')


data_for_cluster = rest_data[c('latitude','longitude')]

#plotting clusters for second best one!
fit2 <- Mclust(data_for_cluster, G = 29, modelNames = "VVV")
plot(fit2, what = "classification")
plot(fit2, what = "uncertainty")

#data_closed = rest_data$text[sentiment_score=0]
clust_class = fit2$classification
clust_uncertainty = fit2$uncertainty
#clust_class
data_for_cluster$cluster_class = clust_class
data_for_cluster$cluster_uncertainty = clust_uncertainty

data_for_cluster$neighborhood = rest_data$neighborhood
data_for_cluster$restraurants = rest_data$business_id

Cluster_records = sqldf('select cluster_class as Cluster ,
                                count(distinct neighborhood) as Number_neighborhoods,
                                count(distinct restraurants) as Number_Restaurants
                                from data_for_cluster
                                group by cluster_class')
Cluster_records
