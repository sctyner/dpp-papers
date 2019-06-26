# base distn for the population: 

means <- c(1.5175,1.5192,1.5220)

library(mixtools)

test <- rnormmix(100, lambda = c(.6, .2, .2), mu = means, sigma = 4e-5)


library(ggplot2)
qplot(x = test, geom = "density")

ggplot(data = glass_pop, aes(x = RI, weight = Count, y = ..density..)) + 
  geom_histogram(binwidth = .0002) + 
  geom_density(data = data.frame(x = test), inherit.aes = F, aes(x = x))


plot(ecdf(test))

library(ggfortify)

mix_fun <- function(x, mus, sigma, ws){
  ws[1]*dnorm(x, mus[1], sigma) + ws[2]*dnorm(x, mus[2], sigma)
   + ws[3]*dnorm(x, mus[3], sigma)
}

ggdistribution(mix_fun, x = seq(1.51,1.53, 0.00005), 
               ws= c(.6, .2, .2), mus = means, sigma = 4e-5)
p <- ggdistribution(dnorm, x= seq(1.51,1.53, 0.00005), mean = means[1], sd = 4e-5 ) 
p  
p2 <- ggdistribution(dnorm, x= seq(1.51,1.53, 0.00005), mean = means[2], sd = 4e-5 , p = p) 
p2
ggdistribution(dnorm, x= seq(1.51,1.53, 0.00005), mean = means[3], sd = 4e-5 , p = p2) 


p <- ggdistribution(pnorm, x = seq(1.51,1.53, 0.001), mean = 1.518458, sd = 4e-5)
p
p2 <- ggdistribution(pnorm, x = seq(1.51,1.53, 0.001), mean = 1.518458-.001, sd = 1e-4, p = p, colour = "red")

ggdistribution(pnorm, x = seq(1.51,1.53, 0.001), mean = 1.518458+.001, sd = 1e-4, p = p2, colour = "blue")


sample_mixnorm <- function(u, mus, sigmas, ws){
  if (u < ws[1]){
    rnorm(1, mean = mus[1], sd = sigmas[1])
  } else if (u < sum(ws[1:2]) & u >= ws[1]){
    rnorm(1, mean = mus[2], sd = sigmas[2])
  } else if (u >= sum(ws[1:2])) {
      rnorm(1, mean = mus[3], sd = sigmas[3])
  }
}


test <- data.frame(u = runif(1000))

test2 <- mutate(test, draws = map_dbl(u, sample_mixnorm, 
                                      mus = means, sigmas = rep(4e-5,3), ws = c(.6,.2,.2)))
test2 <- test2 %>% mutate(draws_bigger = map_dbl(u, sample_mixnorm, 
                         mus = means +.001 , sigmas = rep(4e-4,3), ws = rep(1/3,3)),
                         draws_smaller = map_dbl(u, sample_mixnorm, 
                                                mus = means -.003 , sigmas = rep(4e-4,3), ws = rep(1/3,3)))


test3 <- test2 %>% mutate(draws_ecdf = ecdf(draws)(draws),
                 draws_bigger_ecdf = ecdf(draws_bigger)(draws_bigger), 
                 draws_smaller_ecdf = ecdf(draws_smaller)(draws_smaller))

ggplot(data = test3) + 
  geom_line(aes(x = draws, y = draws_ecdf)) + 
  geom_line(aes(x = draws_bigger, y= draws_bigger_ecdf), color = "blue") + 
  geom_line(aes(x = draws_smaller, y= draws_smaller_ecdf), color = "red")
 

plot(ecdf(test2$draws))

qplot(test2$draws, binwidth = .00002)
means
