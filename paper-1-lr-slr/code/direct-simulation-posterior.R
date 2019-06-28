# simulating from the posterior that is "straight to the posterior" per Terenin & Draper 
# DP(n ,\hat{F}_n)
library(tidyverse)
library(dprocsim)
# population medians: g(theta) where theta = median(RI) in a pane
glass_pop <- read_csv("paper-1-lr-slr/dat/glass_pop_evett.csv")
#head(glass_pop)
glass_pop_all <- rep(glass_pop$RI, glass_pop$Count)
# number of obs in population
n_pop <- length(glass_pop_all)
# empirical cdf function 
ecdf_glass_pop <- ecdf(glass_pop_all)

# Step 1: simulate from posterior for g(theta)
# draws from the posterior for g(theta)
pop_cdf_sims <- dpprior_sim2(M = n_pop , F0 = draw_from_ecdf , sticks = 100, N = 500, ecdf = ecdf_glass_pop)

# plot draws from the posterior for g(theta)
pop_cdf_sims %>% unnest() %>% ggplot(aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(alpha = .05)

# repeated measurements from a single pane: 
data("bennett.df", package = "dafs")
glass_pane <- as.numeric(bennett.df[1,])
median_pane <- median(glass_pane)
glass_pane <- glass_pane - median_pane
n_pane <- length(glass_pane)
ecdf_glass_pane <- ecdf(glass_pane)

# Step 2: simulate from posterior for f(x-theta_0|theta=0)
# draws from the posterior for f(x|theta)
pane_cdf_sims <- dpprior_sim2(M = n_pane , F0 = draw_from_ecdf , sticks = 100, N = 500, ecdf = ecdf_glass_pane)

# plot draws from the posterior for g(theta)
pane_cdf_sims %>% unnest() %>% ggplot(aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(alpha = .05)

# places at which to evaluate f(x|theta) (when evaluating, subtract theta from each of these.)
crime_scene_sample <- read_csv("paper-1-lr-slr/dat/window.csv")
suspect_sample <- c(1.51848, 1.51850, 1.51848, 1.51844, 1.51846)

#### simulating & computing the numerator (integrate f(x,y| theta)*g(theta) wrt theta)

grid_theta <- seq(min(glass_pop_all), max(glass_pop_all), by = .00001)

g_eval_grid <- tibble(theta = grid_theta) %>% 
  mutate(dens_eval_gnum = map(theta, eval_dens, epsilon = .0001, sims = pop_cdf_sims))

head(g_eval_grid) 



g_eval_grid <- g_eval_grid %>% 
  mutate(dens_eval_num =  map_dbl(theta, dnorm, x = xbar, sd = 4e-5, log = FALSE))

g_eval_grid %>% summarize(integral = sum(dens_eval_gnum * dens_eval_num))

xbar_w <- mean(glass_pop_all)
xbar_sus <- mean(suspect_samp)


g_eval_grid <- g_eval_grid %>%
  mutate(dens_eval_den_w =  map_dbl(theta, dnorm, x = xbar_w, sd = 4e-5, log = FALSE),
         dens_eval_den_s =  map_dbl(theta, dnorm, x = xbar_sus, sd = 4e-5, log = FALSE))

g_eval_grid %>% 
  summarize(numerator = sum(dens_eval_gnum * dens_eval_num),
            den1 = sum(dens_eval_gnum * dens_eval_den_w),
            den2 = sum(dens_eval_gnum * dens_eval_den_s)) %>% 
  mutate(den = den1*den2, 
         lr = numerator/den)



g_eval_grid %>% filter(dens_eval_num != max(dens_eval_num)) %>% pull(dens_eval_num) %>% unique()

dnorm(x= xbar, mean = g_eval_grid$theta[1037], sd = 4e-5)

  