# simulating from the posterior that is "straight to the posterior" per Terenin & Draper 
# DP(n ,\hat{F}_n)
library(tidyverse)
glass_pop <- read_csv("paper-1-lr-slr/dat/glass_pop_evett.csv")
head(glass_pop)

glass_pop %>% slice(rep(1:n(), each = Count))

glass_pop_all <- rep(glass_pop$RI, glass_pop$Count)

n <- length(glass_pop_all)

ecdf_glass <- ecdf(glass_pop_all)


draw_from_ecdf <- function(ecdf, n){
  quantile(x = ecdf, probs = runif(n), type = 1)
}


library(dprocsim)
pop_cdf_sims <- dpprior_sim2(M = n , F0 = draw_from_ecdf , sticks = 100, N = 500, ecdf = ecdf_glass)

# draws from the posterior for g(theta)
pop_cdf_sims %>% unnest() %>% ggplot(aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(alpha = .05)

glass_window <- read_csv("paper-1-lr-slr/dat/window.csv")

n <- length(glass_window$RI)

ecdf_window <- ecdf(glass_window$RI)

window_cdf_sims <- dpprior_sim2(M = n, F0 = draw_from_ecdf, sticks = 100, N = 500, ecdf = ecdf_window)


# draws from the posterior for f(y | theta) (in the denominator, left hand side)
window_cdf_sims %>% unnest() %>% ggplot(aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(alpha = .05)

suspect_samp <- c(1.51848, 1.51850, 1.51848, 1.51844, 1.51846)

n <- length(suspect_samp)

ecdf_suspect <- ecdf(suspect_samp)

suspect_sims <- dpprior_sim2(M = n, F0 = draw_from_ecdf, 
                             sticks = 100, N = 500, ecdf = ecdf_suspect)


# draws from the posterior for f(x | theta) (in the denominator, right hand side)
suspect_sims %>% unnest() %>% ggplot(aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(alpha = .05)

# combined 

combined_glass <- c(suspect_samp, glass_window$RI)

n <- length(combined_glass)

ecdf_combined <- ecdf(combined_glass)

combined_sims <- dpprior_sim2(M = n, F0 = draw_from_ecdf, 
                             sticks = 100, N = 500, ecdf = ecdf_combined)

combined_sims %>% unnest() %>% ggplot(aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(alpha = .05)

ggplot(mapping = aes(x = locations, y = weCDF, group=rep))+ 
  geom_line(data = unnest(combined_sims), alpha = .05) + 
  geom_line(data = unnest(suspect_sims), alpha = .05, color = "red") + 
  geom_line(data = unnest(window_cdf_sims), alpha = .05, color = "blue")+ 
  geom_line(data = unnest(pop_cdf_sims), alpha = .05, color = "green")


#### simulating & computing the numerator (integrate f(x,y| theta)*g(theta) wrt theta)

n_tot <- length(combined_glass)
Fnum <- dpprior_sim(M = n_tot, F0 = draw_from_ecdf,sticks = 100, ecdf = ecdf_combined)

n_pop <- length(glass_pop_all)
gnum <- dpprior_sim(M = n_pop, F0 = draw_from_ecdf, sticks = 100, ecdf = ecdf_glass)
gnum <- mutate(gnum, rep = 1)

range(glass_pop_all)
range(combined_glass)

grid_theta <- seq(min(glass_pop_all), max(glass_pop_all), by = .00001)

g_eval_grid <- tibble(theta = grid_theta) %>% 
  mutate(dens_eval_gnum = map_dbl(theta, eval_dens, epsilon = .0001, sims = gnum))

g_eval_grid <- tibble(theta = grid_theta) %>% 
  mutate(dens_eval_gnum = map(theta, eval_dens, epsilon = .0001, sims = pop_cdf_sims))

head(g_eval_grid) 

xbar = mean(combined_glass)

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

sims <- gnum
x = 1.51748
epsilon <- .0001
eval_dens <- function(x, epsilon, sims, ...) {
  passing_sims <- sims %>%
    unnest() %>%
    mutate(in_interval = map_lgl(locations, # column to pass
                                 in_range, # to this function
                                 fixed_value = x, epsilon = epsilon
    )) %>% # other args passed to the function
    filter(in_interval)
  
  if (nrow(passing_sims) == 0){
    return(0)
  } 
  
  dens_ests <- passing_sims %>%
    group_by(rep) %>%
    summarize(
      min_loc = locations[which.min(locations)],
      min_cdf = weCDF[which.min(locations)],
      max_loc = locations[which.max(locations)],
      max_cdf = weCDF[which.max(locations)]
    ) %>%
    mutate(dens_est = max_cdf - min_cdf) %>%
    pull(dens_est)
  
  
}

