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

################################################################################################
###################### Just simulate one draw from the posterior ###############################
################################################################################################
# Step 1: simulate from posterior for g(theta)
# population medians: g(theta) where theta = median(RI) in a pane
glass_pop <- read_csv("paper-1-lr-slr/dat/glass_pop_evett.csv")
#head(glass_pop)
glass_pop_all <- rep(glass_pop$RI, glass_pop$Count)
# number of obs in population
n_pop <- length(glass_pop_all)
# empirical cdf function 
ecdf_glass_pop <- ecdf(glass_pop_all)
# draws from the posterior for g(theta)
pop_cdf_sim <- dpprior_sim(M = n_pop , F0 = draw_from_ecdf , sticks = 100, ecdf = ecdf_glass_pop)

# Step 2: simulate from posterior for f(x-theta_0|theta=0)
# draws from the posterior for f(x|theta)
# repeated measurements from a single pane: 
data("bennett.df", package = "dafs")
glass_pane <- as.numeric(bennett.df[1,])
median_pane <- median(glass_pane)
glass_pane <- glass_pane - median_pane # so just have to subtract theta from the data when evaluating
n_pane <- length(glass_pane)
ecdf_glass_pane <- ecdf(glass_pane)
pane_cdf_sim <- dpprior_sim(M = n_pane , F0 = draw_from_ecdf , sticks = 100, ecdf = ecdf_glass_pane)

# Step 3: Evaluate G(theta + epsilon) - G(theta - epsilon) over a grid of thetas 

grid_theta <- seq(min(glass_pop_all), max(glass_pop_all), by = .00001)

g_eval_grid <- tibble(theta = grid_theta) %>% 
  mutate(dens_eval_gnum = map_dbl(theta, eval_dens, epsilon = .0002, sims = pop_cdf_sim %>% mutate(rep =1)))


# Step 4a: Evaluate F(x - theta + epsilon) - F(x - theta - epsilon) over a grid of thetas 

f_g_eval_grid <- g_eval_grid %>% 
  mutate(dens_eval_fnum.cs1 = map_dbl(theta, function(x) {
    eval_dens(crime_scene_sample$RI[1]- x, epsilon = .0002, sims = pane_cdf_sim %>% mutate(rep = 1)) }
    ))

# Step 4b: repeat 4a for all observations in the joint sample (crime scene + suspect data)

each_window_eval <- function(theta, val) {
  eval_dens(val - theta, epsilon = .0002, sims = pane_cdf_sim %>% mutate(rep = 1))
  }

f_g_eval_grid <- f_g_eval_grid %>% 
  mutate(
   dens_eval_fnum.cs2 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[2]),
   dens_eval_fnum.cs3 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[3]),
   dens_eval_fnum.cs4 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[4]),
   dens_eval_fnum.cs5 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[5]),
   dens_eval_fnum.cs6 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[6]),
   dens_eval_fnum.cs7 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[7]),
   dens_eval_fnum.cs8 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[8]),
   dens_eval_fnum.cs9 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[9]),
   dens_eval_fnum.cs10 = map_dbl(theta, each_window_eval, val = crime_scene_sample$RI[10]), 
   dens_eval_fnum.s1 = map_dbl(theta, each_window_eval, val = suspect_sample[1]),
   dens_eval_fnum.s2 = map_dbl(theta, each_window_eval, val = suspect_sample[2]),
   dens_eval_fnum.s3 = map_dbl(theta, each_window_eval, val = suspect_sample[3]),
   dens_eval_fnum.s4 = map_dbl(theta, each_window_eval, val = suspect_sample[4]),
   dens_eval_fnum.s5 = map_dbl(theta, each_window_eval, val = suspect_sample[5])
)



 f_g_eval_grid %>% nest(dens_eval_fnum.cs1:dens_eval_fnum.s5) %>% 
   select(1,2, combined = data) %>% 
   mutate(num_dat_prod = map_dbl(combined, prod)) %>% # Step 5: Multiply results of 4 for all data combined (for numerator)
   unnest() %>% 
   nest(dens_eval_fnum.cs1:dens_eval_fnum.cs10) %>% 
   rename(crime_scene_dat = data) %>% 
   mutate(den_dat_prod_cs = map_dbl(crime_scene_dat, prod)) %>% # Step 7: Multiply results of 4 for crime scene data only
   unnest() %>% 
   nest(dens_eval_fnum.s1:dens_eval_fnum.s5) %>% 
   rename(suspect_dat = data) %>% 
   mutate(den_dat_prod_sus = map_dbl(suspect_dat, prod)) %>% # Step 8: Multiply results of 4 for suspect data only
   select(theta, dens_eval_gnum, num_dat_prod, den_dat_prod_cs, den_dat_prod_sus) -> summary_values_for_integral



summary_values_for_integral %>% 
  mutate(num_to_sum = dens_eval_gnum * num_dat_prod, # Step 6: Multiply step 3 times step 5 for each value of theta
         den_to_sum1 = dens_eval_gnum * den_dat_prod_cs,  # Step 9: Multiply step 3 times step 7 & step 3 times step 8
         den_to_sum2 = dens_eval_gnum * den_dat_prod_sus) %>% 
  select(num_to_sum, den_to_sum1, den_to_sum2) %>% 
  summarise_all(sum) %>%  # Step 10: Sum over theta then multiply 2 values from step 9
  mutate(LR = num_to_sum / (den_to_sum1 * den_to_sum2)) # Step 11: divide sum of step 6 by  step 10. This is the LR. 


summary_values_for_integral %>%
  ggplot(aes(x = theta)) + 
  geom_line(aes(y = dens_eval_gnum)) +
  geom_line(aes(y = num_dat_prod), color = "red") + 
  geom_line(aes(y = den_dat_prod_cs), color = "blue") + 
  geom_line(aes(y = den_dat_prod_sus), color = "green")
