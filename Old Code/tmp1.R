norm_mod <- lmer(value ~ time + (1 | genotype) + (1 | tank), 
                 data = data)
summary(norm_mod)
car::Anova(norm_mod)

emmeans(norm_mod, ~time)


tst <- lmer(value ~ (1 | time) + (1 | exposure) + (1 | final_disease_state) + (1 | genotype) + (1 | tank), 
     data = data)

summary(tst)

library(rstanarm)
library(tidybayes)

tst <- stan_lmer(value ~ (1 | time) + (1 | exposure) + (1 | final_disease_state) + 
                   (1 | genotype) + (1 | tank), 
            data = data, chains = 2, cores = 2)

select(data, time, final_disease_state) %>%
  distinct %>%
  mutate(exposure = 'tmp',
         tank = 'asdas',
         genotype = 'asfasdf') %>%
  add_epred_draws(tst, re_formula = ~(1 | time) + (1 | final_disease_state),
                  allow_new_levels = TRUE) %>%
  point_interval() %>%
  ggplot(aes(x = time, y = .epred, ymin = .lower, ymax = .upper, colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5))
