
target_data <- cpm(otu_tmm, log = TRUE, prior.count = 2) %>%
  t %>%
  as_tibble(rownames = "sample_id") %>%
  full_join(metadata, by = "sample_id") %>%
  pivot_longer(cols = -any_of(colnames(metadata)), 
               names_to = "asv_names", values_to = "value") %>%
  mutate(across(c(exposure, final_disease_state), factor)) %>%
  filter(time %in% c('T3', 'T7') | (time == "T0" & tank == "HOMO")) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state))


emm_aov <- aov_4(value ~time + (final_disease_state + exposure) + (1 + time | fragment_id), 
                 data = filter(emm_dat, asv_names == 'ASV_3'))


emm_aov2 <- mixed(value ~time + (final_disease_state + exposure) + (0 + dummy(time, c('T3', 'T7')) | fragment_id), 
                 data = filter(target_data, asv_names == 'ASV_3'),
                 method = 'KR')

filter(target_data, asv_names == 'ASV_3') %>%
  count(fragment_id) %>%
  filter(n == 1)
  summarise(sum(n))

#good

emm_aov3 <- mixed(value ~time * (exposure + final_disease_state) +
                    (0 + dummy(time, c('T3', 'T7')) | fragment_id), 
                  data = filter(target_data, asv_names == 'ASV_3') %>%
                    mutate(fragment_id = if_else(time == 'T0', 'homogenate', fragment_id),
                           final_disease_state = factor(as.character(final_disease_state)),
                           exposure = factor(as.character(exposure))),
                  method = 'KR')



tst <- lmerTest::lmer(value ~time * (exposure + final_disease_state) +
       (0 + dummy(time, c('T3', 'T7')) | fragment_id), 
     data = filter(target_data, asv_names == 'ASV_3') %>%
       mutate(fragment_id = if_else(time == 'T0', 'homogenate', fragment_id),
              final_disease_state = factor(as.character(final_disease_state)),
              exposure = factor(as.character(exposure)),
              time = factor(time, levels = c('T0', 'T3', 'T7'), ordered = TRUE)))
summary(tst)

target_data %>%
  select(time, exposure, final_disease_state) %>%
  distinct %>%
  mutate(fit = predict(tst, newdata = ., re.form = NA)) %>%
  ggplot(aes(x = time, y = fit, colour = exposure)) +
  geom_point()


emm_aov3 <- mixed(value ~time * exposure + final_disease_state + time*final_disease_state +
                    (0 + dummy(time, c('T3', 'T7')) | fragment_id), 
                  data = filter(target_data, asv_names == 'ASV_3') %>%
                    mutate(fragment_id = if_else(time == 'T0', 'homogenate', fragment_id),
                           final_disease_state = factor(as.character(final_disease_state)),
                           exposure = factor(as.character(exposure))),
                  method = 'KR')

emmeans(emm_aov3, ~time * exposure + time * final_disease_state) %>%
  as_tibble() %>%
  ggplot(aes(x = time, y = emmean, ymin = lower.CL, ymax = upper.CL, colour = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(~final_disease_state)

emmeans(emm_aov3, ~time * final_disease_state) %>%
  as_tibble() %>%
  ggplot(aes(x = time, y = emmean, ymin = lower.CL, ymax = upper.CL, colour = final_disease_state)) +
  geom_pointrange(position = position_dodge(0.5))

emmeans(emm_aov3, ~exposure) %>%
  as_tibble() %>%
  ggplot(aes(x = exposure, y = emmean, ymin = lower.CL, ymax = upper.CL, colour = exposure)) +
  geom_pointrange(position = position_dodge(0.5))

model.matrix(~0 + time * (final_disease_state + exposure), data = filter(target_data, asv_names == 'ASV_3') %>%
               mutate(fragment_id = if_else(time == 'T0', 'homogenate', fragment_id),
                      final_disease_state = factor(as.character(final_disease_state)),
                      exposure = factor(as.character(exposure))))

emm_aov3$full_model %>% ranef

emm_aov3 %>%
  emmeans(~time * exposure) %>%
  as_tibble %>%
  ggplot(aes(x = time, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange()



tst_data <- mutate(target_data, fragment_id = if_else(time == 'T0', 'homogenate', fragment_id)) %>%
  nest(data = -asv_names) %>%
  sample_n(10) %>%
  unnest(data) %>%
  mutate(exposure = factor(as.character(exposure)),
         final_disease_state = factor(as.character(final_disease_state)))

tst_data %>%
  count(final_disease_state, time, exposure)

with(tst_data, table(time, asv_names))

emm_aov3 <- mixed(value ~time * (exposure + final_disease_state) +
                    (0 + dummy(time, c('T3', 'T7')) | fragment_id), 
                  data = filter(tst_data, asv_names == 'ASV_403'),
                  method = 'KR')


emm_aov3 <- mixed(value ~time * (final_disease_state) * asv_names + exposure * asv_names +
                    (0 + dummy(time, c('T3', 'T7')) + asv_names | fragment_id), 
                  data = tst_data,
                  method = 'KR')

emm_aov3



value ~ asv_names + (0 + asv_names | fragment_id)
value ~ asv_names + (1 | fragment_id:asv_names)
value ~ asv_names + (1 | fragment_id) + (1 | asv_names)