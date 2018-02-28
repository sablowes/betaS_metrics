knitr::opts_knit$set(progress = FALSE)
#---------2: discrete analyses with mobr, and calculate beta_S_rare--------
sim_mob_in = make_mob_in(spdat, plot_attr) 

# compute stats
mob_stats = get_mob_stats(sim_mob_in, 'group', index = c('N', 'S', 'S_n', 'PIE', 'S_PIE'))

# join the community data frame with the plot attributes
comm_dat <- bind_cols(spdat, plot_attr)

# need the vector of individuals for each of the alpha-level samples
# we will rarefy the group-level individual rarefaction curve to these values of N
alpha_N <- mob_stats$samples_stats %>% 
  filter(index == 'N') %>%
  group_by(group) %>%
  mutate(alphaN = value,
         replicate = 1:length(value)) %>%
  select(-value, -index, -effort) %>%
  unite(group, c(group, replicate), sep = '.') %>%
  nest(alphaN)

# we also need the alpha level values of S to calculate beta (gamma/alpha)
alpha_S <- mob_stats$samples_stats %>%
  filter(index=='S') %>%
  group_by(group) %>%
  mutate(alphaS = value,
         replicate = 1:length(value)) %>%
  select(-value, -index, -effort) %>%
  unite(group, c(group, replicate), sep = '.') %>%
  nest(alphaS)

# need a vector of the treatment scale abundances to rarefy
gamma_comms <- comm_dat %>%
  gather(species, abundance, species1:species99) %>%
  as_tibble() %>%
  group_by(group, species) %>%
  summarise(N = sum(abundance)) %>%
  filter(N > 0) %>%
  nest(N) %>% 
  slice(rep(1:n(), each = n_quadrats)) %>%
  group_by(group) %>%
  mutate(Replicate = 1:n()) %>%
  ungroup() %>%
  unite(group, c(group, Replicate), sep = '.')

gamma_comms <- inner_join(gamma_comms, alpha_N, by = 'group')
gamma_comms <- inner_join(gamma_comms, alpha_S, by = 'group')

# do rarefaction and calculate beta_S_rare
gamma_S <- gamma_comms %>%
  mutate(gamma_rareS = purrr::map2(data.x, data.y, ~ rarefaction(.x$N, method = 'indiv', effort = 1:.y$alphaN))) %>%
  mutate(beta_rareS = purrr::map2(data, gamma_rareS, ~ max(.y) / .x$alphaS)) 

# get rarefied gamma-scale rarefactions ready for plotting
gamma_expandS <- gamma_S %>%
  unnest(gamma_rareS) %>%
  group_by(group) %>%
  mutate(N = 1:n(), 
         scale = 'gamma') %>%
  ungroup() %>%
  separate(group, into = c('group', 'Replicate'), sep = '\\.') %>%
  separate(group, into = c('Distribution', 'CV'), sep = '_', remove = FALSE)

# want alpha-scale rarefactions too
alpha_comms <- comm_dat %>%
  gather(species, abundance, species1:species99) %>%
  as_tibble() %>%
  group_by(group, Replicate, species) %>%
  summarise(N = sum(abundance)) %>%
  filter(N > 0) %>%
  nest(N)

alpha_rare <- alpha_comms %>%
  mutate(alpha_rareS = map(data, ~ rarefaction(.x$N, method = 'indiv')))

alpha_expandS <- alpha_rare %>%
  unnest(alpha_rareS) %>%
  group_by(group, Replicate) %>%
  mutate(N = 1:n(),
         scale = 'alpha') %>%
  ungroup() %>%
  separate(group, into = c('Distribution', 'CV'), sep = '_', remove = FALSE)

# check gamma scale: good they are all on top of each other
# ggplot() +
#   facet_wrap(~group) +
#   geom_line(data = gamma_expandS,
#             aes(x = N, y = gamma_rareS, colour = Replicate))

# get the new beta_S_rare ready for plotting  
betaS_rare <- gamma_S %>%
  unnest(beta_rareS) %>%
  ungroup() %>%
  mutate(index = 'beta_S_rare',
         value = beta_rareS) %>%
  select(-beta_rareS) %>%
  separate(group, into = c('group', 'Replicate'), sep = '\\.') %>%
  separate(group, into = c('Distribution', 'CV'), sep = '_', remove = FALSE)

# join new beta metric with others in mobr package
mob_stats$samples_stats <- bind_rows(mob_stats$samples_stats, betaS_rare %>% select(group, index, value))

#--------plot discrete metrics --------
levels <- paste(plot_attr$spatial, plot_attr$SAD_CV, sep='_') %>% unique()
gamma_expandS$group <- factor(gamma_expandS$group, levels = levels)
alpha_expandS$group <- factor(alpha_expandS$group, levels = levels)

# plot the beta-diversity metrics
mob_stats$samples_stats$group <- factor(mob_stats$samples_stats$group, levels = levels)

betaValues <- mob_stats$samples_stats %>%
  filter(index=='beta_S' | index=='beta_S_PIE' | index=='beta_S_rare') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_bw() +
  ggtitle('No species removed')

alphaS <- mob_stats$samples_stats %>%
  filter(index=='S') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab('S') +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

alphaS_PIE <- mob_stats$samples_stats %>%
  filter(index=='S_PIE') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(S[PIE])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

betaS <- mob_stats$samples_stats %>%
  filter(index=='beta_S') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(paste(beta[S]))) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

betaS_PIE <- mob_stats$samples_stats %>%
  filter(index=='beta_S_PIE') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(beta[S[PIE]])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

gammaS <- mob_stats$groups_stats %>%
  filter(index=='S') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol = 1) +
  geom_point(aes(x = group, y = value), pch = 8, size = 3, stroke = 2) +
  scale_y_continuous(breaks = c(20, 30, 40, 50, 100), labels = c(20, 30, 40, 50, 100)) +
  theme_classic() +
  ylab('S') +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

gammaS_PIE <- mob_stats$groups_stats %>%
  filter(index=='S_PIE') %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol = 1) +
  geom_point(aes(x = group, y = value), pch = 8, size = 3, stroke = 2) +
  scale_y_continuous(breaks = c(20, 30, 40, 50, 100), labels = c(20, 30, 40, 50, 100)) +
  theme_classic() +
  ylab(expression(S[PIE])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

plot_grid(alphaS, betaS, gammaS, alphaS_PIE, betaS_PIE, gammaS_PIE, nrow = 2, labels = 'AUTO')
ggsave('sims_existing_metric.png', width=300, height = 200, units = 'mm')

#-----visualise rarefactions for aggregated communities (with and without species removed)--------------
# we already have the new beta_diversity metrics, lets look at the rarefaction curves first...
agg_rarefactions <- ggplot() +
  facet_grid(Replicate ~ group, scales = 'free_x') +
  geom_line(data = filter(gamma_expandS, Distribution == 'agg' & as.numeric(Replicate) <= 9),
            aes(x = N, y = gamma_rareS, group = Replicate, colour = scale)) +
  geom_line(data = filter(alpha_expandS, Distribution == 'agg' & as.numeric(Replicate) <= 9), 
            aes(x = N, y = alpha_rareS, group = Replicate, colour = scale)) +
  ylab('Species richness') +
  # ggtitle('No species removed') +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95), legend.background = element_blank())

# all species
agg_betaS <- mob_stats$samples_stats %>%
  filter(index=='beta_S' & (group=='agg_1' | group=='agg_2' | group=='agg_4')) %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(paste(beta[S]))) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

agg_betaS_PIE <- mob_stats$samples_stats %>%
  filter(index=='beta_S_PIE' & (group=='agg_1' | group=='agg_2' | group=='agg_4')) %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(beta[S[PIE]])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

agg_betaS_rare <- mob_stats$samples_stats %>%
  filter(index=='beta_S_rare' & (group=='agg_1' | group=='agg_2' | group=='agg_4')) %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(beta[S[rare]])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))



agg_ALL_betas <- plot_grid(agg_betaS, agg_betaS_PIE, agg_betaS_rare, ncol = 1, labels = c('B', 'C', 'D'))
plot_grid(agg_rarefactions, agg_ALL_betas, nrow = 1, labels = 'A')
ggsave('agg_rarefactions_betaD.png', width = 290, height = 200, units = 'mm')

#-----visualise rarefactions for poisson communities (with and without species removed)--------------
poi_rarefactions <- ggplot() +
  facet_grid(Replicate ~ group, scales = 'free_x') +
  geom_line(data = filter(gamma_expandS, Distribution == 'poi' & as.numeric(Replicate) <= 9),
            aes(x = N, y = gamma_rareS, group = Replicate, colour = scale)) +
  geom_line(data = filter(alpha_expandS, Distribution == 'poi' & as.numeric(Replicate) <= 9), 
            aes(x = N, y = alpha_rareS, group = Replicate, colour = scale)) +
  ylab('Species richness') +
  # ggtitle('No species removed') +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95), legend.background = element_blank())

poi_betaS <- mob_stats$samples_stats %>%
  filter(index=='beta_S' & (group=='poi_1' | group=='poi_2' | group=='poi_4')) %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(paste(beta[S]))) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

poi_betaS_PIE <- mob_stats$samples_stats %>%
  filter(index=='beta_S_PIE' & (group=='poi_1' | group=='poi_2' | group=='poi_4')) %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(beta[S[PIE]])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

poi_betaS_rare <- mob_stats$samples_stats %>%
  filter(index=='beta_S_rare' & (group=='poi_1' | group=='poi_2' | group=='poi_4')) %>%
  ggplot() +
  facet_wrap(~index, scales = 'free', ncol=1) +
  geom_boxplot(aes(x = group, y = value)) +
  theme_classic() +
  ylab(expression(beta[S[rare]])) +
  xlab('') +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14))

poi_ALL_beta <- plot_grid(poi_betaS, poi_betaS_PIE, poi_betaS_rare, ncol = 1, labels = c('B', 'C', 'D'))
plot_grid(poi_rarefactions, poi_ALL_beta, nrow = 1, labels = 'A')
ggsave('poi_rarefactions_betaD.png', width = 290, height = 200, units = 'mm')

