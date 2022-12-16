Final project: Does benthic diatom diversity respond to expanding
littoral habitat?
================
Vaclava Hazukova
2022-12-16

## Introduction

Lakes in West Greenland have experienced abrupt environmental shifts as
a result of rapid temperature increase during the 1990s. Among the
non-linear, ecological shifts are prolonged open-water period, higher
lakewater temperatures, and increasing water clarity due to decline in
dissolved organic carbon concentrations. Higher transparency of water
also leads to expansion of the littoral area, illuminated shoreline
habitat that benthic diatoms use.

We collected two sets of surface sediments from 13 lakes across the
Kangerlussuaq region in West Greenland; one set was collected in 1996
and characterizes period prior to rapid climate shift, the second set
was collected in 2013 and represents period after the rapid step-shift
occurred. We analyzed diatoms in each set, distinguishing between
planktic and benthic species. Considering that surface sediments
represent accumulation over approximately 10 years prior collection in
these lakes, the diatom assemblages are inter- and intra-annually
integrated.

Paired with the surface sediments, we also have limnological data;
isolated measurements from the late 90s, and continuous monitoring
starting in 2013. Using dissolved organic carbon concentration data, we
are able to calculate shifts in light penetration in each lake.
Conservatively, we use 10% PAR depth to represent benthic diatom
habitat. Using this information, we are able to calculate how lake area
available to benthic diatoms changed between the two periods.

The expectation is that as habitat expands, benthic diatom species
richness will increase. To test the idea, we will use accumulated
species curves and compare the number of species predicted by the curve
based on diatom data from 1996 with data from 2013.

``` r
library(tidyverse); library(lubridate)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.2

    ## Warning: package 'tibble' was built under R version 4.1.2

    ## Warning: package 'tidyr' was built under R version 4.1.2

    ## Warning: package 'readr' was built under R version 4.1.2

    ## Warning: package 'purrr' was built under R version 4.1.2

    ## Warning: package 'dplyr' was built under R version 4.1.2

    ## Warning: package 'stringr' was built under R version 4.1.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
# surface sdiment diatom data

surfsed <- 
  read_csv('data/processed/surface_sediments_combined.csv')
```

    ## Rows: 11046 Columns: 4

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): taxon, lake
    ## dbl (2): year, freq
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# environmental data

# a) water chemistry (layer-level data)

env <- 
  read_csv('data/raw/Greenland_2013-2019_env_MASTER.csv') %>% 
  mutate(year = year(Date)) %>% 
  filter(year == 2013 & Layer == 'epi') %>% 
  select(-ice, -iceout) %>% 
  rename(lake = Lake) %>% 
  group_by(lake) %>% 
  summarise(across(c(chl, DOC, TP, NO3, PAR), mean, na.rm = T))
```

    ## Rows: 638 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr   (4): Lake, Layer, epi, ice
    ## dbl  (16): Depth, chl, DOC, TP, NO3, NH4, PAR, secchi, iceout, a254, a380, S...
    ## lgl   (1): SUVA254
    ## dttm  (1): Date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# b) DOC decline

doc_change <- 
  read_csv('data/processed/doc_calculated.csv') %>% 
  select(c(1:4, 26:28))
```

    ## Rows: 13 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr   (1): lake
    ## dbl  (24): DOC_96, DOC_13, DOC_decline_percent, a380_96, a380_star_96, a380_...
    ## lgl   (2): PAR10, k
    ## date  (1): date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# combine bathymetries and doc change for the 13 lakes where data are available

bathy <- 
  read_csv("data/raw/bathy.csv") %>% 
  inner_join(doc_change %>% select(lake, par10_96, par10_13),
            by = 'lake')
```

    ## Rows: 506 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): lake
    ## dbl (4): lower, upper, volume, area
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# diatom traits

trait <- read_csv('data/raw/diatom_traits_analysis.csv')
```

    ## Rows: 360 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (2): taxon, shape
    ## dbl (26): pl, pl2, L, W, T, biovol, size_class, mobile, teratology, pioneer,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# calculate volume of lake with > 10% PAR

PAR_area <- map_dfr(
  unique(bathy$lake),
  function(x) {
    temp <- 
      bathy %>% 
      filter(lake == x)
    
    a1 = temp %>% filter(upper == 1 ) %>% pull(area)
    
    #1996
    depth_top = temp %>% distinct(trunc(par10_96, 0)) %>% pull()
    depth_bottom = depth_top + 1 
    a_par_round = temp %>% filter(upper == depth_top) %>% pull(area)
    extra = temp %>% summarise(par10_96 - depth_top) %>% distinct() %>% pull()
    diff = temp %>% summarise(area = area[upper == depth_top] - area[upper == depth_bottom]) %>% pull() 
    a_total_at_par = a_par_round - diff*extra
    a_par_96 = a1-a_total_at_par
    
    #2013
    depth_top_13 = temp %>% distinct(trunc(par10_13, 0)) %>% pull()
    depth_bottom_13 = depth_top_13 + 1 
    a_par_round_13 = temp %>% filter(upper == depth_top_13) %>% pull(area)
    extra_13 = temp %>% summarise(par10_13 - depth_top_13) %>% distinct() %>% pull()
    diff_13 = temp %>% summarise(area = area[upper == depth_top_13] - area[upper == depth_bottom_13]) %>% pull() 
    a_total_at_par_13 = a_par_round_13 - diff_13*extra_13
    a_par_13 = a1-a_total_at_par_13
    
    tibble(
      lake = unique(temp$lake),
      PAR_A_1996 = a_par_96,
      PAR_A_2013 = a_par_13
      )
    
  }) %>% 
  mutate(diff_A = PAR_A_2013 - PAR_A_1996,
         diff_A_percent = ((PAR_A_2013 - PAR_A_1996)/PAR_A_2013*100)) %>% 
  left_join(doc_change %>% select(lake, par10_96, par10_13),
            by = 'lake') %>% 
  mutate(diff_par_depth = ((par10_13-par10_96)/par10_13*100)) %>% 
  rename(par10_1996 = par10_96, par10_2013 = par10_13)
  
# join with lake characteristics 

id <- 
  read_csv("data/raw/lake_id.csv") %>% 
  mutate(GR = A^0.25/maxZ,
         CEI = A/avgZ) %>% 
  rename(lake = Lake) %>% 
  left_join(PAR_area,
            by = 'lake')
```

    ## Rows: 38 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (3): Lake, clim, hydro
    ## dbl (5): lat, long, maxZ, avgZ, A
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# assign diatom status based on presence or ansence in the early 96 dataset

status_diatom <- 
  surfsed %>%
  pivot_wider(names_from = year,
              values_from = freq,
              names_prefix = 'year_') %>% 
  filter(!year_2013 == 0 | !year_1996 == 0) %>% 
  mutate(status = case_when(
    year_2013 > 0 & year_1996 == 0 ~ 'immigration',
    year_2013 == 0 & year_1996 > 0 ~ 'extinction',
    year_2013 > 0 & year_1996 > 0 ~ 'stay'))

# total number of species found in each year and for each lake 

numbers_of_species_tot <-
  surfsed %>% 
  filter(!freq == 0) %>% 
  group_by(lake, year) %>%
  summarise(tot_n = n())
```

    ## `summarise()` has grouped output by 'lake'. You can override using the
    ## `.groups` argument.

``` r
# grouped als oby status 

numbers_of_species <-
  surfsed %>% 
  filter(!freq == 0) %>% 
  left_join(status_diatom %>% select(-year_2013, -year_1996),
            by = c('taxon', 'lake')) %>% 
  group_by(lake, year, status) %>%
  summarise(count = n())
```

    ## `summarise()` has grouped output by 'lake', 'year'. You can override using the
    ## `.groups` argument.

``` r
# show count
surfsed %>% 
  filter(!freq == 0) %>% 
  left_join(status_diatom %>% select(-year_2013, -year_1996),
            by = c('taxon', 'lake')) %>% 
  group_by(year, status) %>%
  distinct(taxon) %>% 
  summarise(count = n())
```

    ## `summarise()` has grouped output by 'year'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 4 × 3
    ## # Groups:   year [2]
    ##    year status      count
    ##   <dbl> <chr>       <int>
    ## 1  1996 extinction    184
    ## 2  1996 stay           54
    ## 3  2013 immigration   138
    ## 4  2013 stay           54

``` r
# create the final joint dataset; assign species tratis, filter benthics with more than 0 occurrences

benthos <- 
  surfsed %>%
  left_join(trait %>% select(taxon, pl, pl2, biovol, 10:28),
            by = 'taxon') %>% 
  filter(pl2 == 0) %>% 
  group_by(lake, year) %>% 
  filter(freq > 0) %>% 
  summarise(benthic_count = n()) %>% 
  pivot_wider(names_from = year,
              values_from = benthic_count,
              names_prefix = 'year_') %>% 
  mutate(diff_benthic_diatoms_percent = (year_2013 - year_1996)/year_2013 * 100,
         diff_benthic_diatoms = year_2013 - year_1996) %>% 
  inner_join(PAR_area,
             by = 'lake') %>% 
  inner_join(numbers_of_species %>% 
               filter(status %in% c('immigration')),
             by = 'lake') %>% 
  inner_join(numbers_of_species_tot %>% 
               filter(year == 2013) %>% select(-year),
             by = 'lake')
```

    ## `summarise()` has grouped output by 'lake'. You can override using the
    ## `.groups` argument.

``` r
# SAD 
surfsed %>% 
  group_by(year) %>% 
  ggplot(aes(x = freq)) +
  geom_histogram() +
  facet_wrap(~year)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](final_project_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
benthos %>% 
  ggplot(aes(x = diff_A, y = count)) +
  geom_point() + 
  geom_smooth(method = lm) +
  geom_abline(slope = 1, intercept = 0) +
  ggrepel::geom_text_repel(aes(label = lake)) +
  labs(y = 'Number of new species (2013 - 1996)', x = 'change in 10% PAR surface area (m2)') +
  theme_bw()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](final_project_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
surfsed %>%
  left_join(trait %>% select(taxon, pl, pl2, biovol, 10:28),
            by = 'taxon') %>% 
  filter(freq > 0 &
           !is.na(pl)) %>% 
  group_by(lake, year) %>% 
  summarise(pb = sum(freq[pl == "1"]) / sum(freq[pl == "0"])) %>% 
  inner_join(PAR_area %>% 
               select(c(1:3)) %>%
               pivot_longer(cols = 2:3,
                            names_to = 'year',
                            names_prefix = 'PAR_A_',
                            names_transform = as.double,
                            values_to = 'area'),
             by = c('lake','year')) %>% 
  ggplot(aes(x = area, y = pb, group = lake)) +
  geom_point(aes(color = factor(year)), size = 3) +
  scale_color_manual('Year', values = c('skyblue3', 'orangered')) + 
  geom_line(aes(group = lake)) +
  labs(y = 'P:B ratio', x = '10% PAR surface area (m2)') +
  theme_bw()
```

    ## `summarise()` has grouped output by 'lake'. You can override using the
    ## `.groups` argument.

![](final_project_files/figure-gfm/unnamed-chunk-3-2.png)<!-- --> The
first plot shows the number of ‘new’ species (immigration) that appeared
in each lake in 2013 relative to 1996 against the change in habitat.
Very, very scattered positive relationship.

In the second plot, we look at shifts in plankton vs. benthic species
ratio (the expectation being that as habitat for benthic increase, we
will see a decrease in this metric), and it’s all over the place. Some
lakes seem to have a very similar direction and magnitude of shifts;
however, in others, planktic ratio increases or does not change.

``` r
# rarefy ------------------------------------------------------------------

surf_par <- 
  surfsed %>%
  left_join(trait %>% select(taxon, pl, pl2, biovol),
            by = 'taxon') %>% 
  filter(pl2 == 0) %>% 
  filter(freq > 0) %>% 
  inner_join(PAR_area %>% 
               select(c(1:3)) %>%
               pivot_longer(cols = 2:3,
                            names_to = 'year',
                            names_prefix = 'PAR_A_',
                            names_transform = as.double,
                            values_to = 'area'),
             by = c('lake','year'))

surf_par_96 <- 
  surf_par %>% 
  filter(year == 1996)

cum_sp <- 
  parallel::mclapply(1:100, mc.cores = 8, FUN = function(j) {
    
    names_lake <-  sample(unique(surf_par_96$lake))
    
    x <- 
      lapply(1:length(names_lake), function(i) {
        temp_name <- names_lake[1:i]
        A = surf_par_96 %>% filter(lake %in% temp_name) %>% distinct(lake,area) %>% summarise(area = sum(area)) %>% pull()
        S = surf_par_96 %>% filter(lake %in% temp_name & freq > 0) %>% distinct(taxon) %>% summarise(n = n()) %>% pull()
        run = j
        c(A,S,run)
      }) 
    do.call(rbind, x) })

cum_sp <- do.call(rbind, cum_sp)

surf_par_13 <- 
  surf_par %>% 
  filter(year == 2013)

cum_sp_13 <- 
  parallel::mclapply(1:100, mc.cores = 8, FUN = function(j) {
    
    names_lake <-  sample(unique(surf_par_13$lake))
    
    x <- 
      lapply(1:length(names_lake), function(i) {
        temp_name <- names_lake[1:i]
        A = surf_par_13 %>% filter(lake %in% temp_name) %>% distinct(lake,area) %>% summarise(area = sum(area)) %>% pull()
        S = surf_par_13 %>% filter(lake %in% temp_name & freq > 0) %>% distinct(taxon) %>% summarise(n = n()) %>% pull()
        run = j
        c(A,S,run)
      }) 
    do.call(rbind, x) })

cum_sp_13 <- do.call(rbind, cum_sp_13)

cum_sp <- cum_sp %>% 
  as_tibble() %>% 
  mutate(year = '1996')
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
    ## `.name_repair` is omitted as of tibble 2.0.0.
    ## ℹ Using compatibility `.name_repair`.

``` r
comb <- cum_sp_13 %>% 
  as_tibble() %>% 
  mutate(year = '2013') %>% 
  bind_rows(cum_sp) %>% 
  rename(area = V1, S = V2, run = V3)
   

comb %>% 
  mutate(code = str_c(year, run, sep = '')) %>% 
  ggplot(aes(x = log10(area), y = log10(S), group = code)) +
  geom_point(aes(color = year), size = 0.1)+
  geom_line(aes(group = code, color = year), linewidth = 0.1) +
  geom_smooth(aes(group = year, color = year), method = 'lm') +
  scale_color_manual('Year', values = c('skyblue1', 'orangered')) + 
  labs(y = 'log(S)', x = 'log(Cumulative area (m2))') +
  theme_bw()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](final_project_files/figure-gfm/rarefy-1.png)<!-- -->

``` r
ggsave('output/plots/rarefaction_log_log.pdf', width = 15, height = 15, units = 'cm')
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
library(meteR)

# 1996

diatom_mete <- 
  surfsed %>% 
  inner_join(id %>% 
               filter(!is.na(par10_1996)) %>% 
               select(lake, lat, long, A, PAR_A_1996, PAR_A_2013, diff_A),
             by = 'lake') %>% 
  left_join(trait %>% 
              select(taxon, pl, biovol),
            by = 'taxon') %>% 
  filter(!is.na(biovol) &
           year == 1996)


esf1 <- meteESF(spp = diatom_mete$taxon,
                abund = diatom_mete$freq,
                power = diatom_mete$biovol^(3/4),
                minE = min(diatom_mete$biovol^(3/4)))
esf1  
```

    ## METE object with state variables:
    ##      S0      N0      E0 
    ##    71.0  1112.0 42080.8 
    ## 
    ## with Lagrange multipliers:
    ##         la1         la2 
    ## 0.013348149 0.001733026

``` r
sad1 <- sad(esf1)
sad1
```

    ## Species abundance distribution predicted using raw data 
    ## with parameters: 
    ##      S0      N0      E0 
    ##    71.0  1112.0 42080.8 
    ##     la1     la2 
    ## 0.01335 0.00173

``` r
plot(sad1, ptype = 'rad', log = 'y')
```

    ## Warning in sprintf(xlab, switch(x$type, sad = "Abundance", ipd = "Metabolic
    ## rate", : one argument not used by format 'Rank'

![](final_project_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot(sad1, ptype = 'cdf', log = 'x')
```

    ## Warning in sprintf(ylab, switch(x$type, sad = "Abundance", ipd = "Metabolic
    ## rate", : one argument not used by format 'Cumulative probability'

![](final_project_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
# 2013

diatom_mete <- 
  surfsed %>% 
  inner_join(id %>% 
               filter(!is.na(par10_1996)) %>% 
               select(lake, lat, long, A, PAR_A_1996, PAR_A_2013, diff_A),
             by = 'lake') %>% 
  left_join(trait %>% 
              select(taxon, pl, biovol),
            by = 'taxon') %>% 
  filter(!is.na(biovol) &
           year == 2013)



esf1 <- meteESF(spp = diatom_mete$taxon,
                abund = diatom_mete$freq,
                power = diatom_mete$biovol^(3/4),
                minE = min(diatom_mete$biovol^(3/4)))
esf1  
```

    ## METE object with state variables:
    ##       S0       N0       E0 
    ##    80.00  1125.00 34267.77 
    ## 
    ## with Lagrange multipliers:
    ##         la1         la2 
    ## 0.014938128 0.002413799

``` r
sad1 <- sad(esf1)
sad1
```

    ## Species abundance distribution predicted using raw data 
    ## with parameters: 
    ##       S0       N0       E0 
    ##    80.00  1125.00 34267.78 
    ##     la1     la2 
    ## 0.01494 0.00241

``` r
plot(sad1, ptype = 'rad', log = 'y')
```

    ## Warning in sprintf(xlab, switch(x$type, sad = "Abundance", ipd = "Metabolic
    ## rate", : one argument not used by format 'Rank'

![](final_project_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
plot(sad1, ptype = 'cdf', log = 'x')
```

    ## Warning in sprintf(ylab, switch(x$type, sad = "Abundance", ipd = "Metabolic
    ## rate", : one argument not used by format 'Cumulative probability'

![](final_project_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->
