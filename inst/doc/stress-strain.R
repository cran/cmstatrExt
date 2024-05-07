## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 3
)

## ----setup--------------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(cmstatrExt)
library(tidyverse)

## -----------------------------------------------------------------------------
head(pa12_tension)

## -----------------------------------------------------------------------------
pa12_tension %>%
  ggplot(aes(x = Strain, y = Stress, color = Coupon)) +
  geom_point()

## -----------------------------------------------------------------------------
curve_quadratic <- average_curve_lm(
  pa12_tension, Coupon,
  Stress ~ I(Strain) + I(Strain^2) + 0
)
summary(curve_quadratic)

## -----------------------------------------------------------------------------
curve_quadratic %>%
  augment() %>%
  ggplot(aes(x = Strain)) +
  geom_point(aes(y = Stress, color = Coupon)) +
  geom_line(aes(y = .fit))

## -----------------------------------------------------------------------------
curve_cubic <- average_curve_lm(
  pa12_tension, Coupon,
  Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0
)
summary(curve_cubic)

## -----------------------------------------------------------------------------
curve_cubic %>%
  augment() %>%
  ggplot(aes(x = Strain)) +
  geom_point(aes(y = Stress, color = Coupon)) +
  geom_line(aes(y = .fit))

## -----------------------------------------------------------------------------
average_curve_lm(
  pa12_tension, Coupon,
  Strain ~ I(Stress) + I(Stress^2) + I(Stress^3) + I(Stress^4) + 0
) %>%
  augment() %>%
  ggplot(aes(y = Stress)) +
  geom_point(aes(x = Strain, color = Coupon)) +
  geom_line(aes(x = .fit))

## -----------------------------------------------------------------------------
bilinear_fn <- function(strain, par) {
  c1 <- par[1]
  c2 <- par[2]
  e1 <- par[3]
  if (strain <= e1) {
    return(c1 * strain)
  } else {
    return(c2 * (strain - e1) + c1 * e1)
  }
}

## -----------------------------------------------------------------------------
curve_bilinear <- average_curve_optim(
  pa12_tension,
  Coupon, Strain, Stress,
  bilinear_fn,
  c(1, 1, 0.04) # the initial guess
)
curve_bilinear

## -----------------------------------------------------------------------------
curve_bilinear <- average_curve_optim(
  pa12_tension,
  Coupon, Strain, Stress,
  bilinear_fn,
  c(1, 1, 0.04),
  lower = c(0, 0, 0.025),
  upper = c(2000, 2000, 0.100)
)
curve_bilinear

## -----------------------------------------------------------------------------
curve_bilinear %>%
  augment() %>%
  ggplot(aes(x = Strain)) +
  geom_point(aes(y = Stress, color = Coupon)) +
  geom_line(aes(y = .fit))

## -----------------------------------------------------------------------------
fff_shear %>%
  ggplot(aes(x = Strain, y = Stress, color = Specimen)) +
  geom_point()

## -----------------------------------------------------------------------------
fff_shear %>%
  filter(Stress > 1000 & Stress < 3000) %>%
  group_by(Specimen) %>%
  nest() %>%
  mutate(lm = map(data, ~lm(Strain ~ Stress, data = .))) %>%
  mutate(x_intercept = map(lm, ~predict(.x, data.frame(Stress = 0)))) %>%
  select(-c(lm, data)) %>%
  unnest(x_intercept)

## -----------------------------------------------------------------------------
fff_shear %>%
  filter(Stress > 1000 & Stress < 3000) %>%
  group_by(Specimen) %>%
  nest() %>%
  mutate(lm = map(data, ~lm(Strain ~ Stress, data = .))) %>%
  mutate(x_intercept = map(lm, ~predict(.x, data.frame(Stress = 0)))) %>%
  select(-c(lm, data)) %>%
  unnest(x_intercept) %>%
  inner_join(fff_shear, by = "Specimen") %>%
  head(6)

## -----------------------------------------------------------------------------
fff_shear_offset <- fff_shear %>%
  filter(Stress > 1000 & Stress < 3000) %>%
  group_by(Specimen) %>%
  nest() %>%
  mutate(lm = map(data, ~lm(Strain ~ Stress, data = .))) %>%
  mutate(x_intercept = map(lm, ~predict(.x, data.frame(Stress = 0)))) %>%
  select(-c(lm, data)) %>%
  unnest(x_intercept) %>%
  inner_join(fff_shear, by = "Specimen") %>%
  mutate(Strain = Strain - x_intercept) %>%
  select(-c(x_intercept))

## -----------------------------------------------------------------------------
fff_shear_offset %>%
  ggplot(aes(x = Strain, y = Stress, color = Specimen)) +
  geom_point()

## -----------------------------------------------------------------------------
fff_shear_offset %>%
  group_by(Specimen) %>%
  mutate(Lead_Stress = lead(Stress, 5),
         Lead_Strain = lead(Strain, 5),
         Slope = (Lead_Stress - Stress) / (Lead_Strain - Strain),
         Remove = Slope < -1e5 | is.na(Slope)) %>%
  ggplot(aes(x = Strain, y = Stress, shape = Specimen, color = Remove)) +
  geom_point()

## -----------------------------------------------------------------------------
fff_shear_offset %>%
  group_by(Specimen) %>%
  mutate(Lead_Stress = lead(Stress, 5),
         Lead_Strain = lead(Strain, 5),
         Slope = (Lead_Stress - Stress) / (Lead_Strain - Strain),
         Remove = Slope < -1e5 | is.na(Slope),
         Remove = cumsum(Remove) > 0) %>%
  ggplot(aes(x = Strain, y = Stress, shape = Specimen, color = Remove)) +
  geom_point()

## -----------------------------------------------------------------------------
fff_shear_truncated <- fff_shear_offset %>%
  group_by(Specimen) %>%
  mutate(Lead_Stress = lead(Stress, 5),
         Lead_Strain = lead(Strain, 5),
         Slope = (Lead_Stress - Stress) / (Lead_Strain - Strain),
         Remove = Slope < -1e5 | is.na(Slope),
         Remove = cumsum(Remove) > 0) %>%
  ungroup() %>%
  filter(!Remove) %>%
  select(Specimen, Stress, Strain)

## -----------------------------------------------------------------------------
fff_shear_truncated_no_toe <- fff_shear_truncated %>%
  filter(Stress > 1000)

## -----------------------------------------------------------------------------
fff_shear_truncated_no_toe %>%
  ggplot(aes(x = Strain, y = Stress, color = Specimen)) +
  geom_point() +
  xlim(c(0, NA)) +
  ylim(c(0, NA))

## -----------------------------------------------------------------------------
curve_fff_shear <- fff_shear_truncated_no_toe %>%
  average_curve_lm(
    Specimen,
    Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0
  )
curve_fff_shear

## -----------------------------------------------------------------------------
curve_fff_shear %>%
  augment(fff_shear) %>%
  ggplot(aes(x = Strain)) +
  geom_point(aes(y = Stress, color = Specimen)) +
  geom_line(aes(y = .fit))

## -----------------------------------------------------------------------------
pa12_tension_conditions <-
  bind_rows(
    pa12_tension %>%
      mutate(Condition = "RTA"),
    pa12_tension %>%
      mutate(Condition = "Fake ETA",
             Stress = 0.50 * Stress,
             Strain = 1.25 * Strain)
  )

## -----------------------------------------------------------------------------
curve_cubic_rta <- pa12_tension_conditions %>%
  filter(Condition == "RTA") %>%
  average_curve_lm(
    Coupon,
    Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0
  )
curve_cubic_rta

## -----------------------------------------------------------------------------
curve_cubic_fake_eta <- pa12_tension_conditions %>%
  filter(Condition == "Fake ETA") %>%
  average_curve_lm(
    Coupon,
    Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0
  )
curve_cubic_fake_eta

## -----------------------------------------------------------------------------
bind_rows(
  augment(curve_cubic_rta),
  augment(curve_cubic_fake_eta)
) %>%
  ggplot(aes(x = Strain, y = .fit, color = Condition)) +
  geom_line()

## -----------------------------------------------------------------------------
bind_rows(
  augment(curve_cubic_rta),
  augment(curve_cubic_fake_eta)
) %>%
  group_by(Condition) %>%
  ggplot(aes(x = Strain)) +
  geom_point(aes(y = Stress, color = Condition)) +
  geom_line(aes(y = .fit, group = Condition))

## -----------------------------------------------------------------------------
bind_rows(
  augment(curve_cubic_rta),
  augment(curve_cubic_fake_eta)
) %>%
  ggplot(aes(x = Strain, y = .fit, color = Condition)) +
  geom_line() +
  scale_y_continuous(
    "Stress [MPa]",
    sec.axis = sec_axis(~ . * 0.1450377377, name = "Stress [ksi]")
  )

## -----------------------------------------------------------------------------
bind_rows(
  augment(curve_cubic_rta),
  augment(curve_cubic_fake_eta)
) %>%
  ggplot(aes(x = Strain, y = .fit, color = Condition)) +
  geom_line() +
  scale_y_continuous(
    "Stress [MPa]",
    sec.axis = sec_axis(~ . * 0.1450377377, name = "Stress [ksi]")
  ) +
  theme_bw()

