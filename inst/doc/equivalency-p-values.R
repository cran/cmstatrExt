## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(cmstatrExt)
library(tidyverse)
library(cmstatr)

## -----------------------------------------------------------------------------
dat <- carbon.fabric.2 %>%
  filter(condition == "RTD" & test == "WT")
dat

## -----------------------------------------------------------------------------
qual <- dat %>%
  summarise(n = n(), mean = mean(strength), sd = sd(strength))
qual

## -----------------------------------------------------------------------------
k <- k_equiv_two_sample(0.05, qual$n, 8)
k

## -----------------------------------------------------------------------------
acceptance_limits <- qual$mean - k * qual$sd
acceptance_limits

## -----------------------------------------------------------------------------
p_equiv_two_sample(
  n = qual$n,
  m = 8,
  t1 = (qual$mean - acceptance_limits[1]) / qual$sd,
  t2 = (qual$mean - acceptance_limits[2]) / qual$sd
)

## -----------------------------------------------------------------------------
p_equiv_two_sample(
  n = qual$n,
  m = 8,
  t1 = (qual$mean - 116) / qual$sd,
  t2 = (qual$mean - 138) / qual$sd
)

## -----------------------------------------------------------------------------
curve <- iso_equiv_two_sample(qual$n, 8, 0.05, 4, 1.5, 10)
curve

## -----------------------------------------------------------------------------
curve %>%
  ggplot(aes(x = t1, y = t2)) +
  geom_path() +
  ggtitle("Acceptance criteria for alpha=0.05")

## -----------------------------------------------------------------------------
curve %>%
  ggplot(aes(x = t1, y = t2)) +
  geom_path() +
  geom_hline(yintercept = k[2], color = "red") +
  geom_vline(xintercept = k[1], color = "red") +
  geom_point(data = data.frame(
    t1 = (qual$mean - 116) / qual$sd,
    t2 = (qual$mean - 138) / qual$sd
  ),
  shape = "*", size = 5) +
  ggtitle("Acceptance criteria for alpha=0.05")

## -----------------------------------------------------------------------------
curve %>%
  mutate(x_min = qual$mean - t1 * qual$sd,
         x_mean = qual$mean - t2 * qual$sd) %>%
  ggplot(aes(x = x_min, y = x_mean)) +
  geom_path() +
  geom_hline(yintercept = acceptance_limits[2], color = "red") +
  geom_vline(xintercept = acceptance_limits[1], color = "red") +
  geom_point(data = data.frame(
    x_min = 116,
    x_mean = 138
  ),
  shape = "*", size = 5) +
  ggtitle("Acceptance criteria for alpha=0.05")

