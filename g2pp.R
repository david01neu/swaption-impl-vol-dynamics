# G2++ model: Functions to compute swaption prices and to compute analytical properties
# Reference:
#
#   Brigo, D and F. Mercurio. Interest Rate Models - Theory and
#   Practice. Springer Finance, 2006.

B <- function(z, t, T_) {
  (1 - exp(-z * (T_ - t))) / z
}

V <- function(t, T_) {
  sigma^2 / a^2 * (T_ - t + 2 / a * exp(-a * (T_ - t)) - 1 / (2 * a) * exp(-2 * a * (T_ - t)) - 3 / (2 * a)) +
    eta^2 / b^2 * (T_ - t + 2 / b * exp(-b * (T_ - t)) - 1 / (2 * b) * exp(-2 * b * (T_ - t)) - 3 / (2 * b)) +
    2 * rho * sigma * eta / (a * b) * (T_ - t + (exp(-a * (T_ - t)) - 1) / a +
      (exp(-b * (T_ - t)) - 1) / b - (exp(-(a + b) * (T_ - t)) - 1) / (a + b))
}

PM <- function(end) {
  1 / (1 + nss(end))^end
}

A <- function(t, T_) {
  retVal <- PM(T_) / PM(t) * exp(0.5 * (V(t, T_) - V(0, T_) + V(0, t)))
}

swaption_price_g2f <- function(
    yield_curve_func,
    a, b, sigma, eta, rho,
    N, strike, T_, tenor, tau = 1, payer_swaption = FALSE) {
  # w == -1 if receiver swaption
  # w ==  1 if payer swaption
  w <- 2 * payer_swaption - 1



  n_swaptions <- length(T_)
  swaptionprice <- rep(0, n_swaptions)

  for (swaption_index in n_swaptions) {
    t_i <- seq(T_ + tau, T_ + tenor, by = tau)
    c_i <- rep(strike[swaption_index] * tau, times = length(t_i))
    c_i[length(c_i)] <- c_i[length(c_i)] + 1

    mu_helper <- function(s, t, T_, a, b, sigma, eta, rho) {
      (sigma^2 / a^2 + rho * sigma * eta / (a * b)) * (1 - exp(-a * (t - s)))
      -sigma^2 / (2 * a^2) * (exp(-a * (T_ - t)) - exp(-a * (T_ + t - 2 * s)))
      -rho * sigma * eta / (b * (a + b)) * (exp(-b * (T_ - t)) - exp(-b * T_ - a * t + (a + b) * s))
    }

    mu_x <- -mu_helper(0, T_, T_, a, b, sigma, eta, rho)
    mu_y <- -mu_helper(0, T_, T_, b, a, eta, sigma, rho)

    sigma_x <- sigma * sqrt((1 - exp(-2 * a * T_)) / (2 * a))
    sigma_y <- eta * sqrt((1 - exp(-2 * b * T_)) / (2 * b))
    rho_xy <- rho * sigma * eta / ((a + b) * sigma_x * sigma_y) * (1 - exp(-(a + b) * T_))
    # print(rho_xy)

    x <- seq(mu_x - 10 * sigma_x, mu_x + 10 * sigma_x, length.out = 1001)

    cA <- c_i * A(T_, t_i)
    ybar <- c()
    for (x_val in x) {
      tmp_func <- function(y) {
        sum(cA * exp(-B(a, T_, t_i) * x_val - B(b, T_, t_i) * y)) - 1
      }

      ybar <- tryCatch(
        {
          c(ybar, uniroot(
            function(ybar) {
              sum(cA * exp(-B(a, T_, t_i) * x_val - B(b, T_, t_i) * ybar)) - 1
            },
            lower = -1,
            upper = 1,
            extendInt = "yes"
          )$root)
        },
        error = function(e) {
          return(999)
        }
      )
      if (is.null(ybar)) {
        return(NULL)
      }
    }

    sqrt_rho <- sqrt(1 - rho_xy^2)
    h_1 <- (ybar - mu_y) / (sigma_y * sqrt_rho) - rho_xy * (x - mu_x) / (sigma_x * sqrt_rho)
    h_2 <- matrix(h_1, nrow = 1) |> apply(2, function(f) {
      f + matrix(B(b, T_, t_i) * sigma_y * sqrt_rho, ncol = 1)
    })

    lambda <- matrix(rep(cA, times = length(x)), ncol = length(x), byrow = FALSE) * exp(-matrix(x, nrow = 1) |> apply(2, function(f) {
      f * matrix(B(a, T_, t_i), ncol = 1)
    }))

    Kappa <- -(
      matrix(rho_xy * sigma_y * (x - mu_x) / sigma_x, nrow = 1) |>
        apply(2, function(f) {
          mu_y + f - matrix(0.5 * sqrt_rho^2 * sigma_y^2 * B(b, T_, t_i), ncol = 1)
        })
    ) * matrix(rep(B(b, T_, t_i), times = length(x)), ncol = length(x), byrow = FALSE)


    integrand <- exp(-0.5 * ((x - mu_x) / sigma_x)^2) / (sigma_x * sqrt(2 * pi)) *
      (pnorm(-w[swaption_index] * h_1) - apply(lambda * exp(Kappa) * pnorm(-w[swaption_index] * h_2), 2, sum))

    integral <- trapz(x, integrand)

    return(w[swaption_index] * N[swaption_index] * integral * PM(T_))
  }
}

forward_rate <- function(start, end) {
  return(((1 + nss(end))^end / (1 + nss(start))^start)^(1 / (end - start)) - 1)
}

zero_coupon_bond <- function(T_) {
  return(1 / (1 + nss(T_))^T_)
}

swap_rate <- function(option_period, swap_tenor, payment_times) {
  numerator <- zero_coupon_bond(option_period) - zero_coupon_bond(option_period + swap_tenor)
  denominator <- sum(0.5 * zero_coupon_bond(payment_times))

  return(numerator / denominator)
}

instantaneous_short_rate <- function(T_) {
  nss(T_) +
    sigma^2 / (2 * a^2) * (1 - exp(-a * T_))^2 +
    eta^2 / (2 * b^2) * (1 - exp(-b * T_))^2 +
    rho * sigma * eta / (a * b) * (1 - exp(-a * T_)) * (1 - exp(-b * T_))
}

int_t_T <- function(t, T_) {
  PM(T_) / PM(t) * exp(0.5 * (V(t, T_) - V(0, T_) + V(0, t)))
}

forward_rate_mean <- function(t, tenor) {
  return (instantaneous_short_rate(t))
}

forward_rate_vol <- function(t, tenor, a, b, sigma, eta, rho) {
  beta_a <- (B(a, t, t + tenor) - B(a, t, t)) / tenor

  beta_b <- (B(b, t, t + tenor) - B(b, t, t)) / tenor
  retVal <- sigma^2 * beta_a^2 + eta^2 * beta_b^2 + 2 * rho * sigma * eta * beta_a * beta_b
  return(sqrt(retVal))
}

inst_forward_rate_vol <- function(tenor, a, b, sigma, eta, rho) {
  retVal <- sigma^2 * exp(-2 * a * tenor) + eta^2 * exp(-2 * b * tenor) + 2 * rho * sigma * eta * exp(-(a + b) * tenor)
  return(sqrt(retVal))
}