import numpy as np
import scipy
import scipy.stats


def compute_moments(
    log_likelihood_ratio_func, dens_func, max_order, left=-np.inf, right=np.inf
):
    """
    Compute the moments of the log likelihood ratio up to order n
    """

    moments = [0] * max_order
    errs = [0] * max_order
    for ii in range(max_order):
        order = ii + 1

        def integrand(x):
            return log_likelihood_ratio_func(x) ** order * dens_func(x)

        moments[ii], errs[ii] = scipy.integrate.quad(integrand, left, right)
    return moments, errs


def compute_cumulants(moments):
    """
    Compute the cumulants from the moments up to order 4
    """
    assert len(moments) >= 4, "You must have moments at least up to order 4."

    kappas = [0] * 4
    kappas[0] = moments[0]
    kappas[1] = moments[1] - moments[0] ** 2
    kappas[2] = moments[2] - 3 * moments[1] * moments[0] + 2 * moments[0] ** 3
    kappas[3] = (
        moments[3]
        - 4 * moments[2] * moments[0]
        - 3 * moments[1] ** 2
        + 12 * moments[1] * moments[0] ** 2
        - 6 * moments[0] ** 4
    )
    return kappas


def _approx_delta_Fn_edgeworth(x, sigma_square, kappa_3, kappa_4):
    """
    Compute the approximated value of delta_Fn(x) = Fn(x) - Phi(x), where Phi is CDF
    of the standard Normal distribution.

    Input:
        x - The data point where you want to evaluate delta_Fn.
        sigma_square - An array. It consists of variances of n variables.
        kappa_3 - An array. It consists of 3rd order cumulants of n variables.
        kappa_4 - An array. It consists of 4th order cumulants of n variables.
    """
    inv_sigma_n = 1.0 / np.sqrt(np.sum(sigma_square))
    kap_3 = np.sum(kappa_3)
    kap_4 = np.sum(kappa_4)
    return (
        -1.0 / 6.0 * inv_sigma_n ** 3 * kap_3 * (x ** 2 - 1.0)
        - 1.0 / 24.0 * inv_sigma_n ** 4 * kap_4 * (x ** 3 - 3 * x)
        - 1.0 / 72.0 * inv_sigma_n ** 6 * kap_3 ** 2 * (x ** 5 - 10 * x ** 3 + 15 * x)
    ) * scipy.stats.norm.pdf(x)


def _approx_delta_h_cornish_fisher(alpha, sigma_square, kappa_3, kappa_4):
    """
    Compute the approximated value of delta_h = w - z, where w and z are alpha quantiles
    of the distribution Fn(x) and Phi(x) (standard Normal).

    Input:
        alpha - The data point where you want to evaluate delta_h.
        sigma_square - An array. It consists of variances of n variables.
        kappa_3 - An array. It consists of the 3rd order cumulants of n variables.
        kappa_4 - An array. It consists of the 4th order cumulants of n variables.
    """
    inv_sigma_n = 1.0 / np.sqrt(np.sum(sigma_square))
    kap_3 = np.sum(kappa_3)
    kap_4 = np.sum(kappa_4)
    z = scipy.stats.norm.ppf(alpha)

    return (
        1.0 / 6.0 * inv_sigma_n ** 3 * kap_3 * (z ** 2 - 1.0)
        + 1.0 / 24.0 * inv_sigma_n ** 4 * kap_4 * (z ** 3 - 3 * z)
        - 1.0 / 36.0 * inv_sigma_n ** 6 * kap_3 ** 2 * (2 * z ** 3 - 5 * z)
    )


def _approx_quantile_Fn_edgeworth(alpha, sigma_square, kappa_3, kappa_4):
    """
    Compute the alpha quantile of Fn(x). First, Fn(x) is approximated by the Edgeworth
    expansion. Next, we numerically solve the equation Fn(x) - alpha = 0

    Input:
        alpha - The data point where you want to evaluate delta_h.
        sigma_square - An array. It consists of variances of n variables.
        kappa_3 - An array. It consists of the 3rd order cumulants of n variables.
        kappa_4 - An array. It consists of the 4th order cumulants of n variables.
    """

    def func_to_search(x):
        return (
            scipy.stats.norm.cdf(x)
            + _approx_delta_Fn_edgeworth(x, sigma_square, kappa_3, kappa_4)
            - alpha
        )

    sol = scipy.optimize.fsolve(func_to_search, x0=scipy.stats.norm.ppf(alpha))
    return sol


def approx_f_edgeworth(
    alpha,
    sigma_square_P,
    sigma_square_Q,
    kappa_3_P,
    kappa_3_Q,
    kappa_4_P,
    kappa_4_Q,
    mu,
    use_cornish_fisher=False,
):
    """
    Compute the approximated value of f-DP function at alpha, i.e. f(alpha), using
    the Edgeworth expansion.

    Input:
        alpha - The data point where you want to evaluate f.
        sigma_square_P - An array of the variance of n variables under P
        sigma_square_Q - An array of the variance of n variables under Q
        kappa_3_P  - An array of the 3rd order cumulants of n variables under P.
        kappa_3_Q  - An array of the 3rd order cumulants of n variables under Q.
        kappa_4_P  - An array of the 4th order cumulants of n variables under P.
        kappa_4_Q  - An array of the 4th order cumulants of n variables under Q.
        mu - (Eq[Tn] - Ep[Tn]) /  sqrt(Var_P(Tn))
        use_cornish_fisher - A boolean. If true, use the cornish-fisher expansion to
            compute the 1 - alpha quantile of Fn.
    """

    # Compute quantile of approx Fn(x):
    if use_cornish_fisher:
        z = scipy.stats.norm.ppf(1.0 - alpha)
        delta_h = _approx_delta_h_cornish_fisher(
            1.0 - alpha, sigma_square_P, kappa_3_P, kappa_4_P
        )
        x = z + delta_h - mu
    else:
        h = _approx_quantile_Fn_edgeworth(
            1.0 - alpha, sigma_square_P, kappa_3_P, kappa_4_P
        )
        x = h - mu

    # Approx \tilde{F}n(x):
    if sigma_square_P != sigma_square_Q:
        x = x * np.sqrt(np.sum(sigma_square_P) / np.sum(sigma_square_Q))
    normal_cdf_x = scipy.stats.norm.cdf(x)
    delta_Fn = _approx_delta_Fn_edgeworth(x, sigma_square_Q, kappa_3_Q, kappa_4_Q)
    return np.minimum(1.0, np.maximum(0, normal_cdf_x + delta_Fn))


def compute_gdp_clt(alpha, mu_f):
    assert mu_f >= 0, "mu_f should be non-negative."
    z = scipy.stats.norm.ppf(1.0 - alpha)
    return scipy.stats.norm.cdf(z - mu_f)


def compute_gdp_clt_mixture(alpha, mu_gaussian, p, num_composition):
    mu_clt = p * np.sqrt(num_composition) * np.sqrt(np.exp(mu_gaussian ** 2) - 1)
    return compute_gdp_clt(alpha, mu_clt)


def compute_f_conjugate(y, xval, fx_val):
    """
        Compute the convex conjugate function of f:
            f*(y) := sup_x yx -  f(x)
    """
    return np.max(y * xval - fx_val)


def approx_f_numerical(
    dens_func_Q, log_likelihood_ratio_func, num_composition, epsilon, alpha,
    asymm=False, left=-np.inf, right=np.inf,
):
    """
    Compute the approximated value of f-DP function at alpha: f(alpha) via supporting
    lines defined by a collection of f_{eps, delta} functions.

    Inputs:
        asymm - A boolean. If true, it means the function f is asymmetric and
           the function f_eps_delta should be addressed differently as below.
    """

    def __compute_f_eps_delta(eps, delta, alpha, asymm=False):
        if asymm:
            return np.maximum(0, np.maximum(0, 1.0 - delta - np.exp(eps) * alpha))
        else:
            return np.maximum(
                np.maximum(0, 1.0 - delta - np.exp(eps) * alpha),
                np.exp(-eps) * (1.0 - delta - alpha),
            )

    def __compute_delta(epsilon, current_order):
        delta = np.zeros(epsilon.shape)
        if current_order == 0:
            return np.maximum(1.0 - np.exp(epsilon), 0)
        else:
            delta_prev = __compute_delta(epsilon, current_order - 1)
            for ii in range(epsilon.size):
                def integrand(x):
                    # approximate delta(eps - log_likelihood_ratio_func(x), current_order - 1) * dens_func_Q(x)
                    val = epsilon[ii] - log_likelihood_ratio_func(x)
                    if val < epsilon[0]:
                        dd = 1.0
                    elif val > epsilon[-1]:
                        dd = 0.0
                    else:
                        dd = np.interp(val, epsilon, delta_prev)
                    return dd * dens_func_Q(x)

                delta[ii], _ = scipy.integrate.quad(integrand, left, right)
            return delta

    deltas = __compute_delta(epsilon, num_composition)
    f_eps_delta = np.zeros((epsilon.size, alpha.size))
    ii = 0
    for eps, delta in zip(epsilon, deltas):
        f_eps_delta[ii] = __compute_f_eps_delta(eps, delta, alpha, asymm=asymm)
        ii += 1
    fval = np.amax(f_eps_delta, axis=0)
    return fval, f_eps_delta, deltas
