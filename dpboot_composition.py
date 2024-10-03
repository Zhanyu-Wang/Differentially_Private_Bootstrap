# modified from https://github.com/enosair/gdp-edgeworth/
from scipy.stats import binom
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import time
import pickle
import argparse
import edgeworth

from matplotlib import rc

rc("text", usetex=False)
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
                # print("approx_f_numerical_ii", ii, epsilon.size)
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

def run_expr_dual(
    dens_func_P,
    dens_func_Q,
    log_likelihood_ratio_func,
    num_composition,
    n,
    epsilon_bin=10,
    vertical_line_x=None,
    use_cornish_fisher=False,
    left=-np.inf,
    right=np.inf,
    title=None,
    ax=None,
    log_scale=False,
    asymm=False,
):

    mu_f = np.sqrt(2 - 2 / np.exp(1))
    # print("mu_f", mu_f)

    epsilon = np.linspace(-6, 10, epsilon_bin)
    delta_clt = []
    for eps in epsilon:
        # print(eps)
        delta_clt.append(
            scipy.stats.norm.cdf(-eps / mu_f + mu_f / 2)
            - np.exp(eps) * scipy.stats.norm.cdf(-eps / mu_f - mu_f / 2)
        )
        

    alpha = np.sort(
        np.concatenate(
            (
                np.logspace(-8, -1, 100),
                np.linspace(0.09, 0.99, 500),
                1.0 - np.logspace(-8, -1, 200),
            ),
            axis=None,
        )
    )

    start = time.time()
    f_numerical, f_eps_delta, deltas = approx_f_numerical(
        dens_func_Q, log_likelihood_ratio_func, num_composition, epsilon, alpha=alpha, asymm=asymm
    )
    numerical_time = time.time() - start
    print("Numerical used {} seconds.".format(numerical_time))



    if ax is not None:
        line1, = ax.plot(epsilon, deltas, linewidth=4, color="k", linestyle="-")
        line2, = ax.plot(epsilon, delta_clt, linewidth=4, color="b", linestyle="--")
        # line3, = ax.plot(epsilon, delta_clt2, linewidth=4, color="r")
        # line4, = ax.plot(epsilon, delta_edgeworth, linewidth=4, color="g")
        ax.set_xlabel(r"$\epsilon$", fontsize=30)
        ax.set_ylabel(r"$\delta(\epsilon)$", fontsize=30)
        ax.legend(
            # ["numerical", "CLT", "CLT2", "Edgeworth"],
            ["numerical", "asymptotic"],
            prop={"size": 30},
            loc="upper right",
        )
        ax.xaxis.set_tick_params(size=18, labelsize=30)
        ax.yaxis.set_tick_params(size=18, labelsize=30)
        if log_scale:
            ax.set_yscale("log")
        if title is not None:
            ax.set_title(title, fontdict={"fontsize": 30})

def dpboot_dual_all(num_compositions, n, epsilon_bin, save_fig=False, log_scale=False):
    for ii in range(len(num_compositions)):
        fig, axs = plt.subplots(1, 1, figsize=(10, 10))

        num_composition = num_compositions[ii]

        mu = 1.0 / np.sqrt(num_composition)
        
        seq_1_n = np.array(range(1, n + 1))
        p = 1 / n
        # p_0
        binom_0_n_p = binom.pmf(0, n, p)

        # alpha^*
        P_Q_alpha_eq_beta_prob = np.sum(np.array(list(map(lambda k: norm.cdf((2 * 0 - (k*mu) ** 2) / (2 * (k*mu)), 0, 1), seq_1_n))) * binom.pmf(seq_1_n, n, p)) / (1 - binom_0_n_p)

        # -log(-lambda)
        log_pdf_ratios = np.linspace(-10, 0, epsilon_bin)
        alphas = []
        for log_pdf_ratio in log_pdf_ratios:
            alpha = np.sum(np.array(list(map(lambda k: norm.cdf((2 * log_pdf_ratio - (k*mu) ** 2) / (2 * (k*mu)), 0, 1), seq_1_n))) * binom.pmf(seq_1_n, n, p)) / (1 - binom_0_n_p)
            alphas.append(alpha)
        
        subsampled_neg_slope_left = binom.pmf(0, n, p) + (1 - binom.pmf(0, n, p)) * np.exp(-log_pdf_ratios)
        subsampled_alphas = np.array(alphas)

        alphas = []
        betas = []
        log_pdf_ratios = np.linspace(0, 10, epsilon_bin)
        for log_pdf_ratio in log_pdf_ratios:
            alpha = np.sum(np.array(list(map(lambda k: norm.cdf((2 * log_pdf_ratio - (k*mu) ** 2) / (2 * (k*mu)), 0, 1), seq_1_n))) * binom.pmf(seq_1_n, n, p)) / (1 - binom_0_n_p)
            alphas.append(alpha)
            beta = np.sum(np.array(list(map(lambda k: norm.cdf((-2 * log_pdf_ratio - (k*mu) ** 2) / (2 * (k*mu)), 0, 1), seq_1_n))) * binom.pmf(seq_1_n, n, p)) / (1 - binom_0_n_p)
            betas.append(beta)

        subsampled_neg_slope_right = 1 / (binom.pmf(0, n, p) + np.exp(log_pdf_ratios) * (1 - binom.pmf(0, n, p)))
        subsampled_alphas2 = np.array(alphas) * (1 - binom.pmf(0, n, p)) + binom.pmf(0, n, p) * (1 - np.array(betas))

        err = 1e-10
        def dens_func_Q(x):
            if x < err:
                return 0
            elif x < P_Q_alpha_eq_beta_prob:
                return np.interp(x, subsampled_alphas, subsampled_neg_slope_left)
            elif x < binom_0_n_p + (1 - 2 * binom_0_n_p) * P_Q_alpha_eq_beta_prob:
                return 1
            elif x < 1 - err:
                return np.interp(x, subsampled_alphas2, subsampled_neg_slope_right)
            else:
                return 0
        
        def dens_func_P(x): 
            if x < err:
                return 0
            elif x < 1 - err:
                return 1
            else:
                return 0

        def log_likelihood_ratio_func(x):
            """
            log_likelihood_ratio(x) = log(dens_func_Q(x) / dens_func_P(x))
            """
            if x < err:
                return 0
            elif x < 1 - err:
                return np.log(dens_func_Q(x))
            else:
                return 0
            
        run_expr_dual(
            dens_func_P,
            dens_func_Q,
            log_likelihood_ratio_func,
            num_composition,
            n,
            epsilon_bin,
            title="n = " + str(n) + ", B = " + str(num_composition),
            ax=axs,
            log_scale=log_scale,
            vertical_line_x=num_composition * mu,
            asymm=False,
        )
        plt.show()
        if save_fig:
            fig.savefig("./figures/Figure4/eps_bin="+str(epsilon_bin)+"_n="+str(n)+"_B="+str(num_composition)+".pdf", bbox_inches="tight")

parser = argparse.ArgumentParser(description='PyTorch ImageNet Training')
parser.add_argument('-n', default=3, type=int, metavar='C',
                    help='n (default: 3)')
parser.add_argument('-e', default=10, type=int, metavar='C',
                    help='e (default: 10)')
parser.add_argument('-B', default=1, type=int, metavar='C',
                    help='B (default: 1)')
args = parser.parse_args()                    
print(args)
dpboot_dual_all(num_compositions=[args.B], n=args.n, epsilon_bin=args.e, save_fig=True)
