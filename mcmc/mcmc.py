import functools
import numpy as np
import scipy.stats as sp
from scipy.spatial.distance import mahalanobis
from mcmc.data import Coefficients, CoefficientsChain, Data
from mcmc.utils import multivariate_normal_p
from typing import Any, Callable, Dict, Optional, Union


def _compute_mu_ij(chain: CoefficientsChain, data: Data, k: int) -> np.ndarray:
    mu = np.zeros(data.N)  # size N

    for j in range(data.M):
        idx = data.school == (j + 1)  # pupils belong to school j
        mu[idx] = (
            chain.alpha[k, j, 0]
            + chain.alpha[k, j, 1] * data.LRT[idx]
            + chain.alpha[k, j, 2] * data.VR[idx, 0]
            + chain.beta[k, 0] * data.LRT[idx] ** 2
            + chain.beta[k, 1] * data.VR[idx, 1]
            + chain.beta[k, 2] * data.Gender[idx]
            + chain.beta[k, 3] * data.School_gender[idx, 0]
            + chain.beta[k, 4] * data.School_gender[idx, 1]
            + chain.beta[k, 5] * data.School_denom[idx, 0]
            + chain.beta[k, 6] * data.School_denom[idx, 1]
            + chain.beta[k, 7] * data.School_denom[idx, 2]
        )

    return mu


def _compute_tau_ij(chain: CoefficientsChain, data: Data, k: int) -> np.ndarray:
    return np.exp(chain.theta[k] + chain.phi[k] * data.LRT)


def _sample_alpha(
    chain: CoefficientsChain, data: Data, k: int, prop_sd: float = 1.0
) -> None:
    tau_ij = _compute_tau_ij(chain, data, k)

    for j in range(data.M):
        prop = chain.alpha[k, j] + np.random.normal(loc=0, scale=prop_sd, size=3)
        idx = data.school == (j + 1)

        # compute delta = mu - alpha coefficients
        delta = (
            chain.beta[k, 0] * data.LRT[idx] ** 2
            + chain.beta[k, 1] * data.VR[idx, 1]
            + chain.beta[k, 2] * data.Gender[idx]
            + chain.beta[k, 3] * data.School_gender[idx, 0]
            + chain.beta[k, 4] * data.School_gender[idx, 1]
            + chain.beta[k, 5] * data.School_denom[idx, 0]
            + chain.beta[k, 6] * data.School_denom[idx, 1]
            + chain.beta[k, 7] * data.School_denom[idx, 2]
        )

        def log_pdf(alpha: np.ndarray):
            # p1 = (alpha - chain.gamma[k]).T @ (chain.T[k] @ (alpha - chain.gamma[k]))
            p1 = mahalanobis(alpha, chain.gamma[k], chain.T[k]) ** 2
            tmp = alpha[0] + alpha[1] * data.LRT[idx] + alpha[2] * data.VR[idx, 0]
            p2 = np.sum(tmp**2 - 2 * tmp * (data.Y[idx] - delta) * tau_ij[idx])
            return -(p1 + p2) / 2

        top = log_pdf(prop)
        bottom = log_pdf(chain.alpha[k, j])
        acc = np.exp(top - bottom)

        if np.random.uniform() < acc:
            chain.alpha[k, j] = prop
        # else no change as we have already copied the previous coefficients into the current state


def _sample_beta(
    chain: CoefficientsChain,
    data: Data,
    k: int,
    tau: float = 1e-4,
) -> None:
    # use Gibbs sampler since we know the conditional distribution

    mu_ij = None
    tau_ij = None

    def nu(k: int) -> np.ndarray:
        # get nu
        table: Dict[int, np.ndarray] = {
            0: data.LRT**2,
            1: data.VR[:, 1],
            2: data.Gender,
            3: data.School_gender[:, 0],
            4: data.School_gender[:, 1],
            5: data.School_denom[:, 0],
            6: data.School_denom[:, 1],
            7: data.School_denom[:, 2],
        }
        return table[k]

    for l in range(8):
        # compute mu_ij and tau_ij
        mu_ij = _compute_mu_ij(chain, data, k)
        tau_ij = _compute_tau_ij(chain, data, k)

        # draw beta_l
        nu_ij = nu(l)
        delta_ij = mu_ij - chain.beta[k, l] * nu_ij
        zeta_ij = data.Y - delta_ij

        ## compute sigma
        sigma2 = 1.0 / (tau + np.sum(nu_ij**2 * tau_ij))
        ## compute mu
        mu = np.sum(nu_ij * zeta_ij * tau_ij) * sigma2

        chain.beta[k, l] = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=1)


def _sample_theta(
    chain: CoefficientsChain,
    data: Data,
    k: int,
    tau: float = 1e-4,
    prop_sd: float = 1.0,
) -> None:
    prop = chain.theta[k] + np.random.normal(scale=prop_sd)  # make a proposition
    # compute mu_ij
    mu_ij = _compute_mu_ij(chain, data, k)

    def log_pdf(theta):
        return (
            -(
                theta**2 * tau
                - data.N * theta
                + np.sum(
                    (data.Y - mu_ij) ** 2 * np.exp(theta + chain.phi[k] * data.LRT)
                )
            )
            / 2
        )

    top = log_pdf(prop)
    bottom = log_pdf(chain.theta[k])
    acc = np.exp(top - bottom)
    if np.random.uniform() < acc:
        chain.theta[k] = prop


def _sample_phi(
    chain: CoefficientsChain,
    data: Data,
    k: int,
    tau: float = 1e-4,
    prop_sd: float = 1.0,
) -> None:
    prop = chain.phi[k] + np.random.normal(scale=prop_sd)  # make a proposition
    # compute mu_ij
    mu_ij = _compute_mu_ij(chain, data, k)

    def log_pdf(phi):
        return (
            -(
                phi**2 * tau
                + np.sum(
                    -phi * data.LRT
                    + (data.Y - mu_ij) ** 2 * np.exp(chain.theta[k] + phi * data.LRT)
                )
            )
            / 2
        )

    top = log_pdf(prop)
    bottom = log_pdf(chain.phi[k])
    acc = np.exp(top - bottom)
    if np.random.uniform() < acc:
        chain.phi[k] = prop


def _sample_gamma(
    chain: CoefficientsChain, data: Data, k: int, mu0: np.ndarray, T0: np.ndarray
) -> None:
    prec = T0 + data.M * chain.T[k]
    mu = np.linalg.solve(
        prec, T0.dot(mu0) + chain.T[k].dot(chain.alpha[k].sum(axis=0))
    )  # equivalent to (prec^-1) @ (T0 @ mu0 + sum_j T @ alpha_j)
    chain.gamma[k] = multivariate_normal_p(mu, prec)


def _sample_T(chain: CoefficientsChain, data: Data, k: int) -> None:
    prec = data.R[:]
    for j in range(data.M):
        # compute Gram's matrix and add it to prec
        centered = chain.alpha[k, j] - chain.gamma[k]
        prec += np.outer(centered, centered)

    chain.T[k] = sp.wishart.rvs(3, np.linalg.inv(prec))  # TODO : avoid invert


class MCMC:
    """MCMC sampler."""

    def __init__(self, init: Coefficients, n: int = 10**4) -> None:
        """Create a MCMC sampler.

        :param n: Chain length, defaults to :math:`10^4`.
        """
        self._chain: CoefficientsChain = CoefficientsChain.from_init(init, n)
        self._n = n

    @property
    def n_samples(self) -> int:
        """Returns the number of samples.

        :return: Number of samples
        """
        return self._n

    def fit(
        self,
        data: Data,
        prop_sd: float = 1,
        tau: float = 1e-4,
        mu0: np.ndarray = np.zeros((3,)),
        T0: np.ndarray = np.eye(3),
    ) -> None:
        """Estimate the parameters from the data.

        :param data: Data
        :param init: Initial values of the parameters
        """

        for i in range(self._n):
            # copy from previous state
            self._chain.copy_from_previous(i + 1)
            # compute the coefficients
            ## fixed effects
            _sample_beta(self._chain, data, i + 1, tau)
            _sample_theta(self._chain, data, i + 1, tau, prop_sd)
            _sample_phi(self._chain, data, i + 1, tau, prop_sd)
            ## random coefficients
            _sample_alpha(self._chain, data, i + 1, prop_sd)
            ## hyper priors
            _sample_gamma(self._chain, data, i + 1, mu0, T0)
            _sample_T(self._chain, data, i + 1)

    def predict(self, X: Data, burnin: int = 0) -> np.ndarray:
        """Predict the data using estimated parameters.

        :param X: Data
        :param burnin: Number of elements to drop at the beginning of the chain, defaults to 0
        :return: Predicted data.
        """
        raise NotImplementedError()

    def summary(self, burnin: int = 0) -> Dict[str, Coefficients]:
        """Summarise the results.

        :return: Summary
        """
        # Values:
        ## mean
        ## std
        ## MC_Error -> MC std of the mean
        ## q2.5pc -> 0.025 quantile
        ## median
        ## q97.5pc -> 0.975 quantile
        metrics: Dict[str, Callable] = {
            "mean": functools.partial(np.median, axis=0),
            "std": functools.partial(np.std, axis=0),
            "q2.5pc": functools.partial(np.quantile, q=0.025, axis=0),
            "median": functools.partial(np.median, axis=0),
            "q97.5pc": functools.partial(np.quantile, q=0.975, axis=0),
        }
        results: Dict[str, Coefficients] = {}

        for key, func in metrics.items():
            results[key] = Coefficients(
                alpha=func(self._chain.alpha[burnin:]),
                beta=func(self._chain.beta[burnin:]),
                theta=func(self._chain.theta[burnin:]),
                phi=func(self._chain.phi[burnin:]),
                gamma=func(self._chain.gamma[burnin:]),
                T=func(self._chain.T[burnin:]),
            )

        return results

    def diagnose(self, outfile: Optional[str] = None) -> None:
        """Makes a diagnosis of the chain and display it. You can also save the plot to a file.

        :param outfile: Output filename, defaults to None
        """
        raise NotImplementedError()
