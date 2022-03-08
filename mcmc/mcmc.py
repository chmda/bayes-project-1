from typing import Any, Dict, Optional, Union
import numpy as np
from mcmc.data import Data


def _compute_mu_ij(
    alpha: np.ndarray,
    beta: np.ndarray,
    LRT: np.ndarray,
    VR1: np.ndarray,
    VR2: np.ndarray,
    Girl: np.ndarray,
    GirlsSchool: np.ndarray,
    BoysSchool: np.ndarray,
    CESchool: np.ndarray,
    RCSchool: np.ndarray,
    OtherSchool: np.ndarray,
) -> np.ndarray:
    mu_ij = np.zeros_like(LRT)
    for i in range(mu_ij.shape[0]):
        mu_ij[i] = (
            alpha[0]
            + alpha[1] * LRT[i]
            + alpha[2] * VR1[i]
            + beta[0] * LRT[i] ** 2
            + beta[1] * VR2[i]
            + beta[2] * Girl[i]
            + beta[3] * GirlsSchool
            + beta[4] * BoysSchool
            + beta[5] * CESchool
            + beta[6] * RCSchool
            + beta[7] * OtherSchool
        )

    return mu_ij


def _compute_tau_ij(theta: float, phi: float, LRT: np.ndarray) -> np.ndarray:
    return np.exp(theta + phi * LRT)


def _sample_beta(
    init: np.ndarray,
    Y: np.ndarray,
    x: Data,
    alpha: np.ndarray,
    theta: float,
    phi: float,
    tau: float = 1e-4,
) -> np.ndarray:
    # use Gibbs sampler since we now the conditional distribution

    beta = init[:]  # make a copy of init
    mu_ij = np.zeros_like(Y)
    tau_ij = np.zeros_like(Y)

    def nu(k: int) -> np.ndarray:
        # get nu
        table: Dict[int, np.ndarray] = {
            0: x.LRT,
            1: x["VR2"],
            2: x["Girl"],
            3: x["GirlsSchool"],
            4: x["BoysSchool"],
            5: x["CESchool"],
            6: x["RCSchool"],
            7: x["OtherSchool"],
        }
        return table[k]

    for k in range(8):
        # compute mu_ij and tau_ij
        mu_ij = _compute_mu_ij(
            alpha,
            beta,
            x.LRT,
            x["VR1"],
            x["VR2"],
            x["Girl"],
            x["GirlsSchool"],
            x["BoysSchool"],
            x["CESchool"],
            x["RCSchool"],
            x["OtherSchool"],
        )
        tau_ij = _compute_tau_ij(theta, phi, x.LRT)

        # draw beta_k
        nu_ij = nu(k)
        delta_ij = mu_ij - beta[k] * nu_ij
        zeta_ij = Y - delta_ij

        ## compute sigma
        sigma2 = 1.0 / (tau + np.sum(nu_ij**2 * tau_ij))
        ## compute mu
        mu = np.sum(nu_ij * zeta_ij * tau_ij) * sigma2

        beta[k] = np.random.normal(loc=mu, scale=sigma2, size=1)

    return beta


def _sample_phi(
    phi: float,
    Y: np.ndarray,
    x: Data,
    alpha: np.ndarray,
    beta: np.ndarray,
    theta: float,
    tau: float = 1e-4,
    prop_sd: float = 1.0,
) -> float:
    prop = phi + np.random.normal(scale=prop_sd)  # make a proposition
    # compute mu_ij and tau_ij
    mu_ij = _compute_mu_ij(
        alpha,
        beta,
        x.LRT,
        x["VR1"],
        x["VR2"],
        x["Girl"],
        x["GirlsSchool"],
        x["BoysSchool"],
        x["CESchool"],
        x["RCSchool"],
        x["OtherSchool"],
    )

    def log_pdf(ph: float) -> float:
        # tau_ij is not the same
        tau_ij = _compute_tau_ij(theta, ph, x.LRT)
        return -(ph**2 * tau + np.sum(ph * x.LRT + (Y - mu_ij) ** 2 * tau_ij)) / 2

    top = log_pdf(prop)
    bottom = log_pdf(phi)
    acc = np.exp(top - bottom)
    if np.random.uniform() < acc:
        return prop
    else:
        return phi


def _sample_theta(
    theta: float,
    Y: np.ndarray,
    x: Data,
    alpha: np.ndarray,
    beta: np.ndarray,
    phi: float,
    tau: float = 1e-4,
    prop_sd: float = 1.0,
) -> float:
    prop = theta + np.random.normal(scale=prop_sd)  # make a proposition
    # compute mu_ij and tau_ij
    mu_ij = _compute_mu_ij(
        alpha,
        beta,
        x.LRT,
        x["VR1"],
        x["VR2"],
        x["Girl"],
        x["GirlsSchool"],
        x["BoysSchool"],
        x["CESchool"],
        x["RCSchool"],
        x["OtherSchool"],
    )

    def log_pdf(th: float) -> float:
        # tau_ij is not the same
        tau_ij = _compute_tau_ij(th, phi, x.LRT)
        return -(th**2 * tau + np.sum(th + (Y - mu_ij) ** 2 * tau_ij)) / 2

    top = log_pdf(prop)
    bottom = log_pdf(theta)
    acc = np.exp(top - bottom)
    if np.random.uniform() < acc:
        return prop
    else:
        return theta


class MCMC:
    """MCMC sampler."""

    def __init__(self, n: int = 10**5) -> None:
        """Create a MCMC sampler.

        :param n: Chain size, defaults to :math:`10^5`.
        """
        self._chain: Optional[np.ndarray] = None  # init it later
        self._n = n

    @property
    def n_samples(self) -> int:
        """Returns the number of samples.

        :return: Number of samples
        """
        return self._n

    def fit(
        self, data: Union[Data, Dict[str, Any]], init: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """Estimate the parameters from the data.

        :param data: Data.
        :param init: Initial values of the parameters, defaults to None
        :return: Markov chain.
        """
        raise NotImplementedError()

    def predict(self, X: Dict[str, Any], burnin: int = 0) -> np.ndarray:
        """Predict the data using estimated parameters.

        :param X: Data
        :param burnin: Number of elements to drop at the beginning of the chain, defaults to 0
        :return: Predicted data.
        """
        raise NotImplementedError()

    def summary(self) -> Dict[str, Dict[str, float]]:
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
        raise NotImplementedError()

    def diagnose(self, outfile: Optional[str] = None) -> None:
        """Makes a diagnosis of the chain and display it. You can also save the plot to a file.

        :param outfile: Output filename, defaults to None
        """
        raise NotImplementedError()
