import numpy as np
import json
from dataclasses import dataclass
from typing import Any, Dict


@dataclass
class Data:
    """Dataclass representing the dataset

    :param Gender: A vector of length ``N`` where entries with 1's correspond to girls.
    :param LRT: A vector of length ``N`` containing the standardized London Reading Test (LRT) score
    :param M: Number of schools
    :param N: Number of pupils
    :param R: Scale matrix of shape ``3x3`` representing the prior guess on :math:``\Sigma``
    :param school: A vector of length ``N`` whose entries correspond to the school to which the pupils belong to.
    :param School_denom: A matrix of shape ``Nx3`` where columns are ``CE School``, ``RC School``, ``Other school``
    :param School_gender: A maitrx of shape ``Nx3`` where columns are ``Girls' school``, ``Boys' school``, and if all entries are 0, then ``Mixed``
    :param Y: A vector of length ``N`` containing the standardized mean examination scores
    :param VR: A matrix of shape ``Nx2`` containing the verbal reasoning test category (1, 2 or 3). Note: if all entries are 0, then VR=3
    """

    Gender: np.ndarray
    LRT: np.ndarray
    M: int
    N: int
    R: np.ndarray
    school: np.ndarray
    School_denom: np.ndarray
    School_gender: np.ndarray
    Y: np.ndarray
    VR: np.ndarray

    @staticmethod
    def from_json(data: str) -> "Data":
        """Create a :py:class:`Data` from Json file.

        :param data: Json
        :return: Instance of :py:class:`Data`
        """

        tmp: Dict[str, Any] = json.loads(data)

        # cast List to ndarray
        for k, v in tmp.items():
            if isinstance(v, list):
                tmp[k] = np.asarray(v)

        return Data(**tmp)


@dataclass
class Coefficients:
    """Dataclass representing the coefficients

    :param alpha: A matrix of shape ``Mx3`` representing the coefficients :math:``\alpha_{\cdot j}`` for each school ``j``
    :param beta: A vector of length ``8``
    :param theta: Coefficient defining :math:``\tau_{ij}``
    :param phi: Coefficient defining :math:``\tau{ij}``
    :param gamma: Hyperparameter corresponding to the mean of :math:``\alpha_{\cdot j}``
    :param T: Hyperparameter corresponding to the precision matrix of :math:``\alpha_{\cdot j}``
    """

    alpha: np.ndarray
    beta: np.ndarray
    theta: float
    phi: float
    gamma: np.ndarray
    T: np.ndarray

    @staticmethod
    def from_json(data: str) -> "Coefficients":
        """Create a :py:class:`Coefficients` from Json file.

        :param data: Json
        :return: Instance of :py:class:`Coefficients`
        """

        tmp: Dict[str, Any] = json.loads(data)

        # cast List to ndarray
        for k, v in tmp.items():
            if isinstance(v, list):
                tmp[k] = np.asarray(v)

        return Coefficients(**tmp)


@dataclass
class CoefficientsChain:
    """Dataclass representing the Markov chain

    :param alpha: A tensor of shape ``KxMx3`` representing the coefficients :math:``\alpha_{\cdot j}`` for each school ``j``
    :param beta: A matrix of length ``Kx8``
    :param theta: Coefficient defining :math:``\tau_{ij}``
    :param phi: Coefficient defining :math:``\tau{ij}``
    :param gamma: Hyperparameter corresponding to the mean of :math:``\alpha_{\cdot j}``
    :param T: Hyperparameter corresponding to the precision matrix of :math:``\alpha_{\cdot j}``
    """

    alpha: np.ndarray
    beta: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    gamma: np.ndarray
    T: np.ndarray

    @staticmethod
    def from_init(init: Coefficients, n: int) -> "CoefficientsChain":
        """Create a :py:class:`Coefficients Chain` from the initial coefficients.

        :param init: Initial coefficients
        :param n: Chain size
        :return: Chain
        """
        chain = CoefficientsChain(
            alpha=np.zeros((n + 1, *init.alpha.shape)),
            beta=np.zeros((n + 1, *init.beta.shape)),
            theta=np.zeros((n + 1,)),
            phi=np.zeros((n + 1,)),
            gamma=np.zeros((n + 1, *init.gamma.shape)),
            T=np.zeros((n + 1, *init.T.shape)),
        )

        chain.alpha[0] = init.alpha
        chain.beta[0] = init.beta
        chain.theta[0] = init.theta
        chain.phi[0] = init.phi
        chain.gamma[0] = init.gamma
        chain.T[0] = init.T

        return chain

    def copy_from_previous(self, k: int) -> None:
        """Copy the coefficients from previous state to the current state

        :param k: Current state index
        """
        if k <= 0:
            raise ValueError("'k' must be greater than 0")

        # copy
        self.alpha[k] = self.alpha[k - 1].copy()
        self.beta[k] = self.beta[k - 1].copy()
        self.theta[k] = self.theta[k - 1]
        self.phi[k] = self.phi[k - 1]
        self.gamma[k] = self.gamma[k - 1].copy()
        self.T[k] = self.T[k - 1].copy()
