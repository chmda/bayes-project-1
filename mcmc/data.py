import numpy as np
import json
from dataclasses import dataclass
from typing import Any, Dict


@dataclass
class Data:
    """Dataclass representing the dataset

    :param Gender: A vector of length ``N`` where entries with 1's correspond to boys.
    :param LRT: A vector of length ``N`` containing the standardized London Reading Test (LRT) score
    :param M: Number of schools
    :param N: Number of pupils
    :param R: Scale matrix of shape ``3x3`` representing the prior guess on :math:``\Sigma``
    :param school: A vector of length ``N`` whose entries correspond to the school to which the pupils belong to.
    :param School_denom: A matrix of shape ``Nx3`` where columns are ``CE School``, ``RC School``, ``State school``, and if all entries are 0, then ``Other school``
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

    def copy(self) -> "Coefficients":
        return Coefficients(
            self.alpha[:], self.beta[:], self.theta, self.phi, self.gamma[:], self.T[:]
        )
