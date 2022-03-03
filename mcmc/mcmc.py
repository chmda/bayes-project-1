from typing import Any, Dict, Optional, Union
import numpy as np
from mcmc.data import Data


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
