import numpy as np
import json
from dataclasses import dataclass
from typing import Any, Dict


@dataclass
class Data:
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
