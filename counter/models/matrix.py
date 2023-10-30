from typing import Optional, Dict

from Bio.PDB.Polypeptide import standard_aa_names
import torch


class Matrix:
    def __init__(
        self,
        dictionary: Optional[Dict[str, Dict[str, float]]] = None,
        tensor: Optional[torch.Tensor] = None,
        device: Optional[torch.device] = None,
    ):
        if tensor is not None:
            self._tensor = tensor
        else:
            self._tensor = torch.zeros((20, 20), device=device, dtype=torch.float)

        if dictionary is not None:
            for i, aai in enumerate(standard_aa_names):
                for j, aaj in enumerate(standard_aa_names):
                    self._tensor[i, j] = dictionary[aai][aaj]

    def __getitem__(self, aai: str, aaj: str) -> float:
        i = standard_aa_names.index(aai)
        j = standard_aa_names.index(aaj)

        return self._tensor[i, j].item()

    def __setitem__(self, aai: str, aaj: str, value: float):
        i = standard_aa_names.index(aai)
        j = standard_aa_names.index(aaj)

        self._tensor[i, j] = value

    def sum(self) -> float:
        return self._tensor.sum().item()

    def count_one(self, aai: str, aaj: str):
        i = standard_aa_names.index(aai)
        j = standard_aa_names.index(aaj)

        self._tensor[i, j] += 1

    def __add__(self, other):
        tensor = self._tensor + other._tensor
        return Matrix(tensor=tensor)

    def to_dict(self) -> Dict[str, Dict[str, float]]:
        d = {aai: {aaj: 0.0 for aaj in standard_aa_names} for aai in standard_aa_names}
        for i in range(20):
            aai = standard_aa_names[i]

            for j in range(20):
                aaj = standard_aa_names[j]

                d[aai][aaj] = self._tensor[i, j].item()

        return d
