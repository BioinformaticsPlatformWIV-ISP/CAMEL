from dataclasses import dataclass
from typing import List, Dict

from camel.app.io.tooliofile import ToolIOFile


@dataclass(frozen=True)
class FastqInput:
    read_type: str
    pe: List[ToolIOFile] = None
    se: List[ToolIOFile] = None
    se_fwd: List[ToolIOFile] = None
    se_rev: List[ToolIOFile] = None
    is_trimmed: bool = False
    is_pe: bool = True

    @property
    def is_paired(self) -> bool:
        """
        Returns true if the read input is paired.
        :return: True if paired, False otherwise
        """
        return self.pe is not None

    def to_fq_dict(self) -> Dict[str, List[ToolIOFile]]:
        """
        Converts the FASTQ input to a input dictionary for the workflow.
        :return: FASTQ dictionary
        """
        if self.is_pe:
            fq_dict = {'PE': self.pe}
            if self.se_fwd is not None:
                fq_dict['SE_FWD'] = self.se_fwd
            if self.se_rev is not None:
                fq_dict['SE_REV'] = self.se_rev
            return fq_dict
        else:
            return {'SE': self.se}
