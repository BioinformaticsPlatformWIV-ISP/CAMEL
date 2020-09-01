from pathlib import Path
from typing import Dict, Any, List

from camel.app.tools.deconseq.deconseq import Deconseq

SNAKEFILE_DECONSEQ = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_deconseq = Path('deconseq')
OUTPUT_DECONSEQ_CLEAN_PE = _dir_deconseq / 'fastq_pe_clean.io'
OUTPUT_DECONSEQ_CLEAN_SE = _dir_deconseq / 'fastq_se_clean.io'
OUTPUT_DECONSEQ_INFORMS = _dir_deconseq / 'informs.io'
OUTPUT_DECONSEQ_INFORMS_PE_FWD = _dir_deconseq / 'informs_pe_fwd.io'
OUTPUT_DECONSEQ_INFORMS_PE_REV = _dir_deconseq / 'informs_pe_rev.io'
OUTPUT_DECONSEQ_INFORMS_SE_FWD = _dir_deconseq / 'informs_se_fwd.io'
OUTPUT_DECONSEQ_INFORMS_SE_REV = _dir_deconseq / 'informs_se_rev.io'
OUTPUT_DECONSEQ_REPORT = _dir_deconseq / 'report' / 'html.io'
OUTPUT_DECONSEQ_SUMMARY = _dir_deconseq / 'summary.tsv'


def combine_deconseq_informs(deconseq_pe_fwd: Deconseq, deconseq_pe_rev: Deconseq, deconseq_se_fwd: Deconseq, deconseq_se_rev: Deconseq) -> Dict[str, Dict[str, Any]]:
    return {'PE_FWD': deconseq_pe_fwd.informs,
            'PE_REV': deconseq_pe_rev.informs,
            'SE_FWD': deconseq_se_fwd.informs if deconseq_se_fwd else None,
            'SE_REV': deconseq_se_rev.informs if deconseq_se_rev else None}


def get_processed_dbs(informs: Dict[str, Dict[str, Any]]) -> List[str]:
    dbs = set()
    for read_type in informs.keys():
        if informs[read_type] is not None:
            for db in informs[read_type]['processed_dbs']:
                dbs.add(db)
    return sorted(list(dbs))
