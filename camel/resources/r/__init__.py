import os

_current_dir = os.path.dirname(__file__)

RSCRIPT_QC_COVERAGE = os.path.join(_current_dir, 'qc_visualize_coverage.r')
RSCRIPT_QC_MAPPING_RATE = os.path.join(_current_dir, 'qc_visualize_mapping_rate.r')
RSCRIPT_QC_ST = os.path.join(_current_dir, 'qc_visualize_st.r')
