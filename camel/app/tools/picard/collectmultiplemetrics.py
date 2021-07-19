import logging
import os
import re

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class CollectMultipleMetrics(Picard):

    """
    Class for Picard CollectMultipleMetrics function to calculate various QC metrics for read mapping

    PROGRAM
        Default values:
        - CollectAlignmentSummaryMetrics
        - CollectBaseDistributionByCycle
        - CollectInsertSizeMetrics
        - MeanQualityByCycle
        - QualityScoreDistribution

        This option can be set to 'null' to clear the default list

        Possible values:
        - CollectAlignmentSummaryMetrics
        - CollectInsertSizeMetrics
        - QualityScoreDistribution
        - MeanQualityByCycle
        - CollectBaseDistributionByCycle
        - CollectGcBiasMetrics
        - RnaSeqMetrics
        - CollectSequencingArtifactMetrics
        - CollectQualityYieldMetrics
    """
    # Suffices for all PROGRAM options
    OUTPUT_FILE_SUFFIX = {
        'CollectAlignmentSummaryMetrics': {'AlignmentSummary': '.alignment_summary_metrics'},
        'CollectInsertSizeMetrics': {'InsertSize': '.insert_size_metrics', 'InsertSizeFigure': '.insert_size_histogram.pdf'},
        'QualityScoreDistribution': {'QualityDistribution': '.quality_distribution_metrics', 'QualityDistributionFigure': '.quality_distribution.pdf'},
        'MeanQualityByCycle': {'QualityByCycle': '.quality_by_cycle_metrics', 'QualityByCycleFigure': '.quality_by_cycle.pdf'},
        'CollectBaseDistributionByCycle': {'BaseDistributionByCycle': '.base_distribution_by_cycle_metrics', 'BaseDistributionByCycleFigure': '.base_distribution_by_cycle.pdf'},
        'CollectGcBiasMetrics': {'GcBias': '.gc_bias.detail_metrics', 'GcBiasSummary': '.gc_bias.summary_metrics','GcBiasFigure': '.gc_bias.pdf'},
        'RnaSeqMetrics': {'RnaSeq': '.rna_metrics'},
        'CollectSequencingArtifactMetrics': {'SequencingArtefactDetail':'.bait_bias_detail_metrics', 'SequencingArtefactSummary':'.bait_bias_summary_metrics',
                                             'SequencingArtefactErrorSummary': '.error_summary_metrics', 'SequencingArtefactPreAdapterDetail': '.pre_adapter_detail_metrics',
                                             'SequencingArtefactPreAdapterSummary': '.pre_adapter_summary_metrics'},
        'CollectQualityYieldMetrics': {'QualityYield': '.quality_yield_metrics'}
    }

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CollectMultipleMetrics', '2.23.3', camel)
        self._function_name = 'CollectMultipleMetrics'
        self._required_inputs = ['FASTA_REF']
        self._specific_parameters = ['metric_accumulation_level_multi']
        self._outfile_prefix = None

    def _set_output(self) -> None:
        """
        Set the output for Picard CollectMultipleMetrics function
        :return: None
        """
        self.outfile_prefix = os.path.join(self._folder, self._parameters['output_prefix'].value)

        for key in self._parameters:
            if re.match('metrics_', key):
                # strip 'metrics_' from key
                tool = key.split("_")[1]

                # loop over all possible output suffices for that specific tool
                try:
                    for suffix_key, suffix_value in CollectMultipleMetrics.OUTPUT_FILE_SUFFIX[tool].items():
                        output_file = self.outfile_prefix + suffix_value
                        self._tool_outputs['TXT_' + suffix_key] = [ToolIOFile(output_file)]
                except KeyError:
                    logging.warning(f'Picard CollectMultipleMetrics unsupported metrics {key}, its results will not be analyzed or returned.')

    def _set_informs(self) -> None:
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        # Known Metrics: CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution,
        #   MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics,
        #   CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
        for key in self._parameters:
            if re.match('metrics_', key):
                if key == 'metrics_CollectAlignmentSummaryMetrics':
                    self.__analyze_alignment_summary()
                elif key == 'metrics_CollectInsertSizeMetrics':
                    self.__analyze_insert_size()
                elif key == 'metrics_CollectGcBiasMetrics':
                    self.__analyze_gc_bias()
                elif key == 'metrics_QualityScoreDistribution':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                elif key == 'metrics_CollectBaseDistributionByCycle':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                elif key == 'metrics_MeanQualityByCycle':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                else:  # RnaSeqMetrics, CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
                    logging.warning(
                        f'Picard CollectMultipleMetrics unsupported metrics {key}, its results will not be analyzed or returned.')

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        build_options = self._build_options(excluded_parameters=self._specific_parameters, delimiter='=')

        if 'metric_accumulation_level_multi' in self._parameters:
            attributes = str(self._parameters['metric_accumulation_level_multi'].value).split(",")
            for attribute in attributes:
                build_options.append(f"LEVEL={attribute}")

        option_string = " ".join(build_options)

        self._command.command = " ".join([
            "java", self._java_options, "-jar $PICARD_JAR", self._tool_command, self._java_options_temp_dir,
            self._input_string, self._output_string, option_string, '2>&1'
        ])

    def __analyze_alignment_summary(self) -> None:
        """
        Analyze the Alignment Summary statistics
        :return: None
        """
        # Example content of AlignmentSummaryMetrics
        #
        # ## METRICS CLASS        picard.analysis.AlignmentSummaryMetrics
        # CATEGORY        TOTAL_READS     PF_READS        PCT_PF_READS    PF_NOISE_READS  PF_READS_ALIGNED        PCT_PF_READS_ALIGNED    PF_ALIGNED_BASES        PF_HQ_ALIGNED_READS     PF_HQ_ALIGNED_BASES     PF_HQ_ALIGNED_Q20_BASES PF_HQ_MEDIAN_MISMATCHES PF_MISMATCH_RATE        PF_HQ_ERROR_RATE        PF_INDEL_RATE   MEAN_READ_LENGTH        READS_ALIGNED_IN_PAIRS  PCT_READS_ALIGNED_IN_PAIRS      PF_READS_IMPROPER_PAIRS PCT_PF_READS_IMPROPER_PAIRS     BAD_CYCLES      STRAND_BALANCE  PCT_CHIMERAS    PCT_ADAPTER     SAMPLE  LIBRARY READ_GROUP
        # FIRST_OF_PAIR   268167  268167  1       0       254624  0.949498        57016421        254623  57016401        56709951        2       0.009547        0.009547        0.000105        224.314263      254562  0.999757        444     0.001744        0       0.508267        0.002282        0.000004
        # SECOND_OF_PAIR  268167  268167  1       0       254581  0.949338        50677959        254580  50677938        49922898        2       0.010415        0.010415        0.00011 199.397465      254562  0.999925        401     0.001575        0       0.492111        0.002078        0
        # PAIR    536334  536334  1       0       509205  0.949418        107694380       509203  107694339       106632849       2       0.009955        0.009955        0.000107        211.855864      509124  0.999841        845     0.001659        0       0.50019 0.00218 0.000002
        #
        # NOTE: only PAIR statistics (last line) are extracted
        # Ref: http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics
        align_stats = {}
        with open(self._tool_outputs['TXT_AlignmentSummary'][0].path, 'r') as inf:
            col_nbs = None
            for line in inf:
                if not line.startswith('#') and 'TOTAL_READS' in line and col_nbs is None:  # First line after comments is the header
                    col_nbs = {x: i for i, x in enumerate(line.strip().split('\t'))}
                if re.match('PAIR', line):
                    pair_informs = line.split('\t')
                    align_stats['TOTAL_READS'] = pair_informs[col_nbs['TOTAL_READS']]
                    # PF: Passing (Illumina) Filter
                    align_stats['PassingFilter_READS'] = pair_informs[col_nbs['PF_READS']]
                    align_stats['PF_NOISE_READS'] = pair_informs[col_nbs['PF_NOISE_READS']]
                    align_stats['PF_READS_ALIGNED'] = pair_informs[col_nbs['PF_READS_ALIGNED']]
                    align_stats['PCT_PF_READS_ALIGNED'] = pair_informs[col_nbs['PCT_PF_READS_ALIGNED']]
                    align_stats['PF_ALIGNED_BASES'] = pair_informs[col_nbs['PF_ALIGNED_BASES']]
                    align_stats['PF_MISMATCH_RATE'] = pair_informs[col_nbs['PF_MISMATCH_RATE']]
                    align_stats['PF_INDEL_RATE'] = pair_informs[col_nbs['PF_INDEL_RATE']]
                    align_stats['MEAN_READ_LENGTH'] = pair_informs[col_nbs['MEAN_READ_LENGTH']]
                    align_stats['READS_ALIGNED_IN_PAIRS'] = pair_informs[col_nbs['READS_ALIGNED_IN_PAIRS']]
                    align_stats['PCT_READS_ALIGNED_IN_PAIRS'] = pair_informs[col_nbs['PCT_READS_ALIGNED_IN_PAIRS']]

                    break

        self.informs['AlignmentSummary_stats'] = align_stats

    def __analyze_insert_size(self) -> None:
        """
        Analyze the Insert Size statistics
        :return: None
        """
        # Example content:
        #
        # METRICS CLASS        picard.analysis.InsertSizeMetrics
        # MEDIAN_INSERT_SIZE      MODE_INSERT_SIZE        MEDIAN_ABSOLUTE_DEVIATION       MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE        STANDARD_DEVIATION      READ_PAIRS      PAIR_ORIENTATION        WIDTH_OF_10_PERCENT     WIDTH_OF_20_PERCENT     WIDTH_OF_30_PERCENT     WIDTH_OF_40_PERCENT     WIDTH_OF_50_PERCENT     WIDTH_OF_60_PERCENT     WIDTH_OF_70_PERCENT     WIDTH_OF_80_PERCENT     WIDTH_OF_90_PERCENT     WIDTH_OF_95_PERCENT     WIDTH_OF_99_PERCENT     SAMPLE  LIBRARY READ_GROUP
        # 357     329     105     19      2316    371.815165      166.190341      254312  FR      41      79      121     163     211     267     335     417     535     635     979                     438     22      24      4853839 431.875202      61.702551       497091  FR      9       17      27      35      45      57      69      85      111     763
        #
        data_row = False
        self.informs['InsertSize_stats'] = {}
        with open(self._tool_outputs['TXT_InsertSize'][0].path, 'r') as inf:
            for line in inf:
                if re.match('MEDIAN_INSERT_SIZE', line):
                    data_row = True
                    col_nbs = {x: i for i, x in enumerate(line.strip().split('\t'))}
                    continue
                if data_row:
                    if line.strip() == '':
                        break
                    informs = line.split('\t')
                    # Only take statistics for FR (READ_PAIRS) reads
                    if informs[col_nbs['PAIR_ORIENTATION']] == 'FR':
                        self.informs['InsertSize_stats']['MEDIAN_INSERT_SIZE'] = informs[col_nbs['MEDIAN_INSERT_SIZE']]
                        self.informs['InsertSize_stats']['MEDIAN_ABSOLUTE_DEVIATION'] = informs[col_nbs['MEDIAN_ABSOLUTE_DEVIATION']]
                        self.informs['InsertSize_stats']['MEAN_INSERT_SIZE'] = informs[col_nbs['MEAN_INSERT_SIZE']]
                        self.informs['InsertSize_stats']['STANDARD_DEVIATION'] = informs[col_nbs['STANDARD_DEVIATION']]
                    else:
                        logging.warning(
                            f'Picard CollectMultipleMetrics unsupported READS ORIENTATION {informs[col_nbs["PAIR_ORIENTATION"]]} for InsertSize analysis')


    def __analyze_gc_bias(self) -> None:
        """
        Analyze the GC Bias statistics
        :return: None
        """
        # Only GC Bias summary is analyzed
        self.informs['GCBias_stats'] = {}
        with open(self._tool_outputs['TXT_GcBiasSummary'][0].path, 'r') as inf:
            for line in inf:
                if re.match('All Reads', line):
                    informs = line.split("\t")
                    self.informs['GCBias_stats']['AT_DROPOUT'] = informs[4]
                    self.informs['GCBias_stats']['GC_DROPOUT'] = informs[5]
                    break
