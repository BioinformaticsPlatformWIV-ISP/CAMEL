import logging
import os
import re

from app.io.tooliofile import ToolIOFile
from app.tools.picard.picard import Picard


class CollectMultipleMetrics(Picard):
    """
    Class for Picard CollectMultipleMetrics function to calculate various QC matrics for readmap

    Matrics calculated:
    - CollectAlignmentSummaryMetrics
    - CollectInsertSizeMetrics
    - QualityScoreDistribution
    - CollectGcBiasMetrics

    run_picard.sh CollectMultipleMetrics R=/data/refgenomes/influenza.A/cyril/genomes/A_Alabama_05_2010-H3N2.fa I=bwa_readmap.sorted.bam O=bwa_readmap PROGRAM=CollectInsertSizeMetrics PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=CollectGcBiasMetrics
    """
    OUTPUT_FILE_SUFFIX = {
        'AlignmentSummary': '.alignment_summary_metrics',
        'BaseDistributionByCycle': '.base_distribution_by_cycle_metrics',
        'BaseDistributionByCycleFigure': '.base_distribution_by_cycle.pdf',
        'GcBias': '.gc_bias.detail_metrics',
        'GcBiasSummary': '.gc_bias.summary_metrics',
        'GcBiasFigure': '.gc_bias.pdf',
        'InsertSize': '.insert_size_metrics',
        'InsertSizeFigure': '.insert_size_histogram.pdf',
        'MapQualityDistribution': '.quality_distribution_metrics',
        'MapQualityDistributionFigure': '.quality_distribution.pdf',
        'QualityByCycle': '.quality_by_cycle_metrics',
        'QualityByCycleFigure': '.quality_by_cycle.pdf'
    }

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(CollectMultipleMetrics, self).__init__('Picard CollectMultipleMetrics', '2.6.0', camel)
        self._function_name = 'CollectMultipleMetrics'
        self._required_inputs = ['FASTA_REF']
        self._outfile_prefix = None

    def _set_output(self):
        """
        Set the output for Picard CollectMultipleMetrics function
        :return: None
        """
        self.outfile_prefix = os.path.join(self._folder, self._parameters['output_prefix'].value)

        # Known Matrics: CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution,
        #   MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics,
        #   CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
        for key in self._parameters:
            if re.match('matrics_', key):
                if key == 'matrics_CollectAlignmentSummaryMetrics':
                    self.__set_alignement_summary_output()
                elif key == 'matrics_CollectInsertSizeMetrics':
                    self.__set_insert_size_output()
                elif key == 'matrics_QualityScoreDistribution':
                    self.__set_mapping_quality_output()
                elif key == 'matrics_CollectGcBiasMetrics':
                    self.__set_gc_bias_output()
                else:  # MeanQualityByCycle, CollectBaseDistributionByCycle, RnaSeqMetrics, CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
                    logging.warning(
                        "Picard CollectMultipleMetrics unsupported matrics {!r}, its results will not be analyzed or returned.".format(key))

    def _set_inform(self):
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        # Known Matrics: CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution,
        #   MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics,
        #   CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
        for key in self._parameters:
            if re.match('matrics_', key):
                if key == 'matrics_CollectAlignmentSummaryMetrics':
                    self.__analyze_alignement_summary()
                elif key == 'matrics_CollectInsertSizeMetrics':
                    self.__analyze_insert_size()
                elif key == 'matrics_CollectGcBiasMetrics':
                    self.__analyze_gc_bias()
                elif key == 'matrics_QualityScoreDistribution':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                else:  # MeanQualityByCycle, CollectBaseDistributionByCycle, RnaSeqMetrics, CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
                    logging.warning(
                        "Picard CollectMultipleMetrics unsupported matrics {!r}, its results will not be analyzed or returned.".format(key))

    def __set_alignement_summary_output(self):
        """
        set the Alignment Summary statistics output
        :return: None
        """
        alignment_summary_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['AlignmentSummary']
        self._tool_outputs['TXT_AlignmentSummary'] = [ToolIOFile(alignment_summary_file)]

    def __analyze_alignement_summary(self):
        """
        Analyze the Alignment Summary statistics
        :return: None
        """
        # Example content of AlignmentSummaryMetrics
        #
        # ## METRICS CLASS        picard.analysis.AlignmentSummaryMetrics
        #         CATEGORY        TOTAL_READS     PF_READS        PCT_PF_READS    PF_NOISE_READS  PF_READS_ALIGNED        PCT_PF_READS_ALIGNED    PF_ALIGNED_BASES        PF_HQ_ALIGNED_READS     PF_HQ_ALIGNED_BASES     PF_HQ_ALIGNED_Q20_BASES PF_HQ_MEDIAN_MISMATCHES PF_MISMATCH_RATE        PF_HQ_ERROR_RATE        PF_INDEL_RATE   MEAN_READ_LENGTH        READS_ALIGNED_IN_PAIRS  PCT_READS_ALIGNED_IN_PAIRS      BAD_CYCLES      STRAND_BALANCE  PCT_CHIMERAS    PCT_ADAPTER     SAMPLE  LIBRARY READ_GROUP
        # FIRST_OF_PAIR   545327  545327  1       0       515082  0.944538        64189044        498204  62163596        62093230        0       0.003405        0.003358        0.000026        125.249694      505685  0.981756        0       0.501747        0.002396        0.000013
        # SECOND_OF_PAIR  545327  545327  1       0       509899  0.935033        57174661        493387  55411098        55322270        0       0.000938        0.000885        0.000025        111.532987      505685  0.991736        0       0.498722        0.002396        0
        # PAIR    1090654 1090654 1       0       1024981 0.939786         121363705       991591  117574694       117415500       0       0.002243 0.002193        0.000025        118.39134       1011370 0.986721 0       0.500242   0.002396        0.000006
        #
        # NOTE: only PAIR statistics (last line) are extracted
        # Ref: http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics
        align_stats = {}
        with open(self._tool_outputs['TXT_AlignmentSummary'][0].path, 'r') as inf:
            for l in inf.readlines():
                if re.match('PAIR', l):
                    pair_informs = l.split("\t")
                    align_stats['TOTAL_READS'] = pair_informs[1]
                    # PF: Passing (Illumina) Filter
                    align_stats['PassingFilter_READS'] = pair_informs[2]
                    align_stats['PF_NOISE_READS'] = pair_informs[4]
                    align_stats['PF_READS_ALIGNED'] = pair_informs[5]
                    align_stats['PCT_PF_READS_ALIGNED'] = pair_informs[6]
                    align_stats['PF_ALIGNED_BASES'] = pair_informs[7]
                    align_stats['PF_MISMATCH_RATE'] = pair_informs[12]
                    align_stats['PF_INDEL_RATE'] = pair_informs[14]
                    align_stats['MEAN_READ_LENGTH'] = pair_informs[15]
                    align_stats['READS_ALIGNED_IN_PAIRS'] = pair_informs[16]
                    align_stats['PCT_READS_ALIGNED_IN_PAIRS'] = pair_informs[17]

                    break

        self.informs['AlignmentSummary_stats'] = align_stats

    def __set_insert_size_output(self):
        """
        set the Insert Size statistics output
        :return: None
        """
        metrics_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['InsertSize']
        self._tool_outputs['TXT_InsertSize'] = [ToolIOFile(metrics_file)]
        figure_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['InsertSizeFigure']
        if os.path.exists(figure_file):
            self._tool_outputs['PDF_InsertSizeFigure'] = [ToolIOFile(figure_file)]

    def __analyze_insert_size(self):
        """
        Analyze the Insert Size statistics
        :return: None
        """
        # Example content:
        #
        # ## METRICS CLASS        picard.analysis.InsertSizeMetrics
        # MEDIAN_INSERT_SIZE      MEDIAN_ABSOLUTE_DEVIATION       MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE        STANDARD_DEVIATION      READ_PAIRS      PAIR_ORIENTATION        WIDTH_OF_10_PERCENT     WIDTH_OF_20_PERCENT     WIDTH_OF_30_PERCENT     WIDTH_OF_40_PERCENT     WIDTH_OF_50_PERCENT     WIDTH_OF_60_PERCENT     WIDTH_OF_70_PERCENT     WIDTH_OF_80_PERCENT     WIDTH_OF_90_PERCENT     WIDTH_OF_99_PERCENT     SAMPLE  LIBRARY READ_GROUP
        # 438     22      24      4853839 431.875202      61.702551       497091  FR      9       17      27      35      45      57      69      85      111     763
        #
        data_row = False
        self.informs['InsertSize_stats'] = {}
        with open(self._tool_outputs['TXT_InsertSize'][0].path, 'r') as inf:
            for l in inf.readlines():
                if re.match('MEDIAN_INSERT_SIZE', l):
                    data_row = True
                    continue
                if data_row:
                    if l.strip() == '':
                        break
                    informs = l.split("\t")
                    # Only take statistics for FR reads
                    if informs[7] == 'FR':
                        self.informs['InsertSize_stats']['MEDIAN_INSERT_SIZE'] = informs[0]
                        self.informs['InsertSize_stats']['MEDIAN_ABSOLUTE_DEVIATION'] = informs[1]
                        self.informs['InsertSize_stats']['MEAN_INSERT_SIZE'] = informs[4]
                        self.informs['InsertSize_stats']['STANDARD_DEVIATION'] = informs[5]
                    else:
                        logging.warning(
                            "Picard CollectMultipleMetrics unsupported READS ORIENTATION {} for InsertSize analysis".format(informs[7]))

    def __set_mapping_quality_output(self):
        """
        set the Mapping Quality statistics output
        :return: None
        """
        metrics_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['MapQualityDistribution']
        self._tool_outputs['TXT_MapQualityDistribution'] = [ToolIOFile(metrics_file)]
        figure_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['MapQualityDistributionFigure']
        self._tool_outputs['PDF_MapQualityDistributionFigure'] = [ToolIOFile(figure_file)]

    def __set_gc_bias_output(self):
        """
        set the GC Bias statistics output
        :return: None
        """
        metrics_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['GcBias']
        self._tool_outputs['TXT_GcBiasDetail'] = [ToolIOFile(metrics_file)]
        summary_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['GcBiasSummary']
        self._tool_outputs['TXT_GcBiasSummary'] = [ToolIOFile(summary_file)]
        figure_file = self.outfile_prefix + CollectMultipleMetrics.OUTPUT_FILE_SUFFIX['GcBiasFigure']
        self._tool_outputs['PDF_GcBiasFigure'] = [ToolIOFile(figure_file)]

    def __analyze_gc_bias(self):
        """
        Analyze the GC Bias statistics
        :return: None
        """
        # Only GC Bias summary is analyzed
        self.informs['GCBias_stats'] = {}
        with open(self._tool_outputs['TXT_GcBiasSummary'][0].path, 'r') as inf:
            for l in inf.readlines():
                if re.match('All Reads', l):
                    informs = l.split("\t")
                    self.informs['GCBias_stats']['AT_DROPOUT'] = informs[4]
                    self.informs['GCBias_stats']['GC_DROPOUT'] = informs[5]
                    break
