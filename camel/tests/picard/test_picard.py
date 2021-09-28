import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups
from camel.app.tools.picard.calculatereadgroupchecksum import CalculateReadGroupChecksum
from camel.app.tools.picard.calculatereadgroupchecksum import CollectMultipleMetrics
from camel.app.tools.picard.collectqualityyieldmetrics import CollectQualityYieldMetrics
from camel.app.tools.picard.collectrawwgsmetrics import CollectRawWgsMetrics
from camel.app.tools.picard.collectrawwgsmetrics import CollectWgsMetrics
from camel.app.tools.picard.collectvariantcallingmetrics import CollectVariantCallingMetrics
from camel.app.tools.picard.createsequencedictionary import CreateSequenceDictionary
from camel.app.tools.picard.fastqtosam import FastqToSam
from camel.app.tools.picard.gatherbamfiles import GatherBamFiles
from camel.app.tools.picard.intervallisttools import IntervalListTools
from camel.app.tools.picard.markduplicates import MarkDuplicates
from camel.app.tools.picard.mergebamalignment import MergeBamAlignment
from camel.app.tools.picard.mergevcfs import MergeVCFs
from camel.app.tools.picard.samtofastq import SamToFastq
from camel.app.tools.picard.setnmmdanduqtags import SetNmMdAndUqTags
from camel.app.tools.picard.sortsam import SortSam
from camel.app.tools.picard.validatesamefile import ValidateSamFile

class TestPicard(CamelTestSuite):
    """
    Tests the Picard tool suite.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('picard')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'Homo_sapiens_assembly38_chr22.fasta')
    FILE_VCFdb = ToolIOFile(test_file_dir / 'Homo_sapiens_assembly38.dbsnp138_chr22.vcf.gz')
    FILE_BAM = ToolIOFile(test_file_dir / 'NA12877_chr22_25p.bam')
    FILE_VCF = ToolIOFile(test_file_dir / 'NA12877_chr22.gVCF')
    FILE_VCF1 = ToolIOFile(test_file_dir / 'NA12877_sub1.vcf')
    FILE_VCF2 = ToolIOFile(test_file_dir / 'NA12877_sub2.vcf')
    FILE_VCF3 = ToolIOFile(test_file_dir / 'NA12877_sub3.vcf')
    FILE_FASTQ_R1 = ToolIOFile(test_file_dir / 'NA12877_R1.fastq')
    FILE_FASTQ_R2 = ToolIOFile(test_file_dir / 'NA12877_R2.fastq')
    FILE_BAM1 = ToolIOFile(test_file_dir / 'NA12877_chr22_25p_sub1.bam')
    FILE_BAM2 = ToolIOFile(test_file_dir / 'NA12877_chr22_25p_sub2.bam')
    FILE_BAM3 = ToolIOFile(test_file_dir / 'NA12877_chr22_25p_sub3.bam')
    FILE_uBAM = ToolIOFile(test_file_dir / 'NA12877.unmapped.bam')


    def test_picard_addorreplacereadgroups(self) -> None:
        """
        Test Picard AddOrReplaceReadgroups
        :return: None
        """
        picard_addorreplacereadgroups = AddOrReplaceReadGroups(self.camel)
        picard_addorreplacereadgroups.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_addorreplacereadgroups.run(self.running_dir)
        self.assertTrue('BAM' in picard_addorreplacereadgroups.tool_outputs, "No BAM output generated")
        output_file = Path(picard_addorreplacereadgroups.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_calculatereadgroupchecksum(self) -> None:
        """
        Test Picard CalculateReadgroupChecksum
        :return: None
        """
        picard_calculatereadgroupchecksum = CalculateReadGroupChecksum(self.camel)
        picard_calculatereadgroupchecksum.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_calculatereadgroupchecksum.run(self.running_dir)
        self.assertTrue('TXT_checksum' in picard_calculatereadgroupchecksum.tool_outputs, "No TXT_checksum output generated")
        output_file = Path(picard_calculatereadgroupchecksum.tool_outputs['TXT_checksum'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_collectmultiplemetrics(self) -> None:
        """
        Test Picard CollectMultipleMetrics
        :return: None
        """
        picard_collectmultiplemetrics = CollectMultipleMetrics(self.camel)
        picard_collectmultiplemetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_collectmultiplemetrics.update_parameters(
            assume_sorted = 'true',
            reset_metrics = 'null',
            metrics_CollectAlignmentSummaryMetrics = 'CollectAlignmentSummaryMetrics',
            metrics_CollectGcBiasMetrics = 'CollectGcBiasMetrics',
            metrics_CollectInsertSizeMetrics = 'CollectInsertSizeMetrics',
            metrics_QualityScoreDistribution  = 'QualityScoreDistribution ',
            metrics_MeanQualityByCycle = 'MeanQualityByCycle',
            metrics_CollectBaseDistributionByCycle = 'CollectBaseDistributionByCycle',
            metrics_CollectSequencingArtifactMetrics = 'CollectSequencingArtifactMetrics',
            metrics_CollectQualityYieldMetrics = 'CollectQualityYieldMetrics'
        )
        picard_collectmultiplemetrics.run(self.running_dir)
        self.assertTrue('TXT_CollectAlignmentSummaryMetrics' in picard_collectmultiplemetrics,
                        "No TXT_CollectAlignmentSummaryMetrics output generated")
        self.assertTrue('TXT_CollectInsertSizeMetrics' in picard_collectmultiplemetrics,
                        "No TXT_CollectInsertSizeMetrics output generated")
        self.assertTrue('TXT_QualityScoreDistribution' in picard_collectmultiplemetrics,
                        "No TXT_QualityScoreDistribution output generated")
        self.assertTrue('TXT_MeanQualityByCycle' in picard_collectmultiplemetrics,
                        "No TXT_MeanQualityByCycle output generated")
        self.assertTrue('TXT_CollectBaseDistributionByCycle' in picard_collectmultiplemetrics,
                        "No TXT_CollectBaseDistributionByCycle output generated")
        self.assertTrue('TXT_CollectGcBiasMetrics' in picard_collectmultiplemetrics,
                        "No TXT_CollectGcBiasMetrics output generated")
        self.assertTrue('TXT_CollectSequencingArtifactMetrics' in picard_collectmultiplemetrics,
                        "No TXT_CollectSequencingArtifactMetrics output generated")
        self.assertTrue('TXT_CollectQualityYieldMetrics' in picard_collectmultiplemetrics,
                        "No TXT_CollectQualityYieldMetrics output generated")

        output_file_CollectAlignmentSummaryMetrics = Path(picard_collectmultiplemetrics.tool_outputs['TXT_CollectAlignmentSummaryMetrics'][0].path)
        output_file_CollectInsertSizeMetrics = Path(picard_collectmultiplemetrics.tool_outputs['TXT_CollectInsertSizeMetrics'][0].path)
        output_file_QualityScoreDistribution = Path(picard_collectmultiplemetrics.tool_outputs['TXT_QualityScoreDistribution'][0].path)
        output_file_MeanQualityByCycle = Path(picard_collectmultiplemetrics.tool_outputs['TXT_MeanQualityByCycle'][0].path)
        output_file_CollectBaseDistributionByCycle = Path(picard_collectmultiplemetrics.tool_outputs['TXT_CollectBaseDistributionByCycle'][0].path)
        output_file_CollectGcBiasMetrics = Path(picard_collectmultiplemetrics.tool_outputs['TXT_CollectGcBiasMetrics'][0].path)
        output_file_CollectSequencingArtifactMetrics = Path(picard_collectmultiplemetrics.tool_outputs['TXT_CollectSequencingArtifactMetrics'][0].path)
        output_file_CollectQualityYieldMetrics = Path(picard_collectmultiplemetrics.tool_outputs['TXT_CollectQualityYieldMetrics'][0].path)

        self.assertTrue(output_file_CollectAlignmentSummaryMetrics.exists())
        self.assertTrue(output_file_CollectInsertSizeMetrics.exists())
        self.assertTrue(output_file_QualityScoreDistribution.exists())
        self.assertTrue(output_file_MeanQualityByCycle.exists())
        self.assertTrue(output_file_CollectBaseDistributionByCycle.exists())
        self.assertTrue(output_file_CollectGcBiasMetrics.exists())
        self.assertTrue(output_file_CollectSequencingArtifactMetrics.exists())
        self.assertTrue(output_file_CollectQualityYieldMetrics.exists())

        self.assertGreater(output_file_CollectAlignmentSummaryMetrics.stat().st_size, 0)
        self.assertGreater(output_file_CollectInsertSizeMetrics.stat().st_size, 0)
        self.assertGreater(output_file_QualityScoreDistribution.stat().st_size, 0)
        self.assertGreater(output_file_MeanQualityByCycle.stat().st_size, 0)
        self.assertGreater(output_file_CollectBaseDistributionByCycle.stat().st_size, 0)
        self.assertGreater(output_file_CollectGcBiasMetrics.stat().st_size, 0)
        self.assertGreater(output_file_CollectSequencingArtifactMetrics.stat().st_size, 0)
        self.assertGreater(output_file_CollectQualityYieldMetrics.stat().st_size, 0)

    def test_picard_collectqualityyieldmetrics(self) -> None:
        """
        Test Picard CollectQualityYieldMetrics
        :return: None
        """
        picard_collectqualityyieldmetrics = CollectQualityYieldMetrics(self.camel)
        picard_collectqualityyieldmetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_collectqualityyieldmetrics.run(self.running_dir)
        self.assertTrue('TXT' in picard_collectqualityyieldmetrics.tool_outputs,
                        "No TXT output generated")
        output_file = Path(picard_collectqualityyieldmetrics.tool_outputs['TXT'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_collectrawwgsmetrics(self) -> None:
        """
        Test Picard CollectRawWgsMetrics
        :return: None
        """
        picard_collectrawwgsmetrics = CollectRawWgsMetrics(self.camel)
        picard_collectrawwgsmetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_collectrawwgsmetrics.run(self.running_dir)
        self.assertTrue('TXT_metrics' in picard_collectrawwgsmetrics.tool_outputs,
                        "No TXT_metrics output generated")
        output_file = Path(picard_collectrawwgsmetrics.tool_outputs['TXT_metrics'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_collectwgsmetrics(self) -> None:
        """
        Test Picard CollectWgsMetrics
        :return: None
        """
        picard_collectwgsmetrics = CollectWgsMetrics(self.camel)
        picard_collectwgsmetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_collectwgsmetrics.run(self.running_dir)
        self.assertTrue('TXT_metrics' in picard_collectwgsmetrics.tool_outputs,
                        "No TXT_metrics output generated")
        output_file = Path(picard_collectwgsmetrics.tool_outputs['TXT_metrics'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


    def test_picard_collectvariantcallingmetrics(self) -> None:
        """
        Test Picard CollectVariantCallingMetrics
        :return: None
        """
        picard_collectvariantcallingmetrics = CollectVariantCallingMetrics(self.camel)
        picard_collectvariantcallingmetrics.add_input_files({
            'VCF': [TestPicard.FILE_VCF],
            'VCF_dbsnp': [TestPicard.FILE_VCFdb]
        })
        picard_collectvariantcallingmetrics.run(self.running_dir)
        self.assertTrue('TXT_report' in picard_collectvariantcallingmetrics.tool_outputs,
                        "No TXT_report output generated")
        output_file_summary = Path(picard_collectvariantcallingmetrics.tool_outputs['TXT_report'][0].path)
        self.assertTrue(output_file_summary.exists())
        self.assertGreater(output_file_summary.stat().st_size, 0)

        output_file_details = Path(picard_collectvariantcallingmetrics.tool_outputs['TXT_report'][0].path)
        self.assertTrue(output_file_details.exists())
        self.assertGreater(output_file_details.stat().st_size, 0)


def test_picard_createsequencedictionary(self) -> None:
        """
        Test Picard CreateSequenceDictionary
        :return: None
        """
        picard_createsequencedictionary = CreateSequenceDictionary(self.camel)
        picard_createsequencedictionary.add_input_files({
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_createsequencedictionary.run(self.running_dir)
        self.assertTrue('FASTA_REF' in picard_createsequencedictionary.tool_outputs, "No FASTA_REF output generated")
        output_file = Path(picard_createsequencedictionary.tool_outputs['FASTA_REF'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


    def test_picard_fastqtosam(self) -> None:
        """
        Test Picard FastqToSam
        :return: None
        """
        picard_fastqtosam = FastqToSam(self.camel)
        picard_fastqtosam.add_input_files({
            'FASTQ_PE': [TestPicard.FILE_FASTQ_R1, TestPicard.FILE_FASTQ_R2]
        })
        picard_fastqtosam.run(self.running_dir)
        self.assertTrue('SAM' in picard_fastqtosam.tool_outputs, "No SAM output generated")
        output_file = Path(picard_fastqtosam.tool_outputs['SAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_gatherbamfiles(self) -> None:
        """
        Test Picard GatherBamFiles
        :return: None
        """
        picard_gatherbamfiles = GatherBamFiles(self.camel)
        picard_gatherbamfiles.add_input_files({
            'BAMs': [TestPicard.FILE_BAM1, ToolIOFile.FILE_BAM2, ToolIOFile.FILE_BAM3]
        })
        picard_gatherbamfiles.run(self.running_dir)
        self.assertTrue('BAM' in picard_gatherbamfiles.tool_outputs, "No BAM output generated")
        output_file = Path(picard_gatherbamfiles.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_intervallisttools(self) -> None:
        """
        Test Picard IntervalListTools
        :return: None
        """
        picard_intervallisttools = IntervalListTools(self.camel)
        picard_intervallisttools.add_input_files({
            'VCF': [TestPicard.FILE_VCF]
        })
        picard_intervallisttools.run(self.running_dir)
        self.assertTrue('TXT_intervalLists' in picard_intervallisttools.tool_outputs, "No TXT_intervalLists output generated")
        output_file = Path(picard_intervallisttools.tool_outputs['TXT_intervalLists'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


def test_picard_markduplicates(self) -> None:
        """
        Test Picard MarkDuplicates
        :return: None
        """
        picard_markduplicates = MarkDuplicates(self.camel)
        picard_markduplicates.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_markduplicates.run(self.running_dir)
        self.assertTrue('BAM' in picard_markduplicates.tool_outputs, "No BAM output generated")
        output_file = Path(picard_markduplicates.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

        self.assertTrue('METRICS' in picard_markduplicates.tool_outputs, "No METRICS output generated")
        output_file_metrics = Path(picard_markduplicates.tool_outputs['METRICS'][0].path)
        self.assertTrue(output_file_metrics.exists())
        self.assertGreater(output_file_metrics.stat().st_size, 0)

    def test_picard_mergebamalignment(self) -> None:
        """
        Test Picard MergeBamAlignment
        :return: None
        """
        picard_mergebamalignment = MergeBamAlignment(self.camel)
        picard_mergebamalignment.add_input_files({
            'BAM_UNMAPPED': [TestPicard.FILE_uBAM],
            'BAM_ALIGNED': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_mergebamalignment.run(self.running_dir)
        self.assertTrue('BAM' in picard_mergebamalignment.tool_outputs,"No BAM output generated")
        output_file = Path(picard_mergebamalignment.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_mergevcfs(self) -> None:
        """
        Test Picard MergeVCFs
        :return: None
        """
        picard_mergevcfs = MergeVCFs(self.camel)
        picard_mergevcfs.add_input_files({
            'VCF': [TestPicard.FILE_VCF1, TestPicard.FILE_VCF2, TestPicard.FILE_VCF3],
        })
        picard_mergevcfs.run(self.running_dir)
        self.assertTrue('VCF' in picard_mergevcfs.tool_outputs, "No VCF output generated")
        output_file = Path(picard_mergevcfs.tool_outputs['VCF'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_samtofastq(self) -> None:
        """
        Test Picard SamToFastq
        :return: None
        """
        picard_samtofastq = SamToFastq(self.camel)
        picard_samtofastq.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
        })
        picard_samtofastq.run(self.running_dir)
        self.assertTrue('FASTQ' in picard_samtofastq.tool_outputs, "No FASTQ output generated")
        output_file = Path(picard_samtofastq.tool_outputs['FASTQ'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_setnmmdanduqtags(self) -> None:
        """
        Test Picard SetNmMdAndUqTags
        :return: None
        """
        picard_setnmmdanduqtags = SetNmMdAndUqTags(self.camel)
        picard_setnmmdanduqtags.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_setnmmdanduqtags.run(self.running_dir)
        self.assertTrue('BAM' in picard_setnmmdanduqtags.tool_outputs, "No BAM output generated")
        output_file = Path(picard_setnmmdanduqtags.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_sortsam(self) -> None:
        """
        Test Picard SortSam
        :return: None
        """
        picard_sortsam = SortSam(self.camel)
        picard_sortsam.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_sortsam.run(self.running_dir)
        self.assertTrue('BAM' in picard_sortsam.tool_outputs, "No BAM output generated")
        output_file = Path(picard_sortsam.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_picard_validatesamfile(self) -> None:
        """
        Test Picard ValidateSamFile
        :return: None
        """
        picard_validatesamfile = ValidateSamFile(self.camel)
        picard_validatesamfile.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
        })
        picard_validatesamfile.run(self.running_dir)
        self.assertTrue('TXT_report' in picard_validatesamfile.tool_outputs, "No TXT_report output generated")
        output_file = Path(picard_validatesamfile.tool_outputs['TXT_report'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()