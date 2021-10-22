import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.gatk4.gatk4applybqsr import GATK4ApplyBQSR
from camel.app.tools.gatk4.gatk4baserecalibrator import GATK4BaseRecalibrator
from camel.app.tools.gatk4.gatk4combinegvcfs import GATK4CombineGVCFs
from camel.app.tools.gatk4.gatk4fastaalternatereferencemaker import GATK4FastaAlternateReferenceMaker
from camel.app.tools.gatk4.gatk4gatherbqsrreports import GATK4GatherBQSRReports
from camel.app.tools.gatk4.gatk4genotypegvcfs import GATK4GenotypeGVCFs
from camel.app.tools.gatk4.gatk4haplotypecaller import GATK4HaplotypeCaller
from camel.app.tools.gatk4.gatk4selectvariants import GATK4SelectVariants
from camel.app.tools.gatk4.gatk4validatevariants import GATK4ValidateVariants
from camel.app.tools.gatk4.gatk4variantfiltration import GATK4VariantFiltration


class TestGATK4(CamelTestSuite):
    """
    Tests the GATK4 tool suite.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('gatk4')
    FILE_BAM = ToolIOFile(test_file_dir / 'NA12877_chr22_subRG.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'Homo_sapiens_assembly38_chr22.fasta')
    FILE_BQSR = ToolIOFile(test_file_dir / "recal_data.table")
    FILE_gVCF1 = ToolIOFile(test_file_dir / "NA12877_sub1.g.vcf")
    FILE_gVCF2 = ToolIOFile(test_file_dir / "NA12877_sub2.g.vcf")
    FILE_gVCF3 = ToolIOFile(test_file_dir / "NA12877_sub3.g.vcf")
    FILE_VCF = ToolIOFile(test_file_dir / "NA12877_chr22.VCF.gz")
    FILE_KNOWN_SITES = ToolIOFile(test_file_dir / "Homo_sapiens_assembly38.dbsnp138_chr22.vcf.gz")
    FILE_TXT_INTERVALS = ToolIOFile(test_file_dir / "interval_file.intervals")

    def test_gatk4_applybqsr(self) -> None:
        """
        Test GATK4ApplyBQSR
        :return: None
        """
        apply_bqsr = GATK4ApplyBQSR(self.camel)
        apply_bqsr.add_input_files({
            'BAM': [TestGATK4.FILE_BAM],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
            'BQSR': [TestGATK4.FILE_BQSR]
        })
        apply_bqsr.run(self.running_dir)
        self.assertTrue('BAM' in apply_bqsr.tool_outputs, "No BAM output generated")
        output_file = Path(apply_bqsr.tool_outputs['BAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_baserecalibrator(self) -> None:
        """
        Test GATK4BaseRecalibrator
        :return: None
        """
        baserecalibrator = GATK4BaseRecalibrator(self.camel)
        baserecalibrator.add_input_files({
            'BAM': [TestGATK4.FILE_BAM],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
            'VCF_KNOWN_SNPS': [TestGATK4.FILE_KNOWN_SITES]
        })
        baserecalibrator.run(self.running_dir)
        self.assertTrue('TXT_RecalibrationTable' in baserecalibrator.tool_outputs, "No TXT_RecalibrationTable output generated")
        output_file = Path(baserecalibrator.tool_outputs['TXT_RecalibrationTable'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_combinegvcf(self) -> None:
        """
        Test GATK4CombineGVCFs
        :return: None
        """
        combinegvcf = GATK4CombineGVCFs(self.camel)
        combinegvcf.add_input_files({
            'gVCF': [TestGATK4.FILE_gVCF1, TestGATK4.FILE_gVCF2, TestGATK4.FILE_gVCF3],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        combinegvcf.run(self.running_dir)
        self.assertTrue('gVCF' in combinegvcf.tool_outputs, "No gVCF output generated")
        output_file = Path(combinegvcf.tool_outputs['gVCF'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_fastaalternatereferencemaker(self) -> None:
        """
        Test GATK4FastaAlternateReferenceMaker
        :return: None
        """
        fastaalternatereferencemaker = GATK4FastaAlternateReferenceMaker(self.camel)
        fastaalternatereferencemaker.add_input_files({
            'VCF': [TestGATK4.FILE_VCF],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        fastaalternatereferencemaker.update_parameters(**{'concatenate_sequence_segments': 'true'})
        fastaalternatereferencemaker.run(self.running_dir)
        self.assertTrue('FASTA' in fastaalternatereferencemaker.tool_outputs, "No FASTA output generated")
        output_file = Path(fastaalternatereferencemaker.tool_outputs['FASTA'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_gatherbqsrreports(self) -> None:
        """
        Test GATK4GatherBQSRReports
        :return: None
        """
        gatherbqsrreports = GATK4GatherBQSRReports(self.camel)
        gatherbqsrreports.add_input_files({
            'TXT_intervals': [ToolIOFile(TestGATK4.test_file_dir / "1_recal_data.csv"),
                              ToolIOFile(TestGATK4.test_file_dir / "2_recal_data.csv"),
                              ToolIOFile(TestGATK4.test_file_dir / "3_recal_data.csv")],
        })
        gatherbqsrreports.run(self.running_dir)
        self.assertTrue('TXT_RecalibrationTable' in gatherbqsrreports.tool_outputs, "No TXT_RecalibrationTable output generated")
        output_file = Path(gatherbqsrreports.tool_outputs['TXT_RecalibrationTable'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_genotypegvcfs(self) -> None:
        """
        Test GATK4GenotypeGVCFs
        :return: None
        """
        genotypegvcfs = GATK4GenotypeGVCFs(self.camel)
        genotypegvcfs.add_input_files({
            'gVCF': [TestGATK4.FILE_gVCF1],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        genotypegvcfs.run(self.running_dir)
        self.assertTrue('VCF_MultipleSample' in genotypegvcfs.tool_outputs, "No VCF_MultipleSample output generated")
        output_file = Path(genotypegvcfs.tool_outputs['VCF_MultipleSample'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_haplotypecaller(self) -> None:
        """
        Test GATK4HaplotypeCaller
        :return: None
        """
        haplotypecaller = GATK4HaplotypeCaller(self.camel)
        haplotypecaller.add_input_files({
            'BAM': [TestGATK4.FILE_BAM],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        haplotypecaller.run(self.running_dir)
        self.assertTrue('VCF' in haplotypecaller.tool_outputs, "No VCF output generated")
        output_file = Path(haplotypecaller.tool_outputs['VCF'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_selectvariants(self) -> None:
        """
        Test GATK4SelectVariants
        :return: None
        """
        selectvariants = GATK4SelectVariants(self.camel)
        selectvariants.add_input_files({
            'VCF': [TestGATK4.FILE_VCF],
        })
        selectvariants.run(self.running_dir)
        self.assertTrue('VCF' in selectvariants.tool_outputs, "No VCF output generated")
        output_file = Path(selectvariants.tool_outputs['VCF'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_validatevariants(self) -> None:
        """
        Test GATK4ValidateVariants
        :return: None
        """
        validatevariants = GATK4ValidateVariants(self.camel)
        validatevariants.add_input_files({
            'VCF': [TestGATK4.FILE_VCF],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        validatevariants.run(self.running_dir)
        self.assertTrue('TXT_metrics' in validatevariants.tool_outputs, "No TXT_metrics output generated")
        output_file = Path(validatevariants.tool_outputs['TXT_metrics'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_gatk4_variantfiltration(self) -> None:
        """
        Test GATK4VariantFiltration
        :return: None
        """
        variantfiltration = GATK4VariantFiltration(self.camel)
        variantfiltration.add_input_files({
            'VCF': [TestGATK4.FILE_VCF],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        variantfiltration.update_parameters(**{'filter-names': 'HardQC_filter,lowDP',
                                'filter-expressions': 'QD<2.0||FS>60.0||MQ<40.0||ReadPosRankSum<-8.0,DP<200'})
        variantfiltration.run(self.running_dir)
        self.assertTrue('VCF' in variantfiltration.tool_outputs, "No VCF output generated")
        output_file = Path(variantfiltration.tool_outputs['VCF'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

if __name__ == '__main__':
    unittest.main()