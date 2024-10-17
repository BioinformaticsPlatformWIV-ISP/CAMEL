#!/usr/bin/env python
import argparse
import shutil
from pathlib import Path
from typing import Optional, Sequence, Dict, Any

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.medaka.medakainference import MedakaInference
from camel.app.tools.medaka.medakavcf import MedakaVcf
from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex


class MainCallingMedaka(object):
    """
    This class contains the main script for the Medaka variant calling tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCallingMedaka._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--bam', type=Path, help="Input BAM file", required=True)
        argument_parser.add_argument('--reference', help="The reference fasta file to use", type=Path, required=True)
        argument_parser.add_argument('--reference-name')
        argument_parser.add_argument('--output-dir', type=Path, help='Output directory')
        argument_parser.add_argument('--working-dir', help='Working directory', type=Path, default=Path.cwd())
        argument_parser.add_argument('--threads', type=int, default=4)

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        input_dict = self.__prepare_input()

        # Index BAM file
        logger.info('Using samtools to index input BAM file')
        samtools_index = SamtoolsIndex(Camel.get_instance())
        samtools_index.add_input_files(input_dict)
        samtools_index.run(self._args.working_dir)

        # Run Medaka Inference
        medaka_inference = MedakaInference(Camel.get_instance())
        logger.info(f'Running {medaka_inference.name} to get hdf file')
        medaka_inference.add_input_files(input_dict)
        medaka_inference.update_parameters(threads=self._args.threads)
        hdf_file = Path(self._args.bam).stem + '.hdf'
        medaka_inference.update_parameters(output=hdf_file)
        medaka_inference.run(self._args.working_dir)

        # Run Medaka VCF
        medaka_vcf = MedakaVcf(Camel.get_instance())
        logger.info(f'Running {medaka_vcf.name} to call the variants based on reference fasta and hdf file')
        medaka_vcf.add_input_files({'HDF': [ToolIOFile(Path(hdf_file))], 'FASTA': input_dict['REFERENCE']})
        vcf_file = Path(hdf_file).stem + '.vcf'
        medaka_vcf.update_parameters(output=vcf_file)
        medaka_vcf.run(self._args.working_dir)

        # copy vcf file to output-dir if set
        if self._args.output_dir is not None:
            self._args.output_dir.mkdir(parents=True, exist_ok=True)
            logger.info('Copying output VCF file to output-dir')
            shutil.copyfile(vcf_file, self._args.output_dir / vcf_file)

    def __prepare_input(self) -> Dict[str, Any]:
        """
        Prepares the input for the Medaka variant tool.
        :return: Input dictionary
        """
        self._args.working_dir.mkdir(parents=True, exist_ok=True)
        input_dict = {'BAM': None, 'REFERENCE': None}

        # input BAM file
        bam_file = Path(self._args.bam)
        path_bam = self._args.working_dir / FileSystemHelper.make_valid(bam_file.name)
        path_bam.symlink_to(bam_file)
        input_dict['BAM'] = [ToolIOFile(path_bam)]

        # Reference genome
        reference_file = Path(self._args.reference)
        path_ref = self._args.working_dir / FileSystemHelper.make_valid(reference_file.name)
        path_ref.symlink_to(reference_file)
        input_dict['REFERENCE'] = [ToolIOFile(path_ref)]

        return input_dict


if __name__ == '__main__':
    Camel.get_instance()
    main = MainCallingMedaka()
    main.run()
