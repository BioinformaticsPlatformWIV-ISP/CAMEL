import vcf
from app.error.invalidparametererror import InvalidParameterError


class VCFUtils(object):

    """
    Helper to perform VCF file related functions
    """

    INDEL = 'indel'
    INDEL_INS = 'ins'
    INDEL_DEL = 'del'
    SNP = 'snp'
    SV = 'sv'

    @staticmethod
    def is_multi_sample(vcf_file):
        """
        Function to check whether a VCF file contains multiple sample
        :param vcf_file: the vcf file to be checked (with complete path)
        :return: True if is a multiple sample VCF file
        """
        vcf_reader = vcf.Reader(filename=vcf_file)
        return len(vcf_reader.samples) > 1

    @staticmethod
    def get_reader(vcf_file):
        """
        Function to obtain a vcf file reader
        :param vcf_file: the vcf file to be checked (with complete path)
        :return: vcf_reader that can called to retrieve the record
        """
        return vcf.Reader(filename=vcf_file)

    @staticmethod
    def retrieve_variants(vcf_file, types=None, excluded_types=None):
        """
        Function to retrieve certain types of variants. Parameter
        'types' and 'excluded_types' is mutual exclusive. Either specify
        'types' to retrieved only certain variants or specify
        'excluded_types' to exclude variants.

        KNOWN TYPES: 'tv', 'unknown', 'ts', 'ins', 'del', 'indel', 'snp

        :param vcf_file: [optional] the vcf file to retrieve data
        :param types: [optional] types of variants to be retrieved
        :param excluded_types: types of variants to be excluded
        :return: list of records of types
        """
        types = [] if types is None else types
        excluded_types = [] if excluded_types is None else excluded_types

        vcf_reader = vcf.Reader(filename=vcf_file)
        records = []
        if len(types) > 0 and len(excluded_types) > 0:
            raise InvalidParameterError(
                "Mutually exclusive parameters 'included types' and 'excluded types' are specified. Only one is allowed.")
        elif len(types) > 0:
            # only set types
            for rcd in vcf_reader:
                if any(x in types for x in [rcd.var_type, rcd.var_subtype]):
                    records.append(rcd)
        elif len(excluded_types) > 0:
            # only set excluded_types
            for rcd in vcf_reader:
                if all(x not in excluded_types for x in [rcd.var_type, rcd.var_subtype]):
                    records.append(rcd)
        else:
            # no conditions
            list(vcf_reader)

        return records

    @staticmethod
    def count_variants(vcf_file):
        """
        Counts the number of variants in a VCF file.
        :param vcf_file: VCF file.
        :return: Number of variants
        """
        vcf_reader = vcf.Reader(filename=vcf_file)
        variants = list(vcf_reader)
        return sum(variant.ALT != [None] for variant in variants)
