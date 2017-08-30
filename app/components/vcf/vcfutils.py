import vcf


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
        Function to check whether a VCF file contains mulitple sample
        :param vcf_file: the vcf file to be checked (with complete path)
        :return: True if is a multiple sampel VCF file
        """
        vcf_reader = vcf.Reader(filename=vcf_file)
        if len(vcf_reader.samples) > 1:
            return True
        else:
            return False

    @staticmethod
    def get_reader(vcf_file):
        """
        Function to read in vcf files
        :param vcf_file: the vcf file to be checked (with complete path)
        :return: vcf_reader that can called to retrieve the record
        """
        return vcf.Reader(filename=vcf_file)

    @staticmethod
    def retrieve_variants(vcf_file, types=[], excluded_types=[]):
        """
        Function to retrieve certain types of variants
        :param vcf_file: the vcf file to retrieve data
        :param types: types of variants to be retrieved
        :param excluded_types: types of variants to be excluded
        :return: list of records of types

        KNOWN TYPES: 'tv', 'unknown', 'ts', 'ins', 'del', 'indel', 'snp'
        """
        vcf_reader = vcf.Reader(filename=vcf_file)
        records = []
        if types and excluded_types:
            # both types and excluded_types set
            for re in vcf_reader:
                if any(x in types for x in [re.var_type, re.var_subtype]) and \
                   all(x not in excluded_types for x in [re.var_type, re.var_subtype]):
                    records.append(re)
        elif types:
            # only set types
            for re in vcf_reader:
                if any(x in types for x in [re.var_type, re.var_subtype]):
                    records.append(re)
        elif excluded_types:
            # only set excluded_types
            for re in vcf_reader:
                if all(x not in excluded_types for x in [re.var_type, re.var_subtype]):
                    records.append(re)
        else:
            # no conditions
            for re in vcf_reader:
                records.append(re)

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
