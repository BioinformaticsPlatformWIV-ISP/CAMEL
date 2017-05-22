import vcf


def count_variants(vcf_file):
    """
    Counts the number of variants in a VCF file.
    :param vcf_file: VCF file.
    :return: Number of variants
    """
    vcf_reader = vcf.Reader(open(vcf_file))
    variants = list(vcf_reader)
    return sum(variant.ALT != [None] for variant in variants)
