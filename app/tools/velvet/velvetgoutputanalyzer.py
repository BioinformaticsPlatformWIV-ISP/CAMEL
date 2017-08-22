import re
import logging

import subprocess


class VelvetgOutputAnalyzer(object):
    """Class provides the functions to analyze the velvetg outputs"""

    @staticmethod
    def analyze_log_file(log_file):
        """
        Check Log file (last line) for velvetg run statistics,
        example content:
        Final graph has 2478 nodes and n50 of 97747, max 461560, total 6485226, using 9880515/11198596 reads
        :param log_file: velvetg Log output file (with complete path)
        :return: informs, dictionary containing update velvetg run information
        """
        informs = {}
        info_str = subprocess.check_output(['tail', '-1', log_file])
        res = re.match(
            r'Final graph has (\d+) nodes and n50 of (\d+), max (\d+), total (\d+), using (\d+)/(\d+) reads', info_str)
        if res:
            res_data = res.groups()
            informs['n_nodes'] = res_data[0]
            informs['n50'] = res_data[1]
            informs['max_contig_len'] = res_data[2]
            informs['total_base_in_contigs'] = res_data[3]
            informs['reads_used'] = res_data[4]
            informs['reads_total'] = res_data[5]
        else:
            logging.warning(
                'Fail to extra the basic assembly statistics of velvetg from Log file {}.'.format(log_file))

        return informs

    @staticmethod
    def analyze_stats_file(input_informs, stats_file):
        """
        Analyze stats.txt for contigs stats (idea from VelvetOptimiser): 'number of contigs longer than 1k', 'total bases in
        contigs > 1k'. NOTE that the base count does not consider redundancy (bases mapped to the same position on ref
        genome)!

        :param input_informs: dictionary with the current velvetg run information
        :param stats_file: velvetg stats output file (with complete path)
        :return: informs (dictionary)
        """
        if 'kmer' not in input_informs:
            logging.warning("Require kmer information is missing, Velvetg stats file analysis abort!")
            return {}
        if 'cov_cutoff' not in input_informs:
            logging.warning("Require coverage cutoff information is missing, Velvetg stats file analysis abort!")
            return {}

        ctg_count = 0
        base_count = 0
        output_informs = {}
        with open(stats_file) as f:
            f.readline()             # skip the header

            for l in f.readlines():
                info = l.split("\t")
                # NOTE:
                # - coverage cutoff of velvetg are applied on cov, instead of stricter 0cov
                # - the reported length is kmer-length (refer velvet manual)
                actual_len = int(info[1]) + input_informs['kmer'] - 1
                if float(info[5]) < float(input_informs['cov_cutoff']) or actual_len < 2 * input_informs['kmer']:
                    # cov_cutoff and 2k length filtering
                    continue
                if actual_len >= 1000:
                    ctg_count += 1
                    base_count += actual_len

        output_informs['contigs_>_1k'] = ctg_count
        output_informs['total_bases_in_contigs_>_1k'] = base_count

        return output_informs

    @staticmethod
    def analyze_velvetg_output(informs, output):
        """
        Parse command line output of velvetg for 'Paired Library insert estimation' and 'Coverage stats'
        :param informs: dictionary to store velvetg run information
        :param output: str, velvetg run command line output
        :return: informs, dictionary containing updated velvetg run information

        example output lines:
        - PE insert:
             [0.074875] Paired-end library 1 has length: 442, sample standard deviation: 29
             [0.145957] Paired-end library 1 has length: 443, sample standard deviation: 50
          NOTE: Only the 2nd one (after Scaffolding) is kept!!

        - Coverage
             [0.202295] Estimated Coverage = 11.694646
        """
        for l in output.split("\n"):
            # insert length (max 4 PE libraries)
            res = re.search('Paired-end library (\d) has length: (\d+),', l)
            if res:
                res_match = res.groups()
                informs['ins_length_PE{}'.format(res_match[0])] = res_match[1]

            # coverage
            res = re.search('Estimated Coverage = (.+)', l.strip())
            if res:
                informs['coverage'] = res.groups()[0]
            res = re.search('Estimated Coverage cutoff = (.+)', l.strip())
            if res:
                informs['cov_cutoff'] = res.groups()[0]

        return informs
