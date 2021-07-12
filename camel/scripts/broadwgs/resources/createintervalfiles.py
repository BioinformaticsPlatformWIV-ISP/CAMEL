Path(output.sequence_grouping).mkdir(exist_ok=True)

        dict_genome = SnakemakeUtils.load_object(input.DICT_GENOME)[0]

        with open(f"{dict_genome}", "r") as ref_dict_file:
            sequence_tuple_list = []
            longest_sequence = 0

            for line in ref_dict_file:
                if line.startswith("@SQ"):
                    line_split = line.split("\t")
                    # (Sequence_Name, Sequence_Length)
                    sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

        # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
        # the last element after a :, so we add this as a sacrificial element.
        hg38_protection_tag = ":1+"

        # initialize the tsv string with the first sequence
        tsv_list = [(sequence_tuple_list[0][0] + hg38_protection_tag)]
        temp_size = sequence_tuple_list[0][1]
        tsv_print_list = []

        for sequence_tuple in sequence_tuple_list[1:]:
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                tsv_list.append((sequence_tuple[0] + hg38_protection_tag))
            else:
                tsv_print_list.append(tsv_list)
                tsv_list = [(sequence_tuple[0] + hg38_protection_tag)]
                temp_size = sequence_tuple[1]
            tsv_list.append((sequence_tuple[0] + hg38_protection_tag))

        tsv_print_list.append(tsv_list)

        # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
        tsv_print_list.append(["unmapped"])

        for (n,interval) in enumerate(tsv_print_list):
            with open(Path(output.sequence_grouping) / ("interval_" + str(n) + ".intervals"), "w") as interval_file:
                interval_file.write("\n".join(interval))
                interval_file.close()

            interval_out = Path(output.sequence_grouping) / ("interval_" + str(n) + ".intervals")
            interval_io =  Path(output.sequence_grouping) / ("interval_" + str(n) + ".intervals.io")
            SnakemakeUtils.dump_object([ToolIOFile(interval_out)], str(interval_io))