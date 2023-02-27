#!/usr/bin/env python
import argparse
import ast
import logging
import re
from pathlib import Path
from typing import Sequence, Optional, List, Dict, Any

import humanize
import pandas as pd

from camel.app.camel import Camel


class MainPipelineCombine(object):
    """
    Main script for the pipeline combine tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainPipelineCombine._parse_arguments(args)

    def run(self) -> None:
        """
        Runs this tool.
        """
        # Parse input data
        records_out = []
        for path_input in self._args.input:
            isolate_data = MainPipelineCombine._parse_pipe_out(path_input)
            records_out.append(isolate_data)
        data_out = pd.DataFrame(records_out)
        logging.info(f'Summaries parsed for {len(data_out)} datasets')

        # Check if the pipelines are the same
        if len(data_out['pipeline_name'].unique()) > 1:
            logging.warning('Combining output of multiple pipelines (not recommended!)')

        # Filter columns
        cols_to_keep = data_out.apply(lambda x: MainPipelineCombine._filter_columns(
            x, self._args.skip.split(',')), axis=0)
        data_out_filt = data_out.loc[:, cols_to_keep]
        logging.info(f'Keeping {len(data_out_filt.columns)}/{len(data_out.columns)} columns')

        # Save output file
        data_out_filt.to_csv(self._args.output, sep='\t', index=False)
        logging.info(
            f'Output file created: {self._args.output} ({humanize.naturalsize(self._args.output.stat().st_size)})')

    @staticmethod
    def _parse_pipe_out(path_in: Path) -> Dict[str, Any]:
        """
        Parses a pipeline summary output file.
        :param path_in: Input path
        :return: Parsed info as a dictionary
        """
        isolate_data = {}
        with path_in.open() as handle:
            for line in handle.readlines():
                parts = line.strip().split('\t')
                key = parts[0]
                value = '\t'.join(parts[1:])

                # Gene detection hits
                m = re.match(r'hits_(.*)', key)
                print(key, m)
                if m:
                    hits = ast.literal_eval(value)
                    if len(hits) == 0:
                        continue
                    for hit in hits:
                        isolate_data[f'{m.group(1)}_{hit[0]}'] = f'{hit[1]} (%id={hit[2]}, len={hit[3]})'
                else:
                    isolate_data[key] = value
        return isolate_data

    @staticmethod
    def _filter_columns(column: pd.Index, keys_to_skip: List[str]) -> bool:
        """
        Function to check if the target row should be filtered.
        :param column: Columns to check
        :param keys_to_skip: True if key should be retained, False otherwise
        """
        for key in keys_to_skip:
            # Regular match
            if '*' not in key:
                if column.name == key:
                    return False

            # Regex match
            if '*' in key:
                key = key.replace('*', '.*')
                m = re.match(key, str(column.name))
                if m:
                    return False
        return True

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('input', metavar='N', type=Path, nargs='+', help='Input files')
        parser.add_argument('--output', type=Path, help='Output path')
        parser.add_argument('--skip', type=str, help='Comma separated list of keys to skip')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainPipelineCombine()
    main.run()
