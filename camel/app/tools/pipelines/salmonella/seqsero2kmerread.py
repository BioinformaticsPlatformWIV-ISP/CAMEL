from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2


class SeqSero2KmerRead(SeqSero2):

    def build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        self._informs['_tag'] = 'Kmerread'
        if 'FASTQ' in self._tool_inputs:
            self._command.command = ' '.join([
                self._tool_command, '-t 3 -m k -d {}'.format(self.folder), '-i',
                str(self._tool_inputs['FASTQ'][0].path), " ".join(self._build_options())
            ])
        else:
            self._command.command = ' '.join([
                self._tool_command, '-t 2 -m k -d {}'.format(self.folder), '-i',
                str(self._tool_inputs['FASTQ_PE'][0].path), str(self._tool_inputs['FASTQ_PE'][1].path),
                " ".join(self._build_options())
            ])
