import logging
import subprocess


class Command(object):
    """
    Class meant to handle the execution of commands
    """

    def __init__(self, command=None):
        self._stdout = None
        self._stderr = None
        self._procedure = None
        self._return_code = None
        self._command = command

    @property
    def stderr(self):
        return self._stderr

    @property
    def stdout(self):
        return self._stdout

    @property
    def returncode(self):
        return self._return_code

    @property
    def command(self):
        return self._command

    @command.setter
    def command(self, cmd):
        self._command = cmd

    def run_command(self, folder, stderr_handle=subprocess.PIPE):
        """
        Runs the command given at command initialization
        :param folder: Folder where the command is executed
        :param stderr_handle: Handle for the standard error (e.g. PIPE or STDOUT)
        :return: None
        """
        logging.info('Executing command: {}'.format(self.command))
        if self.command is None:
            raise ValueError("Invalid command 'None'")
        self._procedure = subprocess.run(
            self._command,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            shell=True,
            executable='/bin/bash',
            cwd=folder)
        self._stdout = self._procedure.stdout.decode('utf-8')
        if self._procedure.stderr is not None:
            self._stderr = self._procedure.stderr.decode('utf-8')
        else:
            self._stderr = ''
        self._return_code = self._procedure.returncode
        logging.debug('stdout: {}'.format(self._stdout))
        logging.debug('stderr: {}'.format(self._stderr))
