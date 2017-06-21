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

    def run_command(self, folder):
        """
        Runs the command given at command initialization
        :param folder: Folder where the command is executed
        :return: None
        """
        logging.info('Executing command: {}'.format(self.command))
        if self.command is None:
            raise ValueError("Invalid command 'None'")
        self._procedure = subprocess.run(
            self._command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            executable='/bin/bash',
            cwd=folder)
        self._stdout = self._procedure.stdout.decode('utf-8')
        self._stderr = self._procedure.stderr.decode('utf-8')
        self._return_code = self._procedure.returncode
        logging.debug('stdout: {}'.format(self._stdout))
        logging.debug('stderr: {}'.format(self._stderr))
