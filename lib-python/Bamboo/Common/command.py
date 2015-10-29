#! /usr/local/python/python2.7.3-64bit/bin/python

# Contains Functions for Launching Commands

import os, sys, subprocess, threading, shlex, time

# Code here adapted from code by Sebastian Kelm

class commandManager(object):
    '''
        Enables commands be run on another thread, with various options

        Based on jcollado's solution:
        http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933

        https://gist.github.com/1306188

        Modified by Sebastian Kelm: added timedout parameter
        Modified by Nicholas Pearce: added code inspired by code from the CCP4 dispatcher project
    '''

    def __init__(self, program, cmd_line_args=None, std_inp_lines=None):
        # Name of program
        self.program = shlex.split(program)
        # Command line args and standard input args
        if cmd_line_args: 
            self.cmd_line_args = cmd_line_args
        else:
            self.cmd_line_args = []
        if std_inp_lines:
            self.std_inp_lines = std_inp_lines
        else:
            self.std_inp_lines = []
        # Objects
        self.kwargs = {}
        self.thread = None
        self.process = None
        # Process Control
        self.timeout = None
        self.timedout = False
        # Response
        self.output = ''
        self.error = ''
        # Meta
        self.runtime = None

    def set_timeout(self, timeout):
        """Set various program parameters"""
        self.timeout = timeout

    def add_command_line_arguments(self, *args):
        """Add to the CMD LINE for the program"""
        for arg in args:
            if   isinstance(arg, str):  self.cmd_line_args.extend(shlex.split(arg))
            elif isinstance(arg, list): self.add_command_line_arguments(*arg)

    def add_standard_input(self, args):
        """Add line to the STDIN for the program"""
        for arg in args:
            if   isinstance(arg, str):  self.std_inp_lines.append(arg)
            elif isinstance(arg, list): self.add_standard_input(*arg)

    def print_settings(self):
        """Print out the current settings of the object"""
        print 'Program:\n\t', ' '.join(self.program)
        print 'Command:\n\t', '\n\t'.join(self.cmd_line_args)
        print 'Inputs:\n\t',  '\n\t'.join(self.std_inp_lines)

    def _prepare_inputs_and_outputs(self):
        """Prepare the Input Pipes"""
        if self.std_inp_lines:
            self.__stdin__ = subprocess.PIPE
            self.kwargs['stdin'] = self.__stdin__
        else:
            self.__stdin__ = None
            self.kwargs.pop('stdin',None)
        self.__stdout__ = subprocess.PIPE
        self.kwargs['stdout'] = self.__stdout__
        self.__stderr__ = subprocess.PIPE
        self.kwargs['stderr'] = self.__stderr__

    def run(self):
        """Run command with keyword arguments"""
        self.run_async()
        return self.join_async()

    def run_async(self):
        """Run Command with the Keyword Arguments and then return the thread object"""
        def target():
            self.process = subprocess.Popen(args=self.program+self.cmd_line_args, **self.kwargs)
            self.output, self.error = self.process.communicate('\n'.join(self.std_inp_lines))
        # Prepare pipes
        self._prepare_inputs_and_outputs()
        # Create the new thread
        self.thread = threading.Thread(target=target)
        # Start time
        self.time_start = time.time()
        # Start the Thread
        self.thread.start()

    def join_async(self):
        """Join the thread object returned by self.run_async"""
        # Check for timeout
        self.thread.join(self.timeout)
        # Act accordingly
        if self.thread.is_alive():
            self.process.terminate()
            self.thread.join()
            self.timedout = True
            return 666
        # End time
        self.time_end = time.time()
        # Calculate runtime
        self.runtime = self.time_end-self.time_start
        # Return process exit status
        return self.process.returncode
