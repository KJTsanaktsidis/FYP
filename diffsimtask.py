"""Abstract base class for tasks that can be invoked via the command line interface of diffsim
"""

from abc import ABCMeta, abstractmethod


class DiffSimTask(metaclass=ABCMeta):

    @staticmethod
    @abstractmethod
    def add_arg_options(argparser):
        """This method is overridden for a task to specify some arguments it needs on the argparse module
        """
        pass

    @abstractmethod
    def do_stage1(self):
        """
        This method is the first half of a task; it should do everything that shouldn't be done if --resume
        is not passed into the method
        """
        pass

    @abstractmethod
    def do_stage2(self):
        """
        This method is what should be executed again if --resume is passed as an argument
        """
        pass
