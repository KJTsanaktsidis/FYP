"""Main entry point for revised snazzy method
"""

from argparse import ArgumentParser
from datastore import InputDatastore
from diffsimtask import DiffSimTask
from outstore import OutputDatastore
from singlecomparison import SingleComparisonTask


if __name__ == '__main__':

    aparser = ArgumentParser(description='Main entry point for diffsim simulator')

    aparser.add_argument('--idir', type=str, required=True,
                         help='Directory containing input data files')
    aparser.add_argument('--iprefix', type=str, default='NiCu',
                         help='Prefix for input data files')
    aparser.add_argument('--odir', type=str, default='.',
                         help='Directory to write output files to')
    aparser.add_argument('--oprefix', type=str, default='',
                         help='Prefix of output file names')
    aparser.add_argument('--resume', type=bool, default=False,
                         help='Resume executing from saved data')

    subparsers = aparser.add_subparsers(dest='subparser_invoked')

    SingleComparisonTask.add_arg_options(subparsers)

    args = aparser.parse_args()

    dstore = InputDatastore(args.idir, args.iprefix)
    ostore = OutputDatastore(args.odir, args.oprefix)

    task = None
    if args.subparser_invoked == 'singlecomp':
        task = SingleComparisonTask(dstore, ostore, args)

    assert isinstance(task, DiffSimTask)

    if not args.resume:
        task.do_stage1()
    task.do_stage2()
