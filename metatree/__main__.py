###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import argparse
import logging
import sys
import traceback

from metatree import __description__, __version__
from metatree.common import check_on_path
from metatree.exception import MetaTreeException, MetaTreeExit
from metatree.io import Batchfile
from metatree.io.taxonomy_file import TaxonomyFile
from metatree.logger import logger_setup
from metatree.pipeline import run_pipeline


def print_help():
    lines = [f'metatree v{__version__}',
             'usage: [batchfile] [out_dir] [taxonomy_file] [outgroup] [cpus]']
    print('\n'.join(lines))


def main(args=None):
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('batchfile', type=str,
                        help='First tree must be the reference tree, format: id<tab>path_to_tree')
    parser.add_argument('out_dir', type=str, help='path to the output directory')
    parser.add_argument('taxonomy_file', type=str, help='path to taxonomy file, format: gid<tab>taxonomy')
    parser.add_argument('outgroup', type=str, help='outgroup for rooting')
    parser.add_argument('cpus', type=int, help='number of CPUs to use')

    # Verify that a subparser was selected
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'metatree v{__version__}')
        sys.exit(0)
    else:
        print(f'metatree v{__version__}')
        args = parser.parse_args(args)

        # Setup the logger.
        logger_setup(args.out_dir if hasattr(args, 'out_dir') else None,
                     "metatree.log", "metatree", __version__, False,
                     hasattr(args, 'debug') and args.debug)
        logger = logging.getLogger('timestamp')

        try:
            # Validate the input arguments.
            batchfile = Batchfile(args.batchfile)
            tax_file = TaxonomyFile(args.taxonomy_file)
            if not args.outgroup or args.outgroup[0:3] not in {'d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'}:
                raise MetaTreeExit(f'Invalid outgroup: {args.outgroup}')
            cpus = max(1, args.cpus)

            # Assert that the required programs are on the system path.
            for prog in ('genometreetk', 'phylorank'):
                check_on_path(prog)

            # Run the pipeline.
            run_pipeline(batchfile, args.out_dir, tax_file, args.outgroup, cpus)

        except SystemExit:
            sys.stdout.write('\n')
            sys.stdout.flush()
            logger.error('Controlled exit resulting from early termination.')
            sys.exit(1)
        except KeyboardInterrupt:
            sys.stdout.write('\n')
            sys.stdout.flush()
            logger.error('Controlled exit resulting from interrupt signal.')
            sys.exit(1)
        except MetaTreeExit as e:
            sys.stdout.write('\n')
            sys.stdout.flush()
            if len(str(e)) > 0:
                logger.error('{}'.format(e))
            logger.error('Controlled exit resulting from an unrecoverable error or warning.')
            sys.exit(1)
        except MetaTreeException as e:
            sys.stdout.write('\n')
            sys.stdout.flush()
            msg = 'Controlled exit resulting from an unrecoverable error or warning.\n\n'
            msg += '=' * 80 + '\n'
            msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
            msg += '  MESSAGE: {}\n'.format(e)
            msg += '_' * 80 + '\n\n'
            msg += traceback.format_exc()
            msg += '=' * 80
            logger.error(msg)
            sys.exit(1)
        except Exception as e:
            sys.stdout.write('\n')
            sys.stdout.flush()
            msg = 'Uncontrolled exit resulting from an unexpected error.\n\n'
            msg += '=' * 80 + '\n'
            msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
            msg += '  MESSAGE: {}\n'.format(e)
            msg += '_' * 80 + '\n\n'
            msg += traceback.format_exc()
            msg += '=' * 80
            logger.error(msg)
            sys.exit(1)

    # Done - no errors.
    logger.info('Done.')
    sys.exit(0)


if __name__ == "__main__":
    main()
