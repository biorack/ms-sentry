#primary script to run

import os
import time
import pandas as pd
import itertools
import argparse
import plots
from dataset import RawDataset
from raw_file_validation import get_raw_file_age
from utils import print_progress_bar
from analysis_routines import initialize_data_storage, close_data_storage, read_data_storage, analyze_experiment, export_tables_plots

pos_diagnostic_ions = [111.0808, 129.1045, 140.0793, 176.1135]
pos_ms2_diagnostic = {'pmz':176.1135, 'diagnostic_ions':pos_diagnostic_ions}

neg_diagnostic_ions = [75.01253, 156.07551, 174.09893]
neg_ms2_diagnostic = {'pmz':174.09893, 'diagnostic_ions':neg_diagnostic_ions}

#ms2 diagnostic ions used to plot ms2 data
ms2_diagnostic = {'POS':pos_ms2_diagnostic, 'NEG':neg_ms2_diagnostic}

#current directory of program
current_dir = os.path.dirname(__file__)

class PathExists(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        if not os.path.exists(values):
            msg='experimental directory does not exist'
            raise argparse.ArgumentTypeError(msg)
        setattr(args, self.dest, values)


def arg_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(description='MS-Sentry generates plots and tables useful for performing quality control on Thermo Orbitrap data.')

    parser.add_argument('directory', action=PathExists, help='experimental directory containing Thermo raw files')

    #Analysis options
    analysis_options = parser.add_argument_group()
    analysis_options.add_argument('-num', '--num_files', type=int, action='store', required=False,
                                    help='number of files to analyze in dataset. default is all files currently in directory')
    analysis_options.add_argument('-min','--min_file_age', type=int, action='store', default=5, required=False,
                                    help='minimum file age in minutes. default is 30')
    analysis_options.add_argument('-skip', '--skip_blanks', type=bool, action='store', default=True, required=False,
                                    help='skip blank injections. default is True')
    analysis_options.add_argument('-export', '--export_num', type=int, action='store', default=None, required=False,
                                    help='export partial analysis every n number of files. default is None')

    return parser

def main(args):

    initialize_data_storage(current_dir)
    dataset = RawDataset(args.directory)
    exported_files = 0

    file_name_warnings, file_name_errors = dataset.validate_experiment_filenames()
    file_name_warnings_df = pd.DataFrame(file_name_warnings)
    file_name_errors_df = pd.DataFrame(file_name_errors)
    all_errors = list(itertools.chain.from_iterable(file_name_errors_df['errors']))

    if len(all_errors) > 0:
        print('Filename errors detected!')
        print('Exporting reports...')

        file_name_errors_df.to_csv(os.path.join(args.directory, 'file_name_errors_report.csv'))
        file_name_warnings_df.to_csv(os.path.join(args.directory, 'file_name_warnings_report.csv'))

        close_data_storage(current_dir)
        quit()

    atlas_df = pd.read_csv(os.path.join(current_dir, 'default_atlases\default_{chrom}_atlas.csv').format(chrom=dataset.chromatography.lower()))

    if args.num_files is None:
        args.num_files = len(dataset.files_to_analyze)
        
    print_progress_bar(0, args.num_files, prefix = 'Progress:', suffix = 'Complete', length = 50)

    analysis_complete = False
    while True:

        if args.export_num is not None:

            new_exported_files = len(dataset.analyzed_files) - exported_files

            if new_exported_files == args.export_num:
                print('\nExporting and continuing analysis...')
                exported_files += new_exported_files
                export_tables_plots(current_dir, args, file_name_warnings_df, atlas_df, ms2_diagnostic, exported_files)

        if len(dataset.analyzed_files) >= args.num_files:
            analysis_complete = True

        if not analysis_complete:

            if len(dataset.files_to_analyze) < 1:
                time.sleep(15)
                dataset.set_files_to_analyze()
                continue

            file_ages = [get_raw_file_age(file) for file in dataset.files_to_analyze]
            if max(file_ages) < args.min_file_age:
                wait_time_seconds = (args.min_file_age - max(file_ages)) * 60
                time.sleep(wait_time_seconds+200)

            analaysis_interrupted = analyze_experiment(dataset, atlas_df, ms2_diagnostic, current_dir, args, exported_files, tolerance=0.0015, frag_tolerance=0.005)

            if analaysis_interrupted:
                print('Analysis interrupted. Exiting...')
                close_data_storage(current_dir)
                break

        else:

            print('\nDone! Exporting Results...')

            export_tables_plots(current_dir, args, file_name_warnings_df, atlas_df, ms2_diagnostic, 'full')
            close_data_storage(current_dir)
            break

if __name__ == '__main__':
    parser = arg_parser()
    args = parser.parse_args()

    main(args)








