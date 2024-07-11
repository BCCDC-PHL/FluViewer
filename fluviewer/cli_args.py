import argparse
import os
from pathlib import Path

from . import __version__

def parse_args():
    """
    Parse command line arguments using argparse.

    :return: Dictionary of arguments and their values.
    :rtype: dict
    """
    parser = argparse.ArgumentParser(description='BCCDC-PHL/FluViewer: Influenza A virus consensus sequence generation and variant calling')
    parser.add_argument('-f', '--forward-reads', type=Path, required=True, help='Path to FASTQ file containing forward reads')
    parser.add_argument('-r', '--reverse-reads', type=Path, required=True, help='Path to FASTQ file containing reverse reads')
    parser.add_argument('-d', '--db', type=Path, required=True, help='Path to FASTA file containing FluViewer database')
    parser.add_argument('-o', '--outdir', type=Path, help='Output directory (default=FluViewer_<output-name>)')
    parser.add_argument('-n', '--output-name', type=str, required=True, help='Output name. Includes this name in output files, and in consensus sequence headers')
    parser.add_argument('-i', '--min-identity', type=float, default=90, metavar="[0-100]", help='Minimum percent sequence identity between database reference sequences and contigs (default=90)')
    parser.add_argument('-l', '--min-alignment-length', type=int, default=50, metavar="[32-]", help='Minimum length of alignment between database reference sequences and contigs (default=50)')
    parser.add_argument('-D', '--min-depth', type=int, default=20, metavar="[1-]", help='Minimum read depth for base calling (default=20)')
    parser.add_argument('-q', '--min-mapping-quality', type=int, default=20, metavar="[0-]", help='Minimum PHRED score for mapping quality and base quality during variant calling (default=20)')
    parser.add_argument('-v', '--variant-threshold-calling', type=float, default=0.75, metavar="[0-1]", help='Variant allele fraction threshold for calling variants (default=0.75)')
    parser.add_argument('-V', '--variant-threshold-masking', type=float, default=0.25, metavar="[0-1]", help='Variant allele fraction threshold for masking ambiguous variants (default=0.25)')
    parser.add_argument('-N', '--target-depth', type=int, default=200, metavar="[1-]", help='Target depth for pre-normalization of reads (default=200)')
    parser.add_argument('-L', '--coverage-limit', type=int, default=200, metavar="[1-]", help='Coverage depth limit for variant calling (default=200)')
    parser.add_argument('-t', '--threads', type=int, default=1, metavar="[1-]", help='Threads used for contig/scaffold alignments (default=1)')
    parser.add_argument('-M', '--max-memory', type=int, metavar="[1-]", help='Gigabytes of memory allocated for normalizing reads (default=max)')
    parser.add_argument('-g', '--disable-garbage-collection', action='store_true', help='Disable garbage collection and retain intermediate analysis files')
    parser.add_argument('--skip-depth-normalization', action='store_true', help='Skip read depth normalization (bbnorm) step')
    parser.add_argument('--force', action='store_true', help='Allow overwrite of existing files and directories.')
    parser.add_argument('--log-level', default='info', choices=['info', 'debug'], help='Log level (default=info)')
    parser.add_argument('--version', action='version', version=__version__)

    args = parser.parse_args()

    return args



def validate_args(args):
    """
    Validate the arguments provided by the user.

    :param args: Dictionary of arguments and their values.
    :type args: argparse.Namespace
    """
    independent_validation_rules = {
        'forward_reads': { 'validation_fn': lambda x: os.path.isfile(x),
                           'error_msg': 'Input file does not exist: {0}' },
        'reverse_reads': { 'validation_fn': lambda x: os.path.isfile(x),
                           'error_msg': 'Input file does not exist: {0}' },
        'db':            { 'validation_fn': lambda x: os.path.isfile(x),
                           'error_msg': 'Input file does not exist: {0}' },
        'outdir':        { 'validation_fn': lambda x: True,
                           'error_msg': None },
        'min_identity':         { 'validation_fn': lambda x: 0 <= x <= 100,
                                  'error_msg': 'Minimum percent sequence identity must be between 0 and 100: {0}' },
        'min_alignment_length': { 'validation_fn': lambda x: x >= 32,
                                  'error_msg': 'Minimum alignment length must be at least 32: {0}' },
        'min_depth':            { 'validation_fn': lambda x: x >= 1,
                                  'error_msg': 'Minimum read depth must be at least 1: {0}' },
        'min_mapping_quality':  { 'validation_fn': lambda x: x >= 0,
                                  'error_msg': 'Minimum PHRED score must be at least 0: {0}' },
        'variant_threshold_calling': { 'validation_fn': lambda x: 0 <= x <= 1,
                                       'error_msg': 'Variant allele fraction threshold for calling variants must be between 0 and 1: {0}' },
        'variant_threshold_masking': { 'validation_fn': lambda x: 0 <= x <= 1,
                                       'error_msg': 'Variant allele fraction threshold for masking ambiguous variants must be between 0 and 1: {0}' },
        'target_depth':   { 'validation_fn': lambda x: x >= 1,
                            'error_msg': 'Target depth for pre-normalization of reads must be at least 1: {0}' },
        'coverage_limit': { 'validation_fn': lambda x: x >= 1,
                            'error_msg': 'Coverage depth limit for variant calling must be at least 1: {0}' },
        'threads':        { 'validation_fn': lambda x: x >= 1,
                            'error_msg': 'Threads used for contig/scaffold alignments must be at least 1: {0}' },
        'max_memory':     { 'validation_fn': lambda x: x is None or x >= 1,
                            'error_msg': 'Gigabytes of memory allocated for normalizing reads must be at least 1: {0}' }
    }
    
    for arg_name, rule in independent_validation_rules.items():
        validation_fn = rule['validation_fn']
        arg_value = getattr(args, arg_name)
        valid = validation_fn(arg_value)
        if not valid:
            raise ValueError(rule['error_msg'].format(arg_value))

    interdependent_validation_rules = {
        ('variant_threshold_masking', 'variant_threshold_calling'): {
            'validation_fn': lambda x, y: x < y,
            'error_msg': 'Variant allele fraction for masking ambiguous variants must be lower than variant allele fraction for calling variants. variant_threshold_masking: {0}, variant_threshold_calling: {1}'
        }
    }

    for arg_pair, rule in interdependent_validation_rules.items():
        arg_1_value = getattr(args, arg_pair[0])
        arg_2_value = getattr(args, arg_pair[1])
        validation_fn = rule['validation_fn']
        valid = validation_fn(arg_1_value, arg_2_value)
        if not valid:
            raise ValueError(rule['error_msg'].format(arg_1_value, arg_2_value))
        

    return args
