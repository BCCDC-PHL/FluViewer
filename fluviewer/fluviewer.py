import argparse
import json
import logging
import os
import shutil
import subprocess
import sys

from collections import Counter
from math import ceil
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__

from . import cli_args
from . import parsers
from . import plots
import fluviewer.logging

log = fluviewer.logging.get_logger(__name__, 'info')


def main():
    version = __version__

    args = cli_args.parse_args()
    try:
        args = cli_args.validate_args(args)
    except ValueError as e:
        log.error(e)
        exit(1)

    try:
        level = getattr(logging, args.log_level.upper())
        log.setLevel(level)
    except AttributeError:
        log.error(f"Invalid log level: {level}")
        log.setLevel(logging.INFO)
    

    log.info(f'BCCDC-PHL/FluViewer v{version}')
    version_split = version.split('-')
    log.info(f'Derived from: KevinKuchinski/FluViewer v{version_split[0]}')
    log.info(f'Inputs:')
    log.info(f"Fwd reads: {args.forward_reads}")
    log.info(f"Rev reads: {args.reverse_reads}")
    log.info(f"Reference sequences: {args.database}")

    log.info(f"Outputs:")
    log.info(f"Output directory: {args.outdir}")
    log.info(f"Output name: {args.output_name}")

    log.info(f'Parameters:')
    log.info(f"Minimum percent identity: {args.min_identity}")
    log.info(f"Minimum alignment length: {args.min_alignment_length}")
    log.info(f"Minimum read depth: {args.min_depth}")
    log.info(f"Minimum mapping quality: {args.min_mapping_quality}")
    log.info(f"Variant allele fraction threshold for calling variants: {args.variant_threshold_calling}")
    log.info(f"Variant allele fraction threshold for masking ambiguous variants: {args.variant_threshold_masking}")
    log.info(f"Target depth for pre-normalization of reads: {args.target_depth}")
    log.info(f"Coverage depth limit for variant calling: {args.coverage_limit}")
    

    prepare_outdir(args.outdir, args.output_name)

    check_database(
        args.database,
        args.outdir,
        args.output_name,
        args.disable_garbage_collection
    )    

    analysis_stages = [
        'normalize_depth',
        'assemble_contigs',
        'blast_contigs',
        'scaffolding',
        'blast_scaffolds',
        'read_mapping',
        'variant_calling',
        'consensus_calling',
        'summary_reporting',
    ]
        
        
    log.info('Starting analysis...')

    current_analysis_stage = 'normalize_depth'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    normalize_depth(
        args.outdir,
        args.output_name,
        args.forward_reads,
        args.reverse_reads,
        args.target_depth,
        args.max_memory,
        args.disable_garbage_collection
    )

    log.info(f'Analysis stage complete: {current_analysis_stage}')


    current_analysis_stage = 'assemble_contigs'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    assemble_contigs(
        args.outdir,
        args.output_name,
        args.disable_garbage_collection
    )

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    current_analysis_stage = 'blast_contigs'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    blast_results = blast_contigs(
        args.database,
        args.outdir,
        args.output_name,
        args.disable_garbage_collection,
        args.threads,
        args.min_identity,
        args.min_alignment_length,
    )

    filtered_blast_results = filter_contig_blast_results(
        blast_results,
        args.outdir,
        args.output_name,
        args.disable_garbage_collection
    )

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    current_analysis_stage = 'scaffolding'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    make_scaffold_seqs(
        filtered_blast_results,
        args.outdir,
        args.output_name,
        args.disable_garbage_collection
    )

    scaffold_blast_results = blast_scaffolds(
        args.database,
        args.outdir,
        args.output_name,
        args.disable_garbage_collection,
        args.threads,
    )

    filtered_scaffold_blast_results = filter_scaffold_blast_results(scaffold_blast_results)
    
    log.info(f'Analysis stage complete: {current_analysis_stage}')

    current_analysis_stage = 'read_mapping'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    make_mapping_refs(
        filtered_scaffold_blast_results,
        args.database,
        args.outdir,
        args.output_name,
    )

    map_reads(
        args.outdir,
        args.output_name,
        args.disable_garbage_collection,
        args.min_mapping_quality,
    )

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    current_analysis_stage = 'variant_calling'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')
    
    call_variants(
        args.outdir,
        args.output_name,
        args.min_mapping_quality,
        args.coverage_limit,
        args.disable_garbage_collection,
    )

    mask_ambig_low_cov(
        args.outdir,
        args.output_name,
        args.min_depth,
        args.variant_threshold_calling,
        args.variant_threshold_masking,
        args.min_mapping_quality,
        args.disable_garbage_collection,
    )

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    current_analysis_stage = 'consensus_calling'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    make_consensus_seqs(
        args.outdir,
        args.output_name,
        args.disable_garbage_collection,
    )

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    current_analysis_stage = 'summary_reporting'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')
    
    write_report(
        args.outdir,
        args.output_name,
        args.disable_garbage_collection,
    )

    plots.make_plots(
        args.outdir,
        args.output_name,
        args.disable_garbage_collection,
    )
    
    if not args.disable_garbage_collection:
        garbage_collection(args.outdir, args.output_name)

    log.info('Analysis complete.')
    exit(0)


def prepare_outdir(outdir, out_name):
    """
    Creates the output directory and subdirectories for logs and intermediate
    files.

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name of the output directory.
    :type out_name: str
    :return: None
    :rtype: None
    """
    if not os.path.exists(outdir):
        os.makedirs(os.path.join(outdir, out_name))
        logs_path = os.path.join(outdir, out_name, 'logs')
        os.mkdir(logs_path)
    else:
        log.error(f'Output directory already exists: {outdir}')
        exit(1)
    
    
def check_database(db, outdir, out_name, collect_garbage):
    """
    Checks the contents of the provided reference sequence database to
    ensure proper header formatting and unambiguous sequences.

    :param db: Path to the reference sequence database.
    :type db: Path
    :return: None
    :rtype: None
    """
    log.info('Checking reference sequence database...')
    seqs = {}
    headers = []
    with open(db, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                headers.append(header)
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    headers = Counter(headers)
    for header in headers:
        if headers[header] > 1:
            logging.error(f'The following header is used for multiple sequences:'
                          f'\n{header}\n')
            exit(1)

    num_seqs = len(seqs)
    one_tenth_of_seqs = ceil(num_seqs / 10)
    log.info(f'Found {num_seqs} sequences in database.')

    log.info('Checking database sequences...')

    segment_lengths = {
        'PB2': (2260, 2360),
        'PB1': (2260, 2360),
        'PA': (2120, 2250),
        'HA': (1650, 1800),
        'NP': (1480, 1580),
        'NA': (1250, 1560),
        'M': (975, 1030),
        'NS': (815, 900),
    }

    log.info('Checking sequence headers')
    log.info('Checking sequence lengths within expected ranges for each segment')
    log.info('Checking for ambiguous or lower-case nucleotides')

    num_seqs_checked = 0
    for header, seq in seqs.items():

        #
        # Check header formatting
        if len(header.split('|')) != 4:
            log.error(f'The header for the following database entry does not contain the expected number of |-delimited fields:\n{header}\n')
            exit(1)
        else:
            accession, name, segment, subtype = header.split('|')
            
            # Check that the strain name is formatted correctly
            # (e.g. "A/California/04/2009(H1N1)")
            if any([name.count('(') != 1, name.count(')') != 1,
                    name[-1] != ')', name.index('(') > name.index(')')]):
                log.error(f'The strain_name(strain_subtype) for the following database entry is improperly formatted:\n{header}\n')
                exit(1)
            if segment not in segment_lengths.keys():
                log.error(f'The segment indicated for the following database entry is not recognized:\n{header}\n')
                exit(1)
        if len(seq) != sum(seq.count(base) for base in 'ATGC'):
            log.error(f'The sequence provided for the following database entry '
                      f'contains ambiguous or lower-case nucleotides:\n{header}\n')
            exit(1)

        #
        # Check sequence lengths
        min_length = segment_lengths[segment][0]
        max_length = segment_lengths[segment][1]
        if not (min_length <= len(seq) <= max_length):
            log.error(f'The sequence provided for the following database entry '
                      f'is not within the expected length range for its indicated '
                      f'segment ({segment}: {min_length} to {max_length} bases): '
                      f'\n{header}\n')
            exit(1)

        num_seqs_checked += 1
        if num_seqs_checked % one_tenth_of_seqs == 0:
            percent_seqs_checked = num_seqs_checked / num_seqs * 100
            log.info(f'Checked {num_seqs_checked} ({percent_seqs_checked:.2f}%) database sequences...')

    percent_seqs_checked = num_seqs_checked / num_seqs * 100

    log.info(f'Confirmed {num_seqs} ({percent_seqs_checked}%) database sequence headers are properly formatted.')
    log.info(f'Confirmed {num_seqs} ({percent_seqs_checked}%) database sequences are within expected length ranges for their indicated segments.')
    log.info(f'Confirmed {num_seqs} ({percent_seqs_checked}%) database sequences do not contain ambiguous or lower-case nucleotides.')

    db = os.path.abspath(db)
    db_suffixes = ['nhr', 'nin', 'nsq']
    if any([not os.path.exists(db + '.' + suffix) for suffix in db_suffixes]):
        terminal_command = f'makeblastdb -in {db} -dbtype nucl'
        process_name = 'makeblastdb_contigs'
        error_code = 5
        try:
            run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
        except Exception as e:
            log.error(f'Error creating BLAST database: at {os.path.dirname(db)} {e}')
            exit(1)


def run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage):
    """
    A generalized function for running subprocesses, logging their output, and
    trapping erroneous exit statuses.

    :param terminal_command: The command to be run in the terminal.
    :type terminal_command: str
    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name of the output directory.
    :type out_name: str
    :param process_name: Name of the process being run.
    :type process_name: str
    :param error_code: Exit status if the subprocess fails.
    :type error_code: int
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    :return: None
    :rtype: None
    """
    stdout_file = os.path.join(outdir, out_name, 'logs',
                               f'{process_name}_stdout.txt')
    stderr_file = os.path.join(outdir, out_name, 'logs',
                               f'{process_name}_stderr.txt')
    for file in [stdout_file, stderr_file]:
        if file != None:
            with open(file, 'w') as log_file:
                log_file.write('*' * 80 + '\n')
                log_file.write('Terminal command:' + '\n')
                log_file.write(terminal_command + '\n')
                log_file.write('*' * 80 + '\n')

    complete_process = None
    try:
        with open(stdout_file, 'w') as stdout_file:
            with open(stderr_file, 'w') as stderr_file:
                complete_process = subprocess.run(
                    terminal_command,
                    stdout=stdout_file,
                    stderr=stderr_file,
                    shell=True
                )
    except Exception as e:
        log.error(f'Error running subprocess {process_name}: {e}')
        if collect_garbage:
            garbage_collection(outdir, out_name)
        exit(error_code)

    return_code = complete_process.returncode
    if return_code != 0:
        log.error(f'Subprocess {process_name} failed (Exit status: {return_code})')
        if collect_garbage:
            garbage_collection(outdir, out_name)
        exit(error_code)


def garbage_collection(outdir, out_name):
    """
    Clean up unneccessary intermediate files (many of which occupy
    substantial amounts of storage).

    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name of the output subdirectory.
    :type out_name: str
    :return: None
    """
    if os.path.isdir(os.path.join(outdir, out_name, f'spades_output')):
        shutil.rmtree(os.path.join(outdir, out_name, f'spades_output'))

    files = [
        'R1.fq',
        'R2.fq',
        'contigs_blast.tsv',
        'scaffolds.fa',
        'scaffolds_blast.tsv',
        'alignment.sam',
        'pileup.vcf',
        'variants.bcf',
        'variants.bcf.csi',
        'low_cov.tsv',
        'ambig.tsv',
        'variants.tsv',
        'masked.bed',
        'reads_mapped.tsv',
        'depth_of_cov_samtools.tsv',
        'depth_of_cov_freebayes.tsv',
    ]

    bwa_suffixes = ['amb', 'ann', 'bwt', 'fai', 'pac', 'sa']
    for suffix in bwa_suffixes:
        files.append(f'{out_name}_mapping_refs.fa.{suffix}')

    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    for segment in segments:
        files += [
            f'{segment}_contigs.fa',
            f'{segment}_contigs.afa',
            f'{segment}_contigs.dnd'
        ]

    for file in files:
        file = os.path.join(out_name, file)
        if os.path.isfile(file):
            os.remove(file)
            log.info(f'Removed intermediate file: {file}')


def normalize_depth(
        outdir: Path,
        out_name: str,
        fwd_reads_raw: Path,
        rev_reads_raw: Path,
        depth: int,
        max_memory: int,
        collect_garbage: bool
):
    """
    BBNorm is run on the input reads to downsample regions of deep coverage
    (using a k-mer frequency approach). This balances coverage, increases
    analysis speed, and limits the impacts of artefactual reads.

    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name of the output directory.
    :type out_name: str
    :param fwd_reads_raw: Path to the raw forward reads.
    :type fwd_reads_raw: str
    :param rev_reads_raw: Path to the raw reverse reads.
    :type rev_reads_raw: str
    :param depth: Target depth of coverage.
    :type depth: int
    :param max_memory: Maximum memory to allocate to BBNorm.
    :type max_memory: int
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    """
    log.info('Normalizing depth of coverage and subsampling reads...')
    fwd_reads = os.path.join(outdir, out_name, f'R1.fq')
    rev_reads = os.path.join(outdir, out_name, f'R2.fq')
    terminal_command = (f'bbnorm.sh in={fwd_reads_raw} in2={rev_reads_raw} '
                        f'out={fwd_reads} out2={rev_reads} target={depth}')
    terminal_command = (terminal_command + f' -Xmx{max_memory}g'
                        if max_memory is not None else terminal_command)
    process_name = 'bbnorm'
    error_code = 2
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    bbnorm_stderr_log_path = os.path.join(outdir, out_name, 'logs', 'bbnorm_stderr.txt')
    bbnorm_log = parsers.parse_bbnorm_log(bbnorm_stderr_log_path)

    for pass_num, pass_stats in bbnorm_log.items():
        if not pass_num.startswith('pass_'):
            continue
        pass_num_int = int(pass_num.split('_')[-1])
        log.info(f'Normalization pass {pass_num_int}: Total reads in: {pass_stats["total_reads_in"]}')
        log.info(f'Normalization pass {pass_num_int}: Percent reads kept: {pass_stats["total_reads_kept_percent"]}%')
        log.info(f'Normalization pass {pass_num_int}: Percent unique: {pass_stats["percent_unique"]}%')
        log.info(f'Normalization pass {pass_num_int}: Average depth (unique kmers): {pass_stats["depth_average_unique_kmers"]}')
        log.info(f'Normalization pass {pass_num_int}: Average depth (all kmers): {pass_stats["depth_average_all_kmers"]}')
        log.info(f'Normalization pass {pass_num_int}: Approx. median read depth: {pass_stats["approx_read_depth_median"]}')

    with open(os.path.join(outdir, out_name, 'logs', 'bbnorm_log.json'), 'w') as f:
        json.dump(bbnorm_log, f, indent=4)
        f.write('\n')
    

    
def assemble_contigs(outdir: Path, out_name: str, collect_garbage: bool):
    """
    Normalized, downsampled reads are assembled de novo into contigs
    using SPAdes.

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    """
    log.info('Assembling reads into contigs...')
    spades_output = os.path.join(outdir, out_name, 'spades_output')
    fwd_reads = os.path.join(outdir, out_name, f'R1.fq')
    rev_reads = os.path.join(outdir, out_name, f'R2.fq')
    terminal_command = (f'spades.py --rnaviral --isolate -1 {fwd_reads} '
                        f'-2 {rev_reads} -o {spades_output}')
    process_name = 'spades'
    error_code = 3
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
    if not os.path.isfile(os.path.join(spades_output, 'contigs.fasta')):
        log.error('No contigs assembled! Aborting analysis.')
        if collect_garbage:
            garbage_collection(out_name)
        error_code = 4
        exit(error_code)
    else:
        num_contigs = 0
        with open(os.path.join(spades_output, 'contigs.fasta'), 'r') as f:
            for line in f:
                if line.startswith('>'):
                    num_contigs += 1
        log.info('Contigs assembled successfully.')
        log.info(f'Assembled {num_contigs} contigs.')


def blast_contigs(db, outdir, out_name, collect_garbage, threads, identity, length):
    """
    Contigs are aligned to reference sequences using BLASTn.

    :param db: Path to the reference sequence database.
    :type db: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    :param threads: Number of threads to use for BLASTn.
    :type threads: int
    :param identity: Minimum sequence identity for BLASTn hits.
    :type identity: float
    :param length: Minimum alignment length for BLASTn hits.
    :type length: int
    :return: BLASTn results.
    :rtype: pd.DataFrame
    """
    log.info('Aligning contigs to reference sequences...')
    db = os.path.abspath(db)
    blast_output = os.path.join(outdir, out_name, 'contigs_blast.tsv')
    spades_output = os.path.join(outdir, out_name, 'spades_output')
    contigs = os.path.join(spades_output, 'contigs.fasta')
    cols = [
        'qseqid',
        'sseqid',
        'pident',
        'length',
        'bitscore',
        'sstart',
        'send',
        'qseq',
        'sseq',
        'slen',
    ]
    cols_str = ' '.join(cols)
    terminal_command = (f'blastn -query {contigs} -db {db} '
                        f'-num_threads {threads} -outfmt "6 {cols_str}" '
                        f'> {blast_output}')
    process_name = 'blastn_contigs'
    error_code = 6
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
    blast_results = pd.read_csv(blast_output, names=cols, sep='\t')
    total_num_blast_results = len(blast_results)
    log.info('Contigs aligned to reference sequences.')
    log.info(f'Found {total_num_blast_results} total matches.')
    blast_results = blast_results[blast_results['pident']>=identity]
    blast_results = blast_results[blast_results['length']>=length]
    num_blast_results_after_identity_and_length_filter = len(blast_results)
    log.info(f'Found {num_blast_results_after_identity_and_length_filter} matches with at least {identity}% identity and {length} bp alignment length.')
    percent_blast_results_retained = num_blast_results_after_identity_and_length_filter / total_num_blast_results * 100
    log.info(f'Retained {percent_blast_results_retained:.2f}% matches for further analysis.')

    if len(blast_results) == 0:
        log.error(f'No contigs aligned to reference sequences! Aborting analysis.')
        if collect_garbage:
            garbage_collection(outdir, out_name)
        error_code = 7
        exit(error_code)

    return blast_results


def filter_contig_blast_results(blast_results, outdir, out_name, collect_garbage):
    """
    Contigs alignments are filtered to discard spurious alignments.

    The length and sequence identity of each alignment must exceed certain
    thresholds.

    Afterwards, a single reference sequence is selected for each
    genome segment. The reference sequence with the most positions covered by
    contigs. Ties are broken by:
   
      1. The highest sequence identity
      2. The longest reference sequence length
      3. First alphabetically

    Once a reference sequence has been chosen for each segment, only contig alignments
    to those reference sequences are retained.

    :param blast_results: BLASTn results.
    :type blast_results: pd.DataFrame
    :param out_name: Name used for outputs.
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    :return: Filtered BLASTn results.
    :rtype: pd.DataFrame
    """
    log.info('Filtering contig alignments...')
    # 
    # Annotate each ref seq with its segment and subtype.
    subject_annots = blast_results[['sseqid']].drop_duplicates()

    # sseqid format: accession|strain_name|segment|subtype
    get_segment = lambda row: row['sseqid'].split('|')[2]
    subject_annots['segment'] = subject_annots.apply(get_segment, axis=1)

    get_subtype = lambda row: row['sseqid'].split('|')[3]
    subject_annots['subtype'] = subject_annots.apply(get_subtype, axis=1)

    # Merge segment and subtype annotations back into blast_results
    blast_results = pd.merge(blast_results, subject_annots, on='sseqid')

    # Check for evidence of mixed infections. First, find best alignments(s)
    # for each contig. Next, check if each segment is represented by only one subtype.
    cols = ['qseqid', 'bitscore']
    # Find best alignment(s) for each contig, based on bitscore
    max_bitscores = blast_results[cols].groupby('qseqid').max().reset_index()
    # Use the qseqid and bitscore to find the best alignment(s) for each contig
    best_results = pd.merge(blast_results, max_bitscores, on=cols)

    log.info('Checking for mixed infections...')
    # Check for multiple subtypes for each segment
    segments = blast_results['segment'].unique()
    for segment in segments:
        log.info(f'Checking segment {segment}...')
        segment_results = best_results[best_results['segment']==segment]
        log.info(f'Found {len(segment_results)} contigs for segment {segment}.')
        segment_subtypes = segment_results['subtype'].unique()
        if len(segment_subtypes) > 1:
            segment_subtypes = ', '.join(segment_subtypes)
            log.error(f'Multiple subtypes detected for segment {segment} '
                      f'({segment_subtypes})! Aborting analysis.\n')
            if collect_garbage:
                garbage_collection(out_name)
            error_code = 8
            exit(error_code)
        else:
            if segment_subtypes[0] == 'none':
                log.info(f'No subtype determined for segment {segment}.')
            else:
                log.info(f'Subtype determined for segment {segment}: {segment_subtypes[0]}.')

    #
    # Find ref seq(s) most covered by contigs.
    log.info('Finding reference sequences most covered by contigs...')
    def count_cov_pos(data_frame):
        """
        Count the number of positions covered by contigs.
        Iterates through each alignment, and adds the range of positions
        covered by the alignment to a set. The length of the set is the number
        of covered positions.

        :param data_frame: DataFrame containing BLASTn results.
        :type data_frame: pd.DataFrame
        :return: Number of covered positions.
        :rtype: int
        """
        cov_positions = set()
        for index, row in data_frame.iterrows():
            start = min([row['sstart'], row['send']])
            end = max([row['sstart'], row['send']])
            cov_positions = cov_positions.union(set(range(start, end + 1)))
        return len(cov_positions)

    cols = [
        'sseqid',
        'segment',
        'subtype',
        'sstart',
        'send'
    ]
    group_cols = [
        'sseqid',
        'segment',
        'subtype'
    ]
    cov_pos = blast_results[cols].drop_duplicates()

    # Count the number of positions covered by contigs for each ref seq
    cov_pos = cov_pos.groupby(group_cols).apply(count_cov_pos).reset_index()
    log.info('Found covered positions for each reference sequence.')
    
    cov_pos.columns = [
        'sseqid',
        'segment',
        'subtype',
        'covered_positions'
    ]
    cols = [
        'segment',
        'subtype',
        'covered_positions'
    ]
    group_cols = ['segment', 'subtype']

    # Find the ref seq(s) with the most covered positions for each segment/subtype
    max_cov_pos = cov_pos[cols].drop_duplicates()
    max_cov_pos = max_cov_pos.groupby(group_cols).max().reset_index()
    log.info('Found reference sequence(s) with the most covered positions for each segment/subtype.')
    merge_cols = ['segment', 'subtype', 'covered_positions']
    max_cov_pos = pd.merge(cov_pos, max_cov_pos, on=merge_cols)

    # Log a summary of covered positions for each segment/subtype
    for segment in segments:
        segment_results = max_cov_pos[max_cov_pos['segment']==segment]
        subtype_to_log = '(undetermined)'
        if segment_results['subtype'].values[0] != 'none':
            subtype_to_log = segment_results['subtype'].values[0]
        for subtype in segment_results['subtype'].unique():
            num_covered_positions = segment_results[segment_results['subtype']==subtype]['covered_positions'].values[0]
            log.info(f'Segment: {segment}, Subtype: {subtype_to_log}: {num_covered_positions} covered positions.')

    max_cov_pos = max_cov_pos[['sseqid', 'covered_positions']]
    log.info('Found reference sequence(s) with the most covered positions for each segment/subtype.')
    blast_results = pd.merge(blast_results, max_cov_pos, on='sseqid')
    num_blast_results = len(blast_results)
    log.info(f'Found {num_blast_results} total contig alignments.')

    #
    # Find remaining ref seq(s) with most identical positions.
    def count_id_pos(data_frame):
        """
        Count the number of identical positions between contigs and reference.
        Determined by comparing the query and subject sequences from BLASTn.
        Iterates through each base in the query sequence, and if the base is
        a nucleotide and matches the subject sequence, the subject position is
        added to a set. The length of the set is the number of identical
        positions.

        :param data_frame: DataFrame containing BLASTn results.
        :type data_frame: pd.DataFrame
        :return: Number of identical positions.
        :rtype: int
        """
        identical_positions = set()
        for index, row in data_frame.iterrows():
            start = min([row['sstart'], row['send']])
            increment = 1 if row['sstart'] <= row['send'] else -1
            subject_position = start
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if sbase in 'ATGC' and qbase == sbase:
                    identical_positions.add(subject_position)
                if sbase != '-':
                    subject_position += increment
        return len(identical_positions)

    cols = [
        'sseqid',
        'segment',
        'subtype',
        'sstart',
        'send',
        'qseq',
        'sseq'
    ]
    group_cols = ['sseqid', 'segment', 'subtype']
    ident_pos = blast_results[cols].drop_duplicates()
    # Count the number of identical positions for each ref seq
    ident_pos = ident_pos.groupby(group_cols).apply(count_id_pos).reset_index()
    ident_pos.columns = ['sseqid', 'segment', 'subtype', 'identical_positions']
    cols = ['segment', 'subtype', 'identical_positions']
    group_cols = ['segment', 'subtype']

    # Find the ref seq(s) with the most identical positions for each segment/subtype
    max_ident_pos = ident_pos[cols].drop_duplicates()
    max_ident_pos = max_ident_pos.groupby(group_cols).max().reset_index()
    merge_cols = ['segment', 'subtype', 'identical_positions']
    max_ident_pos = pd.merge(ident_pos, max_ident_pos, on=merge_cols)
    log.info('Found reference sequence(s) with the most identical positions for each segment/subtype.')
    for segment in segments:
        segment_results = max_ident_pos[max_ident_pos['segment']==segment]
        for subtype in segment_results['subtype'].unique():
            subtype_to_log = subtype if subtype != 'none' else 'undetermined'
            num_identical_positions = segment_results[segment_results['subtype']==subtype]['identical_positions'].values[0]
            log.info(f'Segment {segment}, Subtype {subtype_to_log}: {num_identical_positions} identical positions.')

    cols = ['sseqid', 'identical_positions']
    max_ident_pos = max_ident_pos[cols].drop_duplicates()
    blast_results = pd.merge(blast_results, max_ident_pos, on='sseqid')
    
    # Take longest remaining ref seq for each segment/subtype.
    cols = ['segment', 'subtype', 'slen']
    group_cols = ['segment', 'subtype']
    longest_sseq = blast_results[cols].drop_duplicates()
    longest_sseq = longest_sseq.groupby(group_cols).min().reset_index()
    blast_results = pd.merge(blast_results, longest_sseq, on=cols)

    # Take first alphabetical remaining ref seq for each segment/subtype.
    cols = ['segment', 'subtype', 'sseqid']
    group_cols = ['segment', 'subtype']
    first_alpha_sseq = blast_results[cols].drop_duplicates()
    first_alpha_sseq = first_alpha_sseq.groupby(group_cols).min().reset_index()
    blast_results = pd.merge(blast_results, first_alpha_sseq, on=cols)
    for segment in segments:
        segment_results = blast_results[blast_results['segment']==segment]
        for subtype in segment_results['subtype'].unique():
            subtype_to_log = subtype if subtype != 'none' else '(undetermined)'
            ref_seq = segment_results[segment_results['subtype']==subtype]['sseqid'].values[0]
            log.info(f'Segment {segment}, Subtype: {subtype_to_log}. Selected reference sequence: {ref_seq}')
    return blast_results


def make_scaffold_seqs(blast_results, outdir, out_name, collect_garbage):
    """
    A scaffold sequence is created for each genome segment by joining and
    collapsing all the contigs describing that segment.

    Unaligned leading and trailing sequences are trimmed from the
    contigs. Next, leading and trailing Ns are added to the contig so that it is
    properly positioned within the segment (based on the subject-start and
    subject-end coordinates of its alignment to the selected reference sequence).
    Next, clustalW is used to generate a multiple sequence alignment of the
    trimmed, positioned contigs. This multiple sequence alignment is used to
    generate a consensus sequence of the regions of the segment covered by
    contigs.
    """
    log.info('Creating scaffolds...')
    # Make sure contigs are all in the forward orientation.
    
    rev_comp_bases = {'A': 'T',
                      'T': 'A',
                      'G': 'C',
                      'C': 'G',
                      'N': 'N',
                      '-': '-',
                      'W': 'W',
                      'S': 'S',
                      'M': 'K',
                      'K': 'M',
                      'R': 'Y',
                      'Y': 'R',
                      'B': 'V',
                      'D': 'H',
                      'H': 'D',
                      'V': 'B'}
    rev_comp_seq = lambda seq: ''.join(rev_comp_bases[base]
                                       for base in seq[::-1])
    get_start = lambda row: min([row['sstart'], row['send']])

    blast_results['start'] = blast_results.apply(get_start, axis=1)
    get_end = lambda row: max([row['sstart'], row['send']])
    blast_results['end'] = blast_results.apply(get_end, axis=1)

    def flip_qseq(row):
        """
        Flip the query sequence if the alignment is in the reverse orientation.

        :param row: BLASTn result.
        :type row: pd.Series
        :return: Flipped query sequence.
        :rtype: str
        """
        if row['sstart'] > row['send']:
            log.info(f'Flipping seq for contig: {row["qseqid"]}')
            return rev_comp_seq(row['qseq'])
        else:
            return row['qseq']

    blast_results['qseq'] = blast_results.apply(flip_qseq, axis=1)
    def flip_sseq(row):
        """
        Flip the subject sequence if the alignment is in the reverse orientation.

        :param row: BLASTn result.
        :type row: pd.Series
        :return: Flipped subject sequence.
        :rtype: str
        """
        if row['sstart'] > row['send']:
            log.info(f'Flipping seq for ref: {row["sseqid"]}')
            return rev_comp_seq(row['sseq'])
        else:
            return row['sseq']

    blast_results['sseq'] = blast_results.apply(flip_sseq, axis=1)

    # Trim contigs based on their alignments to reference sequences. Also
    # add leading and trailing Ns to contig so that it is properly positioned
    # within the genome segment.
    segments = blast_results['segment'].unique()
    contig_counter = {segment: 0 for segment in segments}
    scaffold_seqs = {}
    for segment in segments:
        contigs = os.path.join(outdir, out_name, f'{segment}_contigs.fa')
        with open(contigs, 'w') as f:
            contig_results = blast_results[blast_results['segment']==segment]
            for index, row in contig_results.iterrows():
                header = f'>{segment}_contig_{contig_counter[segment]}\n'
                f.write(header)
                seq = 'N' * (row['start'] - 1)
                seq += row['qseq'].replace('-', '')
                seq += ('N' * (row['slen'] - row['end']))
                f.write(seq + '\n')
                contig_counter[segment] += 1
        log.info(f'Wrote {contig_counter[segment]} contigs for segment {segment} to {contigs}')
        # Generate multiple sequence alignments of trimmed/positioned contigs.
        log.info(f'Aligning contigs for segment {segment}...')
        aligned_contigs = os.path.join(outdir, out_name, f'{segment}_contigs.afa')
        if contig_counter[segment] > 1:
            log.info(f'Generating multiple sequence alignment for segment {segment}...')
            terminal_command = (f'clustalw -INFILE={contigs} '
                                f'-OUTFILE={aligned_contigs} -OUTPUT=FASTA')
            process_name = f'clustalw_{segment}'
            error_code = 9
            run(terminal_command, outdir, out_name, process_name, error_code,
                collect_garbage)
        else:
            log.info(f'Only one contig for segment {segment}, skipping alignment.')
            shutil.copyfile(contigs, aligned_contigs)
            
        # Replace leading and trailing Ns with dots so that they are ignored
        # when determining consensus bases.
        seqs = {}
        with open(aligned_contigs, 'r') as input_file:
            for line in input_file:
                if line[0] == '>':
                    header = line.strip()
                    seqs[header] = ''
                else:
                    seqs[header] += line.strip()
        clean_seqs = []
        for seq in seqs.values():
            head_len = len(seq) - len(seq.lstrip('N-'))
            tail_len = len(seq) - len(seq.rstrip('N-'))
            seq = seq.strip('N-')
            seq = ('.' * head_len) + seq
            seq += ('.' * tail_len)
            clean_seqs.append(seq)

        # Check that all seqs in the multiple seq alignment are the
        # same length. '''
        alignment_lengths = set(len(seq) for seq in clean_seqs)
        if len(alignment_lengths) > 1:
            log.error(f'Multiple sequence alignment for {segment} '
                      f'generated unequal alignment lengths! Aborting analysis.\n')
            if collect_garbage:
                garbage_collection(out_name)
            error_code = 10
            exit(error_code)

        # Create consensus sequence of multiple seq alignments, i.e. the
        # scaffolds.
        alignment_length = list(alignment_lengths)[0]
        scaffold_seq = ''
        for i in range(alignment_length):
            bases = Counter(seq[i] for seq in clean_seqs if seq[i] not in '.')
            if bases.most_common(1) == []:
                scaffold_seq += 'N'
            elif bases.most_common(1)[0][1] / len(bases) > 0.5:
                scaffold_seq += bases.most_common(1)[0][0]
            else:
                scaffold_seq += 'N'
        scaffold_seqs[segment] = scaffold_seq

    # Write out scaffolds.
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    segments = [segment for segment in segments if segment in scaffold_seqs]
    scaffold_seqs = {segment: scaffold_seqs[segment] for segment in segments}
    scaffolds = os.path.join(outdir, out_name, 'scaffolds.fa')
    with open(scaffolds, 'w') as f:
        for segment, seq in scaffold_seqs.items():
            header = f'>{out_name}|{segment}_scaffold'
            f.write(header + '\n')
            f.write(seq + '\n')

    log.info(f'Wrote {len(scaffold_seqs)} scaffolds to {scaffolds}')

    return scaffold_seqs


def blast_scaffolds(db, outdir, out_name, collect_garbage, threads):
    """
    Scaffold sequences are aligned to reference sequences using BLASTn.

    :param db: Path to the reference sequence database.
    :type db: Path
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    :param threads: Number of threads to use for BLASTn.
    """
    log.info('Aligning scaffolds to reference sequences...')
    num_db_seqs = sum(line[0] == '>' for line in open(db, 'r').readlines())
    missing_db_files = []
    db = os.path.abspath(db)
    for suffix in ['nhr', 'nin', 'nsq']:
        db_file = db + '.' + suffix
        if not os.path.exists(db_file):
            missing_db_files.append(db_file)

    if missing_db_files:
        terminal_command = (f'makeblastdb -in {db} -dbtype nucl')
        process_name = 'makeblastdb_scaffolds'
        error_code = 11
        run(terminal_command, outdir, out_name, process_name, error_code,
            collect_garbage)
    scaffolds = os.path.join(outdir, out_name, 'scaffolds.fa')
    blast_output = os.path.join(outdir, out_name, 'scaffolds_blast.tsv')
    cols = [
        'qseqid',
        'sseqid',
        'bitscore',
        'sstart',
        'send',
        'qseq',
        'sseq'
    ]
    cols_str = ' '.join(cols)
    terminal_command = (f'blastn -query {scaffolds} -db {db} '
                        f'-num_threads {threads} -max_target_seqs {num_db_seqs} '
                        f'-outfmt "6 {cols_str}" > {blast_output}')
    process_name = 'blastn_scaffolds'
    error_code = 12
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
    blast_results = pd.read_csv(blast_output, names=cols, sep='\t')
    num_blast_results = len(blast_results)
    if num_blast_results == 0:
        log.error(f'No scaffolds aligned to reference sequences! '
                  f'Aborting analysis.\n')
        if collect_garbage:
            garbage_collection(out_name)
        error_code = 13
        exit(error_code)
    else:
        log.info('Scaffolds aligned to reference sequences.')
        log.info(f'Found {num_blast_results} total matches.')
        for segment in blast_results['qseqid'].unique():
            segment_results = blast_results[blast_results['qseqid']==segment]
            ref_seq = segment_results['sseqid'].values[0]
            log.info(f'Selected reference sequence for segment {segment}: {ref_seq}')

    return blast_results


def filter_scaffold_blast_results(blast_results):
    """
    A single reference sequence is chosen for each segment scaffold.

    First, the bitscores of all alignments between a scaffold and a reference sequence
    are summed. The bitscore sum is calculated for each scaffold-reference
    sequence pairing. The reference sequence giving the highest bitscore sum is
    selected for each scaffold, with ties being broken by using the first
    alphabetically. Once the best-matching reference sequence has been selected for
    each segment scaffold, all other alignments are discarded.

    :param blast_results: BLASTn results.
    :type blast_results: pd.DataFrame
    :return: Filtered BLASTn results.
    :rtype: pd.DataFrame
    """
    log.info('Filtering scaffold alignments...')

    # Remove reversed alignments (they should be artefactual at this point).
    # check number of reversed alignments
    num_reversed_alignments = len(blast_results[blast_results['send'] < blast_results['sstart']])
    log.info(f'Found {num_reversed_alignments} reversed alignments.')
    log.info('Removing reversed alignments...')
    blast_results = blast_results[blast_results['send'] > blast_results['sstart']]

    #Annotate scaffold seqs with segment.
    query_annots = blast_results[['qseqid']].drop_duplicates()
    get_segment = lambda row: row['qseqid'].split('|')[1].split('_')[0]
    query_annots['segment'] = query_annots.apply(get_segment, axis=1)
    blast_results = pd.merge(blast_results, query_annots, on='qseqid')

    # Find best-matching reference sequence for each segment.
    cols = ['segment', 'sseqid', 'bitscore']
    group_cols = ['segment', 'sseqid']
    combo_scores = blast_results[cols].groupby(group_cols).sum().reset_index()
    cols = ['segment', 'bitscore']
    group_cols = ['segment']
    max_scores = combo_scores[cols].groupby(group_cols).max().reset_index()
    merge_cols = ['segment', 'bitscore']
    max_scores = pd.merge(max_scores, combo_scores, on=merge_cols)
    cols = ['segment', 'sseqid']
    group_cols = ['segment']
    first_alpha = max_scores[cols].groupby(group_cols).min().reset_index()
    merge_cols = ['segment', 'sseqid']
    blast_results = pd.merge(blast_results, first_alpha, on=merge_cols)
    for segment in blast_results['segment'].unique():
        segment_results = blast_results[blast_results['segment']==segment]
        ref_seq = segment_results['sseqid'].values[0]
        log.info(f'Selected reference sequence for segment {segment}: {ref_seq}')

    return blast_results


def make_mapping_refs(blast_results, db, outdir, out_name):
    """
    Mapping references are created for each genome segment. These consist of
    the scaffold for that segment, with all Ns in the scaffold filled-in using
    the corresponding positions from that scaffold's best-matching reference
    sequence.

    :param blast_results: BLASTn results.
    :type blast_results: pd.DataFrame
    :param db: Path to the reference sequence database.
    :type db: Path
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    
    """
    
    log.info('Creating mapping references...')
    # Create dict with best-matching ref seq for each segment.
    sseqids = blast_results['sseqid'].unique()
    best_ref_seqs = {seq_name: '' for seq_name in sseqids}
    with open(db, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip()[1:]
            elif header in best_ref_seqs:
                best_ref_seqs[header] += line.strip()

    # Create mapping ref for each segment.
    def make_map_ref(data_frame):
        data_frame = data_frame.sort_values(by=['sstart', 'send'],
                                            ascending=[True, False])
        ref_seq = best_ref_seqs[data_frame['sseqid'].min()]
        last_position = 0
        seq = ''
        for index, row in data_frame.iterrows():
            if row['sstart'] > last_position:
                seq += ref_seq[last_position:row['sstart'] - 1]
            if row['send'] > last_position:
                qseq = row['qseq'].upper()
                sseq = row['sseq'].upper()
                if row['sstart'] <= last_position:
                    start = (last_position - row['sstart']) + 1
                    qseq = qseq[start:]
                    sseq = sseq[start:]
                for qbase, sbase in zip(qseq, sseq):
                    if qbase in 'ATGC':
                        seq += qbase
                    else:
                        seq += sbase
                last_position = row['send']
        seq += ref_seq[last_position:].upper()
        seq = seq.replace('-', '')
        return seq
    cols = ['sseqid', 'sstart', 'send', 'qseq', 'sseq']
    group_cols = ['sseqid']
    blast_results = blast_results[cols]
    blast_results = blast_results.groupby(group_cols).apply(make_map_ref)
    blast_results = blast_results.reset_index()
    blast_results.columns = ['sseqid', 'mapping_seq']
    # Annotate segment and subtype.
    get_segment = lambda row: row['sseqid'].split('|')[2].split('_')[0]
    blast_results['segment'] = blast_results.apply(get_segment, axis=1)
    get_subtype = lambda row: row['sseqid'].split('|')[3].split('_')[0]
    blast_results['subtype'] = blast_results.apply(get_subtype, axis=1)

    # Write mapping refs to FASTA.
    segment_order = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    get_segment_order = lambda row: segment_order.index(row['segment'])
    blast_results['sort'] = blast_results.apply(get_segment_order, axis=1)
    blast_results = blast_results.sort_values(by='sort')

    mapping_refs = os.path.join(outdir, out_name, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs, 'w') as output_file:
        num_refs = 0
        for index, row in blast_results.iterrows():
            num_refs += 1
            accession, ref_name, segment, subtype = row['sseqid'].split('|')[:4]
            accession = accession.lstrip('>')
            ref_name = ref_name.replace('(', '|').replace(')', '')
            header = f'>{out_name}|{segment}|{subtype}|{accession}|{ref_name}'
            seq = row['mapping_seq']
            output_file.write(header + '\n')
            output_file.write(seq + '\n')

    log.info(f'Wrote {num_refs} mapping references to {mapping_refs}')

    return mapping_refs


def map_reads(outdir, out_name, collect_garbage, min_qual):
    """
    Normalized, downsampled reads (normalize_depth func) are mapped to the
    mapping references (make_mapping_refs func) using BWA mem. The alignment
    is filtered to retain only paired reads, then sorted and indexed.

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    :param min_qual: Minimum mapping quality.
    :type min_qual: int
    """
    log.info('Mapping reads to references...')
    mapping_refs = os.path.join(outdir, out_name, f'{out_name}_mapping_refs.fa')
    terminal_command = (f'bwa index {mapping_refs}')
    process_name = 'bwa_index'
    error_code = 14
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    fwd_reads = os.path.join(outdir, out_name, 'R1.fq')
    rev_reads = os.path.join(outdir, out_name, 'R2.fq')
    alignment = os.path.join(outdir, out_name, 'alignment.sam')
    terminal_command = (f'bwa mem {mapping_refs} {fwd_reads} {rev_reads} '
                        f'> {alignment}')
    process_name = 'bwa_mem'
    error_code = 15
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    filtered_alignment = os.path.join(outdir, out_name, f'{out_name}_alignment.bam')
    samtools_filter_flags = '2828'
    log.info(f'Filtering alignment with samtools flags: {samtools_filter_flags}.')
    log.info('Removing unmapped reads, secondary alignments, and supplementary alignments.')
    log.info(f'Minimum mapping quality: {min_qual}')
    terminal_command = (f'samtools view -f 1 -F {samtools_filter_flags} -q {min_qual} '
                        f'-h {alignment} | samtools sort -o {filtered_alignment}')
    process_name = 'samtools_view'
    error_code = 16
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    log.info(f'Indexing alignment...')
    terminal_command = (f'samtools index {filtered_alignment}')
    process_name = 'samtools_index'
    error_code = 17
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    log.info(f'Wrote alignment to {filtered_alignment}')

    return filtered_alignment
    


def call_variants(outdir, out_name, min_qual, max_depth, collect_garbage):
    """
    FreeBayes is used to create a pileup and call variants from the
    BAM file output (map_reads func).

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param min_qual: Minimum base quality.
    :type min_qual: int
    :param max_depth: Maximum read depth.
    :type max_depth: int
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    """
    log.info('Calling variants...')
    mapping_refs = os.path.join(outdir, out_name, f'{out_name}_mapping_refs.fa')
    filtered_alignment = os.path.join(outdir, out_name, f'{out_name}_alignment.bam')
    pileup = os.path.join(outdir, out_name, 'pileup.vcf')
    log.info(f'Minimum mapping quality: {min_qual}')
    log.info(f'Minimum base quality: {min_qual}')
    log.info(f'Maximum read depth: {max_depth}')
    terminal_command = (f'freebayes -f {mapping_refs} {filtered_alignment} -p 1 '
                        f'--limit-coverage {max_depth} '
                        f'--min-mapping-quality {min_qual} '
                        f'--min-base-quality {min_qual} --pooled-continuous '
                        f'--report-monomorphic --haplotype-length 0 '
                        f'--min-alternate-count 1 --min-alternate-fraction 0 '
                        f'> {pileup}')
    process_name = 'freebayes'
    error_code = 18
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    log.info(f'Wrote pileup to {pileup}')
    


def mask_ambig_low_cov(outdir, out_name, min_depth, vaf_call, vaf_ambig,
                       min_qual, collect_garbage):
    """
    The FreeBayes VCF output is parsed, analyzing total read depth and
    variant read depth at each position. This allows positions to be masked
    for low coverage based on the read depth considered by FreeBayes (which
    could be lower than the read depth in the BAM depending on how FreeBayes
    applies it mapping quality and base quality filters). This also allows
    positions to be masked as ambiguous when the number of reads differing
    from the reference exceeds a threshold, but is not sufficient enough to
    confidently call as a specific variant.

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param min_depth: Minimum read depth.
    :type min_depth: int
    :param vaf_call: Minimum variant allele frequency to call a variant.
    :type vaf_call: float
    :param vaf_ambig: Minimum variant allele frequency to mask as ambiguous.
    :type vaf_ambig: float
    :param min_qual: Minimum base quality.
    :type min_qual: int
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    """
    log.info('Masking ambiguous and low coverage positions...')
    # Open input/output files and initialize dicts.
    pileup = open(os.path.join(outdir, out_name, 'pileup.vcf'), 'r')
    variants = open(os.path.join(outdir, out_name, f'{out_name}_variants.vcf'), 'w')
    depth_of_cov = os.path.join(outdir, out_name, f'depth_of_cov_freebayes.tsv')
    depth_of_cov = open(depth_of_cov, 'w')
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    low_cov_pos = {segment: set() for segment in segments}
    ambig_pos = {segment: set() for segment in segments}
    variant_pos = {segment: set() for segment in segments}
    segment_name, segment_length = dict(), {None: 0}
    segment_length[None] = 0

    # Parse header
    line = pileup.readline()
    while line != '' and line[0] == '#':
        variants.write(line)
        if line[:10] == '##contig=<':
            name = line.strip().split('<ID=')[1].split(',length=')[0]
            segment = name.split('|')[1]
            length = int(line.strip().split(',length=')[1].split('>')[0])
            segment_name[segment] = name
            segment_length[segment] = length
        line = pileup.readline()

    # Parse body
    last_segment = None
    last_position = 0
    while line != '':
        fields = line.strip().split('\t')
        fields = (fields[0], fields[1], fields[3], fields[4], fields[5],
                  fields[8], fields[9])
        name, position, ref, alt, qual, keys, values = fields
        segment = name.split('|')[1]
        if segment != last_segment:
            if last_position < segment_length[last_segment]:
                for p in range(last_position + 1,
                               segment_length[last_segment] + 1):
                    low_cov_pos[last_segment].add(p)
                    depth_of_cov_line = [name, str(p), '0']
                    depth_of_cov_line = '\t'.join(depth_of_cov_line)
                    depth_of_cov.write(depth_of_cov_line + '\n')
            last_position = 0
        last_segment = segment
        position = int(position)
        if position != last_position + 1:
            for p in range(last_position + 1, position):
                low_cov_pos[segment].add(p)
                depth_of_cov_line = [name, str(p), '0']
                depth_of_cov_line = '\t'.join(depth_of_cov_line)
                depth_of_cov.write(depth_of_cov_line + '\n')
        qual = float(qual)
        info = {k: v for k, v in zip(keys.split(':'), values.split(':'))}
        if 'DP' in info and info['DP'].isnumeric():
            total_depth = int(info['DP'])
        else:
            total_depth = 0
        if 'AO' in info:
            alt_depths = tuple(int(i) if i.isnumeric() else 0
                               for i in info['AO'].split(','))
        else:
            alt_depths = (0, )
        max_alt_depth = max(alt_depths)
        total_alt_depth = sum(alt_depths)
        max_vaf = max_alt_depth / total_depth if total_depth > 0 else 0
        total_vaf = total_alt_depth / total_depth if total_depth > 0 else 0
        if all([qual >= min_qual, max_vaf >= vaf_call,
                total_depth >= min_depth]):
            variants.write(line)
            variant_pos[segment].add(position)
        position -= 1
        for p in ref:
            position += 1
            depth_of_cov_line = [name, str(position), str(total_depth)]
            depth_of_cov_line = '\t'.join(depth_of_cov_line)
            depth_of_cov.write(depth_of_cov_line + '\n')
            if total_depth < min_depth:
                low_cov_pos[segment].add(position)
            elif total_vaf >= vaf_ambig and max_vaf < vaf_call:
                ambig_pos[segment].add(position)
        last_position = position
        line = pileup.readline()
    if last_position < segment_length[last_segment]:
        for p in range(last_position + 1, segment_length[last_segment] + 1):
            low_cov_pos[last_segment].add(p)
            depth_of_cov_line = [name, str(p), '0']
            depth_of_cov_line = '\t'.join(depth_of_cov_line)
            depth_of_cov.write(depth_of_cov_line + '\n')

    # Close input/output files
    pileup.close()
    variants.close()
    depth_of_cov.close()

    # Convert sets of low cov positions into tuples representing zero-indexed
    # spans of masked positions (start, end).
    masked_pos = dict()
    for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
        masked_pos[segment] = low_cov_pos[segment].union(ambig_pos[segment])
        masked_pos[segment] = sorted(masked_pos[segment])
    spans = {segment: set() for segment in segments}
    segments = [segment for segment in segments
                if masked_pos[segment] != list()]
    for segment in segments:
        span_start = masked_pos[segment][0]
        for pos_A, pos_B in zip(masked_pos[segment][:-1],
                                masked_pos[segment][1:]):
            if pos_B != pos_A + 1:
                span_end = pos_A
                spans[segment].add((span_start - 1, span_end - 1))
                span_start = pos_B
        span_end = masked_pos[segment][-1]
        spans[segment].add((span_start - 1, span_end - 1))
    spans = {segment: sorted(spans[segment]) for segment in segments}

    # Write spans of low cov positions to TSV file for depth of coverage
    # plots.
    low_cov_path = os.path.join(outdir, out_name, 'low_cov.tsv')
    with open(low_cov_path, 'w') as f:
        for segment, segment_spans in spans.items():
            for (start, end) in segment_spans:
                line = [segment_name[segment], start, end]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote low coverage positions to {low_cov_path}')

    # Write ambiguous positions to TSV file.
    ambig_path = os.path.join(outdir, out_name, 'ambig.tsv')
    with open(ambig_path, 'w') as f:
        for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
            for position in ambig_pos[segment]:
                line = [segment_name[segment], position]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote ambiguous positions to {ambig_path}')

    # Write variant positions to TSV file.
    variant_path = os.path.join(outdir, out_name, 'variants.tsv')
    with open(variant_path, 'w') as f:
        for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
            for position in variant_pos[segment]:
                line = [segment_name[segment], position]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote variant positions to {variant_path}')
              
    # Write spans of masked positions to BED file in BedGraph format.
    masked_path = os.path.join(outdir, out_name, 'masked.bed')          
    with open(masked_path, 'w') as f:
        for segment, segment_spans in spans.items():
            for (start, end) in segment_spans:
                line = [segment_name[segment], start, end + 1, 0]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote masked positions to {masked_path}')


def make_consensus_seqs(outdir, out_name, collect_garbage):
    """
    High quality variants and masked positions (mask_ambig_low_cov func) are
    applied to the mapping references (make_mapping_refs) to generate the final
    consensus sequences for each segment.

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    """
    log.info('Generating consensus sequences...')

    # Zip and index VCF.
    variants = os.path.join(outdir, out_name, f'{out_name}_variants.vcf')
    zipped_variants = os.path.join(outdir, out_name, 'variants.bcf')
    terminal_command = (f'bcftools view {variants} -Ob -o {zipped_variants}')
    process_name = 'bcftools_view'
    error_code = 19
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
    terminal_command = (f'bcftools index {zipped_variants}')
    process_name = 'bcftools_index'
    error_code = 20
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
    # Apply variants to mapping refs.
    log.info('Applying variants to mapping references...')
    mapping_refs = os.path.join(outdir, out_name, f'{out_name}_mapping_refs.fa')
    masked = os.path.join(outdir, out_name, 'masked.bed')
    consensus_seqs = os.path.join(outdir, out_name, f'{out_name}_consensus_seqs.fa')
    terminal_command = (f'cat {mapping_refs} | bcftools consensus -m {masked} '
                        f'{zipped_variants} > {consensus_seqs}')
    process_name = 'bcftools_consensus'
    error_code = 21
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)

    # Reformat FASTA headers and remove whitespace.
    clean_seqs = {}
    with open(consensus_seqs, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip()
                clean_seqs[header] = ''
            else:
                clean_seqs[header] += line.strip().upper()

    with open(consensus_seqs, 'w') as f:
        for header, seq in clean_seqs.items():
            header = '|'.join(header.split('|')[:3]) + '|'
            f.write(header + '\n')
            f.write(seq + '\n')

    log.info(f'Wrote consensus sequences to {consensus_seqs}')

    # Check that consensus seq lenghts are within expected range. '''
    segment_lengths = {'PB2': (2260, 2360), 'PB1': (2260, 2360), 'PA': (2120, 2250),
                       'HA': (1650, 1800), 'NP': (1480, 1580), 'NA': (1250, 1560),
                       'M': (975, 1030), 'NS': (815, 900)}
    for header, seq in clean_seqs.items():
        segment = header.split('|')[1]
        min_length = segment_lengths[segment][0]
        max_length = segment_lengths[segment][1]
        if not (min_length <= len(seq) <= max_length):
            log.error(f'The consensus sequence generated for segment '
                      f'{segment} is not within the expected length range '
                      f'({min_length} to {max_length} bases).\n')
            exit(22)

def write_report(outdir, out_name, collect_garbage):
    """
    Generate a report for each segment.

    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param collect_garbage: Whether to remove intermediate files.
    :type collect_garbage: bool
    """
    log.info('Writing report...')
    
    # Count reads mapped to each segment and add to report.
    filtered_alignment = os.path.join(outdir, out_name, f'{out_name}_alignment.bam')
    reads_mapped = os.path.join(outdir, out_name, 'reads_mapped.tsv')
    terminal_command = (f'samtools idxstats {filtered_alignment} > '
                        f'{reads_mapped}')
    process_name = 'samtools_idxstats'
    error_code = 23
    run(terminal_command, outdir, out_name, process_name, error_code, collect_garbage)
    cols = 'seq_name seq_length reads_mapped reads_unmapped'.split(' ')
    reads_mapped = pd.read_csv(reads_mapped, sep='\t', names=cols)
    reads_mapped = reads_mapped.replace('*', np.nan).dropna()
    get_seq_name = lambda row: '|'.join(row['seq_name'].split('|')[:3])
    reads_mapped['seq_name'] = reads_mapped.apply(get_seq_name, axis=1) 
    cols = ['seq_name', 'reads_mapped', 'seq_length']
    report = reads_mapped[cols].drop_duplicates()
    
    # Annotate segment and subtype.
    get_segment = lambda row: row['seq_name'].split('|')[1]
    report['segment'] = report.apply(get_segment, axis=1)
    get_subtype = lambda row: row['seq_name'].split('|')[2]
    report['subtype'] = report.apply(get_subtype, axis=1)
    
    #Add scaffold completeness to report.
    seqs = {}
    scaffolds = os.path.join(outdir, out_name, 'scaffolds.fa')
    with open(scaffolds, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                seq_name = line[1:].strip()
                seqs[seq_name] = ''
            else:
                seqs[seq_name] += line.strip()
    completeness = {}
    for seq_name, seq in seqs.items():
        segment = seq_name.split('|')[1].split('_')[0]
        perc = sum(seq.count(base) for base in 'ATGC') * 100 / len(seq)
        perc = round(perc, 2)
        completeness[segment] = perc
    report['scaffold_completeness'] = report['segment'].map(completeness)

    # Add consensus completeness to report.
    seqs = {}
    consensus_seqs = os.path.join(outdir, out_name, f'{out_name}_consensus_seqs.fa')
    with open(consensus_seqs, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                seq_name = line[1:].strip()
                seqs[seq_name] = ''
            else:
                seqs[seq_name] += line.strip()
    completeness = {}
    for seq_name, seq in seqs.items():
        segment = seq_name.split('|')[1]
        perc = sum(seq.count(base) for base in 'ATGC') * 100 / len(seq)
        perc = round(perc, 2)
        completeness[segment] = perc
    report['consensus_completeness'] = report['segment'].map(completeness)
    
    # Add best ref seq to report.
    ref_seqs_used = {}
    mapping_refs = os.path.join(outdir, out_name, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                line = line[1:].strip().split('|')
                seq_name, segment, subtype = line[:3]
                accession, ref_name, ref_subtype = line[3:]
                seq_name = f'{seq_name}|{segment}|{subtype}'
                ref_seqs_used[seq_name] = (f'{accession}|{ref_name}'
                                           f'({ref_subtype})')
    report['ref_seq_used'] = report['seq_name'].map(ref_seqs_used)
    
    # Write report to TSV file.
    segment_order ='PB2 PB1 PA HA NP NA M NS'.split(' ')
    get_sort_value = lambda row: segment_order.index(row['segment'])
    report['sort'] = report.apply(get_sort_value, axis=1)
    report = report.sort_values(by='sort')
    cols = ['seq_name', 'segment', 'subtype', 'reads_mapped', 'seq_length',
            'scaffold_completeness', 'consensus_completeness', 'ref_seq_used']
    report = report[cols]
    report.to_csv(os.path.join(outdir, out_name, f'{out_name}_report.tsv'),
                  index=False, sep='\t')


if __name__ == '__main__':
    main()
