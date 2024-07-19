import datetime
import os
import json
import logging
import shutil
import subprocess
import sys

from collections import Counter
from pathlib import Path
from typing import List

import pandas as pd

from fluviewer import parsers


log = logging.getLogger(__name__)


error_messages_by_code = {
    1: 'Error creating output directory.',
    2: 'Error running BBNorm.',
    3: 'Error running SPAdes.',
    4: 'No contigs assembled.',
    5: 'Error running BLASTn.',
    6: 'Error running BLASTn.',
    7: 'No contigs aligned to reference sequences.',
    8: 'Multiple subtypes detected for segment.',
    9: 'Error running ClustalW.',
    10: 'Error generating scaffold sequences.',
    11: 'Error running BLASTn.',
    12: 'No contigs aligned to reference sequences.',
    13: 'Multiple subtypes detected for segment.',
    14: 'Error running ClustalW.',
    15: 'Error generating scaffold sequences.',
    16: 'Error running BLASTn.',
    17: 'No contigs aligned to reference sequences.',
    18: 'Multiple subtypes detected for segment.',
}


def run(terminal_command, outdir, out_name, process_name, error_code):
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
    :return: Exit status of the subprocess.
    :rtype: int
    """
    logs_dir = os.path.join(outdir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    stdout_file = os.path.join(logs_dir, f'{process_name}_stdout.txt')
    stderr_file = os.path.join(logs_dir, f'{process_name}_stderr.txt')

    script_file = os.path.join(outdir, f'{process_name}_script.sh')
    with open(script_file, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write(terminal_command)
        f.write('\n')
    os.chmod(script_file, 0o755)

    complete_process = None
    try:
        with open(stdout_file, 'w') as stdout_file:
            with open(stderr_file, 'w') as stderr_file:
                complete_process = subprocess.run(
                    script_file,
                    stdout=stdout_file,
                    stderr=stderr_file,
                    shell=True
                )
    except Exception as e:
        log.error(f'Error running subprocess {process_name}: {e}')
        return error_code

    return_code = complete_process.returncode

    if return_code != 0:
        log.error(f'Subprocess {process_name} failed (Exit status: {return_code})')
        return error_code

    return return_code


def normalize_depth(
        inputs: dict,
        outdir: Path,
        out_name: str,
        depth: int,
        max_memory: int,
):
    """
    BBNorm is run on the input reads to downsample regions of deep coverage
    (using a k-mer frequency approach). This balances coverage, increases
    analysis speed, and limits the impacts of artefactual reads.

    :param inputs: Dictionary of input files, with keys 'raw_reads_fwd' and 'raw_reads_rev'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name of the output directory.
    :type out_name: str
    :param depth: Target depth of coverage.
    :type depth: int
    :param max_memory: Maximum memory to allocate to BBNorm.
    :type max_memory: int
    :return: Summary of the analysis.
    :rtype: dict
    """
    timestamp_analysis_start = datetime.datetime.now().isoformat()
    log.info('Normalizing depth of coverage and subsampling reads...')
    logs_dir = os.path.join(outdir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)

    input_reads_fwd = inputs.get('input_reads_fwd', None)
    input_reads_rev = inputs.get('input_reads_rev', None)
    input_reads_long = inputs.get('input_reads_long', None)

    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }
    
    normalized_reads_fwd = os.path.join(outdir, f'{out_name}-normalized_R1.fastq')
    normalized_reads_rev = os.path.join(outdir, f'{out_name}-normalized_R2.fastq')
    normalized_reads_long = os.path.join(outdir, f'{out_name}-normalized_RL.fastq')

    if input_reads_fwd and input_reads_rev:
        terminal_command = (f'bbnorm.sh in={input_reads_fwd} in2={input_reads_rev} '
                            f'out={normalized_reads_fwd} out2={normalized_reads_rev} target={depth}')
    elif input_reads_long:
        terminal_command = (f'bbnorm.sh in={input_reads_long} out={normalized_reads_long} target={depth}')
        
    terminal_command = (terminal_command + f' -Xmx{max_memory}g'
                        if max_memory is not None else terminal_command)

    # add gzip compression to output files
    if input_reads_fwd and input_reads_rev:
        terminal_command += f'\n\ngzip -f {normalized_reads_fwd}\n'
        terminal_command += f'\ngzip -f {normalized_reads_rev}\n'
    elif input_reads_long:
        terminal_command += f'\n\ngzip -f {normalized_reads_long}\n'

    process_name = 'bbnorm'
    error_code = 2

    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running BBNorm (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    bbnorm_stderr_log_path = os.path.join(logs_dir, 'bbnorm_stderr.txt')

    parsed_bbnorm_log = parsers.parse_bbnorm_log(bbnorm_stderr_log_path)

    for pass_num, pass_stats in parsed_bbnorm_log.items():
        if not pass_num.startswith('pass_'):
            continue
        pass_num_int = int(pass_num.split('_')[-1])
        log.info(f'Normalization pass {pass_num_int}: Total reads in: {pass_stats["total_reads_in"]}')
        log.info(f'Normalization pass {pass_num_int}: Percent reads kept: {pass_stats["total_reads_kept_percent"]}%')
        log.info(f'Normalization pass {pass_num_int}: Percent unique: {pass_stats["percent_unique"]}%')
        log.info(f'Normalization pass {pass_num_int}: Average depth (unique kmers): {pass_stats["depth_average_unique_kmers"]}')
        log.info(f'Normalization pass {pass_num_int}: Average depth (all kmers): {pass_stats["depth_average_all_kmers"]}')
        log.info(f'Normalization pass {pass_num_int}: Approx. median read depth: {pass_stats["approx_read_depth_median"]}')

    with open(os.path.join(logs_dir, 'bbnorm_log.json'), 'w') as f:
        json.dump(parsed_bbnorm_log, f, indent=4)
        f.write('\n')

    timestamp_analysis_complete = datetime.datetime.now().isoformat()

    if input_reads_fwd and input_reads_rev:
        outputs = {
            'normalized_reads_fwd': os.path.abspath(normalized_reads_fwd) + '.gz',
            'normalized_reads_rev': os.path.abspath(normalized_reads_rev) + '.gz',
        }
    elif input_reads_long:
        outputs = {
            'normalized_reads_long': os.path.abspath(normalized_reads_long) + '.gz',
        }

    analysis_summary = {
        'process_name': process_name,
        'timestamp_analysis_start': timestamp_analysis_start,
        'timestamp_analysis_complete': timestamp_analysis_complete,
        'return_code': return_code,
        'inputs': inputs,
        'outputs': outputs,
    }

    with open(os.path.join(outdir, 'analysis_summary.json'), 'w') as f:
        json.dump(analysis_summary, f, indent=4)
        f.write('\n')

    return analysis_summary


def assemble_contigs(
        inputs: dict,
        outdir: Path,
        out_name: str,
):
    """
    Normalized, downsampled reads are assembled de novo into contigs
    using SPAdes.

    :param inputs: Dictionary of input files, with keys 'normalized_reads_fwd' and 'normalized_reads_rev'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs
    :type out_name: str
    :return: Summary of the analysis.
    :rtype: dict
    """
    log.info('Assembling reads into contigs...')
    timestamp_analysis_start = datetime.datetime.now().isoformat()

    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }

    logs_dir = os.path.join(outdir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    
    spades_output = os.path.join(outdir, 'spades_output')
    fwd_reads = inputs.get('reads_fwd', None)
    rev_reads = inputs.get('reads_rev', None)

    os.makedirs(spades_output, exist_ok=True)

    terminal_command = (f'spades.py --rnaviral --isolate -1 {fwd_reads} '
                        f'-2 {rev_reads} -o {spades_output}')

    process_name = 'spades'
    error_code = 3

    script_file = os.path.join(outdir, f'{process_name}_script.sh')
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    analysis_summary['return_code'] = return_code
    if not os.path.isfile(os.path.join(spades_output, 'contigs.fasta')):
        log.error('No contigs assembled! Aborting analysis.')
        error_code = 4
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        analysis_summary['inputs'] = inputs
        return analysis_summary

    num_contigs = 0
    src_contigs_path = os.path.join(spades_output, 'contigs.fasta')
    dest_contigs_path = os.path.join(outdir, f'{out_name}_contigs.fasta')
    shutil.copy(src_contigs_path, dest_contigs_path)
    with open(dest_contigs_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                num_contigs += 1
    log.info('Contigs assembled successfully.')
    log.info(f'Assembled {num_contigs} contigs.')

    timestamp_analysis_complete = datetime.datetime.now().isoformat()

    outputs = {
        'contigs': os.path.abspath(dest_contigs_path),
    }

    analysis_summary['process_name'] = process_name
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete
    analysis_summary['return_code'] = return_code
    analysis_summary['outputs'] = outputs

    with open(os.path.join(outdir, 'analysis_summary.json'), 'w') as f:
        json.dump(analysis_summary, f, indent=4)
        f.write('\n')

    return analysis_summary


def filter_contig_blast_results(blast_results, outdir, out_name, identity, length):
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
    :param identity: Minimum sequence identity for BLASTn hits.
    :type identity: float
    :param length: Minimum alignment length for BLASTn hits.
    :type length: int
    :return: Filtered BLASTn results.
    :rtype: pd.DataFrame
    """
    log.info('Filtering contig BLAST alignments...')

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
            error_code = 8
            exit(error_code)
        else:
            if len(segment_subtypes) == 0:
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


def blast_contigs(inputs, outdir, out_name, threads, identity, length):
    """
    Contigs are aligned to reference sequences using BLASTn.

    :param inputs: Dictionary of input files, with keys 'contigs'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param threads: Number of threads to use for BLASTn.
    :type threads: int
    :param identity: Minimum sequence identity for BLASTn hits.
    :type identity: float
    :param length: Minimum alignment length for BLASTn hits.
    :type length: int
    :return: BLASTn results.
    :rtype: pd.DataFrame
    """
    log.info('BLASTing contigs to reference sequences...')
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    timestamp_analysis_start = datetime.datetime.now().isoformat()

    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }

    db = os.path.abspath(inputs.get('database', None))
    blast_output_path = os.path.join(outdir, f'{out_name}_contigs_blast.tsv')
    contigs = inputs.get('contigs', None)

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
    with open(blast_output_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')

    cols_str = ' '.join(cols)

    terminal_command = (f'blastn -query {contigs} -db {db} '
                        f'-num_threads {threads} -outfmt "6 {cols_str}" '
                        f'>> {blast_output_path}')

    process_name = 'blastn_contigs'
    error_code = 6
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)

    full_blast_results = pd.read_csv(blast_output_path, sep='\t')

    total_num_blast_results = len(full_blast_results)

    if total_num_blast_results == 0:
        log.error(f'No contigs aligned to reference sequences! Aborting analysis.')
        error_code = 7
        exit(error_code)

    log.info('Contigs BLASTed against reference sequences.')
    log.info(f'Found {total_num_blast_results} total matches.')

    blast_results = full_blast_results[full_blast_results['pident'] >= identity]
    blast_results = blast_results[blast_results['length'] >= length]

    identity_and_length_filtered_blast_output_path = os.path.join(outdir, f'{out_name}_contigs_blast_identity_and_length_filtered.tsv')
    blast_results.to_csv(identity_and_length_filtered_blast_output_path, sep='\t', index=False)
    num_blast_results_after_identity_and_length_filter = len(blast_results)
    log.info(f'Found {num_blast_results_after_identity_and_length_filter} matches with at least {identity}% identity and {length} bp alignment length.')

    percent_blast_results_retained = num_blast_results_after_identity_and_length_filter / total_num_blast_results * 100
    log.info(f'Retained {percent_blast_results_retained:.2f}% matches for further analysis.')

    filtered_blast_results = filter_contig_blast_results(blast_results, outdir, out_name, identity, length)

    timestamp_analysis_complete = datetime.datetime.now().isoformat()

    if len(filtered_blast_results) == 0:
        log.error(f'No contigs aligned to reference sequences! Aborting analysis.')
        error_code = 7
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    num_filtered_blast_results = len(filtered_blast_results)
    log.info(f'Remaining contig alignments after filtering: {num_filtered_blast_results}.')

    filtered_blast_output_path = os.path.join(outdir, f'{out_name}_contigs_blast_filtered.tsv')
    filtered_blast_results.to_csv(filtered_blast_output_path, sep='\t', index=False)

    outputs = {
        'all_contig_blast_results': os.path.abspath(blast_output_path),
        'identity_and_length_filtered_contig_blast_results': os.path.abspath(identity_and_length_filtered_blast_output_path),
        'filtered_contig_blast_results': os.path.abspath(filtered_blast_output_path),
    }

    analysis_summary = {
        'process_name': process_name,
        'timestamp_analysis_start': timestamp_analysis_start,
        'timestamp_analysis_complete': timestamp_analysis_complete,
        'return_code': return_code,
        'inputs': inputs,
        'outputs': outputs,
    }
    with open(os.path.join(outdir, 'analysis_summary.json'), 'w') as f:
        json.dump(analysis_summary, f, indent=4)
        f.write('\n')

    return analysis_summary


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
    log.info('Filtering scaffold BLAST alignments...')

    # Remove reversed alignments (they should be artefactual at this point).
    # check number of reversed alignments
    num_reversed_alignments = len(blast_results[blast_results['send'] < blast_results['sstart']])
    log.info(f'Found {num_reversed_alignments} reversed alignments.')
    log.info('Removing reversed alignments...')
    blast_results = blast_results[blast_results['send'] > blast_results['sstart']]
    num_remaining_alignments = len(blast_results)
    log.info(f'Remaining scaffold alignments after removing reversed alignments: {num_remaining_alignments}.')

    #Annotate scaffold seqs with segment.
    query_annots = blast_results[['qseqid']].drop_duplicates()
    get_segment = lambda row: row['qseqid'].split('|')[1].split('_')[0]
    query_annots['segment'] = query_annots.apply(get_segment, axis=1)
    blast_results = pd.merge(blast_results, query_annots, on='qseqid')

    # Find best-matching reference sequence for each segment.
    # Best match is defined as the reference sequence with the highest bitscore sum.
    # Ties are broken by selecting the first alphabetically.
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


def make_scaffold_seqs(inputs, outdir, out_name):
    """
    A scaffold sequence is created for each genome segment by joining and
    collapsing all the contigs describing that segment.

    Unaligned leading and trailing sequences are trimmed from the
    contigs.

    Next, leading and trailing Ns are added to the contig so that it is
    properly positioned within the segment (based on the subject-start and
    subject-end coordinates of its alignment to the selected reference sequence).

    Next, clustalW is used to generate a multiple sequence alignment of the
    trimmed, positioned contigs. This multiple sequence alignment is used to
    generate a consensus sequence of the regions of the segment covered by
    contigs.

    :param inputs: Dictionary of input files, with keys 'filtered_contig_blast_results'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :return: Summary of the analysis.
    :rtype: dict
    """
    log.info('Creating scaffolds...')
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp_analysis_start = datetime.datetime.now().isoformat()
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }

    outputs = {}
    
    filtered_contig_blast_results_path = inputs.get('filtered_contig_blast_results', None)
    
    # The 'segment' field of the filtered blast results may include 'NA' values, which
    # represent the Neuraminidase (NA) segment.
    # Disable the na_filter option to prevent pandas from converting 'NA' to NaN.
    blast_results = pd.read_csv(filtered_contig_blast_results_path, sep='\t', na_filter=False)

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
    blast_results_flipped_path = os.path.join(outdir, f'{out_name}_contigs_blast_filtered_flipped.tsv')

    blast_results.to_csv(blast_results_flipped_path, sep='\t', index=False)
    log.info(f'Wrote flipped BLASTn results to {blast_results_flipped_path}')
    outputs['flipped_contig_blast_results'] = os.path.abspath(blast_results_flipped_path)

    clustalw_return_code = None
    # Trim contigs based on their alignments to reference sequences. Also
    # add leading and trailing Ns to contig so that it is properly positioned
    # within the genome segment.
    segments = blast_results['segment'].unique()
    contig_counter = {segment: 0 for segment in segments}
    scaffold_seqs = {}
    outputs['segment_contigs'] = {}
    outputs['segment_contigs_alignments'] = {}
    for segment in segments:
        segment_contigs_path = os.path.join(outdir, f'{segment}_contigs.fa')
        with open(segment_contigs_path, 'w') as f:
            contig_results = blast_results[blast_results['segment']==segment]
            for index, row in contig_results.iterrows():
                header = f'>{segment}_contig_{contig_counter[segment]}\n'
                f.write(header)
                seq = 'N' * (row['start'] - 1)
                seq += row['qseq'].replace('-', '')
                seq += ('N' * (row['slen'] - row['end']))
                f.write(seq + '\n')
                contig_counter[segment] += 1
        log.info(f'Wrote {contig_counter[segment]} contigs for segment {segment} to {segment_contigs_path}')
        outputs['segment_contigs'][segment] = os.path.abspath(segment_contigs_path)

        # Generate multiple sequence alignments of trimmed/positioned contigs.
        log.info(f'Aligning contigs for segment {segment} with clustalw...')        
        segment_contigs_alignment_path = os.path.join(outdir, f'{segment}_contigs.afa')
        if contig_counter[segment] > 1:
            log.info(f'Generating multiple sequence alignment for segment {segment}...')
            terminal_command = (f'clustalw -INFILE={segment_contigs_path} '
                                f'-OUTFILE={segment_contigs_alignment_path} -OUTPUT=FASTA')
            process_name = f'clustalw_{segment}'
            error_code = 9
            clustalw_return_code = run(terminal_command, outdir, out_name, process_name, error_code)
            analysis_summary['return_code'] = clustalw_return_code
        else:
            log.info(f'Only one contig for segment {segment}, skipping alignment.')
            shutil.copyfile(segment_contigs_path, segment_contigs_alignment_path)
            analysis_summary['return_code'] = 0

        outputs['segment_contigs_alignments'][segment] = os.path.abspath(segment_contigs_alignment_path)

        # Replace leading and trailing Ns with dots so that they are ignored
        # when determining consensus bases.
        seqs = {}
        with open(segment_contigs_alignment_path, 'r') as input_file:
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
            error_code = 10
            analysis_summary['return_code'] = error_code
            analysis_summary['error_message'] = error_messages_by_code[error_code]
            analysis_summary['outputs'] = outputs
            return analysis_summary
        elif len(alignment_lengths) == 0:
            log.error(f'No sequences in the multiple sequence alignment for {segment}! Aborting analysis.\n')
            error_code = 11
            analysis_summary['return_code'] = error_code
            analysis_summary['error_message'] = error_messages_by_code[error_code]
            analysis_summary['outputs'] = outputs
            return analysis_summary

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
    scaffolds_path = os.path.join(outdir, f'{out_name}_scaffolds.fa')
    with open(scaffolds_path, 'w') as f:
        for segment, seq in scaffold_seqs.items():
            header = f'>{out_name}|{segment}_scaffold'
            f.write(header + '\n')
            f.write(seq + '\n')

    outputs['scaffolds'] = os.path.abspath(scaffolds_path)

    log.info(f'Wrote {len(scaffold_seqs)} scaffolds to {scaffolds_path}')

    timestamp_analysis_complete = datetime.datetime.now().isoformat()

    analysis_summary['process_name'] = 'make_scaffold_seqs'
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete
    analysis_summary['outputs'] = outputs

    
    return analysis_summary


def blast_scaffolds(inputs, outdir, out_name, threads):
    """
    Scaffold sequences are aligned to reference sequences using BLASTn.

    :param inputs: Dictionary of input files, with keys 'scaffolds' and 'database'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param threads: Number of threads to use for BLASTn.
    :type threads: int
    :return: BLASTn results.
    """
    log.info('BLASTing scaffolds against reference sequences...')
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    timestamp_analysis_start = datetime.datetime.now().isoformat()

    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }
    db = inputs.get('database', None)
    
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
        run(terminal_command, outdir, out_name, process_name, error_code)
    else:
        log.info(f'BLAST database already exists for {db} Skipping BLAST db creation.')

    scaffolds_path = inputs.get('scaffolds', None)
    blast_output = os.path.join(outdir, f'{out_name}_scaffolds_blast.tsv')
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
    with open(blast_output, 'w') as f:
        f.write('\t'.join(cols) + '\n')

    num_db_seqs = sum(line[0] == '>' for line in open(db, 'r').readlines())
    terminal_command = (f'blastn -query {scaffolds_path} -db {db} '
                        f'-num_threads {threads} -max_target_seqs {num_db_seqs} '
                        f'-outfmt "6 {cols_str}" >> {blast_output}')
    process_name = 'blastn_scaffolds'
    error_code = 12
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running BLASTn (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    blast_results = pd.read_csv(blast_output, sep='\t')
    num_blast_results = len(blast_results)
    if num_blast_results == 0:
        log.error(f'No scaffolds aligned to reference sequences! '
                  f'Aborting analysis.\n')
        error_code = 13
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary
    
    log.info('Scaffolds aligned to reference sequences.')
    log.info(f'Found {num_blast_results} total scaffold BLAST alignments.')

    filtered_blast_results = filter_scaffold_blast_results(blast_results)

    num_filtered_blast_results = len(filtered_blast_results)
    log.info(f'Remaining scaffold alignments after filtering: {num_filtered_blast_results}.')
    filtered_blast_output = os.path.join(outdir, f'{out_name}_scaffolds_blast_filtered.tsv')
    filtered_blast_results.to_csv(filtered_blast_output, sep='\t', index=False)
    log.info(f'Wrote filtered scaffold alignments to {filtered_blast_output}')
    
    timestamp_analysis_complete = datetime.datetime.now().isoformat()

    outputs = {
        'all_scaffold_blast_results': os.path.abspath(blast_output),
        'filtered_scaffold_blast_results': os.path.abspath(filtered_blast_output),
    }
    analysis_summary['process_name'] = process_name
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete
    analysis_summary['return_code'] = return_code
    analysis_summary['outputs'] = outputs

    with open(os.path.join(outdir, 'analysis_summary.json'), 'w') as f:
        json.dump(analysis_summary, f, indent=4)
        f.write('\n')

    return analysis_summary



def make_mapping_refs(inputs, outdir, out_name):
    """
    Mapping references are created for each genome segment. These consist of
    the scaffold for that segment, with all Ns in the scaffold filled-in using
    the corresponding positions from that scaffold's best-matching reference
    sequence.

    :param inputs: Dictionary of input files, with keys 'filtered_scaffold_blast_results' and 'database'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
  
    """
    log.info('Creating mapping references...')
    
    filtered_scaffold_blast_results_path = inputs.get('filtered_scaffold_blast_results', None)
    blast_results = pd.read_csv(filtered_scaffold_blast_results_path, sep='\t')
    db = inputs.get('database', None)
    db = os.path.abspath(db)
    
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

    mapping_refs_path = os.path.join(outdir, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs_path, 'w') as f:
        num_refs = 0
        for index, row in blast_results.iterrows():
            num_refs += 1
            accession, ref_name, segment, subtype = row['sseqid'].split('|')[:4]
            accession = accession.lstrip('>')
            ref_name = ref_name.replace('(', '|').replace(')', '')
            header = f'>{out_name}|{segment}|{subtype}|{accession}|{ref_name}'
            seq = row['mapping_seq']
            f.write(header + '\n')
            f.write(seq + '\n')

    log.info(f'Wrote {num_refs} mapping references to {mapping_refs_path}')

    return mapping_refs_path


def map_reads(inputs, outdir, out_name, min_qual):
    """
    Reads  are mapped to the mapping references (produced by make_mapping_refs func) using BWA mem.
    The alignment is filtered to retain only paired reads, then sorted and indexed.

    :param inputs: Dictionary of input files, with keys: ['reads_fwd', 'reads_rev' and 'mapping_refs'].
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param min_qual: Minimum mapping quality.
    :type min_qual: int
    """
    log.info('Mapping reads to references...')
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp_analysis_start = datetime.datetime.now().isoformat()

    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }
    
    mapping_refs_path = make_mapping_refs(inputs, outdir, out_name)

    terminal_command = (f'bwa index {mapping_refs_path}')
    process_name = 'bwa_index'
    error_code = 14
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running BWA index (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    fwd_reads = inputs.get('reads_fwd', None)
    rev_reads = inputs.get('reads_rev', None)
    alignment_path = os.path.join(outdir, f'{out_name}_alignment.sam')
    terminal_command = (f'bwa mem {mapping_refs_path} {fwd_reads} {rev_reads} '
                        f'> {alignment_path}')
    process_name = 'bwa_mem'
    error_code = 15
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running BWA mem (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    filtered_alignment_path = os.path.join(outdir, f'{out_name}_alignment.bam')
    
    samtools_require_flags = [
        'PAIRED', # 0x1
    ]
    samtools_require_flags_str = ','.join(samtools_require_flags)
    samtools_exclude_flags = [
        'UNMAP',           # 0x4
        'MUNMAP',          # 0x8
        'SECONDARY',       # 0x100
        'QCFAIL',          # 0x200
        'SUPPLEMENTARY',   # 0x800
    ]
    samtools_exclude_flags_str = ','.join(samtools_exclude_flags)

    log.info(f'Filtering alignment with sam require flags: {samtools_require_flags_str}.')
    log.info(f'Filtering alignment with sam exclude flags: {samtools_exclude_flags_str}.')
    log.info(f'See: http://www.htslib.org/doc/samtools-flags.html for info on sam flags.')
    log.info('Removing unmapped reads, secondary alignments, and supplementary alignments.')
    log.info(f'Applying minimum mapping quality: {min_qual}')
    terminal_command = (f'samtools view '
                        f'--require-flags {samtools_require_flags_str} '
                        f'--exclude-flags {samtools_exclude_flags_str} '
                        f'--min-MQ {min_qual} '
                        f'--with-header {alignment_path} | samtools sort -o {filtered_alignment_path}')
    process_name = 'samtools_view'
    error_code = 16
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running samtools view (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    log.info(f'Indexing alignment...')
    terminal_command = (f'samtools index {filtered_alignment_path}')
    process_name = 'samtools_index'
    error_code = 17
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running samtools index (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    log.info(f'Wrote alignment to {filtered_alignment_path}')

    outputs = {
        'mapping_refs': os.path.abspath(mapping_refs_path),
        'alignment': os.path.abspath(filtered_alignment_path),
        'alignment_index': os.path.abspath(filtered_alignment_path + '.bai'),
    }
    timestamp_analysis_complete = datetime.datetime.now().isoformat()
    analysis_summary['process_name'] = 'map_reads'
    analysis_summary['outputs'] = outputs
    analysis_summary['return_code'] = return_code
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete

    with open(os.path.join(outdir, 'analysis_summary.json'), 'w') as f:
        json.dump(analysis_summary, f, indent=4)
        f.write('\n')

    return analysis_summary


def mask_ambig_low_cov(
        variants_raw_path: Path,
        outdir: Path,
        out_name: str,
        params: dict,
):
    """
    The FreeBayes VCF output is parsed, analyzing total read depth and
    variant read depth at each position.

    This allows positions to be masked for low coverage based on the read
    depth considered by FreeBayes (which could be lower than the read depth
    in the BAM depending on how FreeBayes applies it mapping quality and
    base quality filters).

    This also allows positions to be masked as ambiguous when the number of
    reads differing from the reference exceeds a threshold, but is not sufficient
    enough to confidently call as a specific variant.

    :param variants_raw_path: Path to the raw FreeBayes VCF output.
    :type variants_raw_path: str
    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name used for outputs.
    :type out_name: str
    :param params: Dictionary of parameters, with keys 'min_depth', 'min_qual', 'vaf_call', 'vaf_ambig'.
    :type params: dict
    :return: Output files. Keys are 'filtered_variants' and 'depth_of_cov'.
    :rtype: dict
    """
    log.info('Masking ambiguous and low coverage positions...')
    outputs = {}
    # Open input/output files and initialize dicts.
    variants_raw = open(variants_raw_path, 'r')
    variants_filtered_path = os.path.join(outdir, f'{out_name}_variants.vcf')
    variants_filtered = open(variants_filtered_path, 'w')
    depth_of_cov_path = os.path.join(outdir, f'{out_name}_depth_of_cov_freebayes.tsv')
    depth_of_cov = open(depth_of_cov_path, 'w')
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    low_cov_pos = {segment: set() for segment in segments}
    ambig_pos = {segment: set() for segment in segments}
    variant_pos = {segment: set() for segment in segments}
    segment_name, segment_length = dict(), {None: 0}
    segment_length[None] = 0

    # Parse header
    line = variants_raw.readline()
    while line != '' and line[0] == '#':
        variants_filtered.write(line)
        if line[:10] == '##contig=<':
            name = line.strip().split('<ID=')[1].split(',length=')[0]
            segment = name.split('|')[1]
            length = int(line.strip().split(',length=')[1].split('>')[0])
            segment_name[segment] = name
            segment_length[segment] = length
        line = variants_raw.readline()

    # Parse body
    last_segment = None
    last_position = 0
    min_depth = params['min_depth']
    min_qual = params['min_mapping_quality']
    vaf_call = params['variant_threshold_calling']
    vaf_ambig = params['variant_threshold_masking']
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
            variants_filtered.write(line)
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
        line = variants_raw.readline()
    if last_position < segment_length[last_segment]:
        for p in range(last_position + 1, segment_length[last_segment] + 1):
            low_cov_pos[last_segment].add(p)
            depth_of_cov_line = [name, str(p), '0']
            depth_of_cov_line = '\t'.join(depth_of_cov_line)
            depth_of_cov.write(depth_of_cov_line + '\n')

    # Close input/output files
    variants_raw.close()
    variants_filtered.close()
    depth_of_cov.close()

    outputs['filtered_variants'] = os.path.abspath(variants_filtered_path)
    log.info(f'Wrote filtered variants to {variants_filtered_path}')
    outputs['depth_of_cov'] = os.path.abspath(depth_of_cov_path)
    log.info(f'Wrote depth of coverage to {depth_of_cov_path}')

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
    low_cov_path = os.path.join(outdir, f'{out_name}_low_cov.tsv')
    with open(low_cov_path, 'w') as f:
        for segment, segment_spans in spans.items():
            for (start, end) in segment_spans:
                line = [segment_name[segment], start, end]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote low coverage positions to {low_cov_path}')
    outputs['low_cov'] = os.path.abspath(low_cov_path)

    # Write ambiguous positions to TSV file.
    ambig_path = os.path.join(outdir, f'{out_name}_ambig.tsv')
    with open(ambig_path, 'w') as f:
        for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
            for position in ambig_pos[segment]:
                line = [segment_name[segment], position]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote ambiguous positions to {ambig_path}')
    outputs['ambiguous_positions'] = os.path.abspath(ambig_path)

    # Write variant positions to TSV file.
    variants_tsv_path = os.path.join(outdir, f'{out_name}_variants.tsv')
    with open(variants_tsv_path, 'w') as f:
        for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
            for position in variant_pos[segment]:
                line = [segment_name[segment], position]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote variant positions to {variants_tsv_path}')
    outputs['variant_positions'] = os.path.abspath(variants_tsv_path)

    # Write spans of masked positions to BED file in BedGraph format.
    masked_path = os.path.join(outdir, f'{out_name}_masked.bed')          
    with open(masked_path, 'w') as f:
        for segment, segment_spans in spans.items():
            for (start, end) in segment_spans:
                line = [segment_name[segment], start, end + 1, 0]
                line = '\t'.join(str(i) for i in line)
                f.write(line + '\n')
    log.info(f'Wrote masked positions to {masked_path}')
    outputs['masked_positions'] = os.path.abspath(masked_path)

    return outputs


def call_variants(inputs, outdir, out_name, params):
    """
    FreeBayes is used to create a pileup and call variants from the
    BAM file output (map_reads func).

    :param inputs: Dictionary of input files, with keys 'mapping_refs' and 'alignment'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    :param min_qual: Minimum base quality.
    :type min_qual: int
    :param max_depth: Maximum read depth.
    :type max_depth: int
    :return: Summary of the analysis.
    :rtype: dict
    """
    log.info('Calling variants...')
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    timestamp_analysis_start = datetime.datetime.now().isoformat()
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }

    mapping_refs = inputs.get('mapping_refs', None)
    filtered_alignment = inputs.get('alignment', None)

    variants_raw_path = os.path.join(outdir, f'{out_name}_variants_raw.vcf')
    min_qual = params.get('min_mapping_quality', 20)
    max_depth = params.get('coverage_limit', 1000)

    log.info(f'Minimum mapping quality: {min_qual}')
    log.info(f'Minimum base quality: {min_qual}')
    log.info(f'Maximum read depth: {max_depth}')
    terminal_command = (f'freebayes -f {mapping_refs} {filtered_alignment} -p 1 '
                        f'--limit-coverage {max_depth} '
                        f'--min-mapping-quality {min_qual} '
                        f'--min-base-quality {min_qual} --pooled-continuous '
                        f'--report-monomorphic --haplotype-length 0 '
                        f'--min-alternate-count 1 --min-alternate-fraction 0 '
                        f'> {variants_raw_path}')
    process_name = 'freebayes'
    error_code = 18
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running FreeBayes (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    log.info(f'Wrote raw variants to {variants_raw_path}')

    masking_outputs = mask_ambig_low_cov(
        variants_raw_path,
        outdir,
        out_name,
        params,
    )

    outputs = {
        'variants_raw':           os.path.abspath(variants_raw_path),
        'variants_filtered':      os.path.abspath(masking_outputs['filtered_variants']),
        'masked_positions':       os.path.abspath(masking_outputs['masked_positions']),
        'depth_of_cov_freebayes': os.path.abspath(masking_outputs['depth_of_cov']),
        'low_coverage_positions': os.path.abspath(masking_outputs['low_cov']),
        'ambiguous_positions':    os.path.abspath(masking_outputs['ambiguous_positions']),
        'variant_positions':      os.path.abspath(masking_outputs['variant_positions']),
    }

    analysis_summary['process_name'] = 'call_variants'
    timestamp_analysis_complete = datetime.datetime.now().isoformat()
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete
    analysis_summary['return_code'] = return_code
    analysis_summary['outputs'] = outputs

    with open(os.path.join(outdir, 'analysis_summary.json'), 'w') as f:
        json.dump(analysis_summary, f, indent=4)
        f.write('\n')

    return analysis_summary

    
def make_consensus_seqs(inputs, outdir, out_name):
    """
    High quality variants and masked positions (mask_ambig_low_cov func) are
    applied to the mapping references (make_mapping_refs) to generate the final
    consensus sequences for each segment.

    :param inputs: Dictionary of input files, with keys 'variants_filtered' and 'mapping_refs'.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    """
    log.info('Generating consensus sequences...')
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp_analysis_start = datetime.datetime.now().isoformat()
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }
    outputs = {}

    # Zip and index VCF.
    variants = inputs.get('variants_filtered', None)
    zipped_variants = os.path.join(outdir, f'{out_name}_variants.bcf')
    terminal_command = (f'bcftools view {variants} -Ob -o {zipped_variants}')
    process_name = 'bcftools_view'
    error_code = 19
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running bcftools view (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    terminal_command = (f'bcftools index {zipped_variants}')
    process_name = 'bcftools_index'
    error_code = 20
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running bcftools index (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

    # Apply variants to mapping refs.
    log.info('Applying variants to mapping references...')
    mapping_refs = inputs.get('mapping_refs', None)
    masked = inputs.get('masked_positions', None)
    consensus_seqs = os.path.join(outdir, f'{out_name}_consensus_seqs.fa')
    terminal_command = (f'cat {mapping_refs} | bcftools consensus -m {masked} '
                        f'{zipped_variants} > {consensus_seqs}')
    process_name = 'bcftools_consensus'
    error_code = 21
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error(f'Error running bcftools consensus (Exit status: {return_code})')
        analysis_summary['return_code'] = error_code
        analysis_summary['error_message'] = error_messages_by_code[error_code]
        return analysis_summary

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
    outputs['consensus_seqs'] = os.path.abspath(consensus_seqs)

    # Check that consensus seq lenghts are within expected range. '''
    segment_lengths = {'PB2': (2260, 2360), 'PB1': (2260, 2360), 'PA': (2120, 2250),
                       'HA': (1650, 1800), 'NP': (1480, 1580), 'NA': (1250, 1560),
                       'M': (975, 1030), 'NS': (815, 900)}
    for header, seq in clean_seqs.items():
        segment = header.split('|')[1]
        segment_length = len(seq)
        min_length = segment_lengths[segment][0]
        max_length = segment_lengths[segment][1]
        if not (min_length <= segment_length <= max_length):
            log.error(f'The consensus sequence generated for segment '
                      f'{segment} is {segment_length}, which is not within the expected length range '
                      f'({min_length} to {max_length} bases).\n')
            exit(22)

    timestamp_analysis_complete = datetime.datetime.now().isoformat()
    analysis_summary['process_name'] = 'make_consensus_seqs'
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete
    analysis_summary['return_code'] = return_code
    analysis_summary['outputs'] = outputs

    return analysis_summary
