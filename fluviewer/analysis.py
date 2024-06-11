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
import fluviewer.logging


log = fluviewer.logging.get_logger(__name__, 'info')


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
    11: 'Error generating consensus sequences.',
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
        fwd_reads_raw: Path,
        rev_reads_raw: Path,
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
    :param fwd_reads_raw: Path to the raw forward reads.
    :type fwd_reads_raw: str
    :param rev_reads_raw: Path to the raw reverse reads.
    :type rev_reads_raw: str
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
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }
    
    normalized_reads_fwd = os.path.join(outdir, f'{out_name}-normalized_R1.fastq')
    normalized_reads_rev = os.path.join(outdir, f'{out_name}-normalized_R2.fastq')

    terminal_command = (f'bbnorm.sh in={input_reads_fwd} in2={input_reads_rev} '
                        f'out={normalized_reads_fwd} out2={normalized_reads_rev} target={depth}')
    terminal_command = (terminal_command + f' -Xmx{max_memory}g'
                        if max_memory is not None else terminal_command)

    # add gzip compression to output files
    terminal_command += f'\n\ngzip -f {normalized_reads_fwd}\n'
    terminal_command += f'\ngzip -f {normalized_reads_rev}\n'

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

    outputs = {
        'normalized_reads_fwd': os.path.abspath(normalized_reads_fwd) + '.gz',
        'normalized_reads_rev': os.path.abspath(normalized_reads_rev) + '.gz',
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
    fwd_reads = inputs.get('normalized_reads_fwd', None)
    rev_reads = inputs.get('normalized_reads_rev', None)

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
    log.info('Filtering contig alignments...')
    total_num_blast_results = len(blast_results)

    filtered_blast_results = blast_results[blast_results['pident'] >= identity]
    filtered_blast_results = filtered_blast_results[filtered_blast_results['length'] >= length]
    num_blast_results_after_identity_and_length_filter = len(filtered_blast_results)
    log.info(f'Found {num_blast_results_after_identity_and_length_filter} matches with at least {identity}% identity and {length} bp alignment length.')

    percent_blast_results_retained = num_blast_results_after_identity_and_length_filter / total_num_blast_results * 100
    log.info(f'Retained {percent_blast_results_retained:.2f}% matches for further analysis.')

    if len(filtered_blast_results) == 0:
        log.error(f'No contigs aligned to reference sequences! Aborting analysis.')
        error_code = 7
        exit(error_code)

    # Re-naming the filtered blast results to blast_results
    # For compatibility with the rest of this function
    blast_results = filtered_blast_results

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
    log.info('Aligning contigs to reference sequences...')
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

    blast_results = pd.read_csv(blast_output_path, sep='\t')

    total_num_blast_results = len(blast_results)
    log.info('Contigs aligned to reference sequences.')
    log.info(f'Found {total_num_blast_results} total matches.')

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


def make_scaffold_seqs(inputs, outdir, out_name):
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
    blast_results = pd.read_csv(filtered_contig_blast_results_path, sep='\t')

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
    for segment in segments:
        contigs = os.path.join(outdir, f'{segment}_contigs.fa')
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
        outputs[f'{segment}_contigs'] = os.path.abspath(contigs)
        # Generate multiple sequence alignments of trimmed/positioned contigs.
        log.info(f'Aligning contigs for segment {segment}...')
        aligned_contigs = os.path.join(outdir, f'{segment}_contigs.afa')

        if contig_counter[segment] > 1:
            log.info(f'Generating multiple sequence alignment for segment {segment}...')
            terminal_command = (f'clustalw -INFILE={contigs} '
                                f'-OUTFILE={aligned_contigs} -OUTPUT=FASTA')
            process_name = f'clustalw_{segment}'
            error_code = 9
            clustalw_return_code = run(terminal_command, outdir, out_name, process_name, error_code)
        else:
            log.info(f'Only one contig for segment {segment}, skipping alignment.')
            shutil.copyfile(contigs, aligned_contigs)

        outputs[f'{segment}_aligned_contigs'] = os.path.abspath(aligned_contigs)

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
            error_code = 10
            analysis_summary['return_code'] = error_code
            analysis_summary['error_message'] = error_messages_by_code[error_code]
            return analysis_summary
        elif len(alignment_lengths) == 0:
            log.error(f'No sequences in the multiple sequence alignment for {segment}! Aborting analysis.\n')
            error_code = 11
            analysis_summary['return_code'] = error_code
            analysis_summary['error_message'] = error_messages_by_code[error_code]
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
    scaffolds = os.path.join(outdir, out_name, 'scaffolds.fa')
    with open(scaffolds, 'w') as f:
        for segment, seq in scaffold_seqs.items():
            header = f'>{out_name}|{segment}_scaffold'
            f.write(header + '\n')
            f.write(seq + '\n')

    outputs['scaffolds'] = os.path.abspath(scaffolds)

    log.info(f'Wrote {len(scaffold_seqs)} scaffolds to {scaffolds}')

    timestamp_analysis_complete = datetime.datetime.now().isoformat()

    analysis_summary['process_name'] = 'make_scaffold_seqs'
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete
    analysis_summary['return_code'] = clustalw_return_code
    analysis_summary['outputs'] = outputs
    
    return analysis_summary
