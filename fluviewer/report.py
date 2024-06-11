import logging
import os
import datetime

import pandas as pd
import numpy as np

from fluviewer.analysis import run


log = logging.getLogger(__name__)


def write_report(inputs, outdir, out_name):
    """
    Generate a report for each segment.

    :param inputs: dictionary of input files.
    :type inputs: list
    :param outdir: Path to the output directory.
    :type outdir: Path
    :param out_name: Name used for outputs.
    :type out_name: str
    """
    log.info('Writing report...')
    log_dir = os.path.join(outdir, out_name)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    timestamp_analysis_start = datetime.datetime.now().isoformat()
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }

    outputs = {}

    # Count reads mapped to each segment and add to report.
    filtered_alignment = inputs['alignment']
    reads_mapped = os.path.join(outdir, f'{out_name}_reads_mapped.tsv')
    terminal_command = (f'samtools idxstats {filtered_alignment} > '
                        f'{reads_mapped}')
    process_name = 'samtools_idxstats'
    error_code = 23
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error('Error running samtools idxstats.')
        analysis_summary['error'] = 'Error running samtools idxstats.'
        analysis_summary['return_code'] = error_code
        return analysis_summary

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
    scaffolds = inputs['scaffolds']
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
    consensus_seqs = inputs['consensus_seqs']
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
    mapping_refs = inputs['mapping_refs']
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
    report_path = os.path.join(outdir, f'{out_name}_report.tsv')
    report.to_csv(report_path, index=False, sep='\t')

    outputs['report'] = report_path

    analysis_summary['outputs'] = outputs
    timestamp_analysis_complete = datetime.datetime.now().isoformat()
    analysis_summary['timestamp_analysis_complete'] = timestamp_analysis_complete

    return analysis_summary
    
