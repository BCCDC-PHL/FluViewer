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
from . import database
from . import analysis
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
    

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        log.info(f'Created output directory: {args.outdir}')
    else:
        log.error(f'Output directory already exists: {args.outdir}')
        exit(1)

    database.check_database(
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

    #
    # Stage 0: Normalize depth of reads.
    #
    current_analysis_stage = 'normalize_depth'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')
    log.info(f'Output directory: {current_analysis_stage_outdir}')

    current_analysis_stage_inputs = {
        'input_reads_fwd': os.path.abspath(args.forward_reads),
        'input_reads_rev': os.path.abspath(args.reverse_reads),
    }

    normalize_depth_analysis_summary = analysis.normalize_depth(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
        os.path.abspath(args.forward_reads),
        os.path.abspath(args.reverse_reads),
        args.target_depth,
        args.max_memory,
    )
    if normalize_depth_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {normalize_depth_analysis_summary["return_code"]}')
        exit(normalize_depth_analysis_summary['return_code'])

    outputs_to_publish = {
        'normalized_reads_fwd': os.path.join(args.outdir),
        'normalized_reads_rev': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = normalize_depth_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')


    #
    # Stage 1: Assemble contigs.
    #
    current_analysis_stage = 'assemble_contigs'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    
    current_analysis_stage_inputs = {
        'normalized_reads_fwd': normalize_depth_analysis_summary['outputs']['normalized_reads_fwd'],
        'normalized_reads_rev': normalize_depth_analysis_summary['outputs']['normalized_reads_rev'],
    }
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)

    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    assemble_contigs_analysis_summary = analysis.assemble_contigs(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
    )
    if assemble_contigs_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {assemble_contigs_analysis_summary["return_code"]}')
        exit(assemble_contigs_analysis_summary['return_code'])

    outputs_to_publish = {
        'contigs': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = assemble_contigs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')


    log.info(f'Analysis stage complete: {current_analysis_stage}')

    #
    # Stage 2: Blast contigs.
    #
    current_analysis_stage = 'blast_contigs'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_inputs = {
        'contigs': assemble_contigs_analysis_summary['outputs']['contigs'],
        'database': os.path.abspath(args.database),
    }
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    blast_contigs_analysis_summary = analysis.blast_contigs(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
        args.threads,
        args.min_identity,
        args.min_alignment_length,
    )
    if blast_contigs_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {blast_contigs_analysis_summary["return_code"]}')
        exit(blast_contigs_analysis_summary['return_code'])

    outputs_to_publish = {
        'all_contig_blast_results': os.path.join(args.outdir),
        'filtered_contig_blast_results': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = blast_contigs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    #
    # Stage 3: Scaffolding.
    #
    current_analysis_stage = 'scaffolding'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_inputs = {
        'filtered_contig_blast_results': blast_contigs_analysis_summary['outputs']['filtered_contig_blast_results'],
    }
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
        
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    make_scaffold_seqs_analysis_summary = analysis.make_scaffold_seqs(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
    )
    if make_scaffold_seqs_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {make_scaffold_seqs_analysis_summary["return_code"]}')
        exit(make_scaffold_seqs_analysis_summary['return_code'])

    blast_scaffolds_inputs = {
        'scaffolds': make_scaffold_seqs_analysis_summary['outputs']['scaffolds'],
        'database': os.path.abspath(args.database),
    }
    blast_scaffolds_analysis_summary = analysis.blast_scaffolds(
        blast_scaffolds_inputs,
        args.outdir,
        args.output_name,
        args.threads,
    )

    filtered_scaffold_blast_results = filter_scaffold_blast_results(scaffold_blast_results)

    outputs_to_publish = {
        'scaffolds': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = make_scaffold_seqs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')
    
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
    )

    mask_ambig_low_cov(
        args.outdir,
        args.output_name,
        args.min_depth,
        args.variant_threshold_calling,
        args.variant_threshold_masking,
        args.min_mapping_quality,
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
        for stage in analysis_stages:
            stage_index = analysis_stages.index(stage)
            stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{stage_index:02}_{stage}')
            shutil.rmtree(stage_outdir)
            log.info(f'Removed intermediate files for analysis stage: {stage}')
    else:
        log.info('Garbage collection disabled. Intermediate files were not removed.')

    log.info('Analysis complete.')
    exit(0)







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
