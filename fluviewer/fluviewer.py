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
        current_analysis_stage_outdir,
        args.output_name,
        args.threads,
    )

    outputs_to_publish = {
        'scaffolds': os.path.join(args.outdir),
        'filtered_scaffold_blast_results': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        if output_name in make_scaffold_seqs_analysis_summary['outputs']:
            src_path = make_scaffold_seqs_analysis_summary['outputs'][output_name]
        elif output_name in blast_scaffolds_analysis_summary['outputs']:
            src_path = blast_scaffolds_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')
    
    log.info(f'Analysis stage complete: {current_analysis_stage}')


    #
    # Stage 4: Read mapping.
    #
    current_analysis_stage = 'read_mapping'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
    current_analysis_stage_inputs = {
        'filtered_scaffold_blast_results': blast_scaffolds_analysis_summary['outputs']['filtered_scaffold_blast_results'],
        'database': os.path.abspath(args.database),
        'normalized_reads_fwd': normalize_depth_analysis_summary['outputs']['normalized_reads_fwd'],
        'normalized_reads_rev': normalize_depth_analysis_summary['outputs']['normalized_reads_rev'],
    }
    log.info(f'Beginning analysis stage: {current_analysis_stage}')
        
    map_reads_analysis_summary = analysis.map_reads(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
        args.min_mapping_quality,
    )
    if map_reads_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {map_reads_analysis_summary["return_code"]}')
        exit(map_reads_analysis_summary['return_code'])

    outputs_to_publish = {
        'mapping_refs': os.path.join(args.outdir),
        'alignment': os.path.join(args.outdir),
        'alignment_index': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = map_reads_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')


    #
    # Stage 5: Variant calling.
    #
    current_analysis_stage = 'variant_calling'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
    current_analysis_stage_inputs = {
        'mapping_refs': map_reads_analysis_summary['outputs']['mapping_refs'],
        'alignment': map_reads_analysis_summary['outputs']['alignment'],
        'alignment_index': map_reads_analysis_summary['outputs']['alignment_index'],
    }
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    variant_calling_params = {
        'variant_threshold_calling': args.variant_threshold_calling,
        'variant_threshold_masking': args.variant_threshold_masking,
        'min_depth': args.min_depth,
        'min_mapping_quality': args.min_mapping_quality,
        'coverage_limit': args.coverage_limit,
    }
    
    call_variants_analysis_summary = analysis.call_variants(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
        variant_calling_params,
    )
    if call_variants_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {call_variants_analysis_summary["return_code"]}')
        exit(call_variants_analysis_summary['return_code'])

    outputs_to_publish = {
        'variants_raw': os.path.join(args.outdir),
        'variants_filtered': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = call_variants_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')
    print(json.dumps(call_variants_analysis_summary, indent=4))
    exit(0)


    #
    # Stage 6: Consensus calling.
    #
    current_analysis_stage = 'consensus_calling'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
    current_analysis_stage_inputs = {
        'variants_filtered': call_variants_analysis_summary['outputs']['variants_filtered'],
        'mapping_refs': map_reads_analysis_summary['outputs']['mapping_refs'],
    }
    log.info(f'Beginning analysis stage: {current_analysis_stage}')

    make_consensus_seqs_analysis_summary = analysis.make_consensus_seqs(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
    )
    if make_consensus_seqs_analysis_summary['return_code'] != 0:
        log.error(f'Error in analysis stage: {current_analysis_stage}')
        log.error(f'Error code: {make_consensus_seqs_analysis_summary["return_code"]}')
        exit(make_consensus_seqs_analysis_summary['return_code'])

    outputs_to_publish = {
        'consensus_seqs': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = make_consensus_seqs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

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
