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
from . import logging_config
from . import database
from . import analysis
from . import report
from . import plots


log = logging.getLogger(__name__)

def main():
    version = __version__

    args = cli_args.parse_args()
    try:
        args = cli_args.validate_args(args)
    except ValueError as e:
        print(f'Error: {e}', sys.stderr)
        exit(1)

    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    else:
        print(f'Output directory already exists: {args.outdir}', sys.stderr)
        exit(1)

    log_level = getattr(logging, args.log_level.upper())
    print(log_level)
    if not isinstance(log_level, int):
        raise ValueError(f'Invalid log level: {args.log_level}')

    logs_dir = os.path.join(args.outdir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    log_file = os.path.join(logs_dir, 'fluviewer.log')
    logging_config.configure_logging(log_level, log_file)

    log.info(f'BCCDC-PHL/FluViewer v{version}')
    version_split = version.split('-')
    log.info(f'Derived from: KevinKuchinski/FluViewer v{version_split[0]}')
    log.info(f'Inputs:')
    log.info(f"Fwd reads: {args.forward_reads}")
    log.info(f"Rev reads: {args.reverse_reads}")
    log.info(f"Reference sequences: {args.db}")

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
    
    
    database.check_database(
        args.db,
        args.outdir,
        args.output_name,
    )    

    analysis_stages = [
        'normalize_depth',
        'assemble_contigs',
        'blast_contigs',
        'scaffolding',
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

    #
    # Publish outputs and logs
    outputs_to_publish = {
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = normalize_depth_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

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

    #
    # Publish outputs and logs
    outputs_to_publish = {
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = assemble_contigs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')


    #
    # Stage 2: Blast contigs.
    #
    current_analysis_stage = 'blast_contigs'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_inputs = {
        'contigs': assemble_contigs_analysis_summary['outputs']['contigs'],
        'database': os.path.abspath(args.db),
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

    #
    # Publish outputs and logs
    outputs_to_publish = {
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = blast_contigs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

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
        'database': os.path.abspath(args.db),
    }

    blast_scaffolds_analysis_summary = analysis.blast_scaffolds(
        blast_scaffolds_inputs,
        current_analysis_stage_outdir,
        args.output_name,
        args.threads,
    )

    #
    # Publish outputs and logs
    outputs_to_publish = {
    }
    for output_name, output_dir in outputs_to_publish.items():
        if output_name in make_scaffold_seqs_analysis_summary['outputs']:
            src_path = make_scaffold_seqs_analysis_summary['outputs'][output_name]
        elif output_name in blast_scaffolds_analysis_summary['outputs']:
            src_path = blast_scaffolds_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')
    
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
        'database': os.path.abspath(args.db),
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

    #
    # Publish outputs and logs
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

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

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

    #
    # Publish outputs and logs
    outputs_to_publish = {
        'variants_filtered': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = call_variants_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')

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
        'masked_positions': call_variants_analysis_summary['outputs']['masked_positions'],
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

    #
    # Publish outputs and logs
    outputs_to_publish = {
        'consensus_seqs': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        src_path = make_consensus_seqs_analysis_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')


    #
    # Stage 7: Summary reporting and plotting.
    #
    current_analysis_stage = 'summary_reporting'
    current_analysis_stage_index = analysis_stages.index(current_analysis_stage)
    current_analysis_stage_outdir = os.path.join(args.outdir, 'analysis_by_stage', f'{current_analysis_stage_index:02}_{current_analysis_stage}')
    current_analysis_stage_outdir = os.path.abspath(current_analysis_stage_outdir)
    current_analysis_stage_inputs = {
        'scaffolds': make_scaffold_seqs_analysis_summary['outputs']['scaffolds'],
        'mapping_refs': map_reads_analysis_summary['outputs']['mapping_refs'],
        'alignment': map_reads_analysis_summary['outputs']['alignment'],
        'depth_of_cov_freebayes': call_variants_analysis_summary['outputs']['depth_of_cov_freebayes'],
        'low_coverage_positions': call_variants_analysis_summary['outputs']['low_coverage_positions'],
        'ambiguous_positions': call_variants_analysis_summary['outputs']['ambiguous_positions'],
        'variant_positions': call_variants_analysis_summary['outputs']['variant_positions'],
        'consensus_seqs': make_consensus_seqs_analysis_summary['outputs']['consensus_seqs'],
    }
    log.info(f'Beginning analysis stage: {current_analysis_stage}')
    
    reporting_summary = report.write_report(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
    )

    plotting_summary = plots.make_plots(
        current_analysis_stage_inputs,
        current_analysis_stage_outdir,
        args.output_name,
    )

    #
    # Publish outputs and logs
    outputs_to_publish = {
        'report': os.path.join(args.outdir),
        'depth_of_cov_plots': os.path.join(args.outdir),
    }
    for output_name, output_dir in outputs_to_publish.items():
        if output_name in reporting_summary['outputs']:
            src_path = reporting_summary['outputs'][output_name]
        elif output_name in plotting_summary['outputs']:
            src_path = plotting_summary['outputs'][output_name]
        dest_path = os.path.join(output_dir, os.path.basename(src_path))
        shutil.copy(src_path, dest_path)
        log.info(f'Published output: {output_name} -> {dest_path}')

    analysis_stage_logs_src_dir = os.path.join(current_analysis_stage_outdir, 'logs')
    analysis_stage_logs_dest_dir = os.path.join(logs_dir, os.path.basename(current_analysis_stage_outdir))
    os.makedirs(analysis_stage_logs_dest_dir)
    for log_file in os.listdir(analysis_stage_logs_src_dir):
        src_path = os.path.join(analysis_stage_logs_src_dir, log_file)
        dest_path = os.path.join(analysis_stage_logs_dest_dir, log_file)
        shutil.copy(src_path, dest_path)
        log.info(f'Published log file: {log_file} -> {dest_path}')

    log.info(f'Analysis stage complete: {current_analysis_stage}')

    
    if not args.disable_garbage_collection:
        analysis_by_stage_dir = os.path.join(args.outdir, 'analysis_by_stage')
        for stage in analysis_stages:
            stage_index = analysis_stages.index(stage)
            stage_outdir = os.path.join(analysis_by_stage_dir, f'{stage_index:02}_{stage}')
            shutil.rmtree(stage_outdir)
            log.info(f'Removed intermediate files for analysis stage: {stage}')
        shutil.rmtree(analysis_by_stage_dir)
    else:
        log.info('Garbage collection disabled. Intermediate files were not removed.')

    log.info('Analysis complete.')
    exit(0)


if __name__ == '__main__':
    main()
