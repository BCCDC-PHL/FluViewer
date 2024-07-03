import logging
import os
import datetime

from math import ceil, log10
from pathlib import Path

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

from fluviewer import parsers
from fluviewer.analysis import run

log = logging.getLogger(__name__)


def make_depth_of_coverage_plots(
        alignment: Path,
        depth_of_cov_freebayes: Path,
        low_coverage_positions: Path,
        ambiguous_positions: Path,
        variant_positions: Path,
        outdir: Path,
        out_name: str,
) -> None:
    """
    Generate depth of coverage plots for each segment.

    :param alignment: Path to the alignment (bam) file.
    :type alignment: Path
    :param depth_of_cov_freebayes: Path to the depth of coverage file from FreeBayes.
    :type depth_of_cov_freebayes: Path
    :param low_coverage_positions: Path to the low coverage positions file.
    :type low_coverage_positions: Path
    :param ambiguous_positions: Path to the ambiguous positions file.
    :type ambiguous_positions: Path
    :param variant_positions: Path to the variant positions file.
    :type variant_positions: Path
    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name used for outputs.
    :type out_name: str
    :return: None
    :rtype: None
    """
    log.info('Making depth of coverage plots...')
    log_dir = os.path.join(outdir, out_name)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Get depth of coverage using samtools.
    depth_of_cov_samtools = os.path.join(outdir, f'{out_name}_depth_of_cov_samtools.tsv')
    terminal_command = (f'samtools depth {alignment} > {depth_of_cov_samtools}')
    process_name = 'samtools_depth'
    error_code = 24
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error('Error running samtools depth.')

    cols = ['seq_name', 'position', 'depth_samtools']
    samtools_data = pd.read_csv(depth_of_cov_samtools, sep='\t', names=cols)

    # Get depth of coverage from FreeBayes.
    cols = ['seq_name', 'position', 'depth_freebayes']
    freebayes_data = pd.read_csv(depth_of_cov_freebayes, sep='\t', names=cols)

    # Merge samtools and freebayes data.
    depth_of_cov_data = pd.merge(samtools_data, freebayes_data, on=['seq_name', 'position'], how='left')

    # Annotate with segment.
    get_segment = lambda row: row['seq_name'].split('|')[1]
    depth_of_cov_data['segment'] = depth_of_cov_data.apply(get_segment, axis=1)

    # Make plots.
    sb.set_style('whitegrid')
    segments = depth_of_cov_data['segment'].unique()
    fig_size = (8.5, 2 * len(segments))
    fig, axs = plt.subplots(len(segments), 1, sharex=True, figsize=fig_size)
    max_position = depth_of_cov_data['position'].max()
    x_ticks = [100 * i for i in range(1, ceil(max_position / 100))]
    max_depth = max([depth_of_cov_data['depth_samtools'].max(), depth_of_cov_data['depth_freebayes'].max()])
    y_max = 10 ** ceil(log10(max_depth))
    y_ticks = [10 ** i for i in range(ceil(log10(max_depth)))]
    for segment, ax in zip(segments, axs):
        segment_data = depth_of_cov_data[depth_of_cov_data['segment']==segment]
        sb.lineplot(x='position', y='depth_samtools', ax=ax, data=segment_data,
                    color='grey')
        sb.lineplot(x='position', y='depth_freebayes', ax=ax, data=segment_data,
                    color='black')
        ax.set_xlim(1, max_position)
        ax.set_xlabel('Position')
        ax.set_xticks(x_ticks)
        if ax == axs[-1]:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        else:
            ax.set_xticklabels([])
        ax.set_ylim(1, y_max)
        ax.set_ylabel('Read depth')
        ax.set_yscale('log')
        ax.set_yticks(y_ticks)
        ax.set_title(segment)
        ax.axvspan(segment_data['position'].max(), max_position, color='grey')

        with open(low_coverage_positions, 'r') as input_file:
            for line in input_file:
                seq_name, start, stop = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    for position in range(int(start), int(stop) + 1):
                        ax.axvline(position, color='red', alpha = 0.5)

        with open(ambiguous_positions, 'r') as input_file:
            for line in input_file:
                seq_name, position = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvline(int(position), color='orange', alpha = 0.5)

        with open(variant_positions, 'r') as input_file:
            for line in input_file:
                seq_name, position = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvline(int(position), color='blue', alpha = 0.5)

    plt.tight_layout()
    plots = os.path.join(outdir, f'{out_name}_depth_of_cov.png')
    plt.savefig(plots, dpi=400)
    plt.close()
    
    log.info(f'Depth of coverage plots saved to {plots}')


def make_scaffolding_plots(alignment: list[dict], scaffold: dict, output: Path) -> None:
    """
    Plot the alignment of contigs to the scaffold.

    :param alignment: A list of dictionaries representing the alignment of contigs to the scaffold.
    :type alignment: list[dict]
    :param scaffold: A dictionary representing the scaffold sequence.
    :type scaffold: dict
    :param output: The output plot file.
    :return: None
    :rtype: None
    """
    num_contigs = len(alignment) + 1 # Add 1 for the scaffold
    max_alignment_seq_len = max([len(contig['seq']) for contig in alignment + [scaffold]])
    height_per_contig = 5
    fig, ax = plt.subplots()
    ax.set_ylim(0, num_contigs * height_per_contig)
    ax.set_xlim(0, max_alignment_seq_len)
    contig_ymin = None
    alignment_sorted = sorted(alignment, key=lambda x: x['id'])
    scaffold_plus_alignments_sorted = [scaffold] + alignment_sorted
    for contig in scaffold_plus_alignments_sorted:
        # Need to specify the y-axis value for the contig
        # As an integer, to position the bar.
        # Evenly space the contigs along the y-axis
        if contig_ymin is None:
            contig_ymin = 0
        else:
            contig_ymin += height_per_contig
        contig_id = contig['id']
        matching_regions = contig['matching_regions']
        # Add the contig_id as labels
        center_contig_y = contig_ymin + height_per_contig / 2
        ax.text(0, center_contig_y, contig_id, ha='right', va='center')
        for matching_region in matching_regions:
            start = matching_region['start']
            end = matching_region['end']
            ax.broken_barh([(start, end - start)], (contig_ymin, height_per_contig), facecolors='tab:blue')
            # Draw horizontal lines to separate contigs
            ax.axhline(contig_ymin, color='black', linewidth=0.5)

    ax.set_xlabel('Alignment position')
    # Turn off numerical labels and ticks on y-axis
    ax.set_yticklabels([])
    ax.yaxis.set_ticks([])
    # Only turn on vertical grid lines
    ax.grid(axis='x')

    plt.savefig(output)
    plt.close()
    


def make_plots(inputs, outdir, out_name):
    """
    Generate plots for the analysis.

    :param inputs: The inputs dictionary. Expected keys are:
                   ['depth_of_cov_freebayes', 'alignment', 'depth_of_cov_freebayes',
                    'low_coverage_positions', 'ambiguous_positions', 'variant_positions',
                    'segment_contigs_alignments', 'scaffolds']
    :type inputs: dict
    :param outdir: The output directory.
    :type outdir: str
    :param out_name: The output name.
    :type out_name: str
    :return: The analysis summary. Keys are ['timestamp_analysis_start', 'inputs', 'outputs', 'timestamp_analysis_complete']
    :rtype: dict
    """
    timestamp_analysis_start = datetime.datetime.now().isoformat()
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }
    outputs = {}
    
    depth_of_cov = inputs['depth_of_cov_freebayes']
    alignment = inputs['alignment']
    depth_of_cov_freebayes = inputs['depth_of_cov_freebayes']
    low_coverage_positions = inputs['low_coverage_positions']
    ambiguous_positions = inputs['ambiguous_positions']
    variant_positions = inputs['variant_positions']
    make_depth_of_coverage_plots(
        alignment,
        depth_of_cov_freebayes,
        low_coverage_positions,
        ambiguous_positions,
        variant_positions,
        outdir,
        out_name
    )
    outputs['depth_of_cov_plots'] = os.path.join(outdir, f'{out_name}_depth_of_cov.png')

    segments = [
        'HA',
        'M',
        'NA',
        'NP',
        'NS',
        'PA',
        'PB1',
        'PB2',
    ]
    outputs['scaffold_alignment_plots'] = {}
    for segment in segments:

        if segment not in inputs['segment_contigs_alignments']:
            log.warning(f'No contig alignment file found for segment {segment}. Skipping...')
            continue

        contig_alignment_path  = inputs['segment_contigs_alignments'][segment]
        if not os.path.exists(contig_alignment_path):
            log.warning(f'Contig alignment file {contig_alignment_path} does not exist. Skipping...')
            continue
        contig_alignment = parsers.parse_contig_alignment(contig_alignment_path)
        scaffolds_path = inputs['scaffolds']
        scaffold = parsers.parse_scaffolds(scaffolds_path, segment)
        output = os.path.join(outdir, f'{out_name}_{segment}_alignment.png')
        log.info(f'Making scaffold alignment plot for segment {segment}...')
        make_scaffolding_plots(contig_alignment, scaffold, output)
        outputs['scaffold_alignment_plots'][segment] = output
        log.info(f'Scaffold alignment plot saved to {output}')
    
    analysis_summary['outputs'] = outputs
    analysis_summary['timestamp_analysis_complete'] = datetime.datetime.now().isoformat()

    return analysis_summary
