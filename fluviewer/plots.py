import logging
import os
import datetime

from math import ceil, log10

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

from fluviewer.analysis import run


log = logging.getLogger(__name__)

def make_plots(inputs, outdir, out_name):
    """
    Generate depth of coverage plots for each segment.

    :param inputs: dictionary of input files.
    :type inputs: dict
    :param outdir: Path to the output directory.
    :type outdir: str
    :param out_name: Name used for outputs.
    :type out_name: str
    """
    log.info('Making depth of coverage plots...')
    log_dir = os.path.join(outdir, out_name)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    timestamp_analysis_start = datetime.datetime.now().isoformat()
    analysis_summary = {
        'timestamp_analysis_start': timestamp_analysis_start,
        'inputs': inputs,
    }

    # Get depth of coverage using samtools.
    filtered_alignment = inputs['alignment']
    depth_of_cov = os.path.join(outdir, f'{out_name}_depth_of_cov_samtools.tsv')
    terminal_command = (f'samtools depth {filtered_alignment} > {depth_of_cov}')
    process_name = 'samtools_depth'
    error_code = 24
    return_code = run(terminal_command, outdir, out_name, process_name, error_code)
    if return_code != 0:
        log.error('Error running samtools depth.')

    cols = ['seq_name', 'position', 'depth_samtools']
    samtools_data = pd.read_csv(depth_of_cov, sep='\t', names=cols)

    # Get depth of coverage from FreeBayes.
    depth_of_cov = inputs['depth_of_cov_freebayes']
    cols = ['seq_name', 'position', 'depth_freebayes']
    freebayes_data = pd.read_csv(depth_of_cov, sep='\t', names=cols)

    # Merge samtools and freebayes data.
    data = pd.merge(samtools_data, freebayes_data, on=['seq_name', 'position'], how='left')

    # Annotate with segment.
    get_segment = lambda row: row['seq_name'].split('|')[1]
    data['segment'] = data.apply(get_segment, axis=1)

    # Make plots.
    sb.set_style('whitegrid')
    segments = data['segment'].unique()
    fig_size = (8.5, 2 * len(segments))
    fig, axs = plt.subplots(len(segments), 1, sharex=True, figsize=fig_size)
    max_position = data['position'].max()
    x_ticks = [100 * i for i in range(1, ceil(max_position / 100))]
    max_depth = max([data['depth_samtools'].max(), data['depth_freebayes'].max()])
    y_max = 10 ** ceil(log10(max_depth))
    y_ticks = [10 ** i for i in range(ceil(log10(max_depth)))]
    for segment, ax in zip(segments, axs):
        segment_data = data[data['segment']==segment]
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
        low_cov_path = inputs['low_coverage_positions']
        with open(low_cov_path, 'r') as input_file:
            for line in input_file:
                seq_name, start, stop = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    for position in range(int(start), int(stop) + 1):
                        ax.axvline(position, color='red', alpha = 0.5)
        ambig_path = inputs['ambiguous_positions']
        with open(ambig_path, 'r') as input_file:
            for line in input_file:
                seq_name, position = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvline(int(position), color='orange', alpha = 0.5)
        variant_positions_path = inputs['variant_positions']
        with open(variant_positions_path, 'r') as input_file:
            for line in input_file:
                seq_name, position = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvline(int(position), color='blue', alpha = 0.5)        
    plt.tight_layout()
    plots = os.path.join(outdir, f'{out_name}_depth_of_cov.png')
    plt.savefig(plots, dpi=400)
    plt.close()

    outputs = {'depth_of_cov_plots': plots}
    analysis_summary['outputs'] = outputs
    analysis_summary['timestamp_analysis_complete'] = datetime.datetime.now().isoformat()
    
    log.info(f'Depth of coverage plots saved to {plots}')

    return analysis_summary
