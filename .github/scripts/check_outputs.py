#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import os

from pathlib import Path

def check_expected_files_exist(output_dirs, pipeline_versions, sample_ids, output_file_mapping_by_sample_id):
    """
    Check that the expected files exist in the output directory.

    :param output_dirs: Dictionary with keys ['upstream', 'downstream'] and values as the output directories
    :type output_dirs: Dict[str, Path]
    :param pipeline_versions: Dict with keys ['upstream', 'downstream'] and values as the pipeline versions
    :type pipeline_versions: Dict[str, str]
    :param sample_ids: List of sample IDs
    :type sample_ids: List[str]
    :param output_file_mapping_by_sample_id: Dictionary with keys as sample IDs
                                             and values as dictionaries.
    
    :return: List of dictionaries with keys ['sample_id', 'file_type', 'upstream_path', 'downstream_path', 'upstream_exists', 'downstream_exists']
    :rtype: List[Dict[str, str]]
    """
    expected_file_checks = []
    upstream_pipeline = pipeline_versions['upstream']
    downstream_pipeline = pipeline_versions['downstream']
    for sample_id, output_files in output_file_mapping_by_sample_id.items():
        for file_type, paths_by_pipeline in output_files.items():
            upstream_path = os.path.join(output_dirs['upstream'], paths_by_pipeline[upstream_pipeline])
            downstream_path = os.path.join(output_dirs['downstream'], paths_by_pipeline[downstream_pipeline])
            expected_file_check = {
                'sample_id': sample_id,
                'file_type': file_type,
                'upstream_path': upstream_path,
                'downstream_path': downstream_path,
                'upstream_exists': os.path.exists(upstream_path),
                'downstream_exists': os.path.exists(downstream_path),
                'both_exist': os.path.exists(upstream_path) and os.path.exists(downstream_path)
            }
            expected_file_checks.append(expected_file_check)

    return expected_file_checks


def check_expected_md5sums_match(output_dir, pipeline_version, sample_ids):
    """
    Check that the expected md5sums match the actual md5sums in the output directory.
    """
    pass


def main(args):

    os.makedirs(args.outdir, exist_ok=True)

    # TODO: read this from the 'reads_to_simulate.csv' file
    sample_ids = [
        'MK58361X-H3N2'
    ]
    output_file_mapping_by_sample_id = {}
    for sample_id in sample_ids:
        output_file_mapping = {
            'HA_contigs': {"kevinkuchinski": os.path.join(sample_id, "HA_contigs.fa"),
                           "bccdc-phl":      os.path.join(sample_id, "HA_contigs.fa")},
            'HA_contigs_alignment': {"kevinkuchinski": os.path.join(sample_id, "HA_contigs.afa"),
                                     "bccdc-phl":  os.path.join(sample_id, "HA_contigs.afa")},
            'NA_contigs': {"kevinkuchinski": os.path.join(sample_id, "NA_contigs.fa"),
                           "bccdc-phl":      os.path.join(sample_id, "NA_contigs.fa")},
            'NA_contigs_alignment': {"kevinkuchinski": os.path.join(sample_id, "NA_contigs.afa"),
                                     "bccdc-phl":  os.path.join(sample_id, "NA_contigs.afa")},
            'NP_contigs': {"kevinkuchinski": os.path.join(sample_id, "NP_contigs.fa"),
                           "bccdc-phl":      os.path.join(sample_id, "NP_contigs.fa")},
            'NP_contigs_alignment': {"kevinkuchinski": os.path.join(sample_id, "NP_contigs.afa"),
                                     "bccdc-phl":      os.path.join(sample_id, "NP_contigs.afa")},
            'depth_of_cov_samtools': {"kevinkuchinski": os.path.join(sample_id, "depth_of_cov_samtools.tsv"),
                                      "bccdc-phl":      os.path.join(sample_id, "depth_of_cov_samtools.tsv")},
            'depth_of_cov_freebayes': {"kevinkuchinski": os.path.join(sample_id, "depth_of_cov_freebayes.tsv"),
                                       "bccdc-phl":      os.path.join(sample_id, "depth_of_cov_freebayes.tsv")},
            'normalized_reads_r1': {"kevinkuchinski": os.path.join(sample_id, "R1.fq"),
                                    "bccdc-phl":      os.path.join(sample_id, "R1.fq")},
            'normalized_reads_r2': {"kevinkuchinski": os.path.join(sample_id, "R2.fq"),
                                    "bccdc-phl":      os.path.join(sample_id, "R2.fq")},
        }
        output_file_mapping_by_sample_id[sample_id] = output_file_mapping

    pipeline_versions = {
        "upstream": "kevinkuchinski",
        "downstream": "bccdc-phl"
    }

    pipeline_outdirs = {
        "upstream": args.analysis_outdir_kevinkuchinski,
        "downstream": args.analysis_outdir_bccdc_phl
    }

    expected_files_exist_checks = check_expected_files_exist(
        pipeline_outdirs,
        pipeline_versions,
        sample_ids,
        output_file_mapping_by_sample_id
    )
    expected_outputs_exist_output_path = os.path.join(args.outdir, "check_outputs_exist.csv")
    with open(expected_outputs_exist_output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=expected_files_exist_checks[0].keys(), extrasaction='ignore')
        writer.writeheader()
        for check in expected_files_exist_checks:
            writer.writerow(check)
              
    all_expected_files_exist = all([check['upstream_exists'] and check['downstream_exists'] for check in expected_files_exist_checks])

    # TODO: Add more tests
    tests = [
        {
            "test_name": "all_expected_files_exist",
            "test_passed": all_expected_files_exist,
        },
        {
            "test_name": "all_expected_md5sums_match",
            "test_passed": False,
        },
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    output_path = os.path.join(args.outdir, "check_outputs_summary.csv")
    with open(output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=output_fields, extrasaction='ignore')
        writer.writeheader()
        for test in tests:
            if test["test_passed"]:
                test["test_result"] = "PASS"
            else:
                test["test_result"] = "FAIL"
            writer.writerow(test)

    for test in tests:
        if not test['test_passed']:
            exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check outputs')
    parser.add_argument('--analysis-outdir-kevinkuchinski', type=str, help='Path to the pipeline output directory for the kkuchinski version of fluviewer')
    parser.add_argument('--analysis-outdir-bccdc-phl', type=str, help='Path to the pipeline output directory for the bccdc-phl version of fluviewer')
    parser.add_argument('-o', '--outdir', type=str, help='Path to the directory where the output files will be written')
    args = parser.parse_args()
    main(args)
