#!/usr/bin/env python3

import argparse
import csv
import glob
import hashlib
import json
import os
import sys

from pathlib import Path


def check_expected_files_exist(output_dirs, sample_ids, output_file_mapping_by_sample_id):
    """
    Check that the expected files exist in the output directory.

    :param output_dirs: Dictionary with keys ['upstream', 'downstream'] and values as the output directories
    :type output_dirs: Dict[str, Path]
    :param sample_ids: List of sample IDs
    :type sample_ids: List[str]
    :param output_file_mapping_by_sample_id: Dictionary with keys as sample IDs
                                             and values as dictionaries.
    
    :return: List of dictionaries with keys ['sample_id', 'file_type', 'upstream_path', 'origin_path', 'upstream_exists', 'origin_exists']
    :rtype: List[Dict[str, str]]
    """
    expected_file_checks = []
    for sample_id, output_files in output_file_mapping_by_sample_id.items():
        for file_type, paths_by_pipeline in output_files.items():
            upstream_path = os.path.join(output_dirs['upstream'], paths_by_pipeline['upstream'])
            origin_path = os.path.join(output_dirs['origin'], paths_by_pipeline['origin'])
            expected_file_check = {
                'sample_id': sample_id,
                'file_type': file_type,
                'upstream_path': upstream_path,
                'origin_path': origin_path,
                'upstream_exists': os.path.exists(upstream_path),
                'origin_exists': os.path.exists(origin_path),
                'both_exist': os.path.exists(upstream_path) and os.path.exists(origin_path)
            }
            expected_file_checks.append(expected_file_check)

    return expected_file_checks


def check_expected_md5sums_match(
        output_dirs,
        sample_ids,
        output_file_mapping_by_sample_id,
        files_with_header_added_in_origin
):
    """
    Check that the expected md5sums match the actual md5sums in the output directory.

    :param output_dirs: Dictionary with keys ['upstream', 'downstream'] and values as the output directories
    :type output_dirs: Dict[str, Path]
    :param sample_ids: List of sample IDs
    :type sample_ids: List[str]
    :param output_file_mapping_by_sample_id: Dictionary with keys as sample IDs
                                                and values as dictionaries.
    :type output_file_mapping_by_sample_id: Dict[str, Dict[str, Dict[str, str]]]
    :param files_with_header_added_in_origin: List of file types for which the header is added in the origin version of the file
    :type files_with_header_added_in_origin: List[str]
    :return: List of dictionaries with keys ['sample_id', 'file_type', 'upstream_path', 'origin_path', 'upstream_md5sum', 'origin_md5sum', 'md5sums_match']
    """
    expected_md5sum_checks = []
    for sample_id, output_files in output_file_mapping_by_sample_id.items():
        for file_type, paths_by_pipeline in output_files.items():
            upstream_path = os.path.join(output_dirs['upstream'], paths_by_pipeline['upstream'])
            origin_path = os.path.join(output_dirs['origin'], paths_by_pipeline['origin'])
            upstream_md5sum = None
            origin_md5sum = None
            with open(upstream_path, 'r') as f:
                upstream_md5sum = hashlib.md5(f.read().encode()).hexdigest()
            with open(origin_path, 'r') as f:
                # skip header when calculating checksum for
                # files where header is added in the origin version
                
                if file_type in files_with_header_added_in_origin:
                    f.readline()
                
                origin_md5sum = hashlib.md5(f.read().encode()).hexdigest()
                    
            
            expected_md5sum_check = {
                'sample_id': sample_id,
                'file_type': file_type,
                'upstream_path': upstream_path,
                'origin_path': origin_path,
                'upstream_md5sum': upstream_md5sum,
                'origin_md5sum': origin_md5sum,
                'md5sums_match': upstream_md5sum == origin_md5sum
            }
            expected_md5sum_checks.append(expected_md5sum_check)

    return expected_md5sum_checks
    


def main(args):

    os.makedirs(args.outdir, exist_ok=True)

    # TODO: read this from the 'reads_to_simulate.csv' file
    sample_ids = [
        'MK58361X-H3N2'
    ]
    analysis_stages = [
        '00_normalize_depth',
        '01_assemble_contigs',
        '02_blast_contigs',
        '03_scaffolding',
        '04_read_mapping',
        '05_variant_calling',
        '06_report_variants',
        '07_summary_reporting'
    ]
    output_file_mapping_by_sample_id = {}
    for sample_id in sample_ids:
        output_file_mapping = {
            'HA_contigs':             {"upstream": os.path.join(sample_id, "HA_contigs.fa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "HA_contigs.fa")},
            'HA_contigs_alignment':   {"upstream": os.path.join(sample_id, "HA_contigs.afa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "HA_contigs.afa")},
            'NA_contigs':             {"upstream": os.path.join(sample_id, "NA_contigs.fa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "NA_contigs.fa")},
            'NA_contigs_alignment':   {"upstream": os.path.join(sample_id, "NA_contigs.afa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "NA_contigs.afa")},
            'NP_contigs':             {"upstream": os.path.join(sample_id, "NP_contigs.fa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "NP_contigs.fa")},
            'NP_contigs_alignment':   {"upstream": os.path.join(sample_id, "NP_contigs.afa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "NP_contigs.afa")},
            'PA_contigs':             {"upstream": os.path.join(sample_id, "PA_contigs.fa"),
                                       "origin":  os.path.join(sample_id,
                                                               "analysis_by_stage",
                                                               "03_scaffolding",
                                                               "PA_contigs.fa")},
            'PA_contigs_alignment':   {"upstream": os.path.join(sample_id, "PA_contigs.afa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "PA_contigs.afa")},
            'PB1_contigs':            {"upstream": os.path.join(sample_id, "PB1_contigs.fa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "PB1_contigs.fa")},
            'PB1_contigs_alignment':  {"upstream": os.path.join(sample_id, "PB1_contigs.afa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "PB1_contigs.afa")},
            'PB2_contigs':            {"upstream": os.path.join(sample_id, "PB2_contigs.fa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "PB2_contigs.fa")},
            'PB2_contigs_alignment':  {"upstream": os.path.join(sample_id, "PB2_contigs.afa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                "PB2_contigs.afa")},
            'normalized_reads_r1':    {"upstream": os.path.join(sample_id, "R1.fq"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "00_normalize_depth",
                                                                f"{sample_id}-normalized_R1.fastq.gz")},
            'normalized_reads_r2':    {"upstream": os.path.join(sample_id, "R2.fq"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "00_normalize_depth",
                                                                f"{sample_id}-normalized_R2.fastq.gz")},
            'assembly_contigs':       {"upstream": os.path.join(sample_id, "spades_output", "contigs.fasta"),
                                        "origin":   os.path.join(sample_id,
                                                                 "analysis_by_stage",
                                                                 "01_assemble_contigs",
                                                                 f"{sample_id}_contigs.fasta")},
            'alignment_sam':          {"upstream": os.path.join(sample_id, "alignment.sam"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "04_read_mapping",
                                                                f"{sample_id}_alignment.sam")},
            'ambig_tsv':              {"upstream": os.path.join(sample_id, "ambig.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "05_variant_calling",
                                                                f"{sample_id}_ambig.tsv")},
            'contigs_blast':          {"upstream": os.path.join(sample_id, "contigs_blast.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "02_blast_contigs",
                                                                f"{sample_id}_contigs_blast.tsv")},
            'mapping_refs':           {"upstream": os.path.join(sample_id, f"{sample_id}_mapping_refs.fa"),
                                        "origin":   os.path.join(sample_id,
                                                                 "analysis_by_stage",
                                                                 "04_read_mapping",
                                                                 f"{sample_id}_mapping_refs.fa")},
            'depth_of_cov_freebayes': {"upstream": os.path.join(sample_id, "depth_of_cov_freebayes.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "05_variant_calling",
                                                                f"{sample_id}_depth_of_cov_freebayes.tsv")},
            'depth_of_cov_samtools':  {"upstream": os.path.join(sample_id, "depth_of_cov_samtools.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "07_summary_reporting",
                                                                f"{sample_id}_depth_of_cov_samtools.tsv")},
            'low_cov':                {"upstream": os.path.join(sample_id, "low_cov.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "05_variant_calling",
                                                                f"{sample_id}_low_cov.tsv")},
            'masked_bed':             {"upstream": os.path.join(sample_id, "masked.bed"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "05_variant_calling",
                                                                f"{sample_id}_masked.bed")},
            'pileup_vcf':             {"upstream": os.path.join(sample_id, "pileup.vcf"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "05_variant_calling",
                                                                f"{sample_id}_variants_raw.vcf")},
            'reads_mapped_tsv':       {"upstream": os.path.join(sample_id, "reads_mapped.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "07_summary_reporting",
                                                                f"{sample_id}_reads_mapped.tsv")},
            'scaffolds_fasta':        {"upstream": os.path.join(sample_id, "scaffolds.fa"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                f"{sample_id}_scaffolds.fa")},
            'scaffolds_blast_tsv':    {"upstream": os.path.join(sample_id, "scaffolds_blast.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "03_scaffolding",
                                                                f"{sample_id}_scaffolds_blast.tsv")},
            'variants_tsv':           {"upstream": os.path.join(sample_id, "variants.tsv"),
                                       "origin":   os.path.join(sample_id,
                                                                "analysis_by_stage",
                                                                "05_variant_calling",
                                                                f"{sample_id}_variants.tsv")},
            
        }
        output_file_mapping_by_sample_id[sample_id] = output_file_mapping

    pipeline_outdirs = {
        "upstream": args.analysis_outdir_upstream,
        "origin": args.analysis_outdir_origin
    }

    expected_files_exist_checks = check_expected_files_exist(
        pipeline_outdirs,
        sample_ids,
        output_file_mapping_by_sample_id
    )
    expected_outputs_exist_output_path = os.path.join(args.outdir, "check_outputs_exist.csv")
    with open(expected_outputs_exist_output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=expected_files_exist_checks[0].keys(), extrasaction='ignore')
        writer.writeheader()
        for check in expected_files_exist_checks:
            writer.writerow(check)
              
    all_expected_files_exist = all([check['upstream_exists'] and check['origin_exists'] for check in expected_files_exist_checks])

    files_whose_md5sums_are_not_expected_to_match = [
        'normalized_reads_r1',
        'normalized_reads_r2',
        'alignment_sam',
        'pileup_vcf',
    ]
    files_with_header_added_in_origin = [
        'contigs_blast',
        'scaffolds_blast_tsv',
        
    ]
    output_file_mapping_by_sample_id_for_md5sum_check = {}
    for sample_id, output_files in output_file_mapping_by_sample_id.items():
        for file_type, paths_by_pipeline in output_files.items():
            if file_type in files_whose_md5sums_are_not_expected_to_match:
                continue
            if sample_id not in output_file_mapping_by_sample_id_for_md5sum_check:
                output_file_mapping_by_sample_id_for_md5sum_check[sample_id] = {}
            output_file_mapping_by_sample_id_for_md5sum_check[sample_id][file_type] = paths_by_pipeline

    expected_md5sums_match_checks = check_expected_md5sums_match(
        pipeline_outdirs,
        sample_ids,
        output_file_mapping_by_sample_id_for_md5sum_check,
        files_with_header_added_in_origin
    )
    def get_analysis_stage_number_from_origin_path(origin_path):
        return origin_path.split(os.sep)[-2]
    expected_md5sums_match_checks_sorted = sorted(expected_md5sums_match_checks, key=lambda x: get_analysis_stage_number_from_origin_path(x['origin_path']))
    
    expected_md5sums_match_output_path = os.path.join(args.outdir, "check_md5sums_match.csv")
    with open(expected_md5sums_match_output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=expected_md5sums_match_checks[0].keys(), extrasaction='ignore')
        writer.writeheader()
        for check in expected_md5sums_match_checks_sorted:
            writer.writerow(check)
    all_expected_md5sums_match = all([check['md5sums_match'] for check in expected_md5sums_match_checks])
    
    # TODO: Add more tests
    tests = [
        {
            "test_name": "all_expected_files_exist",
            "test_passed": all_expected_files_exist,
        },
        {
            "test_name": "all_expected_md5sums_match",
            "test_passed": all_expected_md5sums_match,
        },
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    output_path = os.path.join(args.outdir, "check_outputs_summary.csv")
    with open(output_path, 'w') as f:
        file_writer = csv.DictWriter(f, fieldnames=output_fields, extrasaction='ignore')
        stdout_writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, extrasaction='ignore', delimiter='\t')
        stdout_writer.writeheader()
        file_writer.writeheader()
        for test in tests:
            test["test_result"] = "FAIL"
            if test["test_passed"]:
                test["test_result"] = "PASS"

            stdout_writer.writerow(test)
            file_writer.writerow(test)

    for test in tests:
        if not test['test_passed']:
            print(f"\nTest: {test['test_name']} failed.")
            print(f"See {output_path} for more details.")
            exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check outputs')
    parser.add_argument('--analysis-outdir-upstream', type=str, help='Path to the pipeline output directory for the upstream (KevinKuchinski) version of FluViewer')
    parser.add_argument('--analysis-outdir-origin', type=str, help='Path to the pipeline output directory for the origin (BCCDC-PHL) version of FluViewer')
    parser.add_argument('-o', '--outdir', type=str, help='Path to the directory where the output files will be written')
    args = parser.parse_args()
    main(args)
