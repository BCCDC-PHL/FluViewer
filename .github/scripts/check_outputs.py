#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import os

from pathlib import Path

def check_expected_files_exist(output_dir, sample_ids, pipeline_version):
    """
    Check that the expected files exist in the output directory.

    :param output_dir: Path to the output directory
    :type output_dir: Path
    :param sample_ids: List of sample IDs
    :type sample_ids: List[str]
    :param pipeline_version: Version of the pipeline ('kkuchinski' or 'bccdc-phl')
    :type pipeline_version: str
    :return: True if all expected files exist, False otherwise
    :rtype: bool
    """
    for sample_id in sample_ids:
        expected_files = [
        ]

        for expected_file in expected_files:
            expected_file_path = os.path.join(output_dir, expected_file)
            if not os.path.exists(expected_file_path):
                print(f"Expected file {expected_file_path} not found")
                return False

    return True


def main(args):

    output_dir = os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)

    # TODO: Add more tests
    all_expected_files_exist_checks = []
    for pipeline_version in ['kkuchinski', 'bccdc-phl']:
        all_expected_files_exist = check_expected_files_exist(args.pipeline_outdir, sample_ids, pipeline_version)
        all_expected_files_exist_checks.append(all_expected_files_exist)
        
    tests = [
        {
            "test_name": "all_expected_files_exist",
            "test_passed": all(all_expected_files_exist_checks),
        },
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    output_path = args.output
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
    parser.add_argument('--analysis-outdir-kkuchinski', type=str, help='Path to the pipeline output directory for the kkuchinski version of fluviewer')
    parser.add_argument('--analysis-outdir-bccdc-phl', type=str, help='Path to the pipeline output directory for the bccdc-phl version of fluviewer')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(args)
