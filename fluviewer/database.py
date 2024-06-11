import os

from collections import Counter
from math import ceil, log10

import fluviewer.logging

log = fluviewer.logging.get_logger(__name__, 'info')


def check_database(db, outdir, out_name):
    """
    Checks the contents of the provided reference sequence database to
    ensure proper header formatting and unambiguous sequences.

    :param db: Path to the reference sequence database.
    :type db: Path
    :param outdir: Path to the output directory.
    :type outdir: str
    :return: None
    :rtype: None
    """
    log.info('Checking reference sequence database...')
    seqs = {}
    headers = []
    with open(db, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                headers.append(header)
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    headers = Counter(headers)
    for header in headers:
        if headers[header] > 1:
            log.error(f'The following header is used for multiple sequences:'
                      f'\n{header}\n')
            exit(1)

    num_seqs = len(seqs)
    one_tenth_of_seqs = ceil(num_seqs / 10)
    log.info(f'Found {num_seqs} sequences in database.')

    log.info('Checking database sequences...')

    segment_lengths = {
        'PB2': (2260, 2360),
        'PB1': (2260, 2360),
        'PA': (2120, 2250),
        'HA': (1650, 1800),
        'NP': (1480, 1580),
        'NA': (1250, 1560),
        'M': (975, 1030),
        'NS': (815, 900),
    }

    log.info('Checking sequence headers')
    log.info('Checking sequence lengths within expected ranges for each segment')
    log.info('Checking for ambiguous or lower-case nucleotides')

    num_seqs_checked = 0
    for header, seq in seqs.items():

        #
        # Check header formatting
        if len(header.split('|')) != 4:
            log.error(f'The header for the following database entry does not contain the expected number of |-delimited fields:\n{header}\n')
            exit(1)
        else:
            accession, name, segment, subtype = header.split('|')
            
            # Check that the strain name is formatted correctly
            # (e.g. "A/California/04/2009(H1N1)")
            if any([name.count('(') != 1, name.count(')') != 1,
                    name[-1] != ')', name.index('(') > name.index(')')]):
                log.error(f'The strain_name(strain_subtype) for the following database entry is improperly formatted:\n{header}\n')
                exit(1)
            if segment not in segment_lengths.keys():
                log.error(f'The segment indicated for the following database entry is not recognized:\n{header}\n')
                exit(1)
        if len(seq) != sum(seq.count(base) for base in 'ATGC'):
            log.error(f'The sequence provided for the following database entry '
                      f'contains ambiguous or lower-case nucleotides:\n{header}\n')
            exit(1)

        #
        # Check sequence lengths
        min_length = segment_lengths[segment][0]
        max_length = segment_lengths[segment][1]
        if not (min_length <= len(seq) <= max_length):
            log.error(f'The sequence provided for the following database entry '
                      f'is not within the expected length range for its indicated '
                      f'segment ({segment}: {min_length} to {max_length} bases): '
                      f'\n{header}\n')
            exit(1)

        num_seqs_checked += 1
        if num_seqs_checked % one_tenth_of_seqs == 0:
            percent_seqs_checked = num_seqs_checked / num_seqs * 100
            log.info(f'Checked {num_seqs_checked} ({percent_seqs_checked:.2f}%) database sequences...')

    percent_seqs_checked = num_seqs_checked / num_seqs * 100

    log.info(f'Confirmed {num_seqs} ({percent_seqs_checked}%) database sequence headers are properly formatted.')
    log.info(f'Confirmed {num_seqs} ({percent_seqs_checked}%) database sequences are within expected length ranges for their indicated segments.')
    log.info(f'Confirmed {num_seqs} ({percent_seqs_checked}%) database sequences do not contain ambiguous or lower-case nucleotides.')

    db = os.path.abspath(db)
    db_suffixes = ['nhr', 'nin', 'nsq']
    if any([not os.path.exists(db + '.' + suffix) for suffix in db_suffixes]):
        terminal_command = f'makeblastdb -in {db} -dbtype nucl'
        process_name = 'makeblastdb_contigs'
        error_code = 5
        return_code = run(terminal_command, outdir, out_name, process_name, error_code)
        if return_code != 0:
            log.error('Error running makeblastdb.')
            exit(1)
