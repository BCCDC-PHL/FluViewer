import logging
import os
import sys

def configure_logging(log_level, log_file):
    """
    Configure logging to log to stderr and args.outdir/logs/fluviewer.log.
    Should configure in such a way that all other modules can
    call logging.getLogger(__name__) and have the same log settings.

    """
    logger = logging.getLogger()
    logger.setLevel(log_level)

    # Use ISO 8601 format, milliseconds, no timezone
    formatter = logging.Formatter(
        fmt='%(asctime)s.%(msecs)03d - %(levelname)s %(module)s.%(funcName)s L%(lineno)d : %(message)s',
        datefmt='%Y-%m-%dT%H:%M:%S'
    )

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    logs_dir = os.path.join(os.path.dirname(log_file))
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)

    fluviewer_log_file_path = os.path.join(logs_dir, f'fluviewer.log')
    log_file_handler = logging.FileHandler(fluviewer_log_file_path)
    log_file_handler.setFormatter(formatter)
    logger.addHandler(log_file_handler)

