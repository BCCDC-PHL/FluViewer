
import logging




def get_logger(name: str, level_str: str = "info") -> logging.Logger:
    """
    Get a logger with the given name and level.

    Raises a ValueError if the level is invalid.

    :param name: The name of the logger.
    :type name: str
    :param level: The level of the logger.
    :type level: str
    :return: The logger.
    :rtype: logging.Logger
    """
    if level_str not in ["debug", "info", "warning", "error", "critical"]:
        raise ValueError("Invalid level: {level}")

    level = getattr(logging, level_str.upper())
    
    logger = logging.getLogger(name)
    logger.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)
    # Use ISO 8601 format, local timezone
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(levelname)s %(module)s.%(funcName)s L%(lineno)d : %(message)s',
        datefmt='%Y-%m-%dT%H:%M:%S%z'
    )
    
    ch.setFormatter(formatter)

    logger.addHandler(ch)
        

    return logger

    
