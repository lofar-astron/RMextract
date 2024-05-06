#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set up logging for RMextract
"""

import logging
import os


def setup_logging():
    ctrl = os.environ.get("RMEXT_LOGLEVEL", "INFO")
    if ctrl == "ERROR":
        level = logging.ERROR
    elif ctrl == "WARNING":
        level = logging.WARNING
    elif ctrl == "INFO":
        level = logging.INFO

    logger = logging.getLogger("RMextract")
    logger.setLevel(level)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s: %(message)s")
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(level)
    logger.addHandler(ch)
    return logger


logger = setup_logging()
