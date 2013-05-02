"""Converts mySQL tables into csv files"""

import argparse
import csv
import db_toolbox as SQL

class TableWriter:

    def __init__(self, interface):
        self.interface = interface
        

