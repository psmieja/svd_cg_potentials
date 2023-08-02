import sys
import pandas as pd
from my_utilities.data import get_data

if len(sys.argv) != 3:
    sys.exit(f"Usage: python {sys.argv[0]} [dat files dir path] [output filename]")

get_data(sys.argv[1]).to_pickle(sys.argv[2])