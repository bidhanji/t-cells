# Process output from NetCTLpan 1.1 / NetMHCpan 4.1 predictions
import pandas as pd
from utils import load_config

config = load_config()
mhc1_cutoff = config['thresholds']['mhc1_percent_rank']
print(f'Filtering Class I epitopes with %Rank <= {mhc1_cutoff}...')
