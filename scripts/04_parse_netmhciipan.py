# Process output from NetMHCIIpan 4.3 predictions for BoLA-DRB3
import pandas as pd
from utils import load_config

config = load_config()
mhc2_cutoff = config['thresholds']['mhc2_percent_rank']
print(f'Filtering Class II epitopes with %Rank <= {mhc2_cutoff}...')
