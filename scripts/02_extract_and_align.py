# Extract ORFs, run pairwise alignment, and filter for >=90% cross-species identity
from utils import load_config

config = load_config()
cutoff = config['thresholds']['conservation_identity']
print(f'Filtering proteins with cross-species identity >= {cutoff}%...')
