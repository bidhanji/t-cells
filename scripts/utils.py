import yaml
from pathlib import Path

def load_config(config_path='config/config.yaml'):
    '''Load global parameters from YAML configuration.'''
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def parse_fasta(fasta_file):
    '''Yield header and sequence generator.'''
    from Bio import SeqIO
    for record in SeqIO.parse(fasta_file, 'fasta'):
        yield record.id, str(record.seq)
