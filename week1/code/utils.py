import os

def read_fasta(path, name):
    data = []
    with open(os.path.join(path, name), 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line[0] != '>':
                data.append(line)
    return data

def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1

def calculate_n50(contigs):
    """
    Calculate N50 (length of the contig at which 50% of total assembly length is reached)
    
    Args:
        contigs: List of contig sequences
        
    Returns:
        N50 value (length of the contig at which 50% of total assembly length is reached)
    """
    if not contigs:
        return 0
        
    # Get lengths of all contigs
    contig_lengths = [len(contig) for contig in contigs]
    
    # Sort lengths in descending order
    contig_lengths.sort(reverse=True)
    
    # Calculate total assembly length
    total_length = sum(contig_lengths)
    
    # Find N50
    cumulative_length = 0
    for length in contig_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
            
    return 0

def calculate_assembly_stats(contigs):
    """
    Calculate comprehensive assembly statistics
    
    Args:
        contigs: List of contig sequences
        
    Returns:
        Dictionary with various assembly statistics
    """
    if not contigs:
        return {
            'num_contigs': 0,
            'total_length': 0,
            'n50': 0,
            'max_contig': 0,
            'min_contig': 0,
            'mean_contig': 0
        }
        
    contig_lengths = [len(contig) for contig in contigs]
    contig_lengths.sort(reverse=True)
    
    total_length = sum(contig_lengths)
    n50 = calculate_n50(contigs)
    
    stats = {
        'num_contigs': len(contigs),
        'total_length': total_length,
        'n50': n50,
        'max_contig': max(contig_lengths),
        'min_contig': min(contig_lengths),
        'mean_contig': total_length / len(contigs) if contigs else 0
    }
    
    return stats