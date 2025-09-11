from dbg import DBG
from utils import read_data, calculate_nga50
import sys
import os
import time

sys.setrecursionlimit(1000000)

if __name__ == "__main__":
    # Start timing
    start_time = time.time()
    
    argv = sys.argv
    short1, short2, long1 = read_data(os.path.join('./', argv[1]))

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    # dbg.show_count_distribution()
    
    # Store contigs for NGA50 calculation
    contigs = []
    
    with open(os.path.join('./', argv[1], 'contig.fasta'), 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            contigs.append(c)
            # redirect print to file instead of console
            print(f">contig_{i}\n{c}", file=f)
    
    # End timing
    end_time = time.time()
    total_time = end_time - start_time
    
    # Calculate NGA50
    nga50 = calculate_nga50(contigs)
    
    # Write statistics to file
    stats_file = os.path.join('./', argv[1], 'assembly_stats.txt')
    with open(stats_file, 'w') as f:
        print("Assembly Statistics", file=f)
        print("===================", file=f)
        print(f"Total execution time: {total_time:.2f} seconds", file=f)
        print(f"Number of contigs: {len(contigs)}", file=f)
        print(f"NGA50: {nga50}", file=f)
        print(f"Contig lengths: {[len(c) for c in contigs]}", file=f)

    print(f"{argv[1]}  python  {total_time:.2f}  {nga50}")