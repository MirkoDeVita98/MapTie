import sys
from genome_input import reading
from alignment_evaluation import testing
from genome_indexing import indexing
from alignment import processing
from utils.parser import load_config
import multiprocessing as mp
from multiprocessing import Process
import numpy as np
from alignment_output.output import ConcatOutput
from constants import seed_lengths

class MapTie:
    def __init__(self,data):
        self.reference_abs_path = data['input']['reference_abs_path']
        self.reference_name = data['input']['reference_name']
        self.seed_lengths = seed_lengths
        self.fname = reading.path_leaf(data['input']['read1'])[:-4]
        self.read1 = data['input']['read1']
        self.read2 = data['input']['read2']
        self.whole_genome = reading.read_reference_genome(self.reference_abs_path)
        self.suffix_array = indexing.get_indexes(self.whole_genome,self.reference_name, self.seed_lengths)     
        self.lines = reading.get_lines(self.read1,self.read2) 
        # self.lines = testing.get_test_lines(self.whole_genome, length=30000, n=4000)

    def process(self, num_processes=mp.cpu_count()):
        period = 50000
        if num_processes == 1:
            processing.process_lines(
                self.whole_genome, self.suffix_array, self.lines, self.fname,
                self.seed_lengths, 0, period)
        else:
            processes = []
            partitions = np.array_split(self.lines, num_processes)
            for ii in range(num_processes):
                p = Process(target=processing.process_lines, 
                        args=(self.whole_genome, self.suffix_array, partitions[ii], 
                            self.fname, self.seed_lengths, ii, period))
                p.start()
                processes.append(p)
            for p in processes:
                p.join()
        out_writer = ConcatOutput(
            self.whole_genome['id'], len(self.whole_genome['genome']), 
            f'results_{self.fname}.sam')
        out_writer.combine_results(num_processes)

def main():
        data = load_config()
        maptie = MapTie(data)
        if data['execution']['n_jobs'] == -1:
            maptie.process()
        else:
            maptie.process(num_processes= data['execution']['n_jobs'])

if __name__== "__main__" :
    main()

# num_processes = mp.cpu_count()

# processes = []
# partitions = np.array_split(lines, num_processes)

# for ii in range(num_processes):
#     p = Process(target=processing.process_lines, 
#                 args=(whole_genome, suffix_arrays, partitions[ii], 
#                       size, seed_lengths, coverage, ii, 50000))
#     p.start()
#     processes.append(p)

# for p in processes:
#    p.join()

# out_writer = ConcatOutput(
#     whole_genome['id'], len(whole_genome['genome']), 
#     f'results_{coverage}_{size}.sam')
# out_writer.combine_results(mp.cpu_count())