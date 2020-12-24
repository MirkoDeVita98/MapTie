import os
import sys; sys.path.append('..')
import fileinput

import constants
from utils.parser import load_config

class ConcatOutput:

    def __init__(self, rname: str, rlength: int, filename: str):
        self.rname = rname  # REFERENCE NAME
        self.rlength = None  # REFERENCE LENGTH
        self.filename = filename
        os.makedirs(load_config()['output']['output_dir'], exist_ok=True)
        self.file = open(os.path.join(load_config()['output']['output_dir'], filename), "w+")
        self.header = "@HD" + '\t' + 'VN:' + '1.0' + '\t' + 'SO:' + 'unsorted' + '\n' \
                 "@SQ" + '\t' + 'SN:' + rname + '\t' + 'LN:' + str(rlength) + '\n'
        self.file.write(self.header)

    def combine_results(self, number_of_files: int):
        filenames = [
            os.path.join(load_config()['output']['output_dir'], f'process_{ii}_{self.filename}')
            for ii in range(number_of_files)]

        # print(filenames)
        
        with fileinput.input(filenames) as fin:
            for line in fin:
                self.file.write(line)

    def close(self):
        self.file.close()


class TempOutput:

    def __init__(self, rname:str, process_no:int,fname:str):
        data = load_config()
        self.rname = rname

        os.makedirs(data['output']['output_dir'], exist_ok=True)
        filename = '_'.join(
            ['process', str(process_no), f'results_{fname}.sam'])
        self.file = open(os.path.join(data['output']['output_dir'], filename), "w+")
        self.outputs = []

    def append_alignment(self, qname: str, pos: int, cigar: str,
                         alignment: str, genome_lenght: int, read_lenght: int,
                         reverse: bool=False, secondary: bool=False, first: bool=False, last: bool=False) -> int:
        flag = 1 + 2 +  16*reverse + 32 *(not reverse) + 64*first + 128*last + 256*secondary   # default FLAG
        # if reverse:
        #     pos = genome_lenght - pos - read_lenght - 2
        output_line = qname + '\t' + str(flag) + '\t' + self.rname + '\t' + \
                    str(pos) + '\t' + '0' + '\t' + cigar + '\t' + '*' \
                    + '\t' + '0' + '\t' + '0' + '\t' + alignment + '\t' \
                    + '*' + '\n'

        self.outputs.append(output_line)

        cache_size = 50 if read_lenght > 10000 else 1000
        if len(self.outputs) == cache_size:
            self.file.write(''.join(self.outputs))
            self.outputs.clear()
            
        return len(self.outputs) - 1

    def remove_alignment(self, alignment_id: int):
        del self.outputs[alignment_id]

    def flush(self):
        self.file.write(''.join(self.outputs))
        self.outputs = []

    def close(self):
        self.file.write(''.join(self.outputs))
        self.file.close()







