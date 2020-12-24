import re

    # read from our output
def find(filename):
    with open(filename) as f:
        counter = 0
        ll = 0
        output_arguments = []
        for l in f:
            read_arguments = l.strip().split("\t")
            if read_arguments[0][0] != '@':
                output_arguments.append(read_arguments)
            ll += 1
        for o in output_arguments:
            if(len(o[5]) >= 16):
            # if sum(int(s) for s in re.findall(r'\d+', o[5])) >= 130:
            # if '-' in o[5]:
            # if o[5][:2] == '48':
                print("\n", (o[0], o[1], o[3], o[5]))
                counter+=1
            
    print("\n# of wrong cigars: ",counter)
    print("\n% of wrong cigars: ",counter / ll)


find("/home/anej/repos/studies/CBM/CMB_2020_project_1/results/results_output_30xCov.sam")
# find("/home/anej/repos/studies/CBM/CMB_2020_project_1/data_small/output_tiny_30xCov.sam")