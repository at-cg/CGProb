import numpy as np
import sys

class ReadData:
    def __init__(self):
        self.ids = []
        self.start_positions = []
        self.end_positions = []
        self.lengths = []

    def add_read(self, read_id, start, end, length):
        self.ids.append(read_id)
        self.start_positions.append(start)
        self.end_positions.append(end)
        self.lengths.append(length)

    def to_numpy_arrays(self):
        return (np.array(self.ids), np.array(self.start_positions), 
                np.array(self.end_positions), np.array(self.lengths))

def process_reads(read_data1, read_data2, het_locus):
    # Example processing function that takes data from both files
    ids1, start1, end1, lengths1 = read_data1.to_numpy_arrays()
    ids2, start2, end2, lengths2 = read_data2.to_numpy_arrays()

    # Identify x1
    x1 = -1
    for i in range(len(ids1)):
        if start1[i] <= het_locus <= end1[i]:
            if end1[i] > x1:
                x1 = end1[i]

    # Identify x2
    x2 = -1
    for i in range(len(ids2)):
        if start2[i] <= het_locus <= end2[i]:
            if end2[i] > x2:
                x2 = end2[i]
    
    if (x1 == x2) or (x1 == -1) or (x2 == -1):
        return 0
    
    elif x1 > x2:
        # Check for contained reads:
        # Hap 1
        contained_in_hap1 = False
        for i in range(len(ids1)):
            if ((het_locus + 1) <= start1[i] <= x2) and ((1+x2) <= end1[i] <= x1):
                contained_in_hap1 = True
                break
        filling_read_in_hap1 = False
        for i in range(len(ids1)):
            if ((het_locus + 1) <= start1[i] <= x2) and (1+x1 <= end1[i]):
                filling_read_in_hap1 = True
                break
        # Hap 2
        contained_in_hap2 = False
        for i in range(len(ids2)):
            if ((het_locus + 1) <= start2[i] <= x2) and ((1+x2) <= end2[i] <= x1):
                contained_in_hap2 = True
                break
        filling_read_in_hap2 = False
        for i in range(len(ids2)):
            if ((het_locus + 1) <= start2[i] <= x2) and (1+x1 <= end2[i]):
                filling_read_in_hap2 = True
                break
        
        # Test conditions
        if (contained_in_hap1 or contained_in_hap2) and (not filling_read_in_hap1 and not filling_read_in_hap2):
            return 1
        else:
            return 0
    else:
        # Check for contained reads:
        # Hap 1
        contained_in_hap1 = False
        for i in range(len(ids1)):
            if ((het_locus + 1) <= start1[i] <= x1) and ((1+x1) <= end1[i] <= x2):
                contained_in_hap1 = True
                break
        filling_read_in_hap1 = False
        for i in range(len(ids1)):
            if ((het_locus + 1) <= start1[i] <= x1) and (1+x2 <= end1[i]):
                filling_read_in_hap1 = True
                break
        # Hap 2
        contained_in_hap2 = False
        for i in range(len(ids2)):
            if ((het_locus + 1) <= start2[i] <= x1) and ((1+x1) <= end2[i] <= x2):
                contained_in_hap2 = True
                break
        filling_read_in_hap2 = False
        for i in range(len(ids2)):
            if ((het_locus + 1) <= start2[i] <= x1) and (1+x2 <= end2[i]):
                filling_read_in_hap2 = True
                break
        
        # Test conditions
        if (contained_in_hap1 or contained_in_hap2) and (not filling_read_in_hap1 and not filling_read_in_hap2):
            return 1
        else:
            return 0


def main():
    # Check if the user provided the number of iterations
    if len(sys.argv) != 3:
        print("Usage: python countOutputsWithGaps.py <number_of_iterations> <het_locus>")
        return

    n_iterations = int(sys.argv[1])
    het_locus = int(sys.argv[2])

    gap_count = 0
    # Process both files for the specified number of iterations
    for i in range(n_iterations):
        # Read data from in1.txt
        read_data1 = ReadData()
        with open('outputHap1.txt', 'r') as file1:
            lines = file1.readlines()

        # Process lines between !BEGIN and !END
        in_block = False
        for line in lines:
            line = line.strip()
            if line.startswith("!BEGIN"):
                in_block = True
                continue
            elif line.startswith("!END"):
                in_block = False
                continue
            
            if in_block:
                parts = line.split('\t')
                if len(parts) == 4:
                    read_id = int(parts[0])
                    start = int(parts[1])
                    end = int(parts[2])
                    length = int(parts[3])
                    read_data1.add_read(read_id, start, end, length)

        # Read data from in2.txt
        read_data2 = ReadData()
        with open('outputHap2.txt', 'r') as file2:
            lines = file2.readlines()

        # Process lines between !BEGIN and !END
        in_block = False
        for line in lines:
            line = line.strip()
            if line.startswith("!BEGIN"):
                in_block = True
                continue
            elif line.startswith("!END"):
                in_block = False
                continue
            
            if in_block:
                parts = line.split('\t')
                if len(parts) == 4:
                    read_id = int(parts[0])
                    start = int(parts[1])
                    end = int(parts[2])
                    length = int(parts[3])
                    read_data2.add_read(read_id, start, end, length)

        # Process reads from both files
        test = process_reads(read_data1, read_data2, het_locus)
        gap_count += test
    
    print(gap_count)

if __name__ == "__main__":
    main()
