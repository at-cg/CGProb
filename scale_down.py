import sys
from collections import defaultdict

def scale_down(input_filename, output_filename, divisor):
    # Create a dictionary to store and aggregate values based on the scaled key
    scaled_values = defaultdict(int)

    # Read the input file and process each line
    with open(input_filename, 'r') as input_file:
        for line in input_file:
            parts = line.split()
            if len(parts) == 2:
                key, value = int(parts[0]), int(parts[1])
                scaled_key = (key // divisor) + 1
                scaled_values[scaled_key] += value

    # Write the aggregated results to the output file
    with open(output_filename, 'w') as output_file:
        for key, value in sorted(scaled_values.items()):
            output_file.write(f"{key} {value}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 scale_down.py <input_file> <output_file> <divisor>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    divisor = int(sys.argv[3])

    scale_down(input_file, output_file, divisor)