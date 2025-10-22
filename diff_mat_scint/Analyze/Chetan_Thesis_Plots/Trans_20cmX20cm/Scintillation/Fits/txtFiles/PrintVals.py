import sys

# Check if filename is given
if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <file.txt>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, 'r') as f:
        for line in f:
            # Skip empty or whitespace-only lines
            if not line.strip():
                continue

            columns = line.split()  # splits on any whitespace
            if len(columns) >= 2:
                print(columns[1])  # second column (index 1)
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
