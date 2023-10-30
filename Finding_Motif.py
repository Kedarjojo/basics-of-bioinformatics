# Specify the file path
file_path = "your_file.txt"

# Initialize variables to store the lines
sequence = ""
target = ""

# Open the file and read its lines
with open(file_path, 'r') as file:
    sequence = file.readline().strip()  # Read the first line and remove leading/trailing whitespace
    target = file.readline().strip()  # Read the second line and remove leading/trailing whitespace

index=0
start_positions = []
while index < len(sequence):
    index = sequence.find(target,index)
    if index == -1:
        break
    start_positions.append(index)
    index+=1
new_positions = [x+1 for x in start_positions]

print(new_positions)
