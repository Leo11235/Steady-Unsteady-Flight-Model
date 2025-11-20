# Reads json files of the processed data, outputs the max and average thrust values

with open("HOTFIRE4.1.jsonc", 'r') as file:
    thrust = 0
    largest = 0
    total = 0
    i = 0
    avg = 0
    nums = []
    for line in file:
        #print(line)
        if thrust == 1 & (line != "    ]\n") & (line != "}"):
            line = line.strip(",\n")
            nums.append(float(line))
        if line == "    \"thrust\": [\n":
            #print("test")
            thrust = 1
    for num in nums:
        total += num
        i += 1
        if(num > largest):
            largest = num
    avg = total / i
    print(avg)
    print(largest)
