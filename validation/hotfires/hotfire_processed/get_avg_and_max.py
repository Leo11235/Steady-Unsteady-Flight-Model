# Reads json files of the processed data, outputs the max and average thrust values

with open("HOTFIRE4.1.jsonc", 'r') as file:
    MIN = 6630.3
    MAX = 6662
    time = 0
    thrust = 0
    largest = 0
    total = 0
    i = 0
    avg = 0
    times = []
    thrustnums = []
    for line in file:
        #print(line)
        if time == 1 and (line != "{\n") and (line != "    \"seconds\": [\n") and (line != "    \"thrust\": [\n") and (line != "    ],\n"):
            line = line.strip(",\n")
            times.append(float(line))
        if thrust == 1 & (line != "    ]\n") & (line != "}"):
            line = line.strip(",\n")
            thrustnums.append(float(line))
        if line == "    \"seconds\": [\n":
            time = 1
        if line == "    \"thrust\": [\n":
            #print("test")
            thrust = 1
            time = 0
    for num in thrustnums:
        #total += num
        if times[i] >= MIN and times[i] <= MAX:
            total += num
            if(num > largest):
                largest = num
        i += 1
    avg = total / i
    print(avg)
    print(largest)
