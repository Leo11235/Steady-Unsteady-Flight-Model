# Reads json files of the processed data, outputs the max and average thrust values

arr = ["2.2", "2.3", "2.4", "2.6", "2.7", "3.1", "3.4", "3.5", "4.1"]

for n in arr:
    with open(f"validation/hotfires/hotfire_processed/HOTFIRE{n}.jsonc", 'r') as file:
        MIN = [7383.9, 4798, 3622.4, 4322.1, 3473.5, 4697.9, 4339.0, 2205.8, 6630.3]
        MAX = [7390, 4806, 3634, 4333, 3485, 4711, 4349, 2238, 6662]
        time = 0
        thrust = 0
        largest = 0
        total = 0
        i = 0
        avgdivisor = 0
        avg = 0
        times = []
        thrustnums = []
        startinterval = 0
        endinterval = 0

        for line in file:
            #print(line)
            if time == 1 and (line != "{\n") and (line != "    \"seconds\": [\n") and (line != "    \"cc_pressure\": [\n") and (line != "    ],\n"):
                line = line.strip(",\n")
                times.append(float(line))
            if thrust == 1 & (line != "    ]\n") & (line != "}"):
                line = line.strip(",\n")
                thrustnums.append(float(line))
            if line == "    \"seconds\": [\n":
                time = 1
            if line == "    \"cc_pressure\": [\n":
                #print("test")
                thrust = 1
                time = 0
        for num in thrustnums:
            #total += num
            if times[i] >= MIN[arr.index(n)] and times[i] <= MAX[arr.index(n)]:
                #if startinterval == 0:
                #    print(num)
                #    startinterval = 1
                total += num
                if(num > largest):
                    largest = num
                avgdivisor += 1
            i += 1
            #if endinterval == 0 and times[i] >= MAX[arr.index(n)]:
            #    print(num)
            #    endinterval = 1
        avg = total / avgdivisor
        print(f"Avg cc_pressure for hotfire {n}: {avg}")
        print(f"Max cc_pressure for hotfire {n}: {largest}")
