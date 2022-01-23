import sys
from collections import defaultdict

filename = sys.argv[1]

with open(filename, 'r') as f:
    lines = f.readlines()

lines = lines[1:]

# Results is {string => { int => [int] }}
results = defaultdict(lambda: defaultdict(list))
for line in lines:
    line = line.strip()
    strategy, words, guesses = line.split(",")
    results[strategy][int(words)].append(int(guesses))

for strategy, d in results.items():
    points = []

    num_entries = 0
    word_sum = 0
    guess_sum = 0
    for word_count, guess_count_list in sorted(d.items()):
        num_entries += len(guess_count_list)
        word_sum += word_count * len(guess_count_list)
        guess_sum += sum(guess_count_list)
        if num_entries > 100:
            word_avg = word_sum / num_entries
            guess_avg = guess_sum / num_entries
            points.append((word_avg, guess_avg))
            num_entries = 0
            word_sum = 0
            guess_sum = 0

    for x, y in points:
        print("%s, %s" % (x, y))
