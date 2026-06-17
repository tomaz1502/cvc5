import sys

def is_num(c):
    return ord(c) >= ord('0') and ord(c) <= ord('9')

def extract_number(line):
    l = ""
    i = 0
    while not is_num(line[i]):
        i += 1
    while is_num(line[i]):
        l += line[i]
        i += 1
    return int(l)

def main():
    with open(sys.argv[1], 'r') as f:
        projection_time = 0
        total = 0
        for line in f:
            if line[-2:] == ".\n":
                t = extract_number(line)
                if line[:4] == "full":
                    total = t
                else:
                    projection_time += t
        print(f"projection time = {projection_time} microseconds")
        print(f"total time = {total} microseconds")

if __name__ == '__main__':
    main()
