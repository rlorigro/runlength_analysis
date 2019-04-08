from handlers.ShastaRunlengthHandler import ShastaRunlengthHandler
import sys


def main():
    shasta_path = "/home/ryan/code/runlength_analysis/data/ecoli_shasta_pileup_data_NEW.csv"

    shasta_handler = ShastaRunlengthHandler(shasta_path)

    start = 225343 - 1

    name, read = shasta_handler.get_pileup_data(cutoff=start+20)

    with open(shasta_path, "r") as file:
        # print(read.sequence[:10])
        # print(read.lengths[:10])

        for i, (item, line) in enumerate(zip(read.pileup, file)):
            if start < i < start+20:
                print(i)
                print(line.strip())
                print(list(map(str,item)))
                print()


if __name__ == "__main__":
    main()
