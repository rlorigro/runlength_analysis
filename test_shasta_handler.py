from handlers.ShastaRunlengthHandler import ShastaRunlengthHandler


def main():
    shasta_path = "/home/ryan/code/runlength_analysis/data/ecoli_shasta_pileup_data.csv"

    shasta_handler = ShastaRunlengthHandler(shasta_path)

    consensus_sequence, consensus_lengths, pileup_data = shasta_handler.get_pileup_data(cutoff=100)

    print(consensus_sequence[:10])
    print(consensus_lengths[:10])

    for item in pileup_data[:10]:
        print(list(map(str,item)))


if __name__ == "__main__":
    main()