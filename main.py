import math
import argparse
import sys
from DNA.temperatureMelting import MeltingTemp


def check_on_correct_input(nukleotid1, nukleotid2):
    k = 0
    k1 = len(nukleotid1)
    k2 = len(nukleotid2)
    while k < min(k1, k2):
        if not ((nukleotid2[k] == 'A' and nukleotid1[k] == 'T') or (nukleotid2[k] == 'T' and nukleotid1[k] == 'A') or
                (nukleotid2[k] == 'G' and nukleotid1[k] == 'C') or (nukleotid2[k] == 'C' and nukleotid1[k] == 'G')):
            exit('Напишите коректную последовательность')
        k += 1
    if k2 > k1:
        while k < k2:
            if not (nukleotid2[k] == 'T' or nukleotid2[k] == 'A' or nukleotid2[k] == 'C' or nukleotid2[k] == 'G'):
                exit('Напишите коректную последовательность')
            k += 1
    if k1 > k2:
        while k < k1:
            if not (nukleotid1[k] == 'T' or nukleotid1[k] == 'A' or nukleotid1[k] == 'C' or nukleotid1[k] == 'G'):
                exit('Напишите коректную последовательность')
            k += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-seq1", help="Sequence of DNA", type=str)
    parser.add_argument("-seq2", help="Sequence of DNA", type=str)
    parser.add_argument("-Ct", help="the total oligonucleotide strand concentration "
                                    "if the strands are in equal concentration", default=0, type=float)
    parser.add_argument("-Ca", help="the total oligonucleotide 1 strand concentration"
                                    " if the strands are in different concentration", default=0, type=float)
    parser.add_argument("-Cb", help="the total oligonucleotide 2 strand concentration"
                                    " if the strands are in different concentration", default=0, type=float)
    args = parser.parse_args(sys.argv[1:])

    check_on_correct_input(args.seq1, args.seq2)
    temp = MeltingTemp(args.seq1, args.seq2, args.Ct, args.Ca, args.Cb)

    print("%.2fK" % temp.tm())


if __name__ == '__main__':
    main()
