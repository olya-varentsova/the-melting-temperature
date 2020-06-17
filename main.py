import math
import argparse
import sys
from DNA.temperatureMelting import MeltingTemp


def check_on_correct_input(dna1, dna2, dH, dS, Ca, Cb):
    k = 0
    k1 = len(dna1)
    k2 = len(dna2)
    if Ca < 0 or Cb < 0:
        exit('Напишите корректную концентрацию')
    if dH < 0:
        exit('Напишите корректный параметр энтальпии: dH >= 0')
    if dS < 0:
        exit('Напишите корректный параметр энторопии: dS >= 0')
    while k < min(k1, k2):
        if not ((dna2[k] == 'A' and dna1[k] == 'T') or (dna2[k] == 'T' and dna1[k] == 'A') or
                (dna2[k] == 'G' and dna1[k] == 'C') or (dna2[k] == 'C' and dna1[k] == 'G')):
            exit('Напишите корректную последовательность')
        k += 1
    if k2 > k1:
        while k < k2:
            if not (dna2[k] == 'T' or dna2[k] == 'A' or dna2[k] == 'C' or dna2[k] == 'G'):
                exit('Напишите корректную последовательность')
            k += 1
    if k1 > k2:
        while k < k1:
            if not (dna1[k] == 'T' or dna1[k] == 'A' or dna1[k] == 'C' or dna1[k] == 'G'):
                exit('Напишите коректную последовательность')
            k += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-seq1", help="Sequence of DNA", type=str)
    parser.add_argument("-seq2", help="Sequence of DNA", type=str)
    parser.add_argument("-dH", help="Parameter of enthalpy", default=0, type=float)
    parser.add_argument("-dS", help="Parameter of entropy", default=0, type=float)
    parser.add_argument("-Ca", help="the total oligonucleotide 1 strand concentration"
                                    " if the strands are in different or equal concentrations", default=0, type=float)
    parser.add_argument("-Cb", help="the total oligonucleotide 2 strand concentration"
                                    " if the strands are in different concentration", default=0, type=float)
    args = parser.parse_args(sys.argv[1:])

    check_on_correct_input(args.seq1, args.seq2, args.dH, args.dS, args.Ca, args.Cb)
    temp = MeltingTemp(args.seq1, args.seq2, args.dH, args.dS, args.Ca, args.Cb)

    print("%.2fK" % temp.tm())


if __name__ == '__main__':
    main()
