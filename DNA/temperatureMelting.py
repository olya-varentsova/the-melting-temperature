import math


class MeltingTemp:
    def __init__(self, dna1, dna2, dH, dS, Ca, Cb):
        self.dna = [dna1, dna2]
        self.dH = dH  # [kcal/mol]
        self.dS = dS  # [cal/(K*mol)]
        self.Ca = Ca
        self.Cb = Cb

    table_of_nucleotids =  {
            'AA': (-7.9, -22.2),
            'TT': (-7.9, -22.2),
            'AT': (-7.2, -20.4),
            'TA': (-7.2, -21.3),
            'CA': (-8.5, -22.7),
            'GT': (-8.4, -22.4),
            'CT': (-7.8, -21.0),
            'GA': (-8.2, -22.2),
            'CG': (-10.6, -27.2),
            'GC': (-9.8, -24.4),
            'GG': (-8.0, -19.9),
            'CC': (-8.0, -19.9)
        }



    def sum_dH_dS(self):  # count dH and dS
        for strand in self.dna:
            n = self.dna.index(strand)
            k = len(strand)
            i = 0
            while i < k - 1:
                for key in self.table_of_nucleotids:
                    if self.dna[n][i] + self.dna[n][i + 1] == key:
                        self.dH += self.table_of_nucleotids[key][0]
                        self.dS += self.table_of_nucleotids[key][1]
                i = i + 1

    def tm(self):
        R = 1.987  # universal const
        self.sum_dH_dS()

        if self.Ca > 0 and self.Cb == 0:
            return (1000 * self.dH) / (self.dS + R * (math.log(self.Ca)))  # in K
        elif self.Ca > 0 and self.Cb > 0:
            return (1000 * self.dH) / (self.dS + (R * (math.log(self.Ca - (self.Cb / 2)))))  # in K
        else:
            return (1000 * self.dH) / self.dS
