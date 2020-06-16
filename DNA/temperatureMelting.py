import math


class MeltingTemp:
    def __init__(self, nukleotid1, nukleotid2, Ca, Cb):
        self.nukleotid1 = nukleotid1
        self.nukleotid2 = nukleotid2
        self.dH = 0  # [kcal/mol]
        self.dS = 0  # [cal/(K*mol)]
        self.Ca = Ca
        self.Cb = Cb

    table_of_nucleotids = [
            ('AA', 'TT', -7.9, -22.2),
            ('TT', 'AA', -7.9, -22.2),
            ('AT', 'TA', -7.2, -20.4),
            ('TA', 'AT', -7.2, -21.3),
            ('CA', 'GT', -8.5, -22.7),
            ('GT', 'CA', -8.4, -22.4),
            ('CT', 'GA', -7.8, -21.0),
            ('GA', 'CT', -8.2, -22.2),
            ('CG', 'GC', -10.6, -27.2),
            ('GC', 'CG', -9.8, -24.4),
            ('GG', 'CC', -8.0, -19.9),
            ('CC', 'GG', -8.0, -19.9)
    ]

    def sum_dH_dS(self, k):  # count dH and dS
        i = 0
        while i < k - 1:
            j = 0
            while j < 10:
                if self.nukleotid1[i] + self.nukleotid1[i + 1] == self.table_of_nucleotids[j][0]:
                    self.dH += self.table_of_nucleotids[j][2]
                    self.dS += self.table_of_nucleotids[j][3]
                j = j + 1
            i = i + 1

    def tm(self):
        R = 1.987  # universal const
        k1 = len(self.nukleotid1)
        k2 = len(self.nukleotid2)
        self.sum_dH_dS(k1)
        self.sum_dH_dS(k2)

        if self.Ca > 0 and self.Cb == 0:
            return (1000 * self.dH) / (self.dS + R * (math.log(self.Ca)))  # in K
        elif self.Ca > 0 and self.Cb > 0:
            return (1000 * self.dH) / (self.dS + (R * (math.log(self.Ca - (self.Cb / 2)))))  # in K
        else:
            return (1000 * self.dH) / self.dS
