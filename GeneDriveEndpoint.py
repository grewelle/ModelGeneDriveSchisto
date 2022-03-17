import numpy as np
import scipy
import seaborn as sns
import pandas as pd

sns.set(style="white")
#sns.set_palette("Reds")
flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e",
                  "#2ecc71"]
#sns.set_palette(sns.color_palette(flatui))
import matplotlib.pyplot as pl

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (7, 5),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pl.rcParams.update(params)

def main():
    spread = np.random.uniform(0.4, 1, size=100)
    list1 = []
    list2 = []
    list3 = []

    for x in spread:

        numEndpoints = 100
        endpoints = []
        for u in range(numEndpoints+1):
            sig = u/numEndpoints  # selfing rate
            aa = 0.5
            # genotypes must be greater than 0
            A = (sig - np.sqrt(16*aa - 24*sig*aa + sig**2*(1 + 8*aa)))/(4*(-1 + sig))
            B = 1-A
            #sAA = [A**2*100+.0000000000001]
            sAA = [(A**2+A*B*(sig/(2-sig)))*100+.0000000000001]
            #sAB = [2*A*B*100+.0000000000001]
            sAB = [(4*A*B*(1-sig)/(2-sig))*100+.0000000000001]
            #sBB = [B**2*100+.0000000000001]
            sBB = [(B**2+A*B*(sig/(2-sig)))*100+.0000000000001]
            print(sAA, sAB, sBB, sig)

            #sAA = [50]  # set initial population size for wild type AA
            #sAB = [43]  # set initial population size for heterozygous AB
            #sBB = [7]  # set initial population size for homozygous resistant BB
            sABg = [0.00000000000001]  # set initial population size for heterozygous engineered resistant ABg
            sBBg = [0.00000000000001]  # set intial population size for homozygous resistant, heterozygous engineered copy BBg
            #sBgBg = [(u+.001)/(numEndpoints+.001)]  # set initial population size for homozygous engineered resistant BgBg
            #sBgBg = [(2*u/(numEndpoints))+.0000000000001]
            #sBgBg = [1/(u+1)]
            sBgBg = [1]
            resist = [0]
            s = [100]  # set initial total population size
            k = 100  # carrying capacity
            r = 20  # intrinsic growth rate
            inbr = x  # cost of inbreeding
            mu = 0.5  # natural death rate
            eps = 0 # migration parameter (place holder for now)
            rho = 0.15  # force of infection
            nhej = 0
            #mut = u/numEndpoints*10**(-3)
            mut = 0
            beta = (0.2835+200*mut) * r  # reduction in fecundity due to natural resistance
            gamma = 1  # dominance coefficient
            xi = 0.8  # conferred resistance to infection
            g0 = 0.9  # gene drive efficiency
            #g = 0.9*(1-np.exp(-2*u))
            beta_g = 0.1* r  # reduction in fecundity due to engineered resistance
            r_AB = r - gamma * beta  # fecundity of AB genotype
            r_BB = r - beta  # " BB genotype
            r_ABg = r - gamma * beta - beta_g  # " ABg genotype
            r_BBg = r - beta - beta_g  # " BBg genotype
            r_BgBg = r - beta - 2 * beta_g  # " BgBg genotype
            gens = 40  # no. of generations

            # begin simulation
            for i in range(1, gens + 1):
                # add 0 element to arrays as place holder for generation i (aka generation t+1)
                sAA.append(0)
                sAB.append(0)
                sBB.append(0)
                sABg.append(0)
                sBBg.append(0)
                sBgBg.append(0)
                s.append(0)
                g = g0 * (1 - resist[-1])*(1-nhej)

                # k = 100-50*(i%2)
                # mu = 0.2 + 0.1*(i%2)

                # simulate migration in first step
                sAA_migrate = sAA[i - 1] - eps * (sAA[i - 1]  - sAA[0])
                sAB_migrate = sAB[i - 1] - eps * (sAB[i - 1]  - sAB[0])
                sBB_migrate = sBB[i - 1] - eps * (sBB[i - 1]  - sBB[0])
                sABg_migrate = sABg[i - 1] - eps * sABg[i - 1]
                sBBg_migrate = sBBg[i - 1] - eps * sBBg[i - 1]
                sBgBg_migrate = sBgBg[i - 1] - eps * sBgBg[i - 1]

                # simulate mortality and infection in second step
                sAA_death = sAA_migrate - sAA_migrate * (mu + (1 - mu) * rho)
                sAB_death = sAB_migrate - sAB_migrate * (mu + (1 - mu) * rho * (1 - gamma * xi))
                sBB_death = sBB_migrate - sBB_migrate * (mu + (1 - mu) * rho * (1 - xi))
                sABg_death = sABg_migrate - sABg_migrate * (mu + (1 - mu) * rho * (1 - gamma * xi))
                sBBg_death = sBBg_migrate - sBBg_migrate * (mu + (1 - mu) * rho * (1 - xi))
                sBgBg_death = sBgBg_migrate - sBgBg_migrate * (mu + (1 - mu) * rho * (1 - xi))
                s_death = sAA_death + sAB_death + sBB_death + sABg_death + sBBg_death + sBgBg_death

                # calculate genotype frequencies
                p_AA = (sAA_death) / s_death
                p_AB = (sAB_death) / s_death
                p_BB = (sBB_death) / s_death
                p_ABg = sABg_death / s_death
                p_BBg = sBBg_death / s_death
                p_BgBg = sBgBg_death / s_death



                # vector of genotype frequencies
                gen1 = [p_AA, p_AB, p_BB, p_ABg, p_BBg, p_BgBg]

                # outcrossing transition matrix
                genOut = np.matrix([[r * p_AA + ((r + r_AB) * p_AB + (r + r_ABg) * p_ABg) / 4,
                                     ((r + r_BB) / 2) * p_BB + ((r + r_AB) * p_AB + (r + r_BBg) * p_BBg) / 4, 0,
                                     (1 - g) * (((r + r_BgBg) / 2) * p_BgBg + ((r + r_ABg) * p_ABg + (r + r_BBg) * p_BBg) / 4),
                                     0, g * (((r + r_BgBg) / 2) * p_BgBg + ((r + r_ABg) * p_ABg + (r + r_BBg) * p_BBg) / 4)],
                                    [(r + r_AB) * p_AA / 4 + (2 * r_AB * p_AB + (r_AB + r_ABg) * p_ABg) / 8,
                                     ((r + r_AB) * p_AA + 2 * r_AB * p_AB + (r_AB + r_BB) * p_BB) / 4 + (
                                             (r_AB + r_ABg) * p_ABg + (r_AB + r_BBg) * p_BBg) / 8,
                                     (r_AB + r_BB) * p_BB / 4 + (2 * r_AB * p_AB + (r_AB + r_BBg) * p_BBg) / 8, (1 - g) * (
                                             (r_AB + r_BgBg) * p_BgBg / 4 + (
                                             (r_AB + r_ABg) * p_ABg + (r_AB + r_BBg) * p_BBg) / 8), (1 - g) * (
                                             (r_AB + r_BgBg) * p_BgBg / 4 + (
                                             (r_AB + r_ABg) * p_ABg + (r_AB + r_BBg) * p_BBg) / 8), g * (
                                             (r_AB + r_BgBg) * p_BgBg / 2 + (
                                             (r_AB + r_ABg) * p_ABg + (r_AB + r_BBg) * p_BBg) / 4)],
                                    [0, (r + r_BB) * p_AA / 2 + ((r_BB + r_AB) * p_AB + (r_BB + r_ABg) * p_ABg) / 4,
                                     r_BB * p_BB + ((r_BB + r_AB) * p_AB + (r_BB + r_BBg) * p_BBg) / 4, 0, (1 - g) * (
                                             (r_BB + r_BgBg) * p_BgBg / 2 + (
                                             (r_BB + r_ABg) * p_ABg + (r_BB + r_BBg) * p_BBg) / 4), g * (
                                             (r_BB + r_BgBg) * p_BgBg / 2 + (
                                             (r_BB + r_ABg) * p_ABg + (r_BB + r_BBg) * p_BBg) / 4)],
                                    [(r + r_ABg) * p_AA / 4 + ((r_ABg + r_AB) * p_AB + 2 * r_ABg * p_ABg) / 8,
                                     (r_ABg + r_BB) * p_BB / 4 + ((r_ABg + r_AB) * p_AB + (r_ABg + r_BBg) * p_BBg) / 8, 0,
                                     (1 - g) * (((r_ABg + r) * p_AA + 2 * r_ABg * p_ABg + (r_ABg + r_BgBg) * p_BgBg) / 4 + (
                                             (r_ABg + r_AB) * p_AB + (r_ABg + r_BBg) * p_BBg) / 8), (1 - g) * (
                                             (r_ABg + r_BB) * p_BB / 4 + (
                                             (r_ABg + r_AB) * p_AB + (r_ABg + r_BBg) * p_BBg) / 8), g * (
                                             ((r_ABg + r) * p_AA + (r_ABg + r_AB) * p_AB + (r_ABg + r_BB) * p_BB) / 4 + (
                                             2 * r_ABg * p_ABg + (r_ABg + r_BBg) * p_BBg + (
                                             r_ABg + r_BgBg) * p_BgBg) / 4) + (r_ABg + r_BgBg) * p_BgBg / 4 + (
                                             (r_ABg + r_ABg) * p_ABg + (r_BBg + r_ABg) * p_BBg) / 8],
                                    [0, (r + r_BBg) * p_AA / 4 + ((r_BBg + r_AB) * p_AB + (r_BBg + r_ABg) * p_ABg) / 8,
                                     (r_BBg + r_BB) * p_BB / 4 + ((r_BBg + r_AB) * p_AB + 2 * r_BBg * p_BBg) / 8,
                                     (1 - g) * ((r + r_BBg) * p_AA / 4 + ((r_BBg + r_AB) * p_AB + (r_BBg + r_ABg) * p_ABg) / 8),
                                     (1 - g) * (((r_BBg + r_BB) * p_BB + (r_BBg + r_BgBg) * p_BgBg + 2 * r_BBg * p_BBg) / 4 + (
                                             (r_BBg + r_AB) * p_AB + (r_BBg + r_ABg) * p_ABg) / 8), g * (
                                             ((r_BBg + r) * p_AA + (r_BBg + r_AB) * p_AB + (r_BBg + r_BB) * p_BB) / 4 + (
                                             (r_BBg + r_ABg) * p_ABg + (r_BBg + r_BBg) * p_BBg + (
                                             r_BBg + r_BgBg) * p_BgBg) / 4) + (r_BBg + r_BgBg) * p_BgBg / 4 + (
                                             (r_BBg + r_ABg) * p_ABg + (r_BBg + r_BBg) * p_BBg) / 8],
                                    [0, 0, 0, (1 - g) * ((r + r_BgBg) * p_AA / 2 + (
                                            (r_BgBg + r_AB) * p_AB + (r_BgBg + r_ABg) * p_ABg) / 4), (1 - g) * (
                                             (r_BgBg + r_BB) * p_BB / 2 + (
                                             (r_BgBg + r_AB) * p_AB + (r_BgBg + r_BBg) * p_BBg) / 4), (g / 2) * (
                                             (r + r_BgBg) * p_AA + (r_AB + r_BgBg) * p_AB + (r_BB + r_BgBg) * p_BB + (
                                             r_ABg + r_BgBg) * p_ABg / 2 + (r_BBg + r_BgBg) * p_BBg / 2) + (
                                             (r_ABg + r_BgBg) * p_ABg + (r_BBg + r_BgBg) * p_BBg) / 4 + r_BgBg * p_BgBg]])

                # self-fertlization transition matrix
                genIn = np.matrix([[r, 0, 0, 0, 0, 0],
                                   [r_AB / 4, r_AB / 2, r_AB / 4, 0, 0, 0],
                                   [0, 0, r_BB, 0, 0, 0],
                                   [r_ABg / 4, 0, 0, r_ABg * (1 - g) / 2, 0, r_ABg * (g / 2 + 1 / 4)],
                                   [0, 0, r_BBg / 4, 0, r_BBg * (1 - g) / 2, r_BBg * (g / 2 + 1 / 4)],
                                   [0, 0, 0, 0, 0, r_BgBg]])
                #print(genOut[0].sum())

                gen2 = gen1 * ((1 - sig) * genOut + sig * inbr * genIn)  # calculate flux of genes into next generation
                # print(gen2)
                s_r = gen2.sum()  # total growth rate of population
                # print(s_r)
                matrix_R = ((
                                    1 - sig) * genOut + sig * inbr * genIn)  # form whole transition matrix to calculate adjusted growth rates

                # calculate individual adjusted flux
                sAA_AA_r = gen1[0] * matrix_R[0, 0]
                sAB_AA_r = gen1[1] * matrix_R[1, 0]
                sBB_AA_r = gen1[2] * matrix_R[2, 0]
                sABg_AA_r = gen1[3] * matrix_R[3, 0]
                sBBg_AA_r = gen1[4] * matrix_R[4, 0]
                sBgBg_AA_r = gen1[5] * matrix_R[5, 0]
                sAA_AB_r = gen1[0] * matrix_R[0, 1]
                sAB_AB_r = gen1[1] * matrix_R[1, 1]
                sBB_AB_r = gen1[2] * matrix_R[2, 1]
                sABg_AB_r = gen1[3] * matrix_R[3, 1]
                sBBg_AB_r = gen1[4] * matrix_R[4, 1]
                sBgBg_AB_r = gen1[5] * matrix_R[5, 1]
                sAA_BB_r = gen1[0] * matrix_R[0, 2]
                sAB_BB_r = gen1[1] * matrix_R[1, 2]
                sBB_BB_r = gen1[2] * matrix_R[2, 2]
                sABg_BB_r = gen1[3] * matrix_R[3, 2]
                sBBg_BB_r = gen1[4] * matrix_R[4, 2]
                sBgBg_BB_r = gen1[5] * matrix_R[5, 2]
                sAA_ABg_r = gen1[0] * matrix_R[0, 3]
                sAB_ABg_r = gen1[1] * matrix_R[1, 3]
                sBB_ABg_r = gen1[2] * matrix_R[2, 3]
                sABg_ABg_r = gen1[3] * matrix_R[3, 3]
                sBBg_ABg_r = gen1[4] * matrix_R[4, 3]
                sBgBg_ABg_r = gen1[5] * matrix_R[5, 3]
                sAA_BBg_r = gen1[0] * matrix_R[0, 4]
                sAB_BBg_r = gen1[1] * matrix_R[1, 4]
                sBB_BBg_r = gen1[2] * matrix_R[2, 4]
                sABg_BBg_r = gen1[3] * matrix_R[3, 4]
                sBBg_BBg_r = gen1[4] * matrix_R[4, 4]
                sBgBg_BBg_r = gen1[5] * matrix_R[5, 4]
                sAA_BgBg_r = gen1[0] * matrix_R[0, 5]
                sAB_BgBg_r = gen1[1] * matrix_R[1, 5]
                sBB_BgBg_r = gen1[2] * matrix_R[2, 5]
                sABg_BgBg_r = gen1[3] * matrix_R[3, 5]
                sBBg_BgBg_r = gen1[4] * matrix_R[4, 5]
                sBgBg_BgBg_r = gen1[5] * matrix_R[5, 5]

                # form birth-death transition matrix
                transitionMatrix = np.matrix([[1 - (sAA[i - 1] * (mu + (1 - mu) * rho)) / sAA[i - 1] + sAA_AA_r / (
                        sAA[i - 1] * s_r) * (s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sAA_AB_r / (sAA[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               0, sAA_ABg_r / (sAA[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               0, sAA_BgBg_r / (sAA[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death)],
                                              [sAB_AA_r / (sAB[i - 1] * s_r) * (
                                                      s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               1 - (sAB[i - 1] * (mu + (1 - mu) * rho * (1 - gamma * xi))) / sAB[
                                                   i - 1] + sAB_AB_r / (sAB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sAB_BB_r / (sAB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sAB_ABg_r / (sAB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sAB_BBg_r / (sAB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sAB_BgBg_r / (sAB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death)],
                                              [0, sBB_AB_r / (sBB[i - 1] * s_r) * (
                                                      s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               1 - (sBB[i - 1] * (mu + (1 - mu) * rho * (1 - xi))) / sBB[i - 1] + sBB_BB_r / (
                                                       sBB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               0, sBB_BBg_r / (sBB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sBB_BgBg_r / (sBB[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death)],
                                              [sABg_AA_r / (sABg[i - 1] * s_r) * (
                                                      s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sABg_AB_r / (sABg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               0, 1 - (sABg[i - 1] * (mu + (1 - mu) * rho * (1 - gamma * xi))) / sABg[
                                                   i - 1] + sABg_ABg_r / (sABg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sABg_BBg_r / (sABg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sABg_BgBg_r / (sABg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death)],
                                              [0, sBBg_AB_r / (sBBg[i - 1] * s_r) * (
                                                      s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sBBg_BB_r / (sBBg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sBBg_ABg_r / (sBBg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               1 - (sBBg[i - 1] * (mu + (1 - mu) * rho * (1 - xi))) / sBBg[
                                                   i - 1] + sBBg_BBg_r / (sBBg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sBBg_BgBg_r / (sBBg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death)],
                                              [0, 0, 0, sBgBg_ABg_r / (sBgBg[i - 1] * s_r) * (
                                                      s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               sBgBg_BBg_r / (sBgBg[i - 1] * s_r) * (
                                                       s_death * k / (s_death + (k - s_death) * np.exp(-s_r)) - s_death),
                                               1 - (sBgBg[i - 1] * (mu + (1 - mu) * rho * (1 - xi))) / sBgBg[
                                                   i - 1] + sBgBg_BgBg_r / (sBgBg[i - 1] * s_r) * (s_death * k / (
                                                       s_death + (k - s_death) * np.exp(-s_r)) - s_death)]])

                # calculate propensities
                """force = [1, 1, 1, 1, 1, 1]*transitionMatrix
                relativeForce = force/force.sum()
                #print(relativeForce)"""

                # calculate next generation of genotypes as number of individuals (not frequencies)
                nextGen = [sAA[i - 1], sAB[i - 1], sBB[i - 1], sABg[i - 1], sBBg[i - 1], sBgBg[i - 1]] * transitionMatrix

                # calculate cohort sizes in each generation i
                sAA[i] = nextGen[0, 0]
                sAB[i] = nextGen[0, 1]
                sBB[i] = nextGen[0, 2]
                sABg[i] = nextGen[0, 3]
                sBBg[i] = nextGen[0, 4]
                sBgBg[i] = nextGen[0, 5]

                """if 0 <= i%(gens/(u+1)) < 1 :
                    sBgBg[i] += 1/(u+1)"""

                # if sBgBg[i] > 50:
                # print(i)

                # find population size
                s[i] = sAA[i] + sAB[i] + sBB[i] + sABg[i] + sBBg[i] + sBgBg[i]

                # calculate resistance accumulation
                resist.append(resist[-1] + nhej * (1 - resist[-1]) / (1 - g + nhej * (1 - resist[-1])) * (
                        sABg[i] + sBBg[i] - sABg_death - sBBg_death) / s[i])



                #print(v)
                # print(sBgBg[i])

                # print(nextGen)

            t = scipy.linspace(0, gens, gens + 1)
            # transpose vectors for plotting and store as S1, S2, ..., S6
            S1 = scipy.transpose(sAA)
            S2 = scipy.transpose(sAB)
            S3 = scipy.transpose(sBB)
            S4 = scipy.transpose(sABg)
            S5 = scipy.transpose(sBBg)
            S6 = scipy.transpose(sBgBg)

            endpoints.append(S6[-1]/100)

            list1.append(sig)
            list2.append(x)
            list3.append(S6[-1]/100)

    df = pd.DataFrame(list(zip(list1, list2, list3)),
                      columns=['var', 'Selfing', 'GD'])

    projPlot = pl.figure()
    pl.ylim(-0.05, 1.05)
    #pl.xlim(0.8, 5)
    x1 = np.arange(0.0, 0.28, 0.01)
    x2 = np.arange(0.0, 0.42, 0.01)
    x3 = np.arange(0.81, 0.85, 0.01)
    x4 = np.arange(0.71, 0.97, 0.01)
    x5 = np.arange(0.79, 1, 0.01)
    x6 = np.arange(0.0, 0.01, 0.001)
    pl.fill_between(x1, 0.75, 0.85, alpha=0.3, facecolor="blue")
    pl.fill_between(x2, 0.65, 0.75, alpha=0.3, facecolor="#9b59b6")
    pl.fill_between(x3, 0.35, 0.45, alpha=0.3, facecolor="#9b59b6")
    pl.fill_between(x4, 0.55, 0.65, alpha=0.3, facecolor="blue")
    pl.fill_between(x5, 0.45, 0.55, alpha=0.3, facecolor="blue")
    pl.fill_between(x6, 0.85, 0.95, alpha=0.3, facecolor="red")


    pl.axvline(x=0.5, ymin=0, ymax=1.05, linestyle=':')
    # alternative sensitivity plot
    #t = scipy.linspace(0, 1, numEndpoints + 1)
    sns.lineplot(data=df, x="var", y="GD")
    sns.despine()
    #pl.xticks(scipy.linspace(0, 1, 21))
    #pl.grid()
    #pl.legend(loc='best')
    pl.xlabel('Selfing Rate')
    pl.ylabel('Gene Drive Frequency in 10 Yrs')
    #pl.title('Endpoint Analysis')
    projPlot.savefig('C:/Users/Richard/Downloads/GeneDriveFigures/figure_self2.png', bbox_inches = "tight")


    pl.show()


main()