import numpy as np
import scipy
import seaborn as sns; sns.set(style="white", color_codes=True)
import matplotlib.pyplot as pl
from scipy.integrate import quad

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e",
                  "#2ecc71"]
sns.set_palette(sns.color_palette(flatui))
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (7, 5),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pl.rcParams.update(params)


gd = [1,1,.000000000000000001]
ge = [0,.6,.6]
worming = []
lab = ['Gene Drive', 'Gene Drive + MDA', 'MDA']
projPlot = pl.figure()

for y in range(len(gd)):
    #integration function for neg binomial mating
    def integrand(theta, alpha, k):
        return (1-np.cos(theta))/(1+alpha*np.cos(theta))**(1+k)


    def disease(rho, quarter, worm, snail):
        # parameter block
        gens = 1
        print(rho)
        reinfectionConst = 10**3
        treatmentEfficacy = ge[y]


        if (quarter-1)%4 == 0:
            m = (1-treatmentEfficacy)*worm

        else:
            m = worm  # equilibrium mean worm burden

        k = 0.24  # clumping parameter (also known as r in negative binomial)
        b = 1.465 * 10 ** -8 / .48811346
        a = .0142 # 50 cercariae shed per day per snail * 7 days per week * prob of cerc contacting human and maturing
        N = 10 ** 4  # number of total snails
        mu_1 = .004  # + .02666667 #weekly death rate of worms per capita
        mu_2 = 0.25  # weekly death rate of infected snails per capita
        H = 1000  # number of humans
        #rho = 0.5  # fraction of resistant snails [0-1]

        # begin simulation
        X0 = scipy.array([m, snail, 0])  # initials conditions: x0=10  and y0=5
        #print(X0)
        t = scipy.linspace(0, 14 * gens, 14 * gens)

        def dX_dt(X, t=0):  # specify the initial time point t0=0

            alpha = X[0] / (X[0] + k)
            part1 = ((1 - alpha) ** (1 + k)) / (2 * np.pi)
            part2 = quad(integrand, 0, 2 * np.pi, args=(alpha, k))
            phi = 1 - part1 * part2[0]
            w = 0.5 * phi
            beta = b * w
            # beta = 9.5177852*10**(-8)
            #print(beta)
            # y = scipy.array([a * N * X[1] - mu_1 * X[0], beta * H * X[0] * (1 - X[1] - rho) - mu_2 * X[1]])
            #y = scipy.array([a * N * X[1] - mu_1 * X[0], beta * H * X[0] * (1 - X[1] - rho) - mu_2 * X[1], beta * H * X[0] * (1 - X[1] - rho)])
            y = scipy.array([a * N * X[1] - mu_1 * X[0], .0104*(1-np.exp(-reinfectionConst*beta * H * X[0])) * (1 - X[1] - rho) - mu_2 * X[1],
                             .0104 * (1 - np.exp(-reinfectionConst*beta * H * X[0]))])
            return y

        X, infodict = scipy.integrate.odeint(dX_dt, X0, t, full_output=True)

        matedPairs = []
        snailprevalence = []
        humanPrev = []

        for u in range(len(X)):

            alpha_ = X[u, 0] / (X[u, 0] + k)
            part1_ = ((1 - alpha_) ** (1 + k)) / (2 * np.pi)
            part2_ = quad(integrand, 0, 2 * np.pi, args=(alpha_, k))
            phi_ = 1 - part1_ * part2_[0]
            w = 0.5 * phi_
            matedPairs.append(w * X[u, 0])
            sigma = 1 - 2 * (1 + X[u, 0] / (2 * k)) ** -k + (1 + X[u, 0] / k) ** -k
            snailprevalence.append(100 * X[u, 1])
            humanPrev.append(sigma)

        foi = X[-1,2]

        return foi, X[1:,0], X[1:,1], humanPrev[1:]


    def main():


        # genotypes must be greater than 0
        sAA = [50]  # set initial population size for wild type AA
        sAB = [34]  # set initial population size for heterozygous AB
        sBB = [16]  # set initial population size for homozygous resistant BB
        sABg = [0.0000000000000000000000000000001]  # set initial population size for heterozygous engineered resistant ABg
        sBBg = [0.0000000000000000000000000000001]  # set intial population size for homozygous resistant, heterozygous engineered copy BBg
        sBgBg = [gd[y]]  # set initial population size for homozygous engineered resistant BgBg
        s = [100]  # set initial total population size
        sig = 0.5  # selfing rate
        k = 100  # carrying capacity
        r = 20  # intrinsic growth rate
        inbr = .7  # cost of inbreeding
        mu = 0.5  # natural death rate
        eps = 0.0  # migration parameter (place holder for now)
        #rho = 0.15  # force of infection
        beta = 0.2835 * r  # reduction in fecundity due to natural resistance
        gamma = 1  # dominance coefficient
        xi = 0.8  # conferred resistance to infection
        g = 0.9  # gene drive efficiency
        beta_g = 0.1 * r  # reduction in fecundity due to engineered resistance
        r_AB = r - gamma * beta  # fecundity of AB genotype
        r_BB = r - beta  # " BB genotype
        r_ABg = r - gamma * beta - beta_g  # " ABg genotype
        r_BBg = r - beta - beta_g  # " BBg genotype
        r_BgBg = r - beta - 2 * beta_g  # " BgBg genotype
        gens = 4*40  # no. of generations
        worms = [710]
        incidence = [.8]
        snailPrev = [.02]

        # begin simulation
        for i in range(1, gens + 1):
            rho, matedPairs, snailprevalence, humanPrev = disease(1 - sAA[i-1] / s[i-1], i, worms[-1], snailPrev[-1])
            print(snailPrev[-1])
            print(rho)
            worms = [*worms, *matedPairs]
            incidence = [*incidence, *humanPrev]
            snailPrev = [*snailPrev, *snailprevalence]

            # add 0 element to arrays as place holder for generation i (aka generation t+1)
            sAA.append(0)
            sAB.append(0)
            sBB.append(0)
            sABg.append(0)
            sBBg.append(0)
            sBgBg.append(0)
            s.append(0)

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



            # calculate next generation of genotypes as number of individuals (not frequencies)
            nextGen = [sAA[i - 1], sAB[i - 1], sBB[i - 1], sABg[i - 1], sBBg[i - 1], sBgBg[i - 1]] * transitionMatrix

            # calculate cohort sizes in each generation i
            sAA[i] = nextGen[0, 0]
            sAB[i] = nextGen[0, 1]
            sBB[i] = nextGen[0, 2]
            sABg[i] = nextGen[0, 3]
            sBBg[i] = nextGen[0, 4]
            sBgBg[i] = nextGen[0, 5]

            # if sBgBg[i] > 50:
            # print(i)

            # find population size
            s[i] = sAA[i] + sAB[i] + sBB[i] + sABg[i] + sBBg[i] + sBgBg[i]

            """transpMatrix = transitionMatrix.transpose()
            w, v = LA.eig(transitionMatrix)

            print(w)"""
            # print(sBgBg[i])

            # print(nextGen)



        t = scipy.linspace(0, gens, gens + 1)
        t2 = scipy.linspace(0, gens, 13*gens +1)
        # transpose vectors for plotting and store as S1, S2, ..., S6
        S1 = scipy.transpose(sAA)
        S2 = scipy.transpose(sAB)
        S3 = scipy.transpose(sBB)
        S4 = scipy.transpose(sABg)
        S5 = scipy.transpose(sBBg)
        S6 = scipy.transpose(sBgBg)

        worming.append(worms)


        pl.plot(t2, np.array(worming[y])/710, label=lab[y], linewidth = '3')







    main()

pl.legend(loc='best')
pl.xlabel('Years')
pl.ylabel('Fraction of Pre-treatment MWB')
#labels = [0,1,2,3,4,5,6,7,8,9,10]
labels = [0,4,8,12,16,20,24,28,32,36,40]
pl.ylim(-0.05, 1.05)
pl.xticks(scipy.linspace(0, 4*40, 11),labels)  # sets ticks for x-axis
sns.despine()
projPlot.savefig('C:/Users/Richard/Downloads/GeneDriveFigures/figure_MDATreat_long2.png', bbox_inches="tight")
pl.show()