import scipy.stats as ss


class LinearCongruentialGenerator:
    def __init__(self, A=69069, C=858993459, M=4294967296):
        self.A = A
        self.C = C
        self.M = M
        self.previous = 79546

    def generate_next(self):
        self.previous = (self.A * self.previous + self.C) % self.M

    def get_next(self):
        self.generate_next()
        return self.previous

    def get_m(self):
        return self.M

    def reset(self):
        self.previous = 79546


class MultiplicativeCongruentialGenerator:
    def __init__(self, A=283741, M=524287):
        self.A = A
        self.M = M
        self.previous = 58768

    def generate_next(self):
        self.previous = (self.A * self.previous) % self.M

    def get_next(self):
        self.generate_next()
        return self.previous

    def get_m(self):
        return self.M

    def reset(self):
        self.previous = 58768


class CombinedCongruentialGenerators:
    def __init__(self):
        self.LCG = LinearCongruentialGenerator()
        self.MCG = MultiplicativeCongruentialGenerator()
        self.LCG_BUF = [self.LCG.get_next() for i in range(1000)]

    def get_next(self):
        index = self.MCG.get_next() * (len(self.LCG_BUF) / self.MCG.get_m())
        output = self.LCG_BUF[int(index)]
        self.LCG_BUF[int(index)] = self.LCG.get_next()
        return output

    def get_m(self):
        return self.LCG.get_m()

    def reset(self):
        self.LCG.reset()
        self.MCG.reset()
        self.LCG_BUF = [self.LCG.get_next() for i in range(1000)]


def x_square(generator):
    generator.reset()
    buf_size = 1000
    niches_size = 20
    buf = [generator.get_next() for i in range(buf_size)]
    niches = [0 for i in range(niches_size)]

    for i in range(buf_size):
        index = buf[i] / (generator.get_m() / niches_size)
        niches[int(index)] += 1

    expected = buf_size / niches_size
    output = 0

    for niche in niches:
        output += ((niche - expected) * (niche - expected) / expected)

    return output


def serial_correlation(generator):
    generator.reset()
    buf_size = 100000
    buf = [generator.get_next() for i in range(buf_size)]
    sum_m = buf[0] * buf[buf_size - 1]
    sum_of_square = buf[buf_size - 1] * buf[buf_size - 1]
    summa = buf[buf_size - 1]

    for i in range(buf_size - 1):
        sum_m += buf[i] * buf[i + 1]
        sum_of_square += buf[i] * buf[i]
        summa += buf[i]

    summa *= summa

    output = (buf_size * sum_m - summa) / (buf_size * sum_of_square - summa)
    return output

if __name__ == '__main__':
    LCG = LinearCongruentialGenerator()
    print('LinearCongruentialGenerator -> ', [LCG.get_next() for i in range(100)])

    MCG = MultiplicativeCongruentialGenerator()
    print('MultiplicativeCongruentialGenerator -> ', [MCG.get_next() for i in range(100)])

    CCG = CombinedCongruentialGenerators()
    print('CombinedCongruentialGenerators -> ', [CCG.get_next() for i in range(100)])

    print('\n\n')
    Gen = [LCG, MCG, CCG]
    for gen in Gen:
        res = x_square(generator=gen)
        print('X Square for ', gen.__class__.__name__, ' result -> ', res)

    print('\n')
    for gen in Gen:
        res = serial_correlation(generator=gen)
        print('Serial Correlation ', gen.__class__.__name__, ' result -> ', res)

    print('\n')
    for gen in Gen:
        g = [gen.get_next() for i in range(20000)]
        g1 = [gen.get_next() for i in range(20000)]
        res = ss.ks_2samp(g, g1)
        print('K-S ', gen.__class__.__name__, ' result -> ', res)