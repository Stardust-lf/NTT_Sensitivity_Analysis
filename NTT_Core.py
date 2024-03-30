# from sympy import ntt, intt
from sympy.utilities.iterables import ibin, iterable
from sympy.utilities.misc import as_int
from sympy.ntheory import isprime, primitive_root
import random

class RandomFlipper:
    def __init__(self, rand_prob = 0.1):
        self.flips = 0
        self.prob = rand_prob
    def flip_random_bit(self, n):
        binary = bin(n)[2:]
        index = random.randint(0, len(binary) - 1)
        flipped = binary[:index] + ('0' if binary[index] == '1' else '1') + binary[index + 1:]
        #print("{} --> {}".format(binary, flipped))

        self.flips += 1
        return int(flipped, 2)

    def flip_chosen_bit(self, n, bit_pos):
        assert bit_pos>0

        binary = bin(n)[2:].zfill(bit_pos + 1)
        index = -bit_pos
        flipped = binary[:index] + ('0' if binary[index] == '1' else '1') + binary[index + 1:]
        print("{} --> {}".format(binary, flipped))

        self.flips += 1
        return int(flipped, 2)

    def randomize(self, n):
        if random.random() < self.prob:
            #print("---Error Flip Occurs---")
            return self.flip_random_bit(n)
        else:
            return n


def number_theoretic_transform(seq, prime, inverse=False, log=True):
    """Utility function for the Number Theoretic Transform"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of integer coefficients "
                        "for Number Theoretic Transform")

    p = as_int(prime)
    if not isprime(p):
        raise ValueError("Expected prime modulus for "
                         "Number Theoretic Transform")

    a = [as_int(x) % p for x in seq]

    n = len(a)
    if n < 1:
        return a

    b = n.bit_length() - 1
    if n & (n - 1):
        b += 1
        n = 2 ** b

    if (p - 1) % n:
        raise ValueError("Expected prime modulus of the form (m*2**k + 1)")

    a += [0] * (n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        # print("j: ", j)
        if i < j:
            a[i], a[j] = a[j], a[i]

    pr = primitive_root(p)

    rt = pow(pr, (p - 1) // n, p)

    if inverse:
        rt = pow(rt, p - 2, p)

    w = [1] * (n // 2)
    for i in range(1, n // 2):
        w[i] = w[i - 1] * rt % p
    if log:
        print("RT", rt)
        print("Primitive_root:", pr)
        print('After switch a:', a)
        print("W: ", w)
    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                if log: print("a in iteration", a)
                u, v = a[i + j], a[i + j + hf] * w[ut * j]

                a[i + j], a[i + j + hf] = (u + v) % p, (u - v) % p
                if log: print(
                    "Operating a{} = (a{} + a{}w{})mod p = {}".format((i + j), (i + j), (i + j + hf), (ut * j),
                                                                      a[i + j]))
                if log: print(
                    "Operating a{} = (a{} - a{}w{})mod p = {}".format((i + j + hf), (i + j), (i + j + hf), (ut * j),
                                                                      a[i + j + hf]))

        h *= 2

    if inverse:
        rv = pow(n, p - 2, p)
        a = [x * rv % p for x in a]

    return a


def randomized_ntt(seq, prime, inverse=False, rand_prob=0.1, log=True):
    flipper = RandomFlipper(rand_prob=rand_prob)
    """Utility function for the Number Theoretic Transform"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of integer coefficients "
                        "for Number Theoretic Transform")

    p = as_int(prime)
    if not isprime(p):
        raise ValueError("Expected prime modulus for "
                         "Number Theoretic Transform")

    a = [as_int(x) % p for x in seq]

    n = len(a)
    if n < 1:
        return a

    b = n.bit_length() - 1
    if n & (n - 1):
        b += 1
        n = 2 ** b

    if (p - 1) % n:
        raise ValueError("Expected prime modulus of the form (m*2**k + 1)")

    a += [0] * (n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j], a[i]
    pr = primitive_root(p)
    rt = pow(pr, (p - 1) // n, p)
    if inverse:
        rt = pow(rt, p - 2, p)

    w = [1] * (n // 2)
    for i in range(1, n // 2):
        w[i] = w[i - 1] * rt % p
    if log:
        print("PM: ", pr)
        print("RT: ", rt)
        print("W: ", w)
    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                if log: print("a in iteration", a)
                u, v = a[i + j], a[i + j + hf] * w[ut * j]

                a[i + j] = (u + v) % p
                if log: print(
                    "Operating a{} = (a{} + a{}w{})mod p = {}".format((i + j), (i + j), (i + j + hf), (ut * j),
                                                                      a[i + j]))
                a[i + j] = flipper.randomize(a[i + j])
                a[i + j + hf] = (u - v) % p
                if log: print(
                    "Operating a{} = (a{} - a{}w{})mod p = {}".format((i + j + hf), (i + j), (i + j + hf), (ut * j),
                                                                      a[i + j + hf]))
                a[i + j + hf] = flipper.randomize(a[i + j + hf])

        h *= 2

    if inverse:
        rv = pow(n, p - 2, p)
        a = [x * rv % p for x in a]
    if log:
        print("--------------------------------------------------------------------------------")
        print("Flip occurs {} times".format(flipper.flips))
        print("--------------------------------------------------------------------------------")
    return a

def flip_index_ntt(seq, prime, flip_pos, flip_opt_time, inverse=False , log=True):
    assert flip_opt_time >= 1
    flipper = RandomFlipper()
    """Utility function for the Number Theoretic Transform"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of integer coefficients "
                        "for Number Theoretic Transform")

    p = as_int(prime)
    if not isprime(p):
        raise ValueError("Expected prime modulus for "
                         "Number Theoretic Transform")

    a = [as_int(x) % p for x in seq]

    n = len(a)
    if n < 1:
        return a

    b = n.bit_length() - 1
    if n & (n - 1):
        b += 1
        n = 2 ** b

    if (p - 1) % n:
        raise ValueError("Expected prime modulus of the form (m*2**k + 1)")

    a += [0] * (n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j], a[i]
    pr = primitive_root(p)
    rt = pow(pr, (p - 1) // n, p)
    if inverse:
        rt = pow(rt, p - 2, p)

    w = [1] * (n // 2)
    for i in range(1, n // 2):
        w[i] = w[i - 1] * rt % p
    if log:
        print("PM: ", pr)
        print("RT: ", rt)
        print("W: ", w)
    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                if log: print("a in iteration", a)
                u, v = a[i + j], a[i + j + hf] * w[ut * j]

                a[i + j] = (u + v) % p
                if log: print(
                    "Operating a{} = (a{} + a{}w{})mod p = {}".format((i + j), (i + j), (i + j + hf), (ut * j),
                                                                      a[i + j]))
                if flip_opt_time == 1:
                    print("Previous",a[i+j])
                    a[i + j] = flipper.flip_chosen_bit(a[i + j],flip_pos)
                    print("After",a[i+j])

                flip_opt_time -= 1
                a[i + j + hf] = (u - v) % p
                if log: print(
                    "Operating a{} = (a{} - a{}w{})mod p = {}".format((i + j + hf), (i + j), (i + j + hf), (ut * j),
                                                                      a[i + j + hf]))
                if flip_opt_time == 1:
                    a[i + j + hf] = flipper.flip_chosen_bit(a[i + j + hf], flip_pos)
                flip_opt_time -= 1

        h *= 2

    if inverse:
        rv = pow(n, p - 2, p)
        a = [x * rv % p for x in a]
    if log:
        print("--------------------------------------------------------------------------------")
        print("Flip occurs {} times".format(flipper.flips))
        print("--------------------------------------------------------------------------------")
    return a
