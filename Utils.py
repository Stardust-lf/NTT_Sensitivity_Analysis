from NTT_Core import number_theoretic_transform, randomized_ntt, flip_index_ntt
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font',family='Times New Roman')


def hamming_distance_vector(vector1, vector2):
    # print(vector1, vector2)
    total_distance = []
    for num1, num2 in zip(vector1, vector2):
        bin_num1 = bin(int(num1))[2:]
        bin_num2 = bin(int(num2))[2:]

        len1 = len(bin_num1)
        len2 = len(bin_num2)

        if len1 < len2:
            bin_num1 = bin_num1.zfill(len2)
        elif len1 > len2:
            bin_num2 = bin_num2.zfill(len1)

        total_distance.append(sum(c1 != c2 for c1, c2 in zip(bin_num1, bin_num2)))

    return total_distance


def l1_norm_vector(vector1, vector2):
    return np.array(vector1) - np.array(vector2)


def l2_norm_vector(vector1, vector2):
    return np.square(np.array(vector1) - np.array(vector2))


# class RandomizedTest:
#     def __init__(self, prime_num, epochs, length, scale):
#         self.epochs = epochs
#         self.scale = scale
#         self.samples = np.random.randint(0, scale, size=[epochs, length])
#         self.p = prime_num
#
#     def __call__(self, ntt_fun = number_theoretic_transform, dist_fun = hamming_distance_vector, show_figure = False):
#         error = []
#         for i in range(self.epochs):
#             seq = self.samples[i]
#             correct = number_theoretic_transform(seq, self.p, log=False)
#             wrong = ntt_fun(seq, self.p, log=False)
#             error.append(dist_fun(correct, wrong))
#         if show_figure:
#             fig, ax = plt.
#         return error


if __name__ == '__main__':
    a = np.array([1,2,3])
    b = np.array([-3,-2,-1])
    print(n1_norm(a, b))
    print(hamming_distance_vector(a, b))



