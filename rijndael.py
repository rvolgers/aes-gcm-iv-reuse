from gf2poly import *

# irreducible polynomial used to implement GF(2**8)
RD_POLY = (1 << 8) | (1 << 4) | (1 << 3) | (1 << 1) | 1
assert RD_POLY == 0b1_0001_1011

RD_GEN = 0b10

def rd_inverse(x):
    return gf_invmod(x, RD_POLY)

def rd_mul(x, y):
    return gf_mod(gf_mul_noreduce(x, y), RD_POLY)

def rd_pow(x, e):
    return gf_modexp(x, e, RD_POLY)

# linear map and reverse linear map
assert gf_invmod(0b1_1111, 0b1_0000_0001) == 0b100_1010

def sbox_linear_map(x):
    return gf_mod(gf_mul_noreduce(x, 0b1_1111), 0b1_0000_0001)

def sbox_inverse_linear_map(x):
    return gf_mod(gf_mul_noreduce(tmp, 0b100_1010), 0b1_0000_0001)

RD_SBOX = []
for i in range(0, 256):
    tmp = i
    tmp = rd_inverse(tmp) if tmp != 0 else 0
    tmp = sbox_linear_map(tmp)
    tmp ^= 0b110_0011
    RD_SBOX.append(tmp)
    # print(f"{hex(i)}: {hex(tmp)}")

# See Table 4
# https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.197-upd1.pdf
assert RD_SBOX[:16] == [0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76]

# sbox linear map as matrix
M = []
for i in range(8):
    M.append(sbox_linear_map(1 << i))

assert M == [
    0b00011111,
    0b00111110,
    0b01111100,
    0b11111000,
    0b11110001,
    0b11100011,
    0b11000111,
    0b10001111,
]

# use matrix to compute linear map
for x in [0xbb, 0x47, 0x20]: # arbitary test values
    acc = 0
    for i in range(8):
        if x & (1 << i):
            acc ^= M[i]
    assert acc == sbox_linear_map(x)

# these eight values are invariant under sbox_linear_map:
# 0b0, 0b110011, 0b1010101, 0b1100110, 0b10011001, 0b10101010, 0b11001100, 0b11111111
tmp = gfmat_kernel(gfmat_minus_identity(M[:]))
assert tmp == [0b110011, 0b1010101, 0b10011001]
for i in range(2**3):
    acc = 0
    for j in range(3):
        if i & (1 << j):
            acc ^= [0b110011, 0b1010101, 0b10011001][j]
    assert sbox_linear_map(acc) == acc

# find x such that Ax = b, elements of A in the rijndael field
def rdmat_gaussian_elimination(A, b):
    n = len(A)

    # augment matrix
    A = [row + [bi] for row, bi in zip(A, b)]

    # iterate over diagonal
    for i in range(n):
        # ensure A[i][i] != 0, by swapping rows if needed
        for j in range(i, n):
            if A[j][i] != 0:
                A[i], A[j] = A[j], A[i]  # (no effect if i == j)
                break
        else: return None

        # ensure A[i][i] == 1, by dividing by the leading coefficient
        lc_inv = rd_inverse(A[i][i])
        A[i] = [rd_mul(x, lc_inv) for x in A[i]]

        # ensure this column is zero in all other rows,
        # by subtracting multiples of the current row
        for j in range(n):
            if j == i: continue
            factor = A[j][i]
            if factor == 0: continue  # (optimization only)
            for k in range(i, n + 1):
                A[j][k] ^= rd_mul(factor, A[i][k]) 

    return [row[-1] for row in A]

basis = [2**i for i in range(8)]
# note that this is equal to the M we already calculated above
outputs = [sbox_linear_map(x) for x in basis]
moore_matrix = [
    [rd_pow(x, 2**i) for i in range(8)]
    for x in basis
]
coeffs = rdmat_gaussian_elimination(moore_matrix, outputs)

# coefficients of sbox_linear_map as a polynomial with exponents 2**i
# the idea that such a polynomial exists comes from "A simple algebraic representation of Rijndael"
# https://web.archive.org/web/20061215074718/http://www.macfergus.com/pub/rdalgeq.pdf
# (section 2)
# finding the coefficients of this polynomial by using gaussian elimination
# with the moore matrix was suggested by chatgpt.
assert coeffs == [5, 9, 249, 37, 244, 1, 181, 143]
assert coeffs == [0b101, 0b1001, 0b11111001, 0b100101, 0b11110100, 0b1, 0b10110101, 0b10001111]

for x in range(256):
    acc = 0
    for i, c in enumerate(coeffs):
        acc ^= rd_mul(rd_pow(x, 2**i), c)
    assert acc == sbox_linear_map(x)

RD_INV_SBOX = []
for i in range(0, 256):
    tmp = i
    tmp ^= 0b110_0011
    tmp = sbox_inverse_linear_map(tmp)
    tmp = gf_invmod(tmp, RD_POLY) if tmp != 0 else 0
    RD_INV_SBOX.append(tmp)

assert all(RD_SBOX[RD_INV_SBOX[i]] == i for i in range(256))

# AES will use at most 10 of these
ROUND_CONSTANTS = [gf_modexp(RD_GEN, i, RD_POLY) for i in range(29)]

assert ROUND_CONSTANTS[:10] == [0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36]

def xor4(w, v):
    return [wi ^ vi for wi,vi in zip(w,v)]

def wstr(w):
    return f"{w[0]:02x}{w[1]:02x}{w[2]:02x}{w[3]:02x}"

def w4str(w4):
    return '[' + ','.join(wstr(w) for w in w4) + ']'

def key_schedule(k):
    k = [k[i:i+4] for i in range(0, len(k), 4)]

    n = len(k)

    # number of rounds for 128/192/256 bit key
    rounds = {4: 10, 6: 12, 8: 14}[n]

    # key schedule starts with the key itself
    ks = [ki[:] for ki in k]

    # we need 4 words per round, plus 4 more
    for i in range(len(k), 4 * (rounds + 1)):
        # print(str(i))
        tmp = ks[i-1]
        # print(wstr(tmp))
        if i >= n and i % n == 0:
            tmp = tmp[1:] + tmp[:1]  # rotate by 1
            # print(wstr(tmp))
            tmp = [RD_SBOX[b] for b in tmp]
            # print(wstr(tmp))
            tmp = xor4(tmp, [ROUND_CONSTANTS[i // n - 1], 0, 0, 0])
            # print(wstr(tmp))
        elif i >= n and n > 6 and i % n == 4:
            tmp = [RD_SBOX[b] for b in tmp]
            # print(wstr(tmp))
        tmp = xor4(tmp, ks[i-n])
        # print(wstr(tmp))
        ks.append(tmp)

    return (rounds, ks)

# test from the standard, appendix A.1
# https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.197-upd1.pdf
k = [0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0x3c]
assert key_schedule(k)[1][-1] == [0xb6, 0x63, 0x0c, 0xa6]

def sub_bytes(x):
    return [
        [RD_SBOX[b] for b in w]
        for w in x
    ]

def shift_rows(s):
    return [
        [s[0][0], s[1][1], s[2][2], s[3][3]],
        [s[1][0], s[2][1], s[3][2], s[0][3]],
        [s[2][0], s[3][1], s[0][2], s[1][3]],
        [s[3][0], s[0][1], s[1][2], s[2][3]],
    ]

def mix_column(w):
    w2 = [rd_mul(b, 0b10) for b in w]
    w3 = xor4(w, w2)
    return [
        w2[0] ^ w3[1] ^ w[2] ^ w[3],
        w[0] ^ w2[1] ^ w3[2] ^ w[3],
        w[0] ^ w[1] ^ w2[2] ^ w3[3],
        w3[0] ^ w[1] ^ w[2] ^ w2[3],
    ]

def mix_columns(s):
    return [mix_column(w) for w in s]

def add_round_key(s, ks):
    return [xor4(si, ksi) for si, ksi in zip(s, ks)]

def cipher(data, rounds_and_ks):
    (rounds, ks) = rounds_and_ks
    s = [data[i:i+4] for i in range(0, len(data), 4)]
    # print(w4str(s))
    s = add_round_key(s, ks[0:4])
    # print(w4str(s))
    for r in range(1, rounds):
        s = sub_bytes(s)
        # print(w4str(s))
        s = shift_rows(s)
        # print(w4str(s))
        s = mix_columns(s)
        # print(w4str(s))
        s = add_round_key(s, ks[4*r:4*r+4])
        # print(w4str(s))
    s = sub_bytes(s)
    # print(w4str(s))
    s = shift_rows(s)
    # print(w4str(s))
    s = add_round_key(s, ks[4*rounds:4*rounds+4])
    # print(w4str(s))
    return s[0] + s[1] + s[2] + s[3]

data = [0x32,0x43,0xf6,0xa8,0x88,0x5a,0x30,0x8d,0x31,0x31,0x98,0xa2,0xe0,0x37,0x07,0x34]

assert cipher(data, key_schedule(k)) == [0x39,0x25,0x84,0x1d,0x02,0xdc,0x09,0xfb,0xdc,0x11,0x85,0x97,0x19,0x6a,0x0b,0x32]

class Rijndael:
    def __init__(self, key):
        self.rounds_and_ks = key_schedule(key)

    def encrypt_block(self, data):
        assert len(data) == 16
        return bytes(cipher(list(data), self.rounds_and_ks))
