#!/usr/bin/env python3

# to run with sage support:
# PYTHONPATH=/usr/lib/python3/dist-packages ~/src/sage/sage ./gcm.py
# (the PYTHONPATH is needed to pick up modules like cryptography)

# only used for a basic AES-ECB primitive
from collections import Counter
from copy import deepcopy
from math import gcd, isqrt
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

import itertools
from random import getrandbits
from time import time
from functools import reduce
import operator

# count trailing zeroes
def count_trailing_zeros(x):
    # implementations differ in how/if ctz(0) is defined, so avoid it
    assert x != 0

    # (x ^ (x - 1)) isolates the lowest 1 bit
    # 0b110000 - 1 = 0b101111
    # 0b110000 ^ (0b110000 - 1) = 0b10000

    return (x ^ (x - 1)).bit_length() - 1

##############################################
# code for operating on numbers in GF(2^128) #
##############################################


# GF(2^128) with polynomial x^128 + x^7 + x^2 + x + 1
GF_POLY = (1 << 128) | (1 << 7) | (1 << 2) | (1 << 1) | 1

# generator
GF_X = 2  # f(x) = x + 0
GF_GEN = GF_X

# used by _gf_sqrt
MASK128_64 = 0x00000000_00000000_ffffffff_ffffffff
MASK128_32 = 0x00000000_ffffffff_00000000_ffffffff
MASK128_16 = 0x0000ffff_0000ffff_0000ffff_0000ffff
MASK128_8 = 0x00ff00ff_00ff00ff_00ff00ff_00ff00ff

# used by _gf_bitswap and gf_sqrt
MASK128_4 = 0x0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f  # 0b00001111 x 16
MASK128_2 = 0x33333333_33333333_33333333_33333333  # 0b00110011 x 16
MASK128_1 = 0x55555555_55555555_55555555_55555555  # 0b01010101 x 16

# used by gf_from_bytes / gf_to_bytes.
# reverses the order of bits within each byte of a 128 bit value.
# the spec mandates this unusual bit order for some reason.
# it's possible to avoid this swapping by painstakingly adjusting all
# operations to account for this, but that's annoying.
def _gf_bitswap(x):
    x = ((x & MASK128_4) << 4) | ((x >> 4) & MASK128_4)
    x = ((x & MASK128_2) << 2) | ((x >> 2) & MASK128_2)
    x = ((x & MASK128_1) << 1) | ((x >> 1) & MASK128_1)
    return x

def gf_from_bytes(b):
    assert len(b) == 16, "argument must be 16 bytes"
    return _gf_bitswap(int.from_bytes(b, byteorder='little'))

def gf_to_bytes(x):
    assert x < (1 << 128), "element not properly reduced"
    return _gf_bitswap(x).to_bytes(length=16, byteorder='little')

def gf_random():
    return getrandbits(128)

# non-reducing version of gf_mul, used for experimentation
def gf_mul_noreduce(x, y):
    result = 0
    while True:
        if y & 1:
            result ^= x

        y >>= 1
        if y == 0:
            break

        x <<= 1

    return result

ALL64 = (1 << 64) - 1
ALL128 = (1 << 128) - 1

# instrumented version of gf_reduce to show how input bits affect outputs
def gf_reduce_instrumented(x):
    # - only bits in the upper 128 can cause a modular reduction
    # - the upper bit of the modular reduction does not matter,
    #   since it is only updated after it is checked and also
    #   it doesn't end up in the output
    # - 

    # x_sym = [1 << i for i in range(256)]

    # version assuming all odd bits of x are zero, i.e. x is a square
    x_sym = [1 << i if i & 1 == 0 else 0 for i in range(256)]

    sym_shl = lambda x, i: [0] * i + x
    sym_shr = lambda x, i: x[i:]
    sym_xor = lambda a, b: [i ^ j for i,j in itertools.zip_longest(a, b, fillvalue=0)]
    sym_str = lambda x: ' ^ '.join(f"x{i}" for i in range(256) if (x >> i) & 1)

    A = x >> (256 - 1)
    A_sym = sym_shr(x_sym, 256 - 1)
    B = x >> (256 - 2)
    B_sym = sym_shr(x_sym, 256 - 2)
    C = x >> (256 - 7)
    C_sym = sym_shr(x_sym, 256 - 7)

    # bits 128-254 affect bits 0-127 via poly bit 0
    D = x >> (256 - 128)
    D_sym = sym_shr(x_sym, 256 - 128)

    # a polynomial multiplication of two 128 bit values is at most 255 bits
    assert A == 0

    X3D = A ^ B ^ C ^ D
    X3D_sym = sym_xor(sym_xor(A_sym, B_sym), sym_xor(C_sym, D_sym))

    E = X3D << 1
    E_sym = sym_shl(X3D_sym, 1)
    F = X3D << 2
    F_sym = sym_shl(X3D_sym, 2)
    G = X3D << 7
    G_sym = sym_shl(X3D_sym, 7)

    H = X3D ^ E ^ F ^ G
    H_sym = sym_xor(sym_xor(X3D_sym, E_sym), sym_xor(F_sym, G_sym))

    tmp = (x ^ H) & ALL128
    tmp_sym = sym_xor(x_sym, H_sym)[:128]

    print(repr(tmp_sym))
    print("\n".join(sym_str(x) for x in tmp_sym))

    # assert tmp == gf_reduce_64(x)

    return tmp

gf_reduce_instrumented(3213213132132132132)

# further streamlined gf_reduce_128 for python
def gf_reduce(x):
    A = x >> (256 - 1)
    B = x >> (256 - 2)
    C = x >> (256 - 7)
    D = x >> (256 - 128)

    X3D = A ^ B ^ C ^ D

    E = X3D << 1
    F = X3D << 2
    G = X3D << 7
    H = X3D ^ E ^ F ^ G

    tmp = (x ^ H) & ALL128

    # assert tmp == gf_reduce_64(x)

    return tmp

# reduce a less than 256 bit number by the GCM polynomial
# this is just gf_reduce_64 adjusted to operate on 128 bits at a time instead of 64.
# it does *exactly* the same thing.
def gf_reduce_128(x):
    X01 = x & ALL128
    X23 = (x >> 128) & ALL128

    A = X23 >> (63 + 64)
    B = X23 >> (62 + 64)
    C = X23 >> (57 + 64)

    X3D = X23 ^ A ^ B ^ C

    E = (X3D << 1) & ALL128
    F = (X3D << 2) & ALL128
    G = (X3D << 7) & ALL128

    H = X3D ^ E ^ F ^ G

    return X01 ^ H

# algorithm 4 of https://web.archive.org/web/20190806061845/https://software.intel.com/sites/default/files/managed/72/cc/clmul-wp-rev-2.02-2014-04-20.pdf
# also described by https://blog.quarkslab.com/reversing-a-finite-field-multiplication-optimization.html#gk2010
# which credits https://doi.org/10.1016/j.ipl.2010.04.011
# see gf_reduce for a version adapted to use 128 bit ints instead of 64 bit ones.
def gf_reduce_64(x):
    X0 = x & ALL64
    x >>= 64
    X1 = x & ALL64
    x >>= 64
    X2 = x & ALL64
    x >>= 64
    X3 = x & ALL64

    X01 = (X1 << 64) | X0

    A = X3 >> 63
    B = X3 >> 62
    C = X3 >> 57

    D = X2 ^ A ^ B ^ C

    X3D = (X3 << 64) | D

    E = (X3D << 1) & ALL128
    F = (X3D << 2) & ALL128
    G = (X3D << 7) & ALL128

    H = X3D ^ E ^ F ^ G

    return X01 ^ H


# this shows how to perform the modular reduction step in a cpu friendly way
def gf_mul_intrinsic(x, y):

    # Algorithm 2 from # https://web.archive.org/web/20190806061845/https://software.intel.com/sites/default/files/managed/72/cc/clmul-wp-rev-2.02-2014-04-20.pdf

    A0 = x & ALL64
    x >>= 64
    A1 = x & ALL64

    B0 = y & ALL64
    y >>= 64
    B1 = y & ALL64

    if False:
        # one fewer clmul, at the expense of a lot of shuffling
        # these are all 64 x 64 -> 128 carryless multiplications
        C = gf_mul_noreduce(A1, B1)
        D = gf_mul_noreduce(A0, B0)
        E = gf_mul_noreduce(A0 ^ A1, B0 ^ B1)
        E0 = E & ALL64
        E1 = (E >> 64) & ALL64

        tmp = (C << 128) ^ (C << 64) ^ (D << 64) ^ D ^ (E1 << 128) ^ (E0 << 64)
    else:
        # simple implementation with an extra clmul
        # these are all 64 x 64 -> 128 carryless multiplications
        tmp = (
            (gf_mul_noreduce(A1, B1) << 128) ^
            ((gf_mul_noreduce(A1, B0) ^ gf_mul_noreduce(A0, B1)) << 64) ^
            gf_mul_noreduce(A0, B0)
        )

    return gf_reduce(tmp)

gf_mul_count = 0

def gf_mul(x, y):
    global gf_mul_count

    if x == 0 or y == 0:
        return 0

    gf_mul_count += 1

    # tmp = gf_mul_intrinsic(x,y)

    # galois field multiplication aka carryless multiplication.
    # this really is just schoolbook multiplication without carries,
    # though the code is reorganized to be able to efficiently use
    # arithmetic operators and to (almost, except for the overflow
    # check itself) stay within 128 bits.
    result = 0
    while True:
        # note that multiplication in GF2 is the AND operation
        # (the only values are 0 and 1, and only 1*1==1)
        # so the (y&1) is selecting a single bit from y, and then
        # the if is logically multiplying that bit with all the bits
        # in x and then adding that to the result.
        # in case you were wondering why this loop looks linear instead
        # of quadratic complexity, this is the reason: this single line
        # is doing quite a lot of work.
        if y & 1:
            result ^= x

        y >>= 1
        if y == 0:
            break

        x <<= 1
        if x & (1 << 128):
            x ^= GF_POLY

    # assert result == tmp

    return result

def gf_pow(x, e):
    assert e >= 0

    # simple exponentation-by-squaring using gf_mul
    prod = 1
    while True:
        if e & 1:
            prod = gf_mul(prod, x)

        e >>= 1
        if e == 0:
            break

        x = gf_square(x)

    return prod

# for testing only. same as gf_mul but without the reduction step.
def gf_mul_noreduce(x, y):
    result = 0
    while True:
        if y & 1:
            result ^= x

        y >>= 1
        if y == 0:
            break

        x <<= 1

    return result

# used by gf_square
MASK256_64 = 0x00000000_00000000_ffffffff_ffffffff_00000000_00000000_ffffffff_ffffffff
MASK256_32 = 0x00000000_ffffffff_00000000_ffffffff_00000000_ffffffff_00000000_ffffffff
MASK256_16 = 0x0000ffff_0000ffff_0000ffff_0000ffff_0000ffff_0000ffff_0000ffff_0000ffff
MASK256_8 =  0x00ff00ff_00ff00ff_00ff00ff_00ff00ff_00ff00ff_00ff00ff_00ff00ff_00ff00ff
MASK256_4 = 0x0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f
MASK256_2 = 0x33333333_33333333_33333333_33333333_33333333_33333333_33333333_33333333
MASK256_1 = 0x55555555_55555555_55555555_55555555_55555555_55555555_55555555_55555555

# square of a max 64 bit value, requiring no reduction
def gf_square64(x):
    assert x < (1 << 64)
    x = (x | (x << 32)) & MASK128_32
    x = (x | (x << 16)) & MASK128_16
    x = (x | (x << 8)) & MASK128_8
    x = (x | (x << 4)) & MASK128_4
    x = (x | (x << 2)) & MASK128_2
    x = (x | (x << 1)) & MASK128_1

    return x

def gf_square_noreduce(x):
    x = (x | (x << 64)) & MASK256_64
    x = (x | (x << 32)) & MASK256_32
    x = (x | (x << 16)) & MASK256_16
    x = (x | (x << 8)) & MASK256_8
    x = (x | (x << 4)) & MASK256_4
    x = (x | (x << 2)) & MASK256_2
    x = (x | (x << 1)) & MASK256_1
    return x

# equal to gf_mul(x, x), but potentially faster
# note that a gf_mul implementation using a carryless multiply CPU intrinsic
# will definitely beat this. but it should be faster than the naive gf_mul loop.
def gf_square(x):
    orig_x = x

    x = gf_reduce(gf_square_noreduce(x))

    # assert x == gf_mul(orig_x, orig_x)

    return x

GF_SQRT_2 = 0x24924924924924926db6db6db6db6da4
assert GF_SQRT_2 == gf_pow(2, 1 << 127)
assert gf_square(GF_SQRT_2) == 2

# based on https://crypto.stackexchange.com/a/68505
# which cites "Field inversion and point halving revisited" by Fong et al
# at https://cacr.uwaterloo.ca/techreports/2003/corr2003-18.pdf
# and also "Another Look at Square Roots and Traces (and Quadratic Equations) in Fields of Even Characteristic"
# at https://eprint.iacr.org/2007/103.pdf
def gf_sqrt(x):
    # gather even bits into a 64 bit value
    even = x & MASK128_1
    even = (even | (even >> 1)) & MASK128_2
    even = (even | (even >> 2)) & MASK128_4
    even = (even | (even >> 4)) & MASK128_8
    even = (even | (even >> 8)) & MASK128_16
    even = (even | (even >> 16)) & MASK128_32
    even = (even | (even >> 32)) & MASK128_64

    # gather odd bits into a 64 bit value
    odd = (x >> 1) & MASK128_1
    odd = (odd | (odd >> 1)) & MASK128_2
    odd = (odd | (odd >> 2)) & MASK128_4
    odd = (odd | (odd >> 4)) & MASK128_8
    odd = (odd | (odd >> 8)) & MASK128_16
    odd = (odd | (odd >> 16)) & MASK128_32
    odd = (odd | (odd >> 32)) & MASK128_64

    # note that when using CPU intrinsics, the fact that hi is only 64 bits
    # could be used to avoid some of the work.
    result = gf_mul(odd, GF_SQRT_2) ^ even

    # assert result == gf_pow(x, 1 << 127)

    return result

def gf_sqrt_noreduce(x):
    # assert all odd bytes are zero
    assert x < (1 << 128) and x & (MASK128_1 << 1) == 0

    # gather even bits into a 64 bit value
    even = x
    even = (even | (even >> 1)) & MASK128_2
    even = (even | (even >> 2)) & MASK128_4
    even = (even | (even >> 4)) & MASK128_8
    even = (even | (even >> 8)) & MASK128_16
    even = (even | (even >> 16)) & MASK128_32
    even = (even | (even >> 32)) & MASK128_64

    return even

def gf_inverse(x):
    assert x != 0, "zero has no inverse"
    assert x < (1 << 128)

    (d, _, _, inv, _) = gf_extended_euclidean(x, GF_POLY)

    # this is guaranteed if x < GF_POLY since GF_POLY is irreducible
    assert d == 1

    # could've also used fermat's little theorem and gf_pow().
    # easier to understand but a lot more expensive than the above.
    # assert v1 == gf_pow(x, (1<<128) - 2)

    # I think we need to reduce only once max here?
    # it was not necessary under the previous specialized implementation
    # (before we switched to full generic ext eucl.)
    # that version is still in the rust code
    return gf_reduce(inv)

# name conflates two with GF_X but oh well
def gf_inverse_mod_power_of_two(x, e):

    assert x.bit_length() <= e, "x must be smaller than 2**e"
    assert x & 1 == 1, "x must be odd"

    # reduction mod 2**i does not change any of the bits in the remainder,
    # which makes this a simple and satisfying algorithm.

    # loop invariants:
    x_inv = 1 # inverse of x mod (1 << i)
    x_inv_x = x # gf_mul_noreduce(x_inv, x)

    for i in range(1, e):
        # check the loop invariants
        # assert gf_mul_noreduce(x, x_inv) & ((1 << i) - 1) == 1
        # assert x_inv_x == gf_mul_noreduce(x, x_inv)

        if (x_inv_x >> i) & 1:
            # currently digit i of x_inv_x is 1, but we want all digits
            # except digit 0 to be 0. so we flip it by putting a 1 in
            # digit i of of x_inv. see below for more details on how
            # that works; remembering that digit 0 of x is known to be 1.
            x_inv |= 1 << i
            # by adding a digit to x_inv, we are adding another term
            # to the calculation of each digit >=i of x_inv_x.
            # these terms are the digits of x multiplied by the new
            # digit of x_inv (which we just determined to be 1).
            x_inv_x ^= x << i

    # assert gf_mul_noreduce(x, x_inv) & ((1 << e) - 1) == 1

    return x_inv

gf_inverse_mod_power_of_two(1337, (1337).bit_length())
gf_inverse_mod_power_of_two(1337, 20)

def gf_modexp(x, e, m):

    result = 1
    sq = x
    for i in range(e.bit_length()):
        if (e >> i) & 1:
            result = gf_mod(gf_mul_noreduce(result, sq), m)

        sq = gf_mod(gf_square_noreduce(sq), m)
    
    return result

# some small irreducble polynomials over GF2
# https://oeis.org/A014580
GF_IRREDUCIBLES = [2, 3, 7, 11, 13, 19, 25, 31, 37, 41, 47, 55, 59, 61, 67, 73, 87, 91, 97, 103, 109, 115, 117, 131, 137, 143, 145, 157, 167, 171, 185, 191, 193, 203, 211, 213, 229, 239, 241, 247, 253, 283, 285, 299, 301, 313, 319, 333, 351, 355, 357, 361, 369, 375]

# not specific to GF(2**128), works with unreduced inputs
def gf_divmod(x, y):
    q = 0
    r = x

    s = r.bit_length() - y.bit_length()
    while s >= 0:
        r ^= y << s
        q ^= 1 << s
        s = r.bit_length() - y.bit_length()

    return (q, r)

# not specific to GF(2**128), works with unreduced inputs
def gf_mod(x, y):
    return gf_divmod(x, y)[1]

# not specific to GF(2**128), works with unreduced inputs
def gf_div(x, y):
    q, r = gf_divmod(x, y)
    assert r == 0
    return q

# extended euclidean algorithm for integers
# based on https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
# TODO use binary gcd instead, and find a way to use only positive integers
def extended_euclidean(a, b):
    assert a != 0 and b != 0

    (old_r, r) = (a, b)
    (old_s, s) = (1, 0)
    (old_t, t) = (0, 1)

    while r != 0:
        quotient, remainder = divmod(old_r, r)
        # assert remainder == old_r - quotient * r
        (old_r, r) = (r, remainder)
        (old_s, s) = (s, old_s - quotient * s)
        (old_t, t) = (t, old_t - quotient * t)

    # bezout coefficients
    # print(f"{a} * {old_s} + {b} * {old_t} = {a * old_s + b * old_t}, expected {old_r}")
    assert a * old_s + b * old_t == old_r

    # a and b divides by the gcd
    # print(f"{a} // {old_r} = {a // old_r}, expected {abs(t)}")
    # print(f"{b} // {old_r} = {b // old_r}, expected {abs(s)}")
    assert abs(t) == a // old_r and abs(s) == b // old_r

    assert old_r > 0

    return (old_r, abs(t), abs(s), old_s, old_t)

for a in range(1,15):
    for b in range(1,15):
        extended_euclidean(a, b)

# fully-featured extended euclidean algorithm on binary polynomials.
# this is not specific to GF(2**128) and works with unreduced inputs.
# returns tuple (gcd(x, y), x // gcd(x, y), y // gcd(x, y), xc, yc)
# xc and yc are the bezout coefficients such that x * xc + y * yc == gcd(x, y)
# if the gcd is 1, then xc is the inverse of x mod y, and similar for yc
def gf_extended_euclidean(x, y):
    # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
    # https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#B%C3%A9zout's_identity_and_extended_GCD_algorithm
    # https://crypto.stackexchange.com/questions/12956/multiplicative-inverse-in-operatornamegf28/12962#12962
    # https://crypto.stackexchange.com/a/83544

    (u1, u2, u3) = (0, 1, y)
    (v1, v2, v3) = (1, 0, x)

    while v3 != 0:
        # bezout's identity
        assert gf_mul_noreduce(x, u1) ^ gf_mul_noreduce(y, u2) == u3
        assert gf_mul_noreduce(x, v1) ^ gf_mul_noreduce(y, v2) == v3

        (t1, t2, t3) = (u1, u2, u3)
        q = u3.bit_length() - v3.bit_length()
        if q >= 0:
            t1 ^= v1 << q
            t2 ^= v2 << q
            t3 ^= v3 << q
        (u1, u2, u3) = (v1, v2, v3)
        (v1, v2, v3) = (t1, t2, t3)

    # assert gf_mul(gf_mul(gf_div(x, u3), gf_div(y, u3)), gf_square(u3)) == gf_mul(x, y)
    # v2 and v1 are x and y divided by the gcd
    assert v2 == gf_div(x, u3) and v1 == gf_div(y, u3)
    # u1 and u2 are the bezout coefficients
    assert gf_mul_noreduce(x, u1) ^ gf_mul_noreduce(y, u2) == u3
    # if the gcd was 1
    if u3 == 1:
        # then u1 and u2 are the inverses of x mod y and y mod x
        # TODO different way to express this so we don't check for x or y 1
        assert y == 1 or gf_mod(gf_mul_noreduce(x, u1), y) == 1
        assert x == 1 or gf_mod(gf_mul_noreduce(y, u2), x) == 1

    return (u3, v2, v1, u1, u2)

gf_extended_euclidean(gf_mul(GF_IRREDUCIBLES[10], GF_IRREDUCIBLES[20]), gf_mul(GF_IRREDUCIBLES[5], GF_IRREDUCIBLES[10]))
gf_extended_euclidean(gf_mul(GF_IRREDUCIBLES[10], GF_IRREDUCIBLES[20]), gf_mul(GF_IRREDUCIBLES[5], GF_IRREDUCIBLES[15]))

def gf_gcd_split(x, y):
    return tuple(gf_extended_euclidean(x, y)[:3])

def gf_gcd(x, y):
    return gf_extended_euclidean(x, y)[0]

def gf_formal_derivative(x):
    return (x >> 1) & MASK128_1

def mat_mul(a, b):
    assert len(a[0]) == len(b)
    result = []
    for i in range(len(a)):
        row = []
        for j in range(len(b[0])):
            x = 0
            for k in range(len(b)):
                x += a[i][k] * b[k][j]
            row.append(x)
        result.append(row)
    return result

def int_to_bitlist(x, n = None):
    if n is None: n = x.bit_length()
    return [(x >> i) & 1 for i in range(n)]

def bitlist_to_int(x):
    return sum((b & 1) << i for i, b in enumerate(x))

# FIXME untested, unoptimized
def gf_barret_reduction(x, u, u2):
    bits = u.bit_length() - 1

    # u2 should be the inverse of u mod (2**(bits * 2))
    assert gf_mul_noreduce(u, u2) & ((1 << (2*bits)) - 1) == 1

    # the first multiply discards the bottom half of its result,
    # the second one the top half.
    q = gf_mul_noreduce((x >> bits), u2) >> bits
    t = x ^ (gf_mul_noreduce(q, u) & ((1 << bits) - 1))

    return t

def gf_deg(u):
    return u.bit_length() - 1

def gf_factor_canzass(u):

    factors = []

    # cheaply remove powers of x
    power_of_x = count_trailing_zeros(u)
    u >>= power_of_x
    factors.extend([GF_X] * power_of_x)

    # remove repeated factors, making u square free
    (square, u, _) = gf_gcd_split(u, gf_formal_derivative(u))
    if square != 1:
        factors.extend(gf_factor_canzass(gf_sqrt_noreduce(square)) * 2)

    # distinct degree factorization
    todo = {}
    h = GF_X
    deg = 1
    while deg * 2 <= gf_deg(u):
        if deg > 1:
            h = gf_mod(h, u) # in case u became smaller
            h = gf_mod(gf_square_noreduce(h) << 1, u)

        # assert h == gf_modexp(GF_X, (2 ** deg) - 1, u)

        d, _, u = gf_gcd_split(h ^ 1, u)
        if d != 1:
            if gf_deg(d) == deg:
                factors.append(d)
            else:
                todo[deg] = d

        deg += 1

    if u != 1:
        factors.append(u)

    for deg, u in todo.items():
        # TODO same degree factorization
        factors.append(u)

    return factors

# Berlekamp factorization, implemented according to TAOCP vol II 4.6.2 (p. 439)
def gf_factor_berlekamp(u):
    # divide by x as many times as possible, as an optimization.
    # this is extremely cheap compared to letting the regular code handle it.
    x_factors = []
    while u & 1 == 0:
        u >>= 1
        x_factors.append(0b10)

    # make u square free
    (square, u, _) = gf_gcd_split(u, gf_formal_derivative(u))
    square_factors = []
    if square != 1:
        square_factors = gf_factor_berlekamp(gf_sqrt(square))
    
    if u == 1:
        return x_factors + square_factors + square_factors

    # n is the degree of u
    n = u.bit_length() - 1

    # Q is logically an nxn binary matrix.
    # the values in each row are the coefficients of the polynomial
    # x ** (2 * k) mod u, where k is the row number. so we just take
    # the integer representation of that polynomial as the row.
    # note that in the first row, x ** (2 * 0) == x ** 0 == 1
    Q = [1]
    for k in range(1, n):
        Q.append(gf_mod(Q[-1] << 2, u))

    # assert Q == [gf_mod(1 << (2 * k), u) for k in range(n)]

    # index where modular reduction kicks in (divide by 2, rounding up)
    # not used by the algorithm itself. we use this to demonstrate that
    # the first `simple_part` loop iterations are very deterministic
    # and could be skipped (along with the first half of Q) if desired.
    simple_part = n // 2 + (n & 1)
    assert all(Q[k] == 1 << (2 * k) for k in range(simple_part))
    assert all(Q[k] != 1 << (2 * k) for k in range(simple_part, n))

    # subtract identity matrix
    for i in range(n):
        Q[i] ^= 1 << i

    orig_Q = Q[:]

    # bitmask containing bits range(1, simple_part)
    simple_mask = (1 << simple_part) - (1 << 1)

    # original description uses -1 as a sentinel for when c[i] is not set
    # we use a separate c_set bitmap for that.
    # additionally, c itself has key and value reversed to optimize lookups.
    c = [None] * n
    c_set = 0
    # we skip the 1st iteration which sets the first v to 1 and increments r,
    # as suggested in the explanation.
    # we also changed v to be zero-indexed. then, `r` is replaced with len(v).
    v = [1]
    for k in range(1, n):
        Q_k = Q[k]
        # print(f"Q[{k}] = {''.join(map(str, int_to_bitlist(Q_k, n)))}")

        if 1 <= k < simple_part:
            # show that the first simple_part rows are perfectly deterministic
            assert Q_k == (1 << k) | (1 << (2 * k))
        elif k == simple_part:
            # demonstrate we can replicate the effects of the first simple_part
            # iterations on later elements without too much trouble.

            # start with Q[k] without the effects of all previous loops iters
            tmp = orig_Q[k]

            # xor each bit i in range(1, simple_part) with bit i*2**{1..}
            # this will take 7 iterations (or fewer, if lsbs are zero or
            # the polynomial was shorter than 128 bits)
            # note that this is squaring which does not require reduction,
            # due to the mask limiting input bits to the lower half.
            # assuming n is 128, sq is 128 bits and sq & mask is 64 bits.
            # (in fact, if there was reduction, it would mess things up.
            #  we are using squaring just for its nice double-each-bits-
            #  position behavior in GF2)
            sq = gf_square64(tmp & simple_mask)
            for _ in range(6):
                tmp ^= sq
                sq = gf_square64(sq & simple_mask)
                if sq == 0: break
            assert sq == 0

            # show that we produced the correct Q_k
            assert tmp == Q_k

        j = next((j for j in (range(n)) if (Q_k >> j) & 1 != 0 and (c_set >> j) & 1 == 0), None)
        if 1 <= k < simple_part:
            assert j == k

        if j is not None:
            # this has been substantially optimized by lifting the very inner
            # loop to the top and optimizing a lot of binary logic from there.
            # the `i` in the current code is unrelated to the i in the text.

            # bit j is currently 1 in Q_k (we check it when finding j)
            # we don't want to affect it in Q[i], so remove it.
            bit_j = 1 << j
            Q_k_without_bit_j =  Q_k ^ bit_j
            for i in range(k, n):
                if Q[i] & bit_j == 0: continue
                Q[i] ^= Q_k_without_bit_j
            c_set |= bit_j
            c[k] = j
        else:
            # loop iteration k=0 (well actually we hardcode it, but even so)
            # had no bits set in Q_k, and so did not find a j, and so did not
            # set c[0] to anything.
            assert c[0] is None
            # during the simple part, we always pick up j=k because each such
            # row has only two bits set (k and 2*k) and it picks up the first.
            # we could decide to pick the last, which produces a different
            # pattern which seems more difficult to hard-code, so we don't.
            assert c[1:simple_part] == list(range(1, simple_part))

            v_r = (1 << k)
            
            # handle the simple part outside the loop, though it's identical to
            # what the loop would have done
            # technically this should limit to only the bits below k, but it
            # works out (the only time this runs with k<simple_part is k==0,
            # and that Q_k == 0)
            v_r |= (Q_k & simple_mask)

            # handle remaining columns
            for i in range(simple_part, k):
                s = c[i]
                if s is not None:
                    v_r |= ((Q_k >> s) & 1) << i

            v.append(v_r)

        # loop invariant: for every v_i, v_i * Q == 0
        # we can use a generic matmul even though we want binary matmul,
        # because bitlist_to_int discards all but the low bit.
        # for v_i in v:
        #     tmp = mat_mul([int_to_bitlist(v_i, n)], [int_to_bitlist(x, n) for x in Q])
        #     tmp = bitlist_to_int(tmp[0])
        #     assert tmp == 0

    # print("Q:")
    # print('\n'.join(''.join(map(str, int_to_bitlist(r, n))) for r in Q))
    # print("")

    # for i, v_i in enumerate(v):
    #     print(f"v_{i} = {''.join(map(str, int_to_bitlist(v_i, n)))}")

    factors = [u]
    # skip v[0] = 1, which is not useful
    for v_i in v[1:]:
        # ad-hoc optimization: we know u isn't divisible by x so remove those
        # factors from v_i. we could have done this when creating v_i, but that
        # would mess up the nice loop invariant.
        while v_i & 1 == 0:
            v_i >>= 1

        for f in factors[:]:
            # early exit if we've found all the factors
            if len(factors) == len(v): break

            # we are supposed to iterate over v_i - {0,1}
            # however we can skip one of them because:
            # assert u == gf_mul(gf_gcd(v_r ^ 0, u), gf_gcd(v_r ^ 1, u))
            # (see B4 in the cited explanation from TAOCP)
            (d, new_f, new_v_i) = gf_gcd_split(f, v_i)
            if d != 1 and d != f:
                factors.remove(f)
                factors.extend([
                    d,
                    new_f,
                ])
            v_i = new_v_i

    # NOTE if you uncomment, probably remove the early-out in the loop above
    # assert len(factors) == len(v)

    return x_factors + square_factors + square_factors + factors

try:
    from sage.all import GF, Integer, PolynomialRing
    HAS_SAGE = True
    tmp = [Integer((GF_POLY >> i) & 1) for i in range(GF_POLY.bit_length())]
    SAGE_GF = GF(Integer(2)**Integer(128), modulus=tmp, names='b')
    SAGE_POLY = PolynomialRing(SAGE_GF, 'x')

    def gf_to_sage(x): return SAGE_GF.from_integer(x)
    def gf_from_sage(f): return f.to_integer()
    def poly_to_sage(f): return SAGE_POLY([gf_to_sage(x) for x in f])
except ImportError:
    HAS_SAGE = False
    print("No sage support")

for i in range(5):
    x = gf_random()
    
    factors = gf_factor_berlekamp(x)
    factors.sort()
    print(f"factors of {hex(x)}:")
    print(f"berlekamp: {', '.join(hex(f) for f in factors)}")
    tmp = x
    for f in factors:
        q, r = gf_divmod(tmp, f)
        assert r == 0
        tmp = q
    assert q == 1
    factors_cz = gf_factor_canzass(x)
    factors_cz.sort()
    print(f"canzass: {', '.join(hex(f) for f in factors_cz)}")

    if HAS_SAGE:
        tmp = gf_to_sage(x).polynomial().factor()
        assert tmp.unit() == 1
        # polynomials back to field elements
        tmp = [(SAGE_GF(p.list()), e) for p, e in tmp]
        # field elements with multiplicties to flat list of integers
        tmp = [p for p, e in tmp for p in [gf_from_sage(p)] * e]
        tmp.sort()
        print(f"sage: {', '.join(hex(f) for f in tmp)}")
        assert tmp == factors


# factors of the group order (calculated with GNU factor)
GROUP_ORDER = (1 << 128) - 1
GROUP_ORDER_FACTORS = [3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721]
assert GROUP_ORDER == reduce(operator.mul, GROUP_ORDER_FACTORS, 1)

assert 3 * 5 == 0xf
assert 3 * 5 * 17 == 0xff
assert 3 * 5 * 17 * 257 == 0xffff
assert 3 * 5 * 17 * 257 * 65537 == 0xffffffff
assert 3 * 5 * 17 * 257 * 65537 * 641 * 6700417 == 0xffffffffffffffff

# these are the Fermat Numbers, with the first 5 being (the only known) Fermat Primes.
# not to be confused with Mersenne primes (2**e-1), although 3 belongs to both.
assert 3 == 0b11
assert 5 == 0b101
assert 17 == 0b10001
assert 257 == 0b100000001
assert 65537 == 0b10000000000000001
assert 641 * 6700417 == 0b100000000000000000000000000000001
assert 274177 * 67280421310721 == 0b10000000000000000000000000000000000000000000000000000000000000001

# From wikipedia:
# > For any given key and initialization vector value, GCM is limited to encrypting 2**39 − 256 bits
# That translates to a polynomial length of (2**39 - 256) // 8 // 16 + 1 == 0xffffffff
# (Note, this also makes sense, as the 12 bytes of IV + 4 bytes of counter = 16 bytes)
# So that means we only need the "nice" factors to produce roots of unity.

def product(factors):
    return reduce(operator.mul, factors, 1)

# generic euler totient function, though optimized for no duplicate factors
# order_factors is a sorted list of prime factors, such as [2, 2, 5]
def euler_totient(order_factors):
    prod = 1
    prev = None
    for p in order_factors:
        assert prev is None or prev <= p, "list of factors must be sorted"
        if p == prev:
            prod *= p
        else:
            prod *= p - 1
        prev = p
    return prod

# assert that the given element has the given multiplicative order
# order_factors is a sorted list of prime factors of the expected group order, such as [2, 5] for 10
# specialized to the GCM multiplicative group, which means no duplicated factors in the group order are possible
def gf_assert_elem_order(elem, order_factors):
    expected = product(order_factors)

    # establish that the real order at least divides the expected one
    assert gf_pow(elem, expected) == 1

    # assert that all factors are necessary
    prev = None
    for p in order_factors:
        assert prev is None or prev < p, "list of factors must be sorted and not contain duplicates"
        assert gf_pow(elem, expected // p) != 1
        prev = p

# FIXME is this even correct
def gf_find_element_order(elem):
    factors = GROUP_ORDER_FACTORS

    for f in GROUP_ORDER_FACTORS:
        if gf_pow(elem, GROUP_ORDER // f) == 1:
            factors.remove(f)

    return factors

# generator
GF_GEN = 2
gf_assert_elem_order(GF_GEN, GROUP_ORDER_FACTORS)

def num_ordinal(x):
    s = str(x)
    if s.endswith('1'): return s + "st"
    elif s.endswith('2'): return s + "nd"
    elif s.endswith('3'): return s + "rd"
    return s + "th"

FERMAT_PRIMES = [3, 5, 17, 257, 65537]
FERMAT_NUMBERS = [(1 << (1 << b)) | 1 for b in range(7)]
assert all(p in GROUP_ORDER_FACTORS for p in FERMAT_PRIMES)
assert all(p in FERMAT_NUMBERS for p in FERMAT_PRIMES)
assert all(x in FERMAT_NUMBERS for x in [641 * 6700417, 274177 * 67280421310721])

# show multiplicative group structure
prev_order = None
for i in range(1, 128):
    factors = [p for j, p in enumerate(FERMAT_NUMBERS) if (i >> j) & 1]
    if FERMAT_NUMBERS[5] in factors:
        factors.remove(FERMAT_NUMBERS[5])
        factors.extend([641, 6700417])
    if FERMAT_NUMBERS[6] in factors:
        factors.remove(FERMAT_NUMBERS[6])
        factors.extend([274177, 67280421310721])
    factors.sort()
    prod = product(factors)
    tot = euler_totient(factors)
    print(f"cyclic multiplicative subgroup of order {hex(prod)}:")
    if prev_order is not None:
        print(f"    order is {prod / prev_order}x the previous order")
    prev_order = prod
    print(f"    factors of group order: {factors!s}")
    gg = gf_pow(GF_GEN, GROUP_ORDER // prod)
    gf_assert_elem_order(gg, factors)
    print(f"    subgroup has {hex(tot)} generators:")
    print(f"    [gf_pow({hex(gg)}, x) ")
    print(f"        for x in range ({prod})")
    print(f"        if not any(x % p == 0 for p in {factors!r}))]")
    print(f"    each of which has order {hex(prod)}")
    print(f"    which means there are {hex(tot)} primitive {num_ordinal(prod)} roots of unity")

# the gf_square works only because the full group order is a product of these two numbers
print(hex(GROUP_ORDER // 0xffffffffffffffff))
a = gf_pow(GF_GEN, GROUP_ORDER // 0xffffffffffffffff)
print(hex(GROUP_ORDER // 0x10000000000000001))
b = gf_pow(GF_GEN, GROUP_ORDER // 0x10000000000000001)
print(hex(a))
print(hex(gf_mul(b, gf_square(GF_GEN))))
print(hex(b))
print(repr(gf_extended_euclidean(a, b)))

# this has a nasty exponent
print(hex(GROUP_ORDER // 0xffffffff))
a = gf_pow(GF_GEN, GROUP_ORDER // 0xffffffff)
print(hex(GROUP_ORDER // 0x100000001))
b = gf_pow(GF_GEN, GROUP_ORDER // 0x100000001)
print(hex(a))
print(hex(gf_mul(b, gf_pow(GF_GEN, 0x20000000000000002))))
print(hex(b))
print(repr(gf_extended_euclidean(a, b)))

# but it does if we go into a subgroup
gg = gf_pow(GF_GEN, GROUP_ORDER // (0xffffffff * 0x100000001))
print(hex(0xffffffff * 0x100000001 // 0xffffffff))
a = gf_pow(gg, (0xffffffff * 0x100000001) // 0xffffffff)
print(hex(0xffffffff * 0x100000001 // 0x100000001))
print(hex((0xffffffff * 0x100000001) // 0x100000001))
b = gf_pow(gg, (0xffffffff * 0x100000001) // 0x100000001)
print(hex(a))
print(hex(gf_mul(b, gf_pow(gg, 0x2))))
print(hex(b))
print(repr(gf_extended_euclidean(a, b)))
c = gf_pow(gg, 0x100000000)
print(repr(gf_extended_euclidean(a, c)))
print(repr(gf_extended_euclidean(b, c)))

def gfmat_minus_identity(m):
    m = m[:]
    for i in range(len(m)):
        m[i] ^= 1 << i
    return m

def gfmat_kernel(m):

    # most of this was reused from gf_factor_berlekamp, just tidied up

    n = 128
    m = m[:]
    assert len(m) == n

    pivots = []
    pivots_used = 0
    kernel = []
    for i in range(n):
        mi = m[i]

        # print(f"m[{i:3d}] = {''.join(map(str, int_to_bitlist(mi, n)))}")

        # valid pivot columns must be set in mi and not have been used
        available = mi & ~pivots_used

        if available:
            p = count_trailing_zeros(available) # pick the first one.
            for j in range(i, n):
                if m[j] & (1 << p):
                    m[j] ^= mi ^ (1 << p) # xor with mi-without-bit-p
            pivots_used |= 1 << p
            pivots.append(p)
        else:
            tmp = 1 << i
            for j, pj in enumerate(pivots):
                if pj is not None and mi & (1 << pj):
                    tmp |= 1 << j
            kernel.append(tmp)
            pivots.append(None)

        # loop invariant: for every row k in kernel, k * m == 0
        # we can use a generic matmul even though we want binary matmul,
        # because bitlist_to_int discards all but the low bit.
        # for k in kernel:
        #      tmp = mat_mul([int_to_bitlist(k, n)], [int_to_bitlist(x, n) for x in m])
        #      tmp = bitlist_to_int(tmp[0])
        #      assert tmp == 0

    # for j,x in enumerate(m):
    #     print(f"m[{j:3d}] = {''.join(map(str, int_to_bitlist(x, n)))}")

    # for j,x in enumerate(kernel):
    #     print(f"k[{j:3d}] = {''.join(map(str, int_to_bitlist(x, n)))}")

    return kernel

GF_BIT_POWERS = [
    [gf_pow(1 << b, 1 << (1 << e)) for b in range(128)]
    for e in range(8)
]

if HAS_SAGE:
    from sage.all import MatrixSpace, VectorSpace
    SAGE_GFV = VectorSpace(GF(2), 128)
    def gf_to_sage_vector(x): return SAGE_GFV(int_to_bitlist(x, 128))
    def gf_from_sage_vector(v): return bitlist_to_int(map(int, v))
    SAGE_GFM128 = MatrixSpace(GF(2), 128)

    # GF_BIT_POWERS[0] represents a matrix such that multiplying a vector
    # containing the coefficients of a polynomial in the GCM field is equal
    # to squaring and reducing the polynomial directly.
    A = SAGE_GFM128([gf_to_sage_vector(x) for x in GF_BIT_POWERS[0]])

    # print(repr(A.charpoly()))  # x^128 + 1
    # print(repr(A.minpoly().factor()))  # (x + 1)^128

    # verify the the matrix by using it to square a random value
    foo = gf_random()
    assert gf_square(foo) == gf_from_sage_vector(gf_to_sage_vector(foo) * A)

    # 2**32 == 2**2**5
    A32 = A**32
    assert GF_BIT_POWERS[5] == [gf_from_sage_vector(v) for v in A32]

    # GF_KERN = []
    # for i in range(7):
    #     m = SAGE_GFM128([gf_to_sage_vector(x) for x in GF_BIT_POWERS[i]])

    #     # sum with identity matrix and calculate the null space (aka kernel)
    #     sage_kern = (m + SAGE_GFM128.identity_matrix()).kernel()

    #     # get a random element and verify it is invariant under 2**i squarings
    #     foo = gf_from_sage_vector(sage_kern.random_element())
    #     assert foo == gf_pow(foo, 2**2**i)

    #     # use basis() as iterating it directly will give combinations of the basis
    #     GF_KERN.append([gf_from_sage_vector(v) for v in sage_kern.basis()])

# calculate kernel entirely in python
GF_KERN = []
for i in range(7):
    GF_KERN.append(gfmat_kernel(gfmat_minus_identity(GF_BIT_POWERS[i])))

print("GF_KERN = [")
for i in range(len(GF_KERN)):
    print("    [")
    print("        " + '\n        '.join(f"{x:#0130b}," for x in GF_KERN[i]))
    print("    ],")
print("]")

foo = gf_random()
acc = 0
for j,x in enumerate(GF_KERN[5]):
    p = x.bit_length() - 1
    if foo & (1 << p):
        acc ^= x
assert acc == gf_pow(acc, 2**32)


# wouldn't this be nice? but we can't distinguish squares in characteristic 2
# (or rather, everything is a square)
# https://www.ijltemas.in/DigitalLibrary/Vol.5Issue2/06-07.pdf

# this paper suggests x**2 + x as an alternative to quadratic reciprocity in char 2
# note: this operation does not stay within a multiplicative group
# also the property itself is additive, not multiplicative
# https://kconrad.math.uconn.edu/blurbs/ugradnumthy/QRchar2.pdf
# apply x ** 2 + x to g**x:
# g**(2x) + g**x
# (g**x)**2 + g**(sqrt(x))**2
# (g**x)**2 + sqrt(g**x)**2
# (g**x + sqrt(g**x))**2
foo = gf_random() # actually a random exponent < ((1 << 128) - 1)... probably
bar = gf_pow(GF_GEN, foo)
assert gf_square(bar) ^ bar == gf_square(bar ^ gf_sqrt(bar))

# looking at baby-step-giant-step:
# g**x == g**(q*m+r)
# and then setting m = 2**(2**i) so we can use all the cool frobenius identities.
# unfortunately, the multiplications mangle additive structure.

# g**x == (g**q)**m * g**r
# g**x / g**r == (g**q)**m
# (g**x / g**r)**(1/m) == (g**q)
# (g**x)**(1/m) / (g**r)**(1/m) == (g**q)
# (g**x)**(1/m) == (g**q) * (g**r)**(1/m)

# g**x == (g**q)**m * g**r
# g**x / g**r == (g**q)**m
# we can decompose both g**r and g**q into a product-of-squares,
# and we can decompose **m and **(1/m) into a sum-of-products,
# or just see it as a single-term product of squares

# additional option: we can decompose g**r or g**q this way:
# g ** (r[0] * 1 + r[1] * 2 + r[2..4] * 4 + r[4..8] * 8 ...)
# g**(r[0] * 1) * g**(r[1] * 2) * g**(r[2..4] * 4) * g**(r[4..8] * 8) ...
# frob1(g**r[0]) * frob2(g**r[1]) * frob4(g**r[2..4]) * frob8(g**r[4..8]) ...
# or alternatively:
# frob1(g)**r[0] * frob2(g)**r[1] * frob4(g)**r[2..4] * frob8(g)**r[4..8] ...
# remembering that each frobn() represents a sum

p = (1 << 64) | 1;  # not actually prime so just ignore the name
gp = gf_pow(GF_GEN, GROUP_ORDER // p)
m = (1 << 32)
assert gf_pow(gf_pow(gp, m), m) == gf_inverse(gp)
assert gf_pow(gp, (p - 1) // 2) == gf_inverse(gf_sqrt(gp))
assert gf_pow(gp, (p - 1) // 2 + 1) == gf_sqrt(gp)

# do a basic pohlig-hellman discrete logarithm computation as far as we can
# with just the fermat primes, just to see if we can do anything useful with
# that structure.
h = gf_random()

xpp = None
pp = None
hpp = None
gpp = None

for p in FERMAT_PRIMES:
    # the product of all preceding fermat primes is p - 2
    # e.g. for p == 0x10001, pp == 0xffff
    assert pp is None or pp == p - 2

    print(f"prime = {p}")
    print(hex(GROUP_ORDER // (p - 2)))
    assert GROUP_ORDER % (p - 2) == 0
    print(hex(GROUP_ORDER // (p - 1)) + " with remainder " + hex(GROUP_ORDER % (p - 1)))
    print(hex(GROUP_ORDER // p))

    gp = gf_pow(GF_GEN, GROUP_ORDER // p)
    hp = gf_pow(h, GROUP_ORDER // p)

    # find xp such that gf_pow(gp, xp) == hp by dumb brute force
    xp = None
    tmp = 1
    for i in range(0, p):
        if tmp == hp:
            xp = i
            break
        tmp = gf_mul(tmp, gp)

    assert xp is not None
    assert gf_pow(gp, xp) == hp

    if p == 65537:
        # show that hp is part of the vector space we calculated
        acc = 0
        for i,x in enumerate(GF_KERN[5]):
            c = x.bit_length() - 1
            if hp & (1 << c):
                acc ^= x
        assert acc == hp

        # show that hpp is ALSO part of the same vector space
        acc = 0
        for i,x in enumerate(GF_KERN[5]):
            c = x.bit_length() - 1
            if hpp & (1 << c):
                acc ^= x
        assert acc == hpp

        # but it's also in the smaller one, while hp isn't
        # (i.e. pow(hpp, 2**16) == hpp, pow(hp, 2**16) != hp)
        acc = 0
        for i,x in enumerate(GF_KERN[4]):
            c = x.bit_length() - 1
            if hpp & (1 << c):
                acc ^= x
        assert acc == hpp

        # 2**32 % 0xffff == 2**32 % 0x10001 == 1
        assert gf_pow(hpp, 2**32) == hpp
        assert gf_pow(hp, 2**32) == hp
        # 2**16 % 0xffff == 1 but 2**16 % 0x10001 == -1
        assert gf_pow(hpp, 2**16) == hpp
        assert gf_pow(hp, 2**16) == gf_inverse(hp)

    # pretend we found this x in the form of x = (q * m + r) using bsgs
    # both q and r are sort-of in the range 0..m (ignoring the obvious edge case)
    if p != 3:
        assert any(p == 2**2**i + 1 for i in range(8)) # p is of the form 2**2**i + 1
        m = isqrt(p - 1)
        assert any(m == 2**2**i for i in range(8)) # m is of the form 2**2**i
        q, r = divmod(xp, m)

        assert hp == gf_pow(gp, (q * m + r))
        assert hp == gf_mul(gf_pow(gf_pow(gp, q), m), gf_pow(gp, r))

    print(f"new: gf_pow(gf_pow(GF_GEN, GROUP_ORDER // {p}), {xp}) == gf_pow(h, GROUP_ORDER // {p})")

    # because of the group size being of the form 2**2**i + 1, repeated squaring
    # wraps around once, picking up all the inverses, and then cycles.
    # g, g**2, ... g**2**2**i == g**(-1), g**(-2), ... g**2**2**(-i) == g

    initial = hp
    tmp = initial
    elems = []
    for i in range(33):
        elems.append(tmp)
        tmp = gf_square(tmp)
        if tmp == initial:
            print(f"sequence has order {len(elems)}, elems: {repr(elems)}")
            break

    # show where every element's inverse and square root are located
    # (square roots are obvious except the wraparound for i=0)
    for i, e in enumerate(elems):
        assert gf_inverse(e) == elems[(i + len(elems) // 2) % len(elems)]
        assert gf_sqrt(e) == elems[(i - 1) % len(elems)]

    # no element is its own inverse <=> there is no element of order 2
    # so 1 is only in the sequence if it was the start value
    assert elems == [1] or 1 not in elems
    assert elems == [1] or len(elems) == (p.bit_length() - 1) * 2

    # we could take a representative element (say, the maximum) from every such sequence
    # and use that to build an index.

    # failed attempt at another way to construct a cyclic function that returns a subset
    assert gf_pow(gp, p - 1) == gf_inverse(gp)
    foo = (p - 1) // 8
    bar = gf_pow(gp, foo)
    # the problem is we find the inverse of gp here, not the inverse of bar
    assert foo < 2 or gf_pow(bar, 8) == gf_inverse(gp)
    # so instead of being 1, the next item in the sequence is just a lower power of gp
    assert foo < 2 or gf_pow(bar, 9) == gf_pow(gp, foo - 1)

    if pp is None:
        pp = p
        xpp = xp
        hpp = hp
        gpp = gp
    else:
        print(f"prev: gf_pow(gf_pow(GF_GEN, GROUP_ORDER // {pp}), {xpp}) == gf_pow(h, GROUP_ORDER // {pp})")

        assert p - pp == 2
        nexth = gf_pow(h, GROUP_ORDER // (pp * p))
        assert hp == gf_pow(nexth, pp)
        assert hpp == gf_pow(nexth, p)
        assert hpp == gf_mul(hp, gf_square(nexth))
        fakeh = gf_pow(nexth, p - 1)

        nextg = gf_pow(GF_GEN, GROUP_ORDER // (pp * p))
        assert gp == gf_pow(nextg, pp)
        assert gpp == gf_pow(nextg, p)
        assert gpp == gf_mul(gp, gf_square(nextg))
        fakeg = gf_pow(nextg, p - 1)

        assert gf_mul(gf_pow(gp, 10), gf_pow(nextg, 10)) == gf_pow(fakeg, 10)
        assert gf_mul(gf_pow(fakeg, 10), gf_pow(nextg, 10)) == gf_pow(gpp, 10)


        nextx = (xpp * p - xp * pp) * ((pp * p + 1) // 2) % (pp * p)
        assert nextx == (xpp * (pp + 2) - xp * pp) * ((pp * p + 1) // 2) % (pp * p)
        assert gf_pow(nextg, nextx) == nexth

        assert gf_pow(fakeg, nextx) == fakeh

        (d, _, _, a, b) = extended_euclidean(pp, p)
        assert d == 1
        # a is the inverse of pp mod p
        a %= p
        assert (pp * a) % p == 1
        # b is the inverse of p mod pp
        b %= pp
        assert (p * b) % pp == 1
        tmp = (xpp * b * p + xp * a * pp) % (pp * p)
        print(f"xpp * {hex(b)} * p + xp * {hex(a)} * pp")
        assert tmp % pp  == xpp
        assert tmp % p == xp



        # xpp * 0x8 * p + xp * 0x8 * pp
        # 0x8 * (xpp * (t + 1) + xp * (t - 1))
        # 0x8 * (xpp * t + xpp + xp * t - xp)
        # 0x8 * ((xpp + xp) * t + xpp - xp)
        t = p - 1 # this is a power of two
        assert a == b == t // 2 # this is a power of two
        print(repr((t // 2 * ((xpp + xp) * t + xpp - xp))))
        print(repr(pp * p))
        assert tmp == (t // 2 * ((xpp + xp) * t + xpp - xp)) % (pp * p)
        # t//2 * ((xpp + xp) * t + xpp - xp)
        # (xpp + xp) * t**2 // 2 + (xpp - xp) * t // 2
        # ((xpp + xp) * t**2 + (xpp - xp) * t) // 2
        #print(f"ext eucl {t} {pp * p} = " + repr(extended_euclidean(t, pp * p)))
        assert t**2 % (p * pp) == 1 # t is its own inverse mod pp*p
        # ((xpp + xp) * t**2 + (xpp - xp) * t) // 2
        # ((xpp + xp) * 1 + (xpp - xp) * t) // 2
        # ((xpp - xp) * t + xpp + xp) // 2
        #print(f"ext eucl {2} {pp * p} = " + repr(extended_euclidean(2, pp * p)))
        # inverse of 2 is ((pp * p) + 1) // 2 == t**2 // 2
        # this happens to be the t value for the next iteration of the loop
        # this actually applies for all powers of 2: inverse of 2**i is ((pp * p) + 1) // i
        assert tmp == ((xpp - xp) * t + xpp + xp) * (t**2 // 2) % (pp * p)
        # ((xpp - xp) * t + xpp + xp) // 2
        # (xpp*t - xp*t + xpp + xp) // 2
        # (xpp * (t+1) - xp * (t-1)) // 2
        # (xpp * p - xp * pp) // 2
        assert tmp == (xpp * p - xp * pp) * (t**2 // 2) % (pp * p)
        # (xpp * (t+1) - xp * (t-1)) * next_t
        # recurse:
        # (((xpp * (t+1) - xp * (t-1)) * next_t) * (next_t+1) - next_xp * (next_t-1)) // 2
        # (((xpp * (t+1) - xp * (t-1)) * next_t  * (next_t+1)) - next_xp * (next_t-1)) // 2
        # (((xpp * (t+1) - xp * (t-1)) * next_t**2 + (xpp * (t+1) - xp * (t-1)) * next_t)) - next_xp * (next_t-1)) // 2


        xpp = tmp
        pp = pp * p
        hpp = gf_pow(h, GROUP_ORDER // pp)
        gpp = gf_pow(GF_GEN, GROUP_ORDER // pp)

    assert gf_pow(gpp, xpp) == hpp

print(f"gf_pow(gf_pow(GF_GEN, GROUP_ORDER // {pp}), {xpp}) == {hpp}")

# 1 0x80000000 0x80000000
(d, _, _, a, b) = extended_euclidean(0xffffffff, 0xffffffff + 2)
print(f"and then... {d} {a % (0xffffffff + 2)} {b % 0xffffffff}")

# 1 0x8000000000000000 0x8000000000000000
(d, _, _, a, b) = extended_euclidean(0xffffffffffffffff, 0xffffffffffffffff + 2)
print(f"and then... {d} {a % (0xffffffffffffffff + 2)} {b % 0xffffffffffffffff}")



# show how freshman's dream interacts with exponentiation-by-squaring:
# (x+a)**12
# (x+a)**8 * (x+a)**4
# (x**8 + a**8)(x**4 + a**4)
# (x**8 * x**4) + (x**8 * a**4) + (x**4 * a**8) + (a**8 * a**4)
x = gf_random()
a = gf_random()
assert gf_pow(x ^ a, 12) == (
    gf_pow(x, 12)
    ^ gf_mul(gf_pow(x, 8), gf_pow(a, 4))
    ^ gf_mul(gf_pow(x, 4), gf_pow(a, 8))
    ^ gf_pow(a, 12)
)

# separate variable per power of two:
# (x+a)**8 * (x+b)**4
# (x**8 + a**8)(x**4 + b**4)
# (x**8 * x**4) + (x**8 * b**4) + (x**4 * a**8) + (a**8 * b**4)

# for all powers 1,2,4,8 it gets a bit long:
# (x+a)**8 * (x+b)**4 * (x+c)**2 + (x+d)
# (x**8 + a**8)(x**4 + b**4)(x**2 + c**2)(x + d)
#
#   (x**8 * x**4 * x**2 * x)
# + (x**8 * x**4 * x**2 * d)
# + (x**8 * x**4 * c**2 * x)
# + (x**8 * x**4 * c**2 * d)
# + (x**8 * b**4 * x**2 * x)
# + (x**8 * b**4 * x**2 * d)
# + (x**8 * b**4 * c**2 * x)
# + (x**8 * b**4 * c**2 * d)
# + (a**8 * x**4 * x**2 * x)
# + (a**8 * x**4 * x**2 * d)
# + (a**8 * x**4 * c**2 * x)
# + (a**8 * x**4 * c**2 * d)
# + (a**8 * b**4 * x**2 * x)
# + (a**8 * b**4 * x**2 * d)
# + (a**8 * b**4 * c**2 * x)
# + (a**8 * b**4 * c**2 * d)

# reformatted for clarity:
#   (x**15                         )
# + (x**14                      * d)
# + (x**13               * c**2    )
# + (x**12               * c**2 * d)
# + (x**11        * b**4           )
# + (x**10        * b**4        * d)
# + (x**9         * b**4 * c**2    )
# + (x**8         * b**4 * c**2 * d)
# + (x**7  * a**8                  )
# + (x**6  * a**8               * d)
# + (x**5  * a**8        * c**2    )
# + (x**4  * a**8        * c**2 * d)
# + (x**3  * a**8 * b**4           )
# + (x**2  * a**8 * b**4        * d)
# + (x     * a**8 * b**4 * c**2    )
# + (1     * a**8 * b**4 * c**2 * d)



# maybe we can get further by putting even more additive terms in each
# power of two, and then seeing if the sum-of-products form has useful
# structure we can use to rull up a long polynomial into a shorter
# representation? of course we'd then have to still find some way to
# recover roots from that representation...

# TODO write some code to try some of these options

# freshman's dream vs square root
assert gf_square(x) ^ gf_square(a) == gf_square(x ^ a)
assert gf_square(x) ^ a == gf_square(x ^ gf_sqrt(a))

# recalculating the exponent every time is a really bad way to do this,
# but it makes for a nice api. fine as long as it's not used much.
def gf_nth_root(x, e):
    assert all(e % p != 0 for p in GROUP_ORDER_FACTORS)

    e_inv = pow(e, euler_totient(GROUP_ORDER_FACTORS) - 1, GROUP_ORDER)

    return gf_pow(x, e_inv)


# show how to calculate exponentiation by power of two as a sum of bit lookups
# this uses freshman's dream
n = 151057730537251302588469035879534785155
gf_assert_elem_order(n, [65537])

acc = 0
for i in range(128):
    if (n >> i) & 1:
        acc ^= GF_BIT_POWERS[4][i]

assert acc == gf_pow(n, 0x10000)
assert acc == gf_inverse(n)

# flip it around and make a lookup for the value of bit 0 in the output
# one interesting thing to note is that input bit 0 influences NONE of the output
# bits except output bit 0.
# also the lower bits in general tend to have little influence for the smaller groups
bit_lookups = []
for b in range(128):
    acc = 0
    for i in range(128):
        if (GF_BIT_POWERS[3][i] >> b) & 1:
            acc |= 1 << i
    bit_lookups.append(acc)

for i in range(100):
    a = gf_random()
    assert gf_pow(a, 0x100) & 1 == (a & bit_lookups[0]).bit_count() & 1

# assert that the bits in acc are the ONLY ones that influence bit 0 in the output
for i in range(100):
    a = gf_random() & ~bit_lookups[0]
    assert gf_pow(a, 0x100) & 1 == 0

# separate loop so we can generate and print whichever table we want
print("each row shows which input bits are xored to produce that output bit")
m = []
for b in range(128):
    acc = 0
    for i in range(128):
        if (GF_BIT_POWERS[3][i] >> b) & 1:
            acc |= 1 << i
    m.append(acc)
    print(f"{b:3d} {m[b]:0128b}")


if False:

    # use gaussian elimination to find 0x100-th root of bar 
    # yes, this is basically pointless as there are much easier ways

    foo = gf_pow(GF_GEN, GROUP_ORDER // 0x101)
    bar = gf_pow(foo, 0x100)
    assert gf_inverse(foo) == bar

    print(bin(foo))

    # add a result column to complete the equations
    m = [(x << 1) | ((bar >> i) & 1) for i,x in enumerate(bit_lookups)]

    # perform gaussian elimination
    m.sort(reverse=True)
    for i in range(0, 128):
        if m[i].bit_length() < 1: continue
        for j in range(0, 128):
            if i == j: continue
            if m[j] & (1 << (m[i].bit_length() - 1)):
                m[j] ^= m[i]
        m.sort(reverse=True)

    # print matrix
    for i in range(128):
        print(f"{i:3d} {m[i]:0129b}")

    # extract result
    tmp = 0
    for i in range(128):
        # would be 128 - 1 - i as well, but we need to skip the result column
        assert (m[i] >> (128 - i)) & 1
        if m[i] & 1:
            tmp |= 1 << (128 - 1 - i)

    # check it
    assert tmp == foo


# if we didn't already know the order of n, what would this lookup tell us?
# - if the order divides 0xffff, gf_pow(n, 0x10000) == n
# - if the order was 0x10001, gf_pow(n, 0x10000) == gf_inverse(n)

n = 132978334345187347553836768243853666661
print(f"determining order of {n}")
order_factors = gf_find_element_order(n)
print(f"gf_find_element_order: {hex(product(order_factors))} {order_factors}")
n_inv = gf_inverse(n)
x = n
for i in range(129):
    if i > 0:
        order = None
        if x == n:
            order = (1 << i) - 1
        elif x == n_inv:
            order = (1 << i) + 1

        if order is not None:
            factors = [p for p in GROUP_ORDER_FACTORS if order % p == 0]
            remainder = order // product(factors)
            print(f"order of n divides {hex(order)} {factors!r} {hex(remainder) if remainder != 1 else ''}")
        else:
            print(f"order of n does not divide {hex((1 << i) - 1)} or {hex((1 << i) + 1)}")

    x = gf_square(x)

# the exponents needed to calculate the nth-roots for n = 2**2**i
# turn out to also be powers of two
for i in range(0,7):
    e = 1 << (1 << i)
    e_inv = pow(e, euler_totient(GROUP_ORDER_FACTORS) - 1, GROUP_ORDER)
    print(f"e = 2 ** (2 ** {i}) = 2 ** {2 ** i} = {hex(e)}, e_inv = {hex(e_inv)} = 2 ** {e_inv.bit_length() - 1}")
    assert gf_pow(gf_pow(12345, e), e_inv) == 12345



################################################################
# code for operating on polynomials over numbers in GF(2**128) #
################################################################


# the polynomial f(x) = 0
POLY_ZERO = []

# the polynomial f(x) = 1
POLY_ONE = [1]

# the polynomial f(x) = x + 0
POLY_X = [0, 1]

def poly_scalar_mul(f, x):
    if x == 0: return POLY_ZERO
    return [gf_mul(c, x) for c in f]

def poly_trim(f):
    while len(f) > 0 and f[-1] == 0:
        f.pop()
    return f

def poly_eval(f, x):
    result = 0
    xe = 1
    for e, c in enumerate(f):
        # assert xe == gf_pow(x, e)
        result ^= gf_mul(xe, c)
        xe = gf_mul(xe, x)
    return result

def poly_add(f, g):
    if not (len(f) >= len(g)):
        (f,g) = (g,f)

    return [cf ^ cg for cf, cg in zip(f, g)] + f[len(g):]

poly_sub = poly_add

# lc_g_inv means "inverse of the leading coefficient of g".
# if you already have it, passing it in saves some work.
def poly_divmod_simple(f, g, lc_g_inv = None):

    assert g != POLY_ZERO, "cannot divide by zero"

    if f == POLY_ZERO:
        return (POLY_ZERO, POLY_ZERO)

    assert g[-1] != 0, "g is not trimmed"
    assert f[-1] != 0, "f is not trimmed"

    rdigits = len(g) - 1
    qdigits = len(f) - len(g) + 1
    if qdigits <= 0:
        return (POLY_ZERO, f)

    if lc_g_inv is None:
        lc_g_inv = gf_inverse(g[-1])

    q = [0] * qdigits
    r = f[:]

    for i in reversed(range(qdigits)):
        lc = r[i + len(g) - 1]
        # print(f"old: lc {i} = {lc}")
        if lc != 0:
            digit = gf_mul(lc, lc_g_inv)
            q[i] = digit
            r[i + len(g) - 1] = 0
            for j in range(len(g) - 1):
                tmp = gf_mul(digit, g[j])
                # if i + j < rdigits:
                #     print(f"old: r[{i}+{j}] ^= gf_mul(q[{i}], g[{j}])")
                r[i + j] ^= tmp

    # assert r[-qdigits:] == [0] * qdigits

    return (q, poly_trim(r[:-qdigits]))


# poly_divmod_simple but optimized for fewer memory reads and writes
# lc_g_inv means "inverse of the leading coefficient of g".
# if you already have it, passing it in saves some work.
def poly_divmod(f, g, lc_g_inv = None):

    assert g != POLY_ZERO, "cannot divide by zero"

    if f == POLY_ZERO:
        return (POLY_ZERO, POLY_ZERO)

    assert g[-1] != 0, "g is not trimmed"
    assert f[-1] != 0, "f is not trimmed"

    rdigits = len(g) - 1
    qdigits = len(f) - rdigits
    if qdigits <= 0:
        return (POLY_ZERO, f)

    if lc_g_inv is None:
        lc_g_inv = gf_inverse(g[-1])

    qr = f[:]

    for i in reversed(range(0, qdigits + rdigits)):
        tmp = qr[i]
        for j in range(max(0, rdigits - (i + 1)), min(rdigits, qdigits + rdigits - (i + 1))):
            tmp ^= gf_mul(qr[i + 1 + j], g[rdigits - 1 - j])
        if i >= rdigits:
            tmp = gf_mul(tmp, lc_g_inv)
        qr[i] = tmp

    r, q = qr[:rdigits], qr[rdigits:]

    r = poly_trim(r)
    assert len(r) <= rdigits

    # assert (q, r) == poly_divmod_simple(f, g, lc_g_inv)

    return (q, r)

def poly_div(f, g, lc_g_inv = None):
    (q, r) = poly_divmod(f, g, lc_g_inv)
    assert len(r) == 0
    return q

def poly_mod(f, g, lc_g_inv = None):
    (q, r) = poly_divmod(f, g, lc_g_inv)
    return r

def poly_monic(f):
    """divide f by a constant so the leading coefficient becomes 1"""
    if f == POLY_ZERO or f[-1] == 1:
        return f

    return poly_scalar_mul(f, gf_inverse(f[-1]))

def poly_gcd(a, b):
    while b != POLY_ZERO:
        (a, b) = (b, poly_mod(a, b))

    # this seems to be assumed in some algorithms that use gcd,
    # even though it is not strictly in the definition of gcd?
    return poly_monic(a)

def poly_inverse(f, g):
    # basically verbatim from:
    # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Simple_algebraic_field_extensions

    r = g[:]
    newr = f[:]
    t = POLY_ZERO
    newt = POLY_ONE

    while newr != POLY_ZERO:
        (qq, rr) = poly_divmod(r, newr)
        (t, newt) = (newt, poly_sub(t, poly_mul(qq, newt)))
        (r, newr) = (newr, rr)

    assert len(r) == 1

    return poly_scalar_mul(t, gf_inverse(r[0]))

def poly_mul(f, g):
    if len(f) == 0 or len(g) == 0: return []
    result = [0] * (len(f) + len(g) - 1)
    for ef, cf in enumerate(f):
        for eg, cg in enumerate(g):
            result[ef + eg] ^= gf_mul(cf, cg)
    return result

# equivalent to poly_mul(f, g)[:cutoff]
def poly_mul_low(f, g, cutoff, f_is_square=False):
    cutoff = min(cutoff, len(f) + len(g) - 1)
    result = [0] * cutoff
    for ef, cf in enumerate(f):
        if f_is_square and (ef % 2) == 1:
            assert cf == 0
            continue
        for eg, cg in enumerate(g):
            if ef + eg < cutoff:
                result[ef + eg] ^= gf_mul(cf, cg)

    # assert result == poly_mul(f, g)[:cutoff]

    return poly_trim(result)

# equivalent to poly_mul(f, g)[cutoff:]
def poly_mul_high(f, g, cutoff):
    if len(f) + len(g) - 1 <= cutoff:
        return POLY_ZERO[:]
    result = [0] * (len(f) + len(g) - 1 - cutoff)
    for ef, cf in enumerate(f):
        for eg, cg in enumerate(g):
            if ef + eg >= cutoff:
                result[ef + eg - cutoff] ^= gf_mul(cf, cg)

    # assert result == poly_mul(f, g)[cutoff:]

    return poly_trim(result)

# equivalent to poly_mul(f, g)[c]
def poly_mul_coef(f, g, c):

    # we would like to iterate indices 0 <= i <= c in both f and g.
    # however, one or both may be too short to make that possible.
    missing_f = max(0, c + 1 - len(f))
    missing_g = max(0, c + 1 - len(g))

    # note f,g and g,f
    # indices missing from the end of f means we must skip that many from the
    # start of g, and vice versa.
    # NOTE 'start' values may be >= the len, indicating an empty slice.
    #      this will result in the length calculation yielding <= 0.
    f_start = missing_g
    g_start = missing_f

    # normally the last index in both would be c, but some may be missing.
    # we use end in the meaning of slice notation, so one past the last index.
    f_end = c + 1 - missing_f
    g_end = c + 1 - missing_g

    # check that the slice in f is the same length as the slice in g
    len_f = f_end - f_start
    len_g = g_end - g_start
    assert len_f == len_g

    # check that indices sum to c when one of the slices is reversed
    assert f_start + g_end - 1 == c and f_end - 1 + g_start == c

    acc = 0
    for i in range(0, len_f):
        acc ^= gf_mul(f[f_start + i], g[g_end - 1 - i])

    assert acc == poly_mul(f, g)[c] if len_f > 0 else len(poly_mul(f, g)) <= c

    return acc

# calculate first `cutoff` coefficients of f**e
def poly_exp_low(f, e, cutoff):

    if cutoff == 0: return []
    assert cutoff > 0

    f = f[:cutoff]
    prod = None

    while True:
        if e & 1:
            if prod is None:
                prod = f[:cutoff]
            else:
                prod = poly_mul_low(prod, f, cutoff)

        e >>= 1
        if e == 0:
            break

        # TODO thinking about how this progresses, there is some serious
        #      room for optimization here. after a log(n_zeroes) amount of
        #      iterations, all but f[0] will be zero. poly_trim partially
        #      capitalizes on that, but there is more we could do.
        f = poly_trim(poly_square(f)[:cutoff])

    return prod

def poly_square(f):
    if len(f) == 0:
        return []

    # much more efficient than poly_mul(f, f)
    # most of the terms cancel each other out via xor
    result = [0] * (len(f) * 2 - 1)
    for e, c in enumerate(f):
        result[2 * e] = gf_square(c)

    # check it against the unoptimized implementation
    # assert result == poly_mul(f, f)

    return result

def poly_modexp_simple(f, e, g):
    # simple exponentation-by-squaring

    # we need this often, so might as well precalculate it
    lc_g_inv = gf_inverse(g[-1])

    prod = POLY_ONE
    while True:
        if e & 1:
            prod = poly_mod(poly_mul(prod, f), g, lc_g_inv)

        e >>= 1
        if e == 0:
            break

        f = poly_mod(poly_square(f), g, lc_g_inv)

    return prod

poly_modexp = poly_modexp_simple


def poly_modexp_fancy(f, e, g):

    orig_e = e
    orig_f = f[:]

    lc_g_inv = gf_inverse(g[-1])

    table = [None] * (len(g) * 2)
    i = len(g) - 1
    if i & 1 != 0: i += 1
    # all entries could be computed like this but the iterative approach is faster
    table[i] = poly_mod(([0] * i) + [1], g, lc_g_inv)
    # odd entries are unused because those coefficients are always zero (see poly_square)
    i += 2
    while i < len(table):
        table[i] = poly_mod([0, 0] + table[i - 2], g, lc_g_inv)
        i += 2

    prod = POLY_ONE
    while True:
        if e & 1:
            prod = poly_mod(poly_mul(prod, f), g, lc_g_inv)

        e >>= 1
        if e == 0:
            break

        f = poly_square(f)
        tmp = f[:len(g) - 1]
        for i in range(len(g) - 1, len(f)):
            tmp = poly_add(poly_scalar_mul(table[i], f[i]), tmp)
        f = tmp

    assert prod == poly_modexp_simple(orig_f, orig_e, g)

    return prod

test_poly = list(range(3,12))
print(repr(poly_mul(test_poly, test_poly)))
print(repr(poly_mul(test_poly, test_poly[:-1])))
print(repr(poly_mul(test_poly, test_poly[:-2])))
print(repr(poly_mul(test_poly, test_poly[1:])))
print(repr(poly_mul(test_poly, test_poly[2:])))
print(repr(poly_sub(poly_mul(test_poly, test_poly), poly_mul(test_poly, test_poly[1:]))))
print(repr(poly_sub(poly_mul(test_poly, test_poly)[1:], poly_mul(test_poly, test_poly[2:]))))

import re
def human_sort_key(s):
    return [int(t) if i & 1 else t for i, t in enumerate(re.split('([0-9]+)', s))]

def sym_poly_mul(f, g):
    result = [[] for _ in range(len(f) + len(g))]
    for i, ci in enumerate(f):
        for j, cj in enumerate(g):
            prod = tuple(sorted([ci, cj], key=human_sort_key))
            if prod in result[i + j]:
                result[i + j].remove(prod)
            else:
                result[i + j].append(prod)
                result[i + j].sort(key=lambda t: human_sort_key(repr(t)))
    while len(result[-1]) == 0:
        result.pop()
    return result

def sym_poly_str(f):
    tmp = f'total {sum(len(t) for t in f)} terms:\n\t'
    return tmp + '\n\t'.join(f'{i}: ' + ' + '.join(f'{a}*{b}' for a,b in t) for i,t in enumerate(f))

f = [f'f{i}' for i in range(9)]
g = [f'g{i}' for i in range(9)]

print("full:")
print(sym_poly_str(sym_poly_mul(f, g)))

polys = []
for i in range(0, 10, 2):
    print(f"minus {i}:")
    polys.append(sym_poly_mul(g, g[i:]))
    print(sym_poly_str(polys[-1]))

terms = Counter(t for p in polys for x in p for t in x)

print(f"{len(terms)} unique terms")


def mont_reduce(f, g, G, is_square=False):
    m = poly_mul_low(f, G, len(g), f_is_square = is_square)
    # assert poly_mul_low(m, g, len(g)) == poly_trim(f[:len(g)])
    t = poly_add(f[len(g):], poly_mul_high(m, g, len(g)))
    t = poly_trim(t)
    assert len(t) <= len(g)
    return t

def into_mont(f, g, G):
    return poly_mod([0] * len(g) + f, g)

def from_mont(f, g, G):
    return mont_reduce(f, g, G)

def poly_modexp_mont(f, e, g):

    orig_e = e
    orig_f = f[:]
    orig_g = g[:]

    global gf_mul_count
    start_mul_count = gf_mul_count

    # if g is divisible by x (in other words, starts with one or more zeroes)
    # then divide that out for now, we will put it back later via CRT.
    n_zeroes = next(i for i, c in enumerate(g) if c != 0)
    g = g[n_zeroes:]

    R = [0] * len(g) + [1]

    # iteratively build G so that g * G == 1 (i.e. G is the inverse of g mod R)
    # but do it faster than the extended euclidean algorithm by making
    # use of the structure of R by calculating modulo progressively higher
    # powers of f(x) = x.
    # we know G[0] * g[0] == 1 so G[0] is just the inverse of g[0].
    # then, every higher coefficient of gf_mul(g, G) must be zero.
    # adding a new coefficient i to G only adds one new term to the sum
    # constituting gf_mul(g,G)[i], namely gf_mul(G[i], g[0]).
    # this allows us to calculate G[i]:
    # 0 == poly_mul_coef(G[:i], g, i) ^ poly_mul(G[i], g[0])
    # poly_mul_coef(G[:i], g, i) == poly_mul(G[i], g[0])
    # gf_mul(poly_mul_coef(G[:i], g, i), gf_inverse(g[0])) == G[i]
    # gf_mul(poly_mul_coef(G[:i], g, i), G[0]) == G[i]
    G = [gf_inverse(g[0])]
    # if len(g) == 1 then len(G) == 1, so we can skip some work.
    # this is because multiplying two polynomials of length one gives a
    # result of length 1 too, meaning no modular reduction needs to take
    # place at all to find the unique B such that B * g == 1.
    # note that this is purely an optimization, the loop would correctly
    # yield a bunch of zeroes in that case that we then trim off.
    # it's questionable whether this is really worth it as it's quite a
    # niche optimization, but since I went to the trouble of discovering
    # why G sometimes ended in a bunch of zeroes we might as well use it.
    if len(g) > 1:
        while len(G) < len(g):
            G.append(gf_mul(poly_mul_coef(G, g, len(G)), G[0]))
        G = poly_trim(G)

    # assert G == poly_inverse(g, R)
    # assert poly_mod(poly_mul(g, G), R) == POLY_ONE

    fm = into_mont(f, g, G)

    prod = None

    while True:
        if e & 1:
            if prod is None:
                prod = fm
            else:
                prod = mont_reduce(poly_mul(prod, fm), g, G)

        e >>= 1
        if e == 0:
            break

        fm = mont_reduce(poly_square(fm), g, G, is_square=True)

    result = from_mont(prod, g, G)

    if n_zeroes > 0:
        Z = [0] * n_zeroes + [1]

        # perform exponentiation modulo Z
        prod_z = poly_exp_low(orig_f, orig_e, n_zeroes)

        # use chinese remainder theorem to combine result and prod_z

        # compute A, inverse of Z mod g
        A = poly_inverse(poly_mod(Z, g), g)
        result_A = poly_trim([0] * n_zeroes + poly_mul(result,  A))

        print('A: ' + repr(A))

        # assert result_A == poly_trim(poly_mul(result, poly_mul(Z, A)))

        # compute B, inverse of g mod Z
        # basically truncates or adds additional terms to the existing
        # calculation of G as needed, depending on which part of the
        # original g was longer. see comments for G.
        B = G[:n_zeroes]
        if len(g) != 1:
            while len(B) < n_zeroes:
                B.append(gf_mul(poly_mul_coef(B, g, len(B)), B[0]))
            B = poly_trim(B)
        result_B = poly_mul(prod_z, poly_mul(g, B))

        print('B: ' + repr(B))
        # print("n_zeroes = " + str(n_zeroes))
        # assert B == poly_inverse(poly_mod(g, Z), Z)
        # assert poly_mod(poly_mul(B, g), Z) == POLY_ONE

        print('result_A: ' + repr(result_A))
        print('result_B: ' + repr(result_B))
        print("orig_g: " + repr(orig_g))

        result = poly_trim(poly_add(result_A, result_B))

        print(f"len result = {len(result)} len orig_g = {len(orig_g)}")
        result = poly_mod(result, orig_g)

    mont_mul_count = gf_mul_count - start_mul_count

    start_mul_count = gf_mul_count

    tmp = poly_modexp_simple(orig_f, orig_e, orig_g)

    simple_mul_count = gf_mul_count - start_mul_count

    print(f"length {len(orig_g)} simple {simple_mul_count} mont {mont_mul_count}")

    # print('result: ' + repr(tmp))
    # print('correct: ' + repr(result))

    assert result == tmp

    return result

poly_modexp([1, 2], 0xffff, [1, 2, 3, 4])
poly_modexp([3, 4], 0xffff, [0, 5, 6, 7])
poly_modexp([5, 6], 0xffff, [0, 0, 8, 9])
poly_modexp([7, 8], 0xffff, [0, 0, 0, 10])
poly_modexp([7, 8], 0xffff, [0, 0, 0, 1])
poly_modexp([0, 9], 0xffff, [0, 0, 11, 12])

def poly_formal_derivative(f):
    # this is a bit subtle, but the way I understand it, the value to multiply the
    # coefficient by is produced by "multiplying" e with the identity element, but
    # with the usual definition of multiplication (which is "repeated addition"),
    # not the multiplication operation as defined by gf_mul.
    # because we are in characteristic 2, adding the identity element to itself
    # produces 0, so the result of this "multiplication" is e % 2 (== e & 1).
    # at this point we do have a field element so we could use gf_mul to combine
    # it with the existing coefficient, except of course since we are mutiplying
    # by either 0 or 1, we can do better than that.
    # I used the reference given by Wikipedia on the Formal Derivative page:
    # John B. Fraleigh; Victor J. Katz (2002). A First Course in Abstract Algebra. Pearson. p. 443.
    return poly_trim([c if (e & 1) else 0 for e, c in enumerate(f)][1:])

def poly_roots_bta(f):
    if len(f) == 2:
        return [f[0]]
    elif len(f) <= 1:
        return []

    # Compute Tr(ax) mod f
    a = gf_random()
    aXp = [0, a] # 0 + a * x
    Tr = aXp[:]
    for i in range(128-1):
        aXp = poly_mod(poly_square(aXp), f)
        Tr = poly_add(Tr, aXp)

    roots = []
    for c in range(2):
        nf = poly_gcd(f, poly_trim(poly_sub(Tr, [c])))
        if len(nf) >= 2:
            roots.extend(poly_roots_bta(nf))
            f = poly_div(f, nf)
            if len(f) <= 2:
                roots.extend(poly_roots_bta(f))
                break
            Tr = poly_mod(Tr, f)

    return roots

# x**256 divided by GF_POLY is again equal to GF_POLY, though with a remainder.
# the remainder is, as expected, the square of the low part of GF_POLY.
# this is referenced in the paper explaining the gf_reduce optimization.
tmp = [(GF_POLY >> i) & 1 for i in range(GF_POLY.bit_length())]
tmp2 = poly_square(tmp)
q, r = poly_divmod([0] * 256 + [1], tmp)
assert q == tmp and r == [1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

x = 336659828399716795593154996656950128312
x_poly = [(x >> i) & 1 for i in range(x.bit_length())]
# print("foo: " + repr(poly_divmod(tmp2, x_poly)))
# print("bar: " + repr(bin(gf_inverse(x))))

##########################
# AES GCM implementation #
##########################


def split_blocks(data):
    return [data[i:i+16] for i in range(0, len(data), 16)]

def pad_blocksize(data):
    d = len(data) % 16
    if d != 0:
        data += b'\x00' * (16 - d)
    return data

def xor_bytes(x, y):
    assert len(x) == len(y)
    return bytes(i ^ j for i, j in zip(x,y))

def build_ghash_input(ciphertext, auth_data=b''):
    data = b''
    data += pad_blocksize(auth_data)
    data += pad_blocksize(ciphertext)
    data += (8 * len(auth_data)).to_bytes(length=8, byteorder='big')
    data += (8 * len(ciphertext)).to_bytes(length=8, byteorder='big')
    return data

def ghash(auth_key, ghash_input):
    """implements the GHASH function from the AES GCM spec"""

    k = gf_from_bytes(auth_key)

    hash = 0
    for block in split_blocks(ghash_input):
        hash ^= gf_from_bytes(block)
        hash = gf_mul(hash, k)

    print("ghash output: " + hex(hash))

    return gf_to_bytes(hash)

def auth_tag(auth_key, keyblock_0, ciphertext, auth_data=b''):
    """calculate the authentication tag"""

    ghash_input = build_ghash_input(ciphertext, auth_data)

    hash = ghash(auth_key, ghash_input)

    return xor_bytes(hash, keyblock_0)


class AES_GCM:
    """
    AES-GCM implementation intended to be easy to understand and to make
    it easy to do something you should never do: reuse the same IV for
    multiple encryptions.
    """

    def __init__(self, key, iv):
        # basic aes encryption primitive, see _aes_ecb_encrypt
        self._aes_ecb = Cipher(algorithms.AES(key), modes.ECB())

        # the secret used in the ghash function.
        self._auth_key = self._aes_ecb_encrypt(b'\x00' * 16)

        # TODO section 4 mentions a difference here between the NIST spec
        #      and the original GCM spec. look into this.
        #      https://csrc.nist.gov/csrc/media/projects/block-cipher-techniques/documents/bcm/comments/800-38-series-drafts/gcm/joux_comments.pdf

        # if the length of the passed iv is not 12 it is passed through the
        # ghash function and the resulting 128 bit value is *directly* used
        # as the initial value of the counter.
        if len(iv) != 12:
            # *some* minimum length is obviously required for uniqueness.
            # but the exact limits are not consistent between implementations.
            # however most implementations impose limits similar to this.
            assert 8 <= len(iv) <= 128
            hash = ghash(self._auth_key, build_ghash_input(iv))
            init_counter = int.from_bytes(hash, byteorder='big')
        else:
            # if the iv is exactly 12 bytes it is used as the upper 96 bits
            # of the counter. the lower bits are left zero except for the +1.
            # the +1 here is presumably to prevent an all-zero iv from
            # resulting in _auth_key being equal to the auth tag keyblock?
            init_counter = (int.from_bytes(iv, byteorder='big') << 32) + 1

        self._init_counter = init_counter

        print("auth key (decimal): " + str(gf_from_bytes(self._auth_key)))

    def _aes_ecb_encrypt(self, data):
        assert len(data) == 16
        encryptor = self._aes_ecb.encryptor()
        return encryptor.update(data) + encryptor.finalize()

    def _keyblock(self, i):
        counter = self._init_counter + i
        bcounter = counter.to_bytes(length=16, byteorder='big')
        return self._aes_ecb_encrypt(bcounter)

    def _crypt_common(self, data):
        """implements CTR mode encryption and decryption"""

        result = b''
        for i, block in enumerate(split_blocks(data)):
            # +1 because the first keyblock is reserved for the auth tag
            keyblock = self._keyblock(i + 1)
            result += xor_bytes(keyblock[:len(block)], block)

        return result

    def _auth_tag(self, ciphertext, auth_data=b''):
        keyblock_0 = self._keyblock(0)

        return auth_tag(self._auth_key, keyblock_0, ciphertext, auth_data)

    def encrypt(self, plaintext, auth_data=b''):
        ciphertext = self._crypt_common(plaintext)
        assert len(ciphertext) == len(plaintext)

        auth_tag = self._auth_tag(ciphertext, auth_data)

        return ciphertext + auth_tag

    def decrypt(self, ciphertext, auth_data=b''):
        # last 16 bytes of ciphertext are auth tag, split it off
        assert len(ciphertext) >= 16
        ciphertext, ciphertext_auth_tag = ciphertext[:-16], ciphertext[-16:]

        good_auth_tag = self._auth_tag(ciphertext, auth_data)

        if ciphertext_auth_tag != good_auth_tag:
            raise ValueError("invalid tag")

        plaintext = self._crypt_common(ciphertext)

        return plaintext


######################################################################
# code for recovering the auth key from ciphertexts with the same iv #
######################################################################


def ghash_polynomial(data):
    # performs the same operation as ghash(), but since the authentication
    # key is unknown, it returns the result as a polynomial with the
    # authentication key as the variable.

    # rewriting the calculation as done by ghash to standard form:
    # here, '+' is xor, '*' is gf_mul, and '^' means exponentiation.
    # ((((block0 * k) + block1) * k) + block2) * k
    # ((block0 * k) + block1) * k^2 + block2*k
    # block0*k^3 + block1*k^2 + block2*k
    # block0*k^3 + block1*k^2 + block2*k^1 + 0*k^0

    # because everything is multiplied with auth_key at least once,
    # the least significant term is zero.
    poly = [0]

    # we store the polynomial with least significant terms first, because:
    # * it means list indexes are equal to the exponent for that term
    # * we can use len() to determine the order of the polynomial
    for block in reversed(split_blocks(data)):
        poly.append(gf_from_bytes(block))

    return poly

def ciphertext_to_masked_polynomial(ciphertext, auth_data=b''):
    """
    This does not produce a polynomial suitable for solving yet,
    because the tag in the least significant coefficient is still
    xor'ed with a secret masking value.

    But if you xor such a polynomial with another one produced with
    the same iv, this secret masking value will be canceled out and
    the resulting polynomial will equal 0 when evaluated with the
    correct authentication secret.
    """
    ciphertext, tag = ciphertext[:-16], ciphertext[-16:]

    poly = ghash_polynomial(build_ghash_input(ciphertext, auth_data))

    # logically this is xor, but we know poly[0] == 0
    poly[0] = gf_from_bytes(tag)

    return poly

def recover_auth_secret(ciphertexts):

    ciphertexts = set(ciphertexts)

    assert len(ciphertexts) > 1, "need at least two distinct ciphertexts"

    print(f"parsing {len(ciphertexts)} distinct ciphertexts")

    # parse ciphertexts to masked polynomials
    m_polys = [ciphertext_to_masked_polynomial(c) for c in ciphertexts]

    # shorter is better.
    # don't lengthen short polynomials by xoring them with long ones.
    m_polys.sort(key=len)

    print(f"min/max degree: {len(m_polys[0])-1}..{len(m_polys[-1])-1}")

    # unmask each polynomial by xor'ing it with another one
    polys = [poly_sub(f, g) for f, g in zip(m_polys, m_polys[1:])]

    print(f"reducing via gcd")

    # take gcd of all polynomials.
    # makes good use of many ciphertexts, and also helps a lot to reduce
    # the work for long ciphertexts, assuming we have at least three.
    f = poly_monic(polys[0])
    for g in polys[1:]:
        if len(f) <= 2:
            break
        f = poly_gcd(f, g)

    # poly_gcd does this for us
    assert f == poly_monic(f)

    assert f != POLY_ONE, "all ciphertexts should have the same iv"

    print(f"reduced to degree {len(f) - 1}")

    # TODO handle the case where poly is not square-free
    # this does not appear to ever trigger in practice.
    # according to my understanding, we are in a perfect field, and thus this condition is sufficient.
    # https://en.wikipedia.org/wiki/Multiplicity_(mathematics)#Multiplicity_of_a_root_of_a_polynomial
    # "If a {\displaystyle a} is a root of multiplicity k {\displaystyle k} of a polynomial, then it is a root of multiplicity k − 1 {\displaystyle k-1} of the derivative of that polynomial, unless the characteristic of the underlying field is a divisor of k, in which case a {\displaystyle a} is a root of multiplicity at least k {\displaystyle k} of the derivative. "
    c = poly_gcd(f, poly_formal_derivative(f))
    assert c == POLY_ONE, "polynomial is not square-free"

    # from sage.all import GF, Integer, PolynomialRing
    # modulus = [Integer((GF_POLY >> i) & 1) for i in range(GF_POLY.bit_length())]
    # GF128 = GF(Integer(2)**Integer(128), modulus=modulus, names='b')
    # P = PolynomialRing(GF128, 'x')
    # Pf = P([GF128.from_integer(x) for x in f])
    # t = time()
    # tmp = Pf.roots()
    # tmp = [x.to_integer() for x, _ in tmp]
    # print("roots: " + repr(tmp) + f" {time() - t}")
    # from root_find import bta, arm, sra
    # t = time()
    # tmp = bta(Pf)
    # tmp = [x.to_integer() for x in tmp]
    # print("bta: " + repr(tmp) + f" {time() - t}")
    # t = time()
    # tmp = arm(Pf)
    # tmp = [x.to_integer() for x in tmp]
    # print("arm: " + repr(tmp) + f" {time() - t}")
    # t = time()
    # tmp = sra(Pf)
    # tmp = [x.to_integer() for x in tmp]
    # print("sra: " + repr(tmp) + f" {time() - t}")

    t = time();
    tmp = poly_roots_bta(f)
    print("my bta: " + repr(tmp) + f" {time() - t}")

    # A Computational Introduction to Number Theory and Algebra (v2.5)
    # by Victor Shoup
    # https://www.shoup.net/ntb/ntb-v2_5.pdf


    elem = [gf_random() for _ in range(len(f))]
    elem_sq = poly_mod(poly_square(elem), f)

    even = poly_trim([gf_sqrt(c) for i,c in enumerate(elem_sq) if i & 1 == 0])
    odd = poly_trim([gf_sqrt(c) for i,c in enumerate(elem_sq) if i & 1 == 1])

    print("elem " + repr(elem))
    print("elem sq" + repr(elem_sq))
    print("even " + repr(even))
    print("odd " + repr(odd))

    sqrt_poly_x = poly_mod(poly_mul(poly_sub(elem, even), poly_inverse(odd, f)), f)
    assert poly_mod(poly_square(sqrt_poly_x), f) == POLY_X

    print("sqrt x " + repr(sqrt_poly_x))


    t = time()
    print("performing distinct degree factorization")
    # section 20.4.1 distinct degree factorization
    # note that w,p,q are defined at the start of the chapter
    # vastly simplified because we only care about linear factors,
    # which are produced in the first loop iteration.
    h = poly_modexp(POLY_X, 1<<128, f)
    h_minus_x = poly_trim(poly_sub(h, POLY_X))
    f = poly_gcd(h_minus_x, f)

    print(f"reduced to degree {len(f) - 1}")

    print("performing equal-degree factorization")

    # equal degree factorization, specialized for degree 1
    # https://github.com/frereit/frereit.github.io/blob/main/wasm/cantor-zassenhaus/src/factorize.rs
    factors = [f]
    while len(factors) != len(f) - 1:
        # rand = [gf_random() for i in range(len(f) - 1)]
        # suggested by https://arxiv.org/pdf/1012.5322
        rand = [gf_random(), 1]
        g = poly_modexp(rand, (1<<128) // 3, f)

        # (skipped code that does nothing for degree == 1)

        g_plus_one = poly_trim(poly_add(g, POLY_ONE))

        todo = [h for h in factors if len(h) > 2]
        for factor in todo:
            gcd = poly_gcd(factor, g_plus_one)
            if len(gcd) > 1 and len(gcd) < len(factor):
                factors.remove(factor)
                factors.append(poly_div(factor, gcd))
                factors.append(gcd)

    # convert monic linear polynomials into roots
    roots = []
    for x in factors:
        assert len(x) == 2 and x[-1] == 1
        roots.append(x[0])

    print("my roots: "+ repr(roots) + f" {time() - t}")

    # each root is a potential value of the authentication key.
    # for each one, calculate the value used to mask the auth tag.
    results = []
    for root in roots:
        tag_mask = poly_eval(m_polys[0], root)
        results.append((gf_to_bytes(root), gf_to_bytes(tag_mask)))

    return results


############################
# testing and example junk #
############################


if __name__ == '__main__':
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM

    key = AESGCM.generate_key(bit_length=128)
    #iv = b"this iv is a lot longer than 12 bytes"
    #iv = b"short iv"
    #iv = b"A" * 12
    iv = b"\x00" * 12

    orig = b"hello this is a test message"

    gcm = AESGCM(key)
    ciphertext = gcm.encrypt(iv, orig, b"")

    assert orig == gcm.decrypt(iv, ciphertext, b"")

    my_gcm = AES_GCM(key, iv)
    plaintext = my_gcm.decrypt(ciphertext)

    assert plaintext == orig

    if True:
        poly1 = ghash_polynomial(build_ghash_input(ciphertext[:-16]))
        k = gf_from_bytes(my_gcm._auth_key)
        print("should equal ghash output: " + hex(poly_eval(poly1, k)))

    ciphertext2 = gcm.encrypt(iv, b"hi this is another message which is longer but not too long", b"")

    poly = poly_sub(
        ciphertext_to_masked_polynomial(ciphertext),
        ciphertext_to_masked_polynomial(ciphertext2),
    )

    if True:
        k = gf_from_bytes(my_gcm._auth_key)
        print("should equal 0: " + hex(poly_eval(poly, k)))


    ciphertext3 = gcm.encrypt(iv, b"hey it's a third ciphertext, neat", b"")

    poly2 = poly_sub(
        ciphertext_to_masked_polynomial(ciphertext),
        ciphertext_to_masked_polynomial(ciphertext3),
    )

    poly = poly_gcd(poly, poly2)

    from random import randbytes

    recovered = recover_auth_secret([
        gcm.encrypt(iv, randbytes(316), b""),
        gcm.encrypt(iv, randbytes(316), b""),
        # gcm.encrypt(iv, randbytes(100), b""),
    ])

    if True:
        k = my_gcm._auth_key
        tag_mask = my_gcm._keyblock(0)
        assert (k, tag_mask) in recovered

    from binascii import hexlify

    for auth_key, keyblock_0 in recovered:
        ciphertext = b"fake ciphertext"
        ciphertext += auth_tag(auth_key, keyblock_0, ciphertext)
        try:
            my_gcm.decrypt(ciphertext)
            print(f"auth_key == {hexlify(auth_key)}")
            print(f"keyblock_0 == {hexlify(keyblock_0)}")
        except ValueError:
            print(f"auth_key != {hexlify(auth_key)}")
            print(f"keyblock_0 != {hexlify(keyblock_0)}")
