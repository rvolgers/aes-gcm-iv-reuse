#!/usr/bin/env python3

# only used for a basic AES-ECB primitive
from collections import Counter
from copy import deepcopy
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

from random import getrandbits
from time import time
from functools import reduce
import operator

##############################################
# code for operating on numbers in GF(2^128) #
##############################################


# GF(2^128) with polynomial x^128 + x^7 + x^2 + x + 1
GF_POLY = (1 << 128) | (1 << 7) | (1 << 2) | (1 << 1) | 1

# used by _gf_bitswap.
MASK4 = 0x0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f  # 0b00001111 x 16
MASK2 = 0x33333333_33333333_33333333_33333333  # 0b00110011 x 16
MASK1 = 0x55555555_55555555_55555555_55555555  # 0b01010101 x 16

# used by gf_from_bytes / gf_to_bytes.
# reverses the order of bits within each byte of a 128 bit value.
# the spec mandates this unusual bit order for some reason.
# it's possible to avoid this swapping by painstakingly adjusting all
# operations to account for this, but that's annoying.
def _gf_bitswap(x):
    x = ((x & MASK4) << 4) | ((x >> 4) & MASK4)
    x = ((x & MASK2) << 2) | ((x >> 2) & MASK2)
    x = ((x & MASK1) << 1) | ((x >> 1) & MASK1)
    return x

def gf_from_bytes(b):
    assert len(b) == 16, "argument must be 16 bytes"
    return _gf_bitswap(int.from_bytes(b, byteorder='little'))

def gf_to_bytes(x):
    assert x < (1 << 128), "element not properly reduced"
    return _gf_bitswap(x).to_bytes(length=16, byteorder='little')

def gf_random():
    return gf_from_bytes(randbytes(16))

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

# https://web.archive.org/web/20190806061845/https://software.intel.com/sites/default/files/managed/72/cc/clmul-wp-rev-2.02-2014-04-20.pdf
# this shows how to perform the modular reduction step in a cpu friendly way
def gf_mul_intrinsic(x, y):

    # Algorithm 2

    A0 = x & ALL64
    A1 = (x >> 64) & ALL64

    B0 = y & ALL64
    B1 = (y >> 64) & ALL64

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

    # Algorithm 4

    if False:
        X0 = tmp & ALL64
        tmp >>= 64
        X1 = tmp & ALL64
        tmp >>= 64
        X2 = tmp & ALL64
        tmp >>= 64
        X3 = tmp & ALL64

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
    else:
        X01 = tmp & ALL128
        X23 = (tmp >> 128) & ALL128

        A = X23 >> (63 + 64)
        B = X23 >> (62 + 64)
        C = X23 >> (57 + 64)

        X3D = X23 ^ A ^ B ^ C

        E = (X3D << 1) & ALL128
        F = (X3D << 2) & ALL128
        G = (X3D << 7) & ALL128

        H = X3D ^ E ^ F ^ G

        return X01 ^ H

def gf_mul(x, y):
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
    # simple exponentation-by-squaring using gf_mul
    prod = 1
    while True:
        if e & 1:
            prod = gf_mul(prod, x)

        e >>= 1
        if e == 0:
            break

        x = gf_mul(x, x)

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

def gf_inverse(x):
    assert x != 0, "zero has no inverse"

    # inverse by extended euclidean algorithm.
    # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
    # https://crypto.stackexchange.com/questions/12956/multiplicative-inverse-in-operatornamegf28/12962#12962
    # https://crypto.stackexchange.com/a/83544
    (u1, u2, u3) = (0, 1, GF_POLY)
    (v1, v2, v3) = (1, 0, x)

    # notice that the GF_POLY term of the invariant is zero modulo GF_POLY.
    # so we exit once we have v1 such that gf_mul(x, v1) == 1,
    # which means v1 is the inverse of x modulo GF_POLY.
    # note that v3 is actually the gcd of GF_POLY and x.
    # we know it will always be 1, because GF_POLY is irreducible.
    while v3 != 1:
        # loop invariant
        # assert gf_mul_noreduce(x, u1) ^ gf_mul_noreduce(GF_POLY, u2) == u3
        # assert gf_mul_noreduce(x, v1) ^ gf_mul_noreduce(GF_POLY, v2) == v3

        (t1, t2, t3) = (u1, u2, u3)
        q = u3.bit_length() - v3.bit_length()
        if q >= 0:
            t1 ^= v1 << q
            t2 ^= v2 << q
            t3 ^= v3 << q
        (u1, u2, u3) = (v1, v2, v3)
        (v1, v2, v3) = (t1, t2, t3)

        # after the first loop iteration we should be within 128 bits again
        # assert u1 < (1<<128) and u2 < (1<<128) and u3 < (1<<128)
        # assert v1 < (1<<128) and v2 < (1<<128) and v3 < (1<<128)

        # so now we can express the loop invariant modulo GF_POLY
        # assert gf_mul(x, u1) == u3
        # assert gf_mul(x, v1) == v3

    # could've also used fermat's little theorem and gf_pow().
    # easier to understand but a lot more expensive than the above.
    # assert v1 == gf_pow(x, (1<<128) - 2)

    return v1


# factors of the group order
GROUP_ORDER = (1 << 128) - 1
GROUP_ORDER_FACTORS = [3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721]
assert GROUP_ORDER == reduce(operator.mul, GROUP_ORDER_FACTORS, 1)

# generator
GF_GEN = 2
assert not any(gf_pow(GF_GEN, x) == (1 << 128) - 1 for x in GROUP_ORDER_FACTORS)

# non-trivial roots of unity (note that 1 is also technically a root of unity)
ROOTS_OF_UNITY = [gf_pow(GF_GEN, GROUP_ORDER // x) for x in GROUP_ORDER_FACTORS]
assert all(gf_pow(x, e) == 1 for x, e in zip(ROOTS_OF_UNITY, GROUP_ORDER_FACTORS))


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

if False:
    import re
    def human_sort_key(s):
        return [int(t) if i & 1 else t for i, t in enumerate(re.split('([0-9]+)', s))]

    # all values are in sum-of-products form
    def prod_to_str(p):
        # return ' * '.join(f if e == 1 else f"{f}**{e}" for f,e in sorted(Counter(p).items()))
        return ' * '.join(f for f in p)
    def sym_to_str(s):
        terms = [prod_to_str(p) for p in s if not ('0' in p)]
        terms.sort(key = human_sort_key)
        # terms = [f"{t} * {m}" if m > 1 else t for t, m in sorted(Counter(terms).items())]
        return ' + '.join(terms)
    def sym_mul_const(val, const):
        return [p + [const] for p in val]
    f = [[[f"0"]] for i in range(14)]
    f[-1] = [["lc(f)"]]
    g = [[[f"g{i}"]] for i in range(4)]
    g[-1] = [["lc(g)"]]
    qdigits = len(f) - len(g) + 1
    r = f[:]
    n = 0
    digits = []
    for i in reversed(range(len(f))):
        lc = r[i]
        digit_without_lc_g_inv = lc
        digit = sym_mul_const(lc, f'inv(g[{len(g) - 1}])')
        digits.append(digit)
        # if n < 3:
        #     digit = [[f'digit_{n}']]
        #     digit_without_lc_g_inv = [[f'digit_{n} * g[{len(g) - 1}]']]
        # preserve r[i] at the moment they would otherwise have been zeroed
        #r[i] = [['0']]
        for j in reversed(range(len(g) - 1)):
            if i - len(g) + 1 + j < 0: break
            #r[i - len(g) + 1 + j] = deepcopy(r[i - len(g) + 1 + j]) + sym_mul_const(deepcopy(digit), f'g[{j}]')
            # use lc instead of digit and gg instead of g. this moves the lc_g_inv factor from digit to gg.
            r[i - len(g) + 1 + j] = deepcopy(r[i - len(g) + 1 + j]) + sym_mul_const(deepcopy(digit_without_lc_g_inv), f'gg[{j}]')
        n += 1

    # for i, t in enumerate(digits):
    #     if i > 6: break
    #     prefix = f"digit_{i} = ("
    #     indent = " " * len(prefix)
    #     print(prefix + sym_to_str(t).replace(' + ', ' +\n' + indent) + ')')

    for i, t in enumerate(r):
        prefix = f"r_{i} = ("
        indent = " " * len(prefix)
        print(prefix + sym_to_str(t).replace(' + ', ' +\n' + indent) + ')')

if False:
    import re
    def human_sort_key(s):
        return [int(t) if i & 1 else t for i, t in enumerate(re.split('([0-9]+)', s))]

    # all values are in sum-of-products form
    def prod_to_str(p):
        # return ' * '.join(f if e == 1 else f"{f}**{e}" for f,e in sorted(Counter(p).items()))
        return ' * '.join(f for f in p)
    def sym_to_str(s):
        terms = [prod_to_str(p) for p in s if not ('0' in p)]
        terms.sort(key = human_sort_key)
        # terms = [f"{t} * {m}" if m > 1 else t for t, m in sorted(Counter(terms).items())]
        return ' + '.join(terms)
    def sym_mul_const(val, const):
        return [p + [const] for p in val]
    f = [[[f"0"]] for i in range(14)]
    f[-1] = [["lc(f)"]]
    g = [[[f"g{i}"]] for i in range(4)]
    g[-1] = [["lc(g)"]]
    qdigits = len(f) - len(g) + 1
    r = f[:]
    n = 0
    digits = []
    for i in reversed(range(len(f))):
        lc = r[i]
        digit_without_lc_g_inv = lc
        digit = sym_mul_const(lc, f'inv(g[{len(g) - 1}])')
        digits.append(digit)
        if n < 3:
            digit = [[f'digit_{n}']]
            digit_without_lc_g_inv = [[f'digit_{n} * g[{len(g) - 1}]']]
        # preserve r[i] at the moment they would otherwise have been zeroed
        #r[i] = [['0']]
        for j in reversed(range(len(g))):
            if i - len(g) + 1 + j < 0: break
            #r[i - len(g) + 1 + j] = deepcopy(r[i - len(g) + 1 + j]) + sym_mul_const(deepcopy(digit), f'g[{j}]')
            # use lc instead of digit and gg instead of g. this moves the lc_g_inv factor from digit to gg.
            r[i - len(g) + 1 + j] = deepcopy(r[i - len(g) + 1 + j]) + sym_mul_const(deepcopy(digit_without_lc_g_inv), f'gg[{j}]')
        n += 1

    for i, t in enumerate(digits):
        if i > 6: break
        prefix = f"digit_{i} = ("
        indent = " " * len(prefix)
        print(prefix + sym_to_str(t).replace(' + ', ' +\n' + indent) + ')')

    for i, t in enumerate(r):
        prefix = f"r_{i} = ("
        indent = " " * len(prefix)
        print(prefix + sym_to_str(t).replace(' + ', ' +\n' + indent) + ')')
    '''
    # to make ide shut up
    lc_g_inv = None
    gg = []

    # gg[x] = g[x] * lc_g_inv
    digit_0 = (lc(f) * lc_g_inv)
    digit_1 = (digit_0 * gg[2])
    digit_2 = (digit_0 * gg[1] +
               digit_1 * gg[2])
    digit_3 = (digit_0 * gg[0] +
               digit_1 * gg[1] +
               digit_2 * gg[2])
    digit_4 = (digit_1 * gg[0] +
               digit_2 * gg[1] +
               digit_3 * gg[2])
    
    digit_3 = (digit_0 * gg[0] +
               (digit_0 * gg[2]) * gg[1] +
               (digit_0 * gg[1] + (digit_0 * gg[2]) * gg[2]) * gg[2])
    digit_3 = (digit_0 * gg[0] +
               digit_0 * gg[2] * gg[1] +
               digit_0 * gg[1] * gg[2] + digit_0 * gg[2] * gg[2] * gg[2])
    digit_3 = (gg[0] +
               gg[2] * gg[1] +
               gg[1] * gg[2] +
               gg[2] * gg[2] * gg[2]) * digit_0
    digit_3 = (gg[0] + gg[2]*gg[2]*gg[2]) * digit_0


    # recursively substitute to depend only on digit_0
    digit_4 = ((digit_0 * gg[2]) * gg[0] +
               (digit_0 * gg[1] + (digit_0 * gg[2]) * gg[2]) * gg[1] +
               (digit_0 * gg[0] + (digit_0 * gg[2]) * gg[1] + (digit_0 * gg[1] + (digit_0 * gg[2]) * gg[2]) * gg[2]) * gg[2])
    # flatten a bit
    digit_4 = ((digit_0 * gg[2]) * gg[0] +
               (digit_0 * gg[1] + digit_0 * gg[2] * gg[2]) * gg[1] +
               (digit_0 * gg[0] + digit_0 * gg[2] * gg[1] + digit_0 * gg[1] * gg[2] + digit_0 * gg[2] * gg[2] * gg[2]) * gg[2])
    # pull out digit_0
    digit_4 = (gg[2] * gg[0] +
               (gg[1] + gg[2] * gg[2]) * gg[1] +
               (gg[0] + gg[2] * gg[1] + gg[1] * gg[2] + gg[2] * gg[2] * gg[2]) * gg[2]) * digit_0
    # gg[2]*gg[1] + gg[1]*gg[2] == 0
    digit_4 = (gg[2]*gg[0] +
               (gg[1] + gg[2]*gg[2]) * gg[1] +
               (gg[0] + gg[2]*gg[2]*gg[2]) * gg[2]) * digit_0
    # flatten
    digit_4 = (gg[0]*gg[2] +
               gg[1]*gg[1] + gg[1]*gg[2]*gg[2] +
               gg[0]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]) * digit_0
    # gg[0]*gg[2] + gg[0]*gg[2] == 0
    digit_4 = (gg[1]*gg[1] + gg[1]*gg[2]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]) * digit_0

    # next!
    digit_5 = ((digit_0 * gg[1] + (digit_0 * gg[2]) * gg[2]) * gg[0] +
               ((gg[0] + gg[2]*gg[2]*gg[2]) * digit_0) * gg[1] +
               ((gg[1]*gg[1] + gg[1]*gg[2]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]) * digit_0) * gg[2])
    digit_5 = (gg[1]*gg[0] + gg[2]*gg[2]*gg[0] +
               gg[0]*gg[1] + gg[2]*gg[2]*gg[2]*gg[1] +
               gg[1]*gg[1]*gg[2] + gg[1]*gg[2]*gg[2]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0
    digit_5 = (gg[2]*gg[2]*gg[0] +
               gg[2]*gg[2]*gg[2]*gg[1] +
               gg[1]*gg[1]*gg[2] + gg[1]*gg[2]*gg[2]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0
    digit_5 = (gg[2]*gg[2]*gg[0] + gg[1]*gg[1]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0

    # next!
    digit_6 = (((gg[0] + gg[2]*gg[2]*gg[2])) * gg[0] +
           ((gg[1]*gg[1] + gg[1]*gg[2]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2])) * gg[1] +
           ((gg[2]*gg[2]*gg[0] + gg[1]*gg[1]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2])) * gg[2]) * digit_0
    digit_6 = (gg[0]*gg[0] + gg[2]*gg[2]*gg[2]*gg[0] +
           gg[1]*gg[1]*gg[1] + gg[1]*gg[2]*gg[2]*gg[1] + gg[2]*gg[2]*gg[2]*gg[2]*gg[1] +
           gg[2]*gg[2]*gg[2]*gg[0] + gg[1]*gg[1]*gg[2]*gg[2] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0
    digit_6 = (gg[0]*gg[0] + gg[2]*gg[2]*gg[2]*gg[0] +
           gg[1]*gg[1]*gg[1] + gg[2]*gg[2]*gg[2]*gg[2]*gg[1] +
           gg[2]*gg[2]*gg[2]*gg[0] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0
    digit_6 = (gg[0]*gg[0] + 
           gg[1]*gg[1]*gg[1] + gg[2]*gg[2]*gg[2]*gg[2]*gg[1] +
           gg[2]*gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0
    digit_6 = (gg[0]*gg[0] + gg[1]*gg[1]*gg[1] + gg[2]*gg[2]*gg[2]*gg[2]*gg[1] + gg[2]*gg[2]*gg[2]*gg[2]*gg[2]*gg[2]) * digit_0


    # new idea, instead of producing a single new value, we produce three new values from three inputs
    digit_3 = (digit_0 * gg[0] +
               digit_1 * gg[1] +
               digit_2 * gg[2])
    # becomes:
    f(digit_0, digit_1, digit_2) = (digit_1, digit_2, (digit_0 * gg[0] + digit_1 * gg[1] + digit_2 * gg[2]))
    # written as a matrix:
    [
        [0, 0, gg[0]],
        [1, 0, gg[1]],
        [0, 1, gg[2]],
    ]
    # hm that doesn't seem useful, forget it

    # initial 3 values
    digit_0 = 1 * digit_0
    digit_1 = gg[2] * digit_0
    digit_2 = (gg[1] + (gg[2] * gg[2])) * digit_0

    # next 3 values as a function of the first three
    digit_3 = (digit_0 * gg[0] +
               digit_1 * gg[1] +
               digit_2 * gg[2])

    digit_4 = (digit_1 * gg[0] +
               digit_2 * gg[1] +
               (digit_0 * gg[0] + digit_1 * gg[1] + digit_2 * gg[2]) * gg[2])
    digit_4 = (digit_1 * gg[0] +
               digit_2 * gg[1] +
               digit_0 * gg[0] * gg[2] +
               digit_1 * gg[1] * gg[2] +
               digit_2 * gg[2] * gg[2])
    digit_4 = (digit_0 * gg[0] * gg[2] +
               digit_1 * (gg[0] + gg[1] * gg[2]) +
               digit_2 * (gg[1] + gg[2] * gg[2]) +

    digit_5 = (digit_2 * gg[0] +
               (digit_0 * gg[0] + digit_1 * gg[1] + digit_2 * gg[2]) * gg[1] +
               (digit_1 * gg[0] + digit_2 * gg[1] + digit_0 * gg[0] * gg[2] + digit_1 * gg[1] * gg[2] + digit_2 * gg[2] * gg[2]) * gg[2])
    digit_5 = (digit_2 * gg[0] +
               digit_0 * gg[0] * gg[1] +
               digit_1 * gg[1] * gg[1] +
               digit_2 * gg[2] * gg[1] +
               digit_1 * gg[0] * gg[2] +
               digit_2 * gg[1] * gg[2] +
               digit_0 * gg[0] * gg[2] * gg[2] +
               digit_1 * gg[1] * gg[2] * gg[2] +
               digit_2 * gg[2] * gg[2] * gg[2])
    digit_5 = (digit_0 * gg[0] * gg[1] +
               digit_0 * gg[0] * gg[2] * gg[2] +
               digit_1 * gg[1] * gg[1] +
               digit_1 * gg[0] * gg[2] +
               digit_1 * gg[1] * gg[2] * gg[2] +
               digit_2 * gg[2] * gg[2] * gg[2]
               digit_2 * gg[1] * gg[2] +
               digit_2 * gg[0] +
               digit_2 * gg[2] * gg[1])
    digit_5 = (digit_0 * (gg[0] * gg[1] + gg[0] * gg[2] * gg[2]) +
               digit_1 * (gg[1] * gg[1] + gg[0] * gg[2] + gg[1] * gg[2] * gg[2]) +
               digit_2 * (gg[2] * gg[2] * gg[2] + gg[1] * gg[2] + gg[0] + gg[2] * gg[1]))

    # okay so this might work, but it's fundamentally quadratic so kinda expensive.
    # but lets see how it goes, maybe we can optimize it somehow...

    

    '''

    def mat_mul(a, b):
        assert len(a[0]) == len(b)
        result = []
        for i in range(len(a)):
            row = []
            for j in range(len(b[0])):
                x = 0
                for k in range(len(b)):
                    x ^= gf_mul(a[i][k], b[k][j])
                row.append(x)
            result.append(row)
        return result

    def sym_mat_mul(a, b):
        assert len(a[0]) == len(b)
        result = []
        for i in range(len(a)):
            row = []
            for j in range(len(b[0])):
                x = '0'
                for k in range(len(b)):
                    if a[i][k] != '0' and b[k][j] != '0': 
                        if x == '0':
                            x = ''
                        else:
                            x += ' + '
                        x += f"{a[i][k]} * {b[k][j]}"
                row.append(x)
            result.append(row)
        return result

    L = [[f'l{i}{j}' if i >= j else '0' for j in range(3)] for i in range(3)]
    U = [[f'u{i}{j}' if i <= j else '0' for j in range(3)] for i in range(3)]

    print(repr(sym_mat_mul(L, U)))

    # [
    #     ['l00 * u00', 'l00 * u01',             'l00 * u02'],
    #     ['l10 * u00', 'l10 * u01 + l11 * u11', 'l10 * u02 + l11 * u12'],
    #     ['l20 * u00', 'l20 * u01 + l21 * u11', 'l20 * u02 + l21 * u12 + l22 * u22']
    # ]
    # comparing this with sym_m:
    # by the first column, u00 can't be zero, so l00 must be zero.
    # but by the first row, l00 can't be zero.
    
    sym_m = [
        ['0', '0', 'gg[0]'],
        ['1', '0', 'gg[1]'],
        ['0', '1', 'gg[2]'],
    ]

    assert(sym_mat_mul([['digit_0', 'digit_1', 'digit_2']], sym_m) == [[
        'digit_1 * 1',
        'digit_2 * 1',
        'digit_0 * gg[0] + digit_1 * gg[1] + digit_2 * gg[2]',
    ]])




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
    result = [0] * (len(f) + len(g) - 1)
    for ef, cf in enumerate(f):
        for eg, cg in enumerate(g):
            result[ef + eg] ^= gf_mul(cf, cg)
    return result

def poly_square(f):
    # much more efficient than poly_mul(f, f)
    # most of the terms cancel each other out via xor
    result = [0] * (len(f) * 2 - 1)
    for e, c in enumerate(f):
        result[2 * e] = gf_mul(c, c)

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


def poly_modexp(f, e, g):

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

def mont_reduce(f, g, G):
    # to see speed wins, the multiplications should be optimized to only
    # compute the terms that are actually used. even then, there doesn't
    # seem to be an obvious advantage.
    m = poly_mul(f[:len(g)], G)[:len(g)]
    t = poly_add(f, poly_mul(m, g))[len(g):]
    t = poly_trim(t)
    assert len(t) <= len(g)
    return poly_mod(t, g)

def into_mont(f, g, G):
    return poly_mod([0] * len(g) + f, g)

def from_mont(f, g, G):
    return mont_reduce(f, g, G)

def poly_modexp_montgomery(f, e, g):

    orig_e = e
    orig_f = f[:]

    R = [0] * len(g) + [1]
    G = poly_inverse(g, R)

    # assert poly_mod(poly_mul(g, G), R) == POLY_ONE

    fm = into_mont(f, g, G)

    prod = into_mont(POLY_ONE, g, G)

    while True:
        if e & 1:
            prod = mont_reduce(poly_mul(prod, fm), g, G)

        e >>= 1
        if e == 0:
            break

        fm = mont_reduce(poly_square(fm), g, G)

    result = from_mont(prod, g, G)

    # assert result == poly_modexp_simple(orig_f, orig_e, g)

    return result

def poly_formal_derivative(f):
    return [gf_mul(c, e) for e, c in enumerate(f)][1:]


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
    c = poly_gcd(f, poly_formal_derivative(f))
    assert c == POLY_ONE, "polynomial is not square-free"

    from sage.all import GF, Integer, PolynomialRing
    modulus = [Integer((GF_POLY >> i) & 1) for i in range(GF_POLY.bit_length())]
    GF128 = GF(Integer(2)**Integer(128), modulus=modulus, names='b')
    P = PolynomialRing(GF128, 'x')
    Pf = P([GF128.from_integer(x) for x in f])
    t = time()
    tmp = Pf.roots()
    tmp = [x.to_integer() for x, _ in tmp]
    print("roots: " + repr(tmp) + f" {time() - t}")
    from root_find import bta, arm, sra
    t = time()
    tmp = bta(Pf)
    tmp = [x.to_integer() for x in tmp]
    print("bta: " + repr(tmp) + f" {time() - t}")
    t = time()
    tmp = arm(Pf)
    tmp = [x.to_integer() for x in tmp]
    print("arm: " + repr(tmp) + f" {time() - t}")
    t = time()
    tmp = sra(Pf)
    tmp = [x.to_integer() for x in tmp]
    print("sra: " + repr(tmp) + f" {time() - t}")

    t = time();
    tmp = poly_roots_bta(f)
    print("my bta: " + repr(tmp) + f" {time() - t}")

    # A Computational Introduction to Number Theory and Algebra (v2.5)
    # by Victor Shoup
    # https://www.shoup.net/ntb/ntb-v2_5.pdf

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
        rand = [getrandbits(128) for i in range(len(f) - 1)]
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
        gcm.encrypt(iv, randbytes(100+16), b""),
        gcm.encrypt(iv, randbytes(100+16), b""),
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
