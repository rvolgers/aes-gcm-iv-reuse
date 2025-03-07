#!/usr/bin/env python3

# only used for a basic AES-ECB primitive
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

from random import getrandbits
from time import time

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

def gf_mul(x, y):
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

def gf_inverse(x):
    assert x != 0, "zero has no inverse"

    # inverse by extended euclidean algorithm.
    # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
    # https://crypto.stackexchange.com/questions/12956/multiplicative-inverse-in-operatornamegf28/12962#12962
    (u1, u3) = (0, GF_POLY)
    (v1, v3) = (1, x)

    while v3 != 0:
        (t1, t3) = (u1, u3)
        q = u3.bit_length() - v3.bit_length()
        if q >= 0:
            t1 ^= v1 << q
            t3 ^= v3 << q
        (u1, u3) = (v1, v3)
        (v1, v3) = (t1, t3)

    if u1 & (1 << 128):
        u1 ^= GF_POLY

    # could've also used fermat's little theorem and gf_pow().
    # easier to understand but a lot more expensive than the above.
    # assert u1 == gf_pow(x, (1<<128) - 2)

    return u1


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
def poly_divmod(f, g, lc_g_inv = None):

    assert g != POLY_ZERO, "cannot divide by zero"

    if f == POLY_ZERO:
        return (POLY_ZERO, POLY_ZERO)

    assert g[-1] != 0, "g is not trimmed"
    assert f[-1] != 0, "f is not trimmed"

    qdigits = len(f) - len(g) + 1
    if qdigits <= 0:
        return (POLY_ZERO, f)

    if lc_g_inv is None:
        lc_g_inv = gf_inverse(g[-1])

    q = [0] * qdigits
    r = f[:]

    for i in reversed(range(qdigits)):
        lc = r[i + len(g) - 1]
        if lc != 0:
            digit = gf_mul(lc, lc_g_inv)
            q[i] = digit
            r[i + len(g) - 1] = 0
            for j in range(len(g) - 1):
                r[i + j] ^= gf_mul(digit, g[j])

    # assert r[-qdigits:] == [0] * qdigits

    return (q, poly_trim(r[:-qdigits]))

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

def poly_modexp(f, e, g):
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

def poly_formal_derivative(f):
    return [gf_mul(c, e) for e, c in enumerate(f)][1:]


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

class AES_GCM:
    def __init__(self, key):
        # basic aes encryption primitive, see _aes_ecb_encrypt
        self._aes_ecb = Cipher(algorithms.AES(key), modes.ECB())

        # the secret used in the ghash function.
        self._auth_key = self._aes_ecb_encrypt(b'\x00' * 16)

        print("auth key: " + hex(gf_from_bytes(self._auth_key)))
        print("auth key (decimal): " + str(gf_from_bytes(self._auth_key)))

    def _aes_ecb_encrypt(self, data):
        encryptor = self._aes_ecb.encryptor()
        return encryptor.update(data) + encryptor.finalize()

    def _counter_from_iv(self, iv):
        # if the length of the passed iv is not 12 it is passed through the
        # ghash function and the resulting 128 bit value is *directly* used
        # as the initial value of the counter.
        if len(iv) != 12:
            # *some* minimum length is obviously required for uniqueness.
            # but the exact limits are not consistent between implementations.
            # however most implementations impose limits similar to this.
            assert 8 <= len(iv) <= 128
            hash = ghash(self._auth_key, build_ghash_input(iv))
            return int.from_bytes(hash, byteorder='big')

        # if the iv is exactly 12 bytes it is used as the upper 96 bits
        # of the counter. the lower bits are left zero except for the +1.
        # the +1 here is presumably to prevent an all-zero iv from resulting
        # in self._auth_key being equal to the auth tag keyblock?
        return (int.from_bytes(iv, byteorder='big') << 32) + 1

    def _keyblock_from_counter(self, counter):
        bcounter = counter.to_bytes(length=16, byteorder='big')
        return self._aes_ecb_encrypt(bcounter)

    def _crypt_common(self, init_counter, data):
        """implements CTR mode encryption and decryption"""

        # +1 because the first one is reserved for masking the auth tag
        counter = init_counter + 1

        result = b''
        for block in split_blocks(data):
            keyblock = self._keyblock_from_counter(counter)
            result += xor_bytes(keyblock[:len(block)], block)
            counter += 1

        return result

    def _auth_tag(self, init_counter, auth_data, ciphertext):
        """calculate the authentication tag"""

        ghash_input = build_ghash_input(ciphertext, auth_data)

        hash = ghash(self._auth_key, ghash_input)

        keyblock = self._keyblock_from_counter(init_counter)
        return xor_bytes(hash, keyblock)

    def encrypt(self, iv, plaintext, auth_data=b''):
        init_counter = self._counter_from_iv(iv)

        ciphertext = self._crypt_common(init_counter, plaintext)
        assert len(ciphertext) == len(plaintext)

        auth_tag = self._auth_tag(init_counter, auth_data, ciphertext)

        return ciphertext + auth_tag

    def decrypt(self, iv, ciphertext, auth_data=b''):
        init_counter = self._counter_from_iv(iv)

        # last 16 bytes of ciphertext are auth tag, split it off
        assert len(ciphertext) >= 16
        ciphertext, ciphertext_auth_tag = ciphertext[:-16], ciphertext[-16:]

        good_auth_tag = self._auth_tag(init_counter, auth_data, ciphertext)

        if ciphertext_auth_tag != good_auth_tag:
            raise ValueError("invalid tag")

        plaintext = self._crypt_common(init_counter, ciphertext)

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

    if len(f) == 2:
        print("early out because gcd produced a linear factor")
        return [f[0]]

    assert len(f) > 2

    print(f"reduced to degree {len(f) - 1}")

    # TODO handle the case where poly is not square-free
    # this does not appear to ever trigger in practice.
    c = poly_gcd(f, poly_formal_derivative(f))
    assert c == POLY_ONE, "polynomial is not square-free"

    # A Computational Introduction to Number Theory and Algebra (v2.5)
    # by Victor Shoup
    # https://www.shoup.net/ntb/ntb-v2_5.pdf

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

    for x in factors:
        assert len(x) == 2 and x[-1] == 1

    return [x[0] for x in factors]


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

    my_gcm = AES_GCM(key)
    plaintext = my_gcm.decrypt(iv, ciphertext)

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
        gcm.encrypt(iv, randbytes(100), b""),
        gcm.encrypt(iv, randbytes(100), b""),
        #gcm.encrypt(iv, randbytes(1000), b""),
    ])

    assert gf_from_bytes(my_gcm._auth_key) in recovered

    print(f"possible auth secret values: {recovered}")
