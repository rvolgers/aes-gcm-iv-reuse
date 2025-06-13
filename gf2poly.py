

# count trailing zeroes
def count_trailing_zeros(x):
    # implementations differ in how/if ctz(0) is defined, so avoid it
    assert x != 0

    # (x ^ (x - 1)) isolates the lowest 1 bit
    # 0b110000 - 1 = 0b101111
    # 0b110000 ^ (0b110000 - 1) = 0b10000

    return (x ^ (x - 1)).bit_length() - 1

def gf_deg(u):
    return u.bit_length() - 1

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

def gf_divmod(x, y):
    q = 0
    r = x

    s = r.bit_length() - y.bit_length()
    while s >= 0:
        r ^= y << s
        q ^= 1 << s
        s = r.bit_length() - y.bit_length()

    return (q, r)

def gf_mod(x, y):
    return gf_divmod(x, y)[1]

def gf_div(x, y):
    q, r = gf_divmod(x, y)
    assert r == 0
    return q

# fully-featured extended euclidean algorithm on binary polynomials.
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

def gf_gcd_split(x, y):
    return tuple(gf_extended_euclidean(x, y)[:3])

def gf_gcd(x, y):
    return gf_extended_euclidean(x, y)[0]

def gf_invmod(x, y):
    assert x != 0, "zero has no inverse"

    (d, _, _, inv, _) = gf_extended_euclidean(x, y)

    assert d == 1, "x and y must be coprime"

    if inv.bit_length() == y.bit_length():
        inv ^= y

    assert inv.bit_length() < y.bit_length()

    return inv

def gf_modexp(x, e, m):

    result = 1
    sq = x
    for i in range(e.bit_length()):
        if (e >> i) & 1:
            result = gf_mod(gf_mul_noreduce(result, sq), m)

        # TODO use gf_square_noreduce
        sq = gf_mod(gf_mul_noreduce(sq, sq), m)
    
    return result

def gfmat_minus_identity(m):
    m = m[:]
    for i in range(len(m)):
        m[i] ^= 1 << i
    return m

def gfmat_kernel(m):

    # most of this was reused from gf_factor_berlekamp, just tidied up.
    # could have rewritten that to use this function, but that would lose
    # all the comments/asserts/exposition in that function only relevant
    # to that specific use case.

    # see "Algorithm N: Null space algorithm" TAOCP vol II 4.6.2 (p. 439)
    # although variable names are different and a lot of stuff has been
    # specialized / optimized for GF(2).

    # we assume the matrix is n*n square. we can't really tell because python
    # ints have no length independent of their value.
    n = len(m)
    m = m[:]

    pivots = []
    pivots_used = 0
    kernel = []
    for r in range(n):
        mr = m[r]

        # print(f"m[{r:3d}] = {''.join(map(str, int_to_bitlist(mr, n)))}")

        # valid pivot columns must be set in mr and not have been used
        available = mr & ~pivots_used

        if available:
            # just pick the first available one.
            # I believe we have freedom to choose, except that the rationale
            # for skipping the first (r - 1) rows in the next loop requires
            # the criteria for selection to be consistent between rows.
            # messing with the forward progress of the algorithm by changing
            # previous rows like that is probably also a bad thing.
            pc = count_trailing_zeros(available)

            # zero out all the other columns in the current row, by adding
            # (i.e. xor'ing) column pc to other columns which are 1 in mr.
            # note that previous rows are 0 in column pc. so the operation
            # would do nothing, and we can skip them.
            for j in range(r, n):
                if m[j] & (1 << pc):
                    m[j] ^= mr ^ (1 << pc) # xor with mr-without-bit-pc

            pivots_used |= 1 << pc
            pivots.append(pc)
        else:
            k = 1 << r
            assert len(pivots) == r
            for pr, pc in enumerate(pivots):
                if pc is not None and mr & (1 << pc):
                    k |= 1 << pr
            kernel.append(k)
            pivots.append(None)

        # loop invariant: for every row k in kernel, k * m == 0
        # we can use a generic matmul even though we want binary matmul,
        # because bitlist_to_int discards all but the low bit.
        # for k in kernel:
        #      tmp = mat_mul([int_to_bitlist(k, n)], [int_to_bitlist(mr, n) for mr in m])
        #      tmp = bitlist_to_int(tmp[0])
        #      assert tmp == 0

    # for r, mr in enumerate(m):
    #     print(f"m[{r:3d}] = {''.join(map(str, int_to_bitlist(mr, n)))}")

    # for i, k in enumerate(kernel):
    #     print(f"k[{i:3d}] = {''.join(map(str, int_to_bitlist(k, n)))}")

    return kernel