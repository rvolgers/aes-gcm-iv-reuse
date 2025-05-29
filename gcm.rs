// rustc -C target-cpu=native -C opt-level=3 -g gcm.rs

// rustc -C target-cpu=native -C opt-level=3 -g gcm.rs -L. -lstatic:+verbatim=find_roots_ntl.o -lntl -lstdc++


#![feature(random)]
#![feature(array_chunks)]
#![feature(let_chains)]
#![feature(bigint_helper_methods)]

use std::iter;
use std::convert::TryInto;
use std::convert::TryFrom;
use std::random::random;
use std::random::DefaultRandomSource;
use std::random::RandomSource;
use std::ffi::c_int;
use std::time::Instant;

unsafe extern "C" {
    // int find_roots_ntl(unsigned char *poly, int coeff_count, unsigned char* roots)
    fn find_roots_ntl(poly: *const u8, coeff_count: c_int, roots: *mut u8) -> c_int;
}

fn find_roots_ntl_wrapper(poly: &[u128]) -> Vec<u128> {
    let poly_bytes: Vec<u8> = poly.iter().copied().flat_map(u128::to_le_bytes).collect();

    let mut roots_bytes: Vec<u8> = Vec::with_capacity(poly_bytes.len());

    unsafe {
        let n = find_roots_ntl(poly_bytes.as_ptr(), poly.len() as c_int, roots_bytes.spare_capacity_mut().as_mut_ptr() as *mut u8);

        roots_bytes.set_len(n as usize * 16);
    }

    roots_bytes.array_chunks().copied().map(u128::from_le_bytes).collect()
}

// high bit 128 is implicit
const GF_POLY: u128 = (1 << 7) | (1 << 2) | (1 << 1) | 1;


const GF_X: u128 = (1 << 1);
const GF_GEN: u128 = GF_X;

fn gf_from_bytes(b: [u8; 16]) -> u128 {
    u128::from_le_bytes(b).reverse_bits()
}

fn gf_to_bytes(x: u128) -> [u8; 16] {
    x.reverse_bits().to_le_bytes()
}

fn gf_random() -> u128 {
    random()
}

fn gf_reduce(lo: u128, hi: u128) -> u128 {
    let a = hi >> (128 - 1);
    let b = hi >> (128 - 2);
    let c = hi >> (128 - 7);

    let d = hi ^ a ^ b ^ c;

    let e = d << 1;
    let f = d << 2;
    let g = d << 7;

    let h = d ^ e ^ f ^ g;

    return lo ^ h;
}

#[cfg(all(
    target_feature = "pclmulqdq",
    any(target_arch = "x86", target_arch = "x86_64")
))]
fn gf_square(x: u128) -> u128 {
    // if we have a clmul intrinsic, we can't do better than that
    gf_mul(x, x)
}

// intersperses the bits of x with zero bits
fn spread64(x: u64) -> u128 {
    let x = x as u128;
    let x = (x | (x << 32)) & 0x00000000_ffffffff_00000000_ffffffff;
    let x = (x | (x << 16)) & 0x0000ffff_0000ffff_0000ffff_0000ffff;
    let x = (x | (x << 8)) & 0x00ff00ff_00ff00ff_00ff00ff_00ff00ff;
    let x = (x | (x << 4)) & 0x0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f;
    let x = (x | (x << 2)) & 0x33333333_33333333_33333333_33333333;
    let x = (x | (x << 1)) & 0x55555555_55555555_55555555_55555555;
    x
}

// reverse of spread64
fn unspread64(x: u128) -> u64 {
    let x = x & 0x55555555_55555555_55555555_55555555;
    let x = (x | (x >> 1)) & 0x33333333_33333333_33333333_33333333;
    let x = (x | (x >> 2)) & 0x0f0f0f0f_0f0f0f0f_0f0f0f0f_0f0f0f0f;
    let x = (x | (x >> 4)) & 0x00ff00ff_00ff00ff_00ff00ff_00ff00ff;
    let x = (x | (x >> 8)) & 0x0000ffff_0000ffff_0000ffff_0000ffff;
    let x = (x | (x >> 16)) & 0x00000000_ffffffff_00000000_ffffffff;
    let x = (x | (x >> 32)) as u64;
    x
}

#[cfg(not(all(
    target_feature = "pclmulqdq",
    any(target_arch = "x86", target_arch = "x86_64")
)))]
fn gf_square(x: u128) -> u128 {
    let lo = spread64(x as u64);
    let hi = spread64((x >> 64) as u64);
    gf_reduce(lo, hi)
}

#[cfg(not(all(
    target_feature = "pclmulqdq",
    any(target_arch = "x86", target_arch = "x86_64")
)))]
fn gf_mul(mut x: u128, mut y: u128) -> u128 {

    let mut result: u128 = 0;

    for _ in 0..128 {
        if y & 1 == 1 {
            result ^= x;
        }
        y >>= 1;

        let carry;
        (x, carry) = x.carrying_add(x, false);

        if carry {
            x ^= GF_POLY;
        }
    }

    result
}

#[cfg(all(
    target_feature = "pclmulqdq",
    any(target_arch = "x86", target_arch = "x86_64")
))]
fn gf_mul(a: u128, b: u128) -> u128 {
    use std::arch::x86_64::*;

    fn u128_to_m128(x: u128) -> __m128i {
        unsafe {
            _mm_set_epi64x((x >> 64) as u64 as i64, x as u64 as i64)
        }
    }
    
    fn m128_to_u128(x: __m128i) -> u128 {
        unsafe {
            ((_mm_extract_epi64(x, 1) as u64 as u128) << 64) |
            (_mm_cvtsi128_si64(x) as u64 as u128)
        }
    }

    // https://web.archive.org/web/20190806061845/https://software.intel.com/sites/default/files/managed/72/cc/clmul-wp-rev-2.02-2014-04-20.pdf
    unsafe {
        let a = u128_to_m128(a);
        let b = u128_to_m128(b);
        let XMMMASK = _mm_setr_epi32(-1, 0x0, 0x0, 0x0);
        let tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
        let tmp6 = _mm_clmulepi64_si128(a, b, 0x11);
        let tmp4 = _mm_shuffle_epi32(a,78);
        let tmp5 = _mm_shuffle_epi32(b,78);
        let tmp4 = _mm_xor_si128(tmp4, a);
        let tmp5 = _mm_xor_si128(tmp5, b);
        let tmp4 = _mm_clmulepi64_si128(tmp4, tmp5, 0x00);
        let tmp4 = _mm_xor_si128(tmp4, tmp3);
        let tmp4 = _mm_xor_si128(tmp4, tmp6);
        let tmp5 = _mm_slli_si128(tmp4, 8);
        let tmp4 = _mm_srli_si128(tmp4, 8);
        let tmp3 = _mm_xor_si128(tmp3, tmp5);
        let tmp6 = _mm_xor_si128(tmp6, tmp4);
        let tmp7 = _mm_srli_epi32(tmp6, 31);
        let tmp8 = _mm_srli_epi32(tmp6, 30);
        let tmp9 = _mm_srli_epi32(tmp6, 25);
        let tmp7 = _mm_xor_si128(tmp7, tmp8);
        let tmp7 = _mm_xor_si128(tmp7, tmp9);
        let tmp8 = _mm_shuffle_epi32(tmp7, 147);
        let tmp7 = _mm_and_si128(XMMMASK, tmp8);
        let tmp8 = _mm_andnot_si128(XMMMASK, tmp8);
        let tmp3 = _mm_xor_si128(tmp3, tmp8);
        let tmp6 = _mm_xor_si128(tmp6, tmp7);
        let tmp10 = _mm_slli_epi32(tmp6, 1);
        let tmp3 = _mm_xor_si128(tmp3, tmp10);
        let tmp11 = _mm_slli_epi32(tmp6, 2);
        let tmp3 = _mm_xor_si128(tmp3, tmp11);
        let tmp12 = _mm_slli_epi32(tmp6, 7);
        let tmp3 = _mm_xor_si128(tmp3, tmp12);
        m128_to_u128(_mm_xor_si128(tmp3, tmp6))
    }
}

#[cfg(all(
    target_feature = "pclmulqdq",
    any(target_arch = "x86", target_arch = "x86_64")
))]
fn gf_square64(x: u64) -> u128 {
    use std::arch::x86_64::*;
    unsafe {
        let a = _mm_set1_epi64x(x as i64);
        let b = _mm_clmulepi64_si128(a, a, 0x00);
        ((_mm_extract_epi64(b, 1) as u64 as u128) << 64) |
            (_mm_cvtsi128_si64(b) as u64 as u128)
    }
}

#[cfg(not(all(
    target_feature = "pclmulqdq",
    any(target_arch = "x86", target_arch = "x86_64")
)))]
fn gf_square64(x: u64) -> u128 {
    spread64(x)
}


// GF_SQRT_X == gf_pow(GF_X, 1 << 127)
// gf_square(GF_SQRT_X) == GF_X
const GF_SQRT_X: u128 = 0x24924924924924926db6db6db6db6da4;

fn gf_sqrt(x: u128) -> u128 {
    let even = unspread64(x) as u128;
    let odd = unspread64(x >> 1) as u128;
    gf_mul(odd, GF_SQRT_X) ^ even
}

fn gf_pow(mut x: u128, mut e: u128) -> u128 {
    let mut result: u128 = 1;

    while e != 0 {
        if e & 1 == 1 {
            result = gf_mul(result, x);
        }
        e >>= 1;

        x = gf_mul(x, x);
    }

    result
}

trait BitLen {
    fn bit_len(self) -> u32;
}

impl BitLen for u128 {
    fn bit_len(self) -> u32 {
        128 - self.leading_zeros()
    }
}

fn gf_inverse(x: u128) -> u128 {

    assert!(x != 0);

    // inverse by extended euclidean algorithm.
    // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures

    // u2/v2 are unused but useful for understanding the algorithm. compiler will optimize them out.
    let (mut u1, mut u2, mut u3) = (0u128, 1u128, GF_POLY);
    let (mut v1, mut v2, mut v3) = (1u128, 0u128, x);

    // hardcode u3.bit_len() for the first iteration to 129 to account for
    // the implicit 1 bit at position 128.
    // since the algorithm proceeds by iteratively cancelling the highest bit of u3,
    // after the first loop iteration everything will be within 128 bits again.
    // inside the loop v3 >= 2, so v3.bit_len() >= 2, so q <= 127, so v1 << q remains within 128 bits.
    let mut q = 129 - v3.bit_len();

    // notice that the GF_POLY term of the invariant is zero modulo GF_POLY.
    // so we exit once we have v1 such that gf_mul(x, v1) == 1,
    // which means v1 is the inverse of x modulo GF_POLY.
    // note that the final v3 is actually the gcd of GF_POLY and x.
    // we know it will always be 1, because GF_POLY is irreducible.
    // this also cleanly handles the x == 1 case.
    while v3 != 1 {
        // loop invariants (imagine GF_POLY includes its implicit 128th bit here as well)
        // q (= u3.bit_len() - v3.bit_len()) >= 0
        // clmul(x, u1) ^ clmul(GF_POLY, u2) == u3
        // clmul(x, v1) ^ clmul(GF_POLY, v2) == v3
        // (these last two are Bezout's Identity)
        // note that clmul does not reduce modulo a polynomial like gf_mul does.
        // after a single loop iteration the invariant holds with gf_mul as well,
        // although of course the GF_POLY term becomes zero.

        // this is a little different from usual extended euclidean, but relies on the same
        // rule that the normal one does: for any m, gcd(a + m * b, b) == gcd(a, b)
        // we can apply this repeatedly to make u3 and v3 smaller.
        u1 ^= v1 << q;
        u2 ^= v2 << q;
        u3 ^= v3 << q;

        let mut u3_len = u3.bit_len();
        let mut v3_len = v3.bit_len();
        if u3_len < v3_len {
            // ensure q >= 0 by swapping u and v
            (u1, u2, u3, u3_len, v1, v2, v3, v3_len) = (v1, v2, v3, v3_len, u1, u2, u3, u3_len);
        }
        q = u3_len - v3_len;
    }

    v1
}

fn gf_divmod(x: u128, y: u128) -> (u128, u128) {
    let mut q = 0;
    let mut r = x;

    while let Some(s) = r.bit_len().checked_sub(y.bit_len()) {
        r ^= y << s;
        q ^= 1 << s;
    }

    (q, r)
}

fn gf_div_exact(x: u128, y: u128) -> u128 {
    let (q, r) = gf_divmod(x, y);

    assert!(r == 0);

    q
}

fn gf_mod(x: u128, y: u128) -> u128 {
    let (q, r) = gf_divmod(x, y);
    r
}

fn gf_extended_euclidean(x: u128, y: u128) -> (u128, u128, u128, u128, u128) {
    let (mut u1, mut u2, mut u3) = (0u128, 1u128, y);
    let (mut v1, mut v2, mut v3) = (1u128, 0u128, x);

    while v3 != 0 {
        // assert!(gf_mul(x, u1) ^ gf_mul(y, u2) == u3);
        // assert!(gf_mul(x, v1) ^ gf_mul(y, v2) == v3);

        let (mut t1, mut t2, mut t3) = (u1, u2, u3);
        if let Some(q) = u3.bit_len().checked_sub(v3.bit_len()) {
            t1 ^= v1 << q;
            t2 ^= v2 << q;
            t3 ^= v3 << q;
        }
        (u1, u2, u3) = (v1, v2, v3);
        (v1, v2, v3) = (t1, t2, t3);
    }

    // assert!(v2 == gf_div_exact(x, u3));
    // assert!(v1 == gf_div_exact(y, u3));
    // println!("x={x:x} y={y:x} v1={v1:x} v2={v2:x} u1={u1:x} u2={u2:x} u3={u3:x}");
    // dbg!(gf_mul(x, v1) ^ gf_mul(y, v2));
    // assert!(gf_mul(x, v1) ^ gf_mul(y, v2) == 0);
    // dbg!(gf_mul(x, u1) ^ gf_mul(y, u2));
    // assert!(gf_mul(x, u1) ^ gf_mul(y, u2) == u3);

    (u3, v2, v1, u1, u2)
}

fn gf_formal_derivative(x: u128) -> u128 {
    (x >> 1) & 0x55555555_55555555_55555555_55555555
}

// this implements berlekamp factorization, which is not a great
// choice really, it just has a lot of fun moving parts and I
// wanted to play with it and see how much room there was for
// optimization. the answer is "more than you'd think, but not
// enough to make it competitive with nondeterministic methods"
// as always a lot of prototyping was done in the python file.
#[inline(never)]
fn gf_factor(u: u128) -> Vec<u128> {

    let power_of_x = u.trailing_zeros();
    let u = u >> power_of_x;

    // there is probably a lot to be gained by adding some good
    // trial division code to weed out more small factors here.
    // however, that's kind of boring and it makes optimizations
    // to the rest of the code less impactful.
    // it does actually help worst case performance though, since
    // many smaller factors means more work left to the last phase.

    let mut factors = vec![GF_X; power_of_x as usize];

    // FIXME think this through a bit more, are we handling all cases correctly?
    let (square, u_div_gcd, ..) = gf_extended_euclidean(u, gf_formal_derivative(u));
    let square_factors = if square != 1 {
        gf_factor(gf_sqrt(square))
    } else {
        vec![]
    };
    factors.extend_from_slice(&square_factors);
    factors.extend_from_slice(&square_factors);
    let u = u_div_gcd;

    if u == 1 {
        return factors;
    }

    let n = u.bit_len() as usize - 1;

    let mut q = [0u128; 128];
    let high_bit = 1u128 << n;

    let simple_part = n.div_ceil(2);

    // bits 1..simple_part are set
    let simple_mask = (1u128 << simple_part) - (1u128 << 1);

    // this implementation completely omits doing anything
    // for the first `simple_part` entries of q, instead
    // opting to do a little extra work to emulate the effect
    // of those loop iterations on each later row in a
    // single pass.
    // I've not bothered to actually shorten q (and c)
    // because it's unlikely to have much effect and
    // rearranging the indices is a hassle.
    let mut prev = 1 << (2 * (simple_part - 1));
    for i in simple_part..n {
        let mut tmp = prev;
        tmp <<= 1;
        if tmp & high_bit != 0 {
            tmp ^= u;
        }
        tmp <<= 1;
        if tmp & high_bit != 0 {
            tmp ^= u;
        }
        prev = tmp;

        tmp ^= 1 << i;

        let mut sq = gf_square64((tmp & simple_mask) as u64);
        for _ in 0..6 {
            tmp ^= sq;
            sq = gf_square64((sq & simple_mask) as u64);
        }
        assert!(sq == 0);

        q[i] = tmp;
    }

    let mut c = [-1i8; 128];
    for i in 1..simple_part {
        c[i] = i as i8;
    }
    let mut c_set = simple_mask;

    let mut v = vec![];

    // println!("n={n} simple_part={simple_part}");

    let byte_aligned_mask = 0x01010101_01010101_01010101_01010101;

    for k in simple_part..n {
        let q_k = q[k];
        let tmp = (q_k & !c_set);
        let mut j = (tmp & (byte_aligned_mask << (k & 7))).trailing_zeros();
        if j == 128 {
            j = tmp.trailing_zeros();
        }
        // assert!(j == 128 || (q_k >> j) & 1 != 0);
        // assert!(j == 128 || (c_set >> j) & 1 == 0);

        if j < 128 {
            let bit_j = 1 << j;
            let q_k_without_bit_j = q_k ^ bit_j;
            for i in k..n {
                if q[i] & bit_j != 0 {
                    q[i] ^= q_k_without_bit_j;
                }
            }
            c_set |= bit_j;
            c[k] = j as i8;
        } else {
            let mut v_r = (1 << k) | (q_k & simple_mask);

            for i in simple_part..k {
                if let Ok(s) = u32::try_from(c[i]) {
                    v_r |= ((q_k >> s) & 1) << i;
                }
            }

            v.push(v_r);
        }

        let q_k = q[k];
        let good = if (j as usize - (k & 7)) & 7 == 0 { 1 } else { 0 };
        // println!("q_{k} = {q_k:#0128b} {good}");
    }

    let mut tmp = vec![u,];
    for mut v_i in v.iter().copied() {
        v_i >>= v_i.trailing_zeros();

        let origlen = tmp.len();
        for j in 0..origlen {
            if tmp.len() == v.len() {
                break;
            }

            let gcd;
            (gcd, tmp[j], ..) = gf_extended_euclidean(tmp[j], v_i);
            if gcd != 1 {
                tmp.push(gcd);
            }
        }
    }

    factors.extend(tmp);
    factors
}

const POLY_ZERO: &[u128] = &[];

const POLY_ONE: &[u128] = &[1];

const POLY_X: &[u128] = &[0, 1];

fn poly_trim(mut f: Vec<u128>) -> Vec<u128> {
    while let Some(0) = f.last() {
        f.pop();
    }
    f
}

fn poly_eval(f: &[u128], x: u128) -> u128 {
    let mut result = 0;
    let mut xe = 1;
    for c in f.iter().copied() {
        result ^= gf_mul(xe, c);
        xe = gf_mul(xe, x);
    }
    result
}

fn poly_add(mut f: Vec<u128>, g: &[u128]) -> Vec<u128> {
    for (x, y) in f.iter_mut().zip(g.iter()) {
        *x ^= y
    }
    if g.len() > f.len() {
        f.extend_from_slice(&g[f.len()..]);
    }
    f
}

fn poly_divmod_simple(f: Vec<u128>, g: &[u128]) -> (Vec<u128>, Vec<u128>) {

    let lc_g = g.last().copied().expect("division by 0");

    let Some(lc_f) = f.last().copied() else {
        return (POLY_ZERO.to_owned(), POLY_ZERO.to_owned());
    };

    assert!(lc_f != 0 && lc_g != 0, "polynomials not trimmed");

    if f.len() < g.len() {
        return (POLY_ZERO.to_owned(), f.to_owned());
    }

    let lc_g_inv = gf_inverse(lc_g);

    let rdigits = g.len() - 1;
    let qdigits = f.len() - rdigits;

    let mut q = vec![0; qdigits];
    let mut r = f;

    for i in (0..qdigits).rev() {
        let lc = r[i + rdigits];
        if lc != 0 {
            let digit = gf_mul(lc, lc_g_inv);
            q[i] = digit;
            r[i + rdigits] = 0;
            for j in 0..rdigits {
                r[i + j] ^= gf_mul(digit, g[j]);
            }
        }
    }

    (q, poly_trim(r))
}

// in-place divmod of f by g.
// returns (qr, rdigits), where the lower rdigits of qr are the remainder and
// the upper digits are the quotient.
// in-place is possible because it is always possible to divide by lc(g) without remainder.
fn poly_divmod(f: Vec<u128>, g: &[u128]) -> (Vec<u128>, usize) {
    // return poly_divmod_simple(f, g);

    let lc_g = g.last().copied().expect("division by 0");

    let f_len = f.len();
    if f_len < g.len() {
        return (f, f_len);
    }

    let lc_g_inv = gf_inverse(lc_g);

    let rdigits = g.len() - 1;
    let qdigits = f.len() - rdigits;

    let mut qr = f;

    for i in (0..qdigits+rdigits).rev() {
        let mut tmp = qr[i];
        // add the product of the quotient digits and coefficients of g that affect this digit.
        // qr[i + 1] is the most recent digit of q, which pairs with the second-highest coefficient of
        // g, which is g[rdigits - 1]. from there, indices in qr proceed upwards and g downwards.
        // loop conditions ensure we stay within all of the following slices:
        // * g[..rdigits] (note that leading coefficient (at g[rdigits]) is only used as lc_g_inv)
        // * qr[i+1..] (digits of q, most recent first, but note next item)
        // * qr[rdigits..] (indices below that are part of r, not q)
        for j in rdigits.saturating_sub(i + 1)..rdigits.min(qdigits + rdigits - (i + 1)) {
            tmp ^= gf_mul(qr[i + 1 + j], g[rdigits - 1 - j]);
        }
        // if it's part of the quotient, divide it by leading coefficient of g,
        // producing a new quotient digit.
        // we do not need to store a remainder digit because it will always be perfectly
        // cancelled out.
        if i >= rdigits {
            tmp = gf_mul(tmp, lc_g_inv);
        }
        qr[i] = tmp;
    }

    (qr, rdigits)
}

fn poly_div(f: Vec<u128>, g: &[u128]) -> Vec<u128> {
    let (mut qr, rdigits) = poly_divmod(f, g);
    qr.copy_within(rdigits.., 0);
    qr.truncate(qr.len() - rdigits);
    qr
}

fn poly_mod(f: Vec<u128>, g: &[u128]) -> Vec<u128> {
    let (mut qr, rdigits) = poly_divmod(f, g);
    qr.truncate(rdigits);
    poly_trim(qr)
}

fn poly_monic(mut f: Vec<u128>) -> Vec<u128> {
    let Some(lc) = f.last().copied() else {
        return f;
    };

    if lc == 1 {
        return f;
    }

    let lc_inv = gf_inverse(lc);

    for x in f.iter_mut() {
        *x = gf_mul(*x, lc_inv);
    }

    f
}

fn poly_gcd(mut f: Vec<u128>, g: &[u128]) -> Vec<u128> {
    // assert!(f.len() > 0 && g.len() > 0);

    let mut g = g.to_owned();

    while g.len() != 0 {
        let tmp = poly_mod(f, &g);
        (f, g) = (g, tmp);
    }

    poly_monic(f)
}

fn poly_mul(f: &[u128], g: &[u128]) -> Vec<u128> {
    let mut result = vec![0; f.len() + g.len() - 1];
    for (ef, cf) in f.iter().copied().enumerate() {
        for (eg, cg) in g.iter().copied().enumerate() {
            result[ef + eg] ^= gf_mul(cf, cg);
        }
    }

    result
}

// much more efficient than poly_mul(f, f) because almost all terms cancel in characteristic 2
fn poly_square(mut f: Vec<u128>) -> Vec<u128> {
    let orig_len = f.len();
    f.resize(orig_len * 2 - 1, 0);
    for i in (0..orig_len).rev() {
        let c = std::mem::take(&mut f[i]);
        f[i * 2] = gf_square(c);
    }

    f
}

fn poly_modexp(mut f: Vec<u128>, mut e: u128, g: &[u128]) -> Vec<u128> {

    // assert!(f.len() > 0);

    // if f.len() > g.len() {
    //     f = poly_mod(f, g);
    // }

    assert!(f.len() > 0 && e > 0);

    let mut prod = POLY_ONE.to_owned();

    loop {
        if e & 1 == 1 {
            prod = poly_mod(poly_mul(&prod, &f), g);
        }

        e >>= 1;
        if e == 0 {
            break;
        }

        f = poly_mod(poly_square(f), g);

        assert!(f.len() > 0);
    }

    prod
}

fn poly_formal_derivative(g: &[u128]) -> Vec<u128> {
    g.iter().copied().enumerate().skip(1).map(|(e, c)| gf_mul(c, e as u128)).collect()
}

fn poly_roots_bta(f: &Vec<u128>) -> Vec<u128> {
    if f.len() == 2 {
        return vec![f[0]];
    } else if f.len() <= 1 {
        return vec![];
    }

    let a = gf_random();
    let mut axp = vec![0, a];
    let mut tr = axp.clone();
    for i in 0..127 {
        axp = poly_mod(poly_square(axp), f);
        tr = poly_add(tr, &axp);
    }

    let mut f = f.clone();

    let mut roots = vec![];
    for c in 0..2 {
        let nf = poly_gcd(f.clone(), &poly_trim(poly_add(tr.clone(), &[c])));
        if nf.len() >= 2 {
            roots.extend_from_slice(&poly_roots_bta(&nf));
            f = poly_div(f, &nf);
            if f.len() <= 2 {
                roots.extend_from_slice(&poly_roots_bta(&f));
                break;
            }
            tr = poly_mod(tr, &f);
        }
    }

    return roots;
}

fn pad_blocksize(data: &mut Vec<u8>) {
    if data.len() % 16 != 0 {
        data.extend(std::iter::repeat(0).take(16 - (data.len() % 16)));
    }
}

fn build_ghash_input(ciphertext: &[u8], auth_data: &[u8]) -> Vec<u8> {
    let mut data = vec![];
    data.extend_from_slice(ciphertext);
    pad_blocksize(&mut data);
    data.extend_from_slice(auth_data);
    pad_blocksize(&mut data);
    data.extend_from_slice(&(auth_data.len() as u64 * 8).to_be_bytes());
    data.extend_from_slice(&(ciphertext.len() as u64 * 8).to_be_bytes());
    data
}

fn ghash(auth_key: [u8; 16], ghash_input: &[u8]) -> [u8; 16] {
    let k = gf_from_bytes(auth_key);

    let mut hash = 0;
    for block in ghash_input.array_chunks().copied() {
        hash ^= gf_from_bytes(block);
        hash = gf_mul(hash, k);
    }

    gf_to_bytes(hash)
}

fn auth_tag(auth_key: [u8; 16], keyblock_0: [u8; 16], ciphertext: &[u8], auth_data: &[u8]) -> [u8; 16] {
    let ghash_input = build_ghash_input(ciphertext, auth_data);

    let mut hash = ghash(auth_key, &ghash_input);

    for i in 0..16 {
        hash[i] ^= keyblock_0[i];
    }

    hash
}

fn ghash_polynomial(data: &[u8]) -> Vec<u128> {
    assert!(data.len() % 16 == 0);

    let mut poly = Vec::with_capacity(data.len() / 16 + 1);
    poly.push(0);
    poly.extend(data.array_chunks().rev().copied().map(gf_from_bytes));

    poly
}

fn ciphertext_to_masked_polynomial(ciphertext: &[u8], auth_data: &[u8]) -> Vec<u128> {

    let (ciphertext, tag) = ciphertext.split_at(ciphertext.len() - 16);
    let tag: [u8; 16] = tag.try_into().unwrap();

    let mut poly = ghash_polynomial(&build_ghash_input(ciphertext, auth_data));

    poly[0] = gf_from_bytes(tag);

    poly
}

fn recover_auth_secret(mut ciphertexts: Vec<Vec<u8>>) -> Vec<([u8; 16], [u8; 16])> {

    println!("preprocessing");

    // sorting keeps short ciphertexts grouped together.
    // also makes it easy to spot and filter out duplicates.
    ciphertexts.sort_by(|a, b| a.cmp(b).reverse());

    let mut polys: Vec<_> = ciphertexts.iter().map(|c| ciphertext_to_masked_polynomial(c, &[])).collect();

    // xor each polynomial with another one
    for i in 0..(polys.len() - 1) {
        let tmp = std::mem::take(&mut polys[i]);
        polys[i] = poly_trim(poly_add(tmp, &polys[i+1]));
    }

    // we combined masked polynomials pairwise into unmasked ones,
    // leaving one unmasked polynomial at the end of the array.
    // remove it & keep it, since we need a single masked polynomial later
    // to produce a matching keyblock_0 for each auth secret candidate.
    let test_poly = polys.pop().unwrap();

    // zeroes occur if multiple identical ciphertexts were passed, remove them
    polys.retain(|f| f.len() > 0);

    assert!(polys.len() > 0);

    println!("gcd");

    // combine all polynomials using gcd
    let mut it = polys.iter().rev();
    let mut f = it.next().unwrap().clone();
    while f.len() > 2 && let Some(g) = it.next() {
        f = poly_gcd(f, g);
    }

    // println!("invoking NTL");
    // let roots = find_roots_ntl_wrapper(&f);
    // println!("NTL roots = {:?}", roots);

    assert!(f.len() >= 2);

    // check resulting polynomial is square-free
    let c = poly_gcd(poly_formal_derivative(&f), &f);
    assert!(c == POLY_ONE);

    let now = Instant::now();
    let roots = poly_roots_bta(&f);
    println!("BTA roots = {:?} {:?}", roots, now.elapsed());

    let now = Instant::now();

    println!("distinct degree factorization");

    // obtain product of linear factors ("distinct degree factorization")
    let h = poly_modexp(POLY_X.to_owned(), 1<<127, &f);
    let h = poly_modexp(h, 2, &f);
    let h = poly_trim(poly_add(h, POLY_X));
    f = poly_gcd(f, &h);

    println!("equal degree factorization");

    // break it into linear factors ("equal degree factorization")
    let mut factors = vec![f.clone()];
    while factors.len() != f.len() - 1 {
        let rand: Vec<u128> = (0..(f.len() - 1)).map(|_| random()).collect();
        // (2**128-1) / 3, but without overflowing
        let g = poly_modexp(rand, u128::MAX / 3, &f);

        let g = poly_trim(poly_add(g, POLY_ONE));
        factors = factors.into_iter().flat_map(|factor| {
            if factor.len() > 2 {
                let gcd = poly_gcd(factor.clone(), &g);
                if gcd.len() > 1 && gcd.len() < factor.len() {
                    return vec![poly_div(factor, &gcd), gcd];
                }
            }
            vec![factor]
        }).collect();
    }

    println!("factoring took {:?}", now.elapsed());

    factors.into_iter().map(|f| {
        assert!(f.len() == 2 && f[1] == 1, "factor not linear and monic");
        let auth_key = dbg!(f[0]);
        let keyblock_0 = poly_eval(&test_poly, auth_key);
        (gf_to_bytes(auth_key), gf_to_bytes(keyblock_0))
    }).collect()
}

fn main() {

    let mut keyblock_0: [u8; 16] = [0; 16];
    (&mut DefaultRandomSource).fill_bytes(&mut keyblock_0);

    let mut auth_key: [u8; 16] = [0; 16];
    (&mut DefaultRandomSource).fill_bytes(&mut auth_key);

    let mut iv: [u8; 12] = [0; 12];
    (&mut DefaultRandomSource).fill_bytes(&mut iv);

    let mut ciphertext1 = vec![0; 14 * 1024];
    (&mut DefaultRandomSource).fill_bytes(&mut ciphertext1);
    ciphertext1.extend_from_slice(&auth_tag(auth_key, keyblock_0, &ciphertext1, &[]));

    let mut ciphertext2 = vec![0; 14 * 1024];
    (&mut DefaultRandomSource).fill_bytes(&mut ciphertext2);
    ciphertext2.extend_from_slice(&auth_tag(auth_key, keyblock_0, &ciphertext2, &[]));

    let mut ciphertext3 = vec![0; 4 * 1024];
    (&mut DefaultRandomSource).fill_bytes(&mut ciphertext3);
    ciphertext3.extend_from_slice(&auth_tag(auth_key, keyblock_0, &ciphertext3, &[]));

    let ciphertexts = vec![ciphertext1, ciphertext2];

    let candidates = recover_auth_secret(ciphertexts);

    assert!(candidates.contains(&(auth_key, keyblock_0)));


    let mut ciphertext1 = vec![0; 16 * 100 * 1000];
    (&mut DefaultRandomSource).fill_bytes(&mut ciphertext1);
    let factor_inputs = ciphertext1.array_chunks().copied().map(gf_from_bytes);

    let now = Instant::now();
    for x in factor_inputs {
        // println!("factoring {x}");
        let factors = gf_factor(x);
        let mut tmp = 1;
        for f in factors {
            // println!("    factor {f}");
            tmp = gf_mul(tmp, f);
        }
        assert_eq!(tmp, x);
    }

    println!("gf_factor took {:?}", now.elapsed());

}