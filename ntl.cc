// sudo apt install libntl-dev

// g++ -c -o find_roots_ntl.o ntl.cc

#include <NTL/pair_GF2EX_long.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2EXFactoring.h>
#include <assert.h>

using namespace std;
using namespace NTL;

// NTL defines import/export of bytes in terms of unsigned chars
static_assert(sizeof(unsigned char) == 1, "wide char not supported");

#define ELEM_SIZE (128 / 8)

extern "C" int find_roots_ntl(unsigned char *poly, int coeff_count, unsigned char* roots)
{
    GF2X P(INIT_MONO, 128);
    SetCoeff(P, 0);
    SetCoeff(P, 1);
    SetCoeff(P, 2);
    SetCoeff(P, 7);

    GF2E::init(P);

    GF2EX f(INIT_SIZE, coeff_count);

    for (int i = 0; i < coeff_count; i++) {
        // note that GF2XFromBytes normalizes, so no need to worry about that.
        // as for the polynomial itself, caller should ensure it's normalized.
        SetCoeff(f, i, to_GF2E(GF2XFromBytes(poly + i * ELEM_SIZE, ELEM_SIZE)));
    }

    MakeMonic(f);

    // cin >> f;

    vec_pair_GF2EX_long result;
    result = CanZass(f, 1);
    // result = berlekamp(f, 1);

    // cout << result << "\n";

    int roots_count = 0;

    for (long i = 0; i < result.length(); i++) {

        // this does nothing, but double-check range of roots_count just to be safe
        if (roots_count == coeff_count) break;

        if (deg(result[i].a) == 1) {
            MakeMonic(result[i].a);
            BytesFromGF2X(roots + roots_count * ELEM_SIZE, rep(ConstTerm(result[i].a)), ELEM_SIZE);
            roots_count++;
        }
    }

    return roots_count;
}