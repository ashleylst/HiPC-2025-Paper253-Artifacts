#include "sme_helpers.hpp"

#include <cassert>
#include <complex>

/*
 * Good: - small N allowed
 *
 * Neutral: - VLA for vlen >=6
 *
 * Bad : - only 1 subtile in use -> raw/waw dependency for each fmopa
 *       - only 6/VLEN of compute utilized
 *
 */

extern "C" void mm3x3t3xn_nega_1xtile(
        const std::complex<double>* a,
        const std::complex<double>* b,
        std::complex<double> *c,
        std::size_t N)
{
    std::size_t vlen = get_sme_vlen(); 
    // TODO: handle tail (whilelo p0.....)
    assert(0 == N % vlen);

    // Matrix columns should fit into single vector registers
    assert(vlen > 6);

    asm(
        "smstart\n\t"
        // predicate for N (b loads, opa N dim)
        "ptrue p0.d\n\t"
        // predicate for negating every other element
        "pfalse p2.b\n\t"
        "zip1 p2.d,p0.d,p2.d\n\t"

        "mov x6, #6\n\t"
        // predicate for M (a loads, opa M dim)
        "whilelo p1.d, xzr, x6\n\t"
        "lsl x7, x6, #1\n\t" // 12
        
        // z0.d = {Re(a00), Im(a00), Re(a10), Im(a10), Re(a20), Im(a20), 0, 0, ...}
        "ld1d {z0.d}, p1/z, [%[a]]\n\t"
        // z1.d = {Im(a00), Re(a00), Im(a10), Re(a10), Im(a20), Re(a20), 0, 0, ...}
        "revd z1.q, p1/m, z0.q\n\t"
        // z1.d = {-Im(a00), Re(a00), -Im(a10), Re(a10), -Im(a20), Re(a20), 0, 0, ...}
        "fneg z1.d, p2/m, z1.d\n\t"

        // z2.d = {Re(a01), Im(a01), Re(a11), Im(a11), Re(a21), Im(a21), 0, 0, ...}
        "ld1d {z2.d}, p1/z, [%[a], x6, LSL #3]\n\t"
        // z3.d = {Im(a01), Re(a01), Im(a11), Re(a11), Im(a21), Re(a21), 0, 0, ...}
        "revd z3.q, p1/m, z2.q\n\t"
        // z3.d = {-Im(a01), Re(a01), -Im(a11), Re(a11), -Im(a21), Re(a21), 0, 0, ...}
        "fneg z3.d, p2/m, z3.d\n\t"

        // z4.d = {Re(a02), Im(a02), Re(a12), Im(a12), Re(a22), Im(a22), 0, 0, ...}
        "ld1d {z4.d}, p1/z, [%[a], x7, LSL #3]\n\t"
        // z5.d = {Im(a02), Re(a02), Im(a12), Re(a12), Im(a22), Re(a22), 0, 0, ...}
        "revd z5.q, p1/m, z4.q\n\t"
        // z5.d = {-Im(a02), Re(a02), -Im(a12), Re(a12), -Im(a22), Re(a22), 0, 0, ...}
        "fneg z5.d, p2/m, z5.d\n\t"
        
        "mov x6, 0\n\t"
        "lsl x6, %[N], 1\n\t" // Nxcomplex
        "lsl x7,   x6, 1\n\t" // 2xN complex

        // array select regs for storing result
        "mov x12, 0\n\t"
        "mov x13, 2\n\t"
        "mov x14, 4\n\t"

        "mov x8,  %[b]\n\t"
        "mov x10, %[c]\n\t"
        "add x9, x8, x6, LSL #3\n\t" // x9 is after end of first row and can be used as stop condition 

        ".matmulloop%=:\n\t"
        "zero {za0.d}\n\t"
        // z16.d = {Re(b00), Re(b10), Re(b20), Re(b30), Re(b40), ...}
        // Z17.d = {Im(b00), Im(b10), Im(b20), Im(b30), Im(b40), ...}
        "ld2d {z16.d,z17.d}, p0/z, [x8]\n\t"

        /*
         *          | Re(a00)*Re(b00) Re(a00)*Re(b10) Re(a00)*Re(b20) Re(a00)*Re(b30) ...
         *          |
         * za0.d  = | Im(a00)*Re(b00) Im(a00)*Re(b10) Im(a00)*Re(b20) Im(a00)*Re(b30) ...
         *          |
         *          | Re(a10)*Re(b00) Re(a10)*Re(b10) Re(a10)*Re(b20) Re(a10)*Re(b30) ...
         *          |
         *          | Im(a10)*Re(b00) Im(a10)*Re(b10) Im(a10)*Re(b20) Im(a10)*Re(b30) ...
         *          |
         *          | Re(a20)*Re(b00) Re(a20)*Re(b10) Re(a20)*Re(b20) Re(a20)*Re(b30) ...
         *          |
         *          | Im(a20)*Re(b00) Im(a20)*Re(b10) Im(a20)*Re(b20) Im(a20)*Re(b30) ...
         *          |
         *          |       0                0              0               0        
         *          |       ⋮                ⋮              ⋮               ⋮          0
         */
        "fmopa za0.d, p1/m,p0/m, z0.d, z16.d\n\t"

        /*
         *          |  Re(a00)*Re(b00)  Re(a00)*Re(b10)  Re(a00)*Re(b20)  Re(a00)*Re(b30) 
         *          | -Im(a00)*Im(b00) -Im(a00)*Im(b10) -Im(a00)*Im(b20) -Im(a00)*Im(b30) ...
         *          |
         * za0.d  = |  Im(a00)*Re(b00)  Im(a00)*Re(b10)  Im(a00)*Re(b20)  Im(a00)*Re(b30)
         *          | +Re(a00)*Im(b00) +Re(a00)*Im(b10) +Re(a00)*Im(b20) +Re(a00)*Im(b30) ...
         *          |
         *          |  Re(a10)*Re(b00)  Re(a10)*Re(b10)  Re(a10)*Re(b20)  Re(a10)*Re(b30) 
         *          | -Im(a10)*Im(b00) -Im(a10)*Im(b10) -Im(a10)*Im(b20) -Im(a10)*Im(b30) ...
         *          |
         *          |  Im(a10)*Re(b00)  Im(a10)*Re(b10)  Im(a10)*Re(b20)  Im(a10)*Re(b30)
         *          | +Re(a10)*Im(b00) +Re(a10)*Im(b10) +Re(a10)*Im(b20) +Re(a10)*Im(b30) ...
         *          |
         *          |  Re(a20)*Re(b00)  Re(a20)*Re(b10)  Re(a20)*Re(b20)  Re(a20)*Re(b30) 
         *          | -Im(a20)*Im(b00) -Im(a20)*Im(b10) -Im(a20)*Im(b20) -Im(a20)*Im(b30) ...
         *          |
         *          |  Im(a20)*Re(b00)  Im(a20)*Re(b10)  Im(a20)*Re(b20)  Im(a20)*Re(b30)
         *          | +Re(a20)*Im(b00) +Re(a20)*Im(b10) +Re(a20)*Im(b20) +Re(a20)*Im(b30) ...
         *          |
         *          |       0                0              0               0        
         *          |       ⋮                ⋮              ⋮               ⋮          0
         */
        "fmopa za0.d, p1/m,p0/m, z1.d, z17.d\n\t"

        // z18.d = {Re(b01), Re(b11), Re(b21), 0, 0, ...}
        // Z19.d = {Im(b01), Im(b11), Im(b21), 0, 0, ...}
        "ld2d {z18.d,z19.d}, p0/z, [x8, x6, LSL #3]\n\t"
        // Add ax1 o bx1 outer product
        "fmopa za0.d, p1/m,p0/m, z2.d, z18.d\n\t"
        "fmopa za0.d, p1/m,p0/m, z3.d, z19.d\n\t"

        // z20.d = {Re(b02), Re(b12), Re(b22), 0, 0, ...}
        // Z21.d = {Im(b02), Im(b12), Im(b22), 0, 0, ...}
        "ld2d {z20.d,z21.d}, p0/z, [x8, x7, LSL #3]\n\t"
        // Add ax2 o bx2 outer product
        "fmopa za0.d, p1/m,p0/m, z4.d, z20.d\n\t"
        "fmopa za0.d, p1/m,p0/m, z5.d, z21.d\n\t"

        




        // Store results without a loop, since we know exactly 6 rows have valid values
        "mova z16.d, p0/m, za0h.d[w12,0]\n\t"
        "mova z17.d, p0/m, za0h.d[w12,1]\n\t"
        "mova z18.d, p0/m, za0h.d[w13,0]\n\t"
        "mova z19.d, p0/m, za0h.d[w13,1]\n\t"
        "mova z20.d, p0/m, za0h.d[w14,0]\n\t"
        "mova z21.d, p0/m, za0h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z16.d,z17.d}, p0, [x10]\n\t"
        "st2d {z18.d,z19.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z20.d,z21.d}, p0, [x10, x7, LSL #3]\n\t"


        // next 2xVLEN-size batch
        "incb x8, ALL, MUL #2\n\t"
        "incb x10, ALL, MUL #2\n\t"
        "cmp x8, x9\n\t"
        "b.ne .matmulloop%=\n\t"

        "smstop\n\t"
        : [dummy_c] "+m"(*(std::complex<double>(*)[])c)
        : [a] "r" (a), [b] "r" (b), [c] "r" (c), [N] "r" (N)
        : "z0", "z1", "z2", "z3", "z4", "z5",
          "z16", "z17", "z18", "z19", "z20", "z21",
          "za", "x6", "x7", "x8", "x9", "x10", "x12", "x13", "x14"
    );
}
