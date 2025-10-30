#include "sme_helpers.hpp"

#include <cassert>
#include <complex>

/*
 * Good: - VLA for vlen >= 3 (separate implementation required for vlen < 3)
 *       - small N allowed
 *
 * Bad : - 1xrevd + 1xfneg for each ld1d in the loop
 *       - only 1 subtile in use -> raw/waw dependency for each fmopa
 *       - only 3/VLEN of compute utilized
 *
 */
extern "C" void mm3x3t3xn_negb_1xtile(
        const std::complex<double>* a,
        const std::complex<double>* b,
        std::complex<double> *c,
        std::size_t N)
{
    /*  Layouts:
     *           _       _
     *    | a00 | | a01 | | a02
     * A: v a10 | v a11 | v a12
     *    | a20 | | a21 | | a22
     *    |_____| |_____|
     *
     * 
     *     <------------------N---------------->
     * B:  ------->--------------->----------...
     *     b00 b10 b20 b30 b40 b50 b60 b70 ...  |
     *    ,----------------------------------...'
     *    '------->--------------->----------...
     *     b01 b11 b21 b31 b41 b51 b61 b71 ...  |
     *    ,----------------------------------...'
     *    '------->--------------->----------...
     *     b02 b12 b22 b32 b42 b52 b62 b72 ...  
     *
     * Possibly better layout for B:
     *     <---- VLEN ---->
     * B:  ------->--------, ,--------->--------, ,---- ...
     *     b00 b10 b20 b30 | |  b40 b50 b60 b70 | |
     *    ,----------------' | ,----------------' |
     *    '------->--------, | '------->--------, |
     *     b01 b11 b21 b31 | |  b41 b51 b61 b71 | |
     *    ,----------------' | ,----------------' |
     *    '------->----------' '------->----------'
     *     b02 b12 b22 b32      b42 b52 b62 b72 
     *
     */
    std::size_t vlen = get_sme_vlen();
    // TODO: handle tail (whilelo p0.....)
    assert(0 == N % vlen);

    // Matrix columns should fit into single vector registers
    assert(vlen > 3);

    asm(
        "smstart\n\t"
        // predicate for N (b loads, opa N dim)
        "ptrue p0.d\n\t"
        // predicate for negating every other element
        "pfalse p2.b\n\t"
        "zip1 p2.d,p0.d,p2.d\n\t"

        "mov x6, #3\n\t"
        // predicate for M (a loads, opa M dim)
        "whilelo p1.d, xzr, x6\n\t"
        "lsl x6, x6, #1\n\t" // 6
        "lsl x7, x6, #1\n\t" // 12
        // z0.d = {Re(a00), Re(a10), Re(a20), 0, 0, ...}
        // Z1.d = {Im(a00), Im(a10), Im(a10), 0, 0, ...}
        "ld2d {z0.d,z1.d}, p1/z, [%[a]]\n\t"
        // z2.d = {Re(a01), Re(a11), Re(a21), 0, 0, ...}
        // Z3.d = {Im(a01), Im(a11), Im(a21), 0, 0, ...}
        "ld2d {z2.d,z3.d}, p1/z, [%[a], x6, LSL #3]\n\t"
        // z4.d = {Re(a02), Re(a12), Re(a22), 0, 0, ...}
        // Z5.d = {Im(a02), Im(a12), Im(a22), 0, 0, ...}
        "ld2d {z4.d,z5.d}, p1/z, [%[a], x7, LSL #3]\n\t"
        
        "mov x6, 0\n\t"
        "lsl x6, %[N], 1\n\t" // Nxcomplex
        "lsl x7,   x6, 1\n\t" // 2xN complex

        // array select regs for storing result
        "mov x12, 0\n\t"
        "mov x13, 2\n\t"

        "mov x8,  %[b]\n\t"
        "mov x10, %[c]\n\t"
        "add x9, x8, x6, LSL #3\n\t"
        ".matmulloop%=:\n\t"
        "zero {za0.d}\n\t"
        // z16.d = {Re(b00), Im(b00), Re(b10), Im(b10), ... }
        "ld1d {z16.d}, p0/z, [x8]\n\t"
        // z17.d = {Im(b00), Re(b00), Im(b10), Re(b10), ... }
        "revd z17.q, p0/m, z16.q\n\t"
        // z17.d = {-Im(b00), Re(b00), -Im(b10), Re(b10), ... }
        "fneg z17.d, p2/m, z17.d\n\t"

        /*
         *          | Re(a00)*Re(b00) Re(a00)*Im(b00) Re(a00)*Re(b10) Re(a00)*Im(b10) ...
         *          |
         * za0.d  = | Re(a10)*Re(b00) Re(a10)*Im(b00) Re(a10)*Re(b10) Re(a10)*Im(b10) ...
         *          |
         *          | Re(a20)*Re(b00) Re(a20)*Im(b00) Re(a20)*Re(b10) Re(a20)*Im(b10) ...
         *          |                                                              
         *          |       0                0              0               0        
         *          |       ⋮                ⋮              ⋮               ⋮          0
         */
        "fmopa za0.d, p1/m,p0/m, z0.d, z16.d\n\t"
        /*
         *          |  Re(a00)*Re(b00)  Re(a00)*Im(b00)  Re(a00)*Re(b10)  Re(a00)*Im(b10) 
         *          | -Im(a00)*Im(b00) +Im(a00)*Re(b00) -Im(a00)*Im(b10) +Im(a00)*Re(b10) ...
         *          |
         * za0.d  = |  Re(a10)*Re(b00)  Re(a10)*Im(b00)  Re(a10)*Re(b10)  Re(a10)*Im(b10) 
         *          | -Im(a10)*Im(b00) +Im(a10)*Re(b00) -Im(a10)*Im(b10) +Im(a10)*Re(b10) ...
         *          |
         *          |  Re(a20)*Re(b00)  Re(a20)*Im(b00)  Re(a20)*Re(b10)  Re(a20)*Im(b10) 
         *          | -Im(a20)*Im(b00) +Im(a20)*Re(b00) -Im(a20)*Im(b10) +Im(a20)*Re(b10) ...
         *          |                                                              
         *          |       0                0              0               0        
         *          |       ⋮                ⋮              ⋮               ⋮          0
         */
        "fmopa za0.d, p1/m,p0/m, z1.d, z17.d\n\t"

        // z18.d = {Re(b01), Im(b01), Re(b11), Im(b11), ... }
        "ld1d {z18.d}, p0/z, [x8, x6, LSL #3]\n\t"
        // z19.d = {Im(b01), Re(b01), Im(b11), Re(b11), ... }
        "revd z19.q, p0/m, z18.q\n\t"
        // z19.d = {-Im(b01), Re(b01), -Im(b11), Re(b11), ... }
        "fneg z19.d, p2/m, z19.d\n\t"

        // Add ax1 o bx1 outer product
        "fmopa za0.d, p1/m,p0/m, z2.d, z18.d\n\t"
        "fmopa za0.d, p1/m,p0/m, z3.d, z19.d\n\t"
        

        // z20.d = {Re(b02), Im(b02), Re(b12), Im(b12), ... }
        "ld1d {z20.d}, p0/z, [x8, x7, LSL #3]\n\t"
        // z21.d = {Im(b02), Re(b02), Im(b12), Re(b12), ... }
        "revd z21.q, p0/m, z20.q\n\t"
        // z21.d = {-Im(b02), Re(b02), -Im(b12), Re(b12), ... }
        "fneg z21.d, p2/m, z21.d\n\t"

        // Add ax2 o bx2 outer product
        "fmopa za0.d, p1/m,p0/m, z4.d, z20.d\n\t"
        "fmopa za0.d, p1/m,p0/m, z5.d, z21.d\n\t"


        // Store results without a loop, since we know exactly 3 rows have valid values
        "mova z16.d, p0/m, za0h.d[w12,0]\n\t"
        "mova z17.d, p0/m, za0h.d[w12,1]\n\t"
        "mova z18.d, p0/m, za0h.d[w13,0]\n\t"

        "st1d z16.d, p0,[x10]\n\t"
        "st1d z17.d, p0,[x10, x6, LSL #3]\n\t"
        "st1d z18.d, p0,[x10, x7, LSL #3]\n\t"


        // next VLEN-size batch
        "incb x8\n\t"
        "incb x10\n\t"
        "cmp x8, x9\n\t"
        "b.ne .matmulloop%=\n\t"

        "smstop\n\t"
        : [dummy_c] "+m"(*(std::complex<double>(*)[])c)
        : [a] "r" (a), [b] "r" (b), [c] "r" (c), [N] "r" (N)
        : "z0", "z1", "z2", "z3", "z4", "z5",
          "z16", "z17", "z18", "z19", "z20", "z21",
          "za", "x6", "x7", "x8", "x9", "x10", "x12", "x13"
    );
}
