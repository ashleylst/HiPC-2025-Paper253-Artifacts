#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include <cblas.h>

inline std::size_t get_sme_vlen()
{
    std::size_t vlen = 0;
    asm(
        "smstart\n\t"
        "mov %[vlen],0\n\t"
        "incd %[vlen]\n\t"
        "smstop\n\t"
        : [vlen] "=r" (vlen)
        :
        :
    );  

    return vlen;
}

/*
 * Good: - VLA for vlen >= 3 (separate implementation required for vlen < 3)
 *       - small N allowed
 *
 * Bad : - 1xrevd + 1xfneg for each ld1d in the loop
 *       - only 1 subtile in use -> raw/waw dependency for each fmopa
 *       - only 3/VLEN of compute utilized
 *
 */
extern "C" void mm3x3t3xn_negb_1xtile(const double* a, const double* b, double *c, std::size_t N)
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
        : [dummy_c] "+m"(*(double(*)[])c)
        : [a] "r" (a), [b] "r" (b), [c] "r" (c), [N] "r" (N)
        : "z0", "z1", "z2", "z3", "z4", "z5",
          "z16", "z17", "z18", "z19", "z20", "z21",
          "za", "x6", "x7", "x8", "x9", "x10", "x12", "x13"
    );
}

/*
 * Good: - small N allowed
 *
 * Neutral: - VLA for vlen >=6
 *
 * Bad : - only 1 subtile in use -> raw/waw dependency for each fmopa
 *       - only 6/VLEN of compute utilized
 *
 */
extern "C" void mm3x3t3xn_nega_1xtile(const double* a, const double* b, double *c, std::size_t N)
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
        : [dummy_c] "+m"(*(double(*)[])c)
        : [a] "r" (a), [b] "r" (b), [c] "r" (c), [N] "r" (N)
        : "z0", "z1", "z2", "z3", "z4", "z5",
          "z16", "z17", "z18", "z19", "z20", "z21",
          "za", "x6", "x7", "x8", "x9", "x10", "x12", "x13", "x14"
    );
}

/*
 * Good: - maximized za tile usage (8)
 *
 * Neutral: - VLA for vlen >=6
 *
 * Bad : - N must be multiple of 8*VLEN
 *       - only 6/VLEN of compute utilized
 *
 * Possible further optimizations:
 *  - preload/increase distance between ld2d and fmopa
 *  - better memory layout, so that all loads can be made
 *    contiguous
 *  - better interleaving of ld/st/fmopa
 *
 *  Better memory layout:
 *     <----- 8xVLEN ----->    <----- 8xVLEN -----> 
 * B:  ------->--------..., ,--------->--------..., ,---- ...
 *     b00 b10 b20 b30... | |  bX0 bY0 bZ0 bW0... | |
 *    ,----------------...' | ,----------------...' |
 *    '------->--------..., | '------->--------..., |
 *     b01 b11 b21 b31... | |  bX1 bY1 bZ1 bW1... | |
 *    ,----------------...' | ,----------------...' |
 *    '------->--------...--' '------->--------...--'
 *     b02 b12 b22 b32...      bX2 bY2 bZ2 bW2... 
 *
 *
 *
 */
extern "C" void mm3x3t3xn_nega_8xtile(const double* a, const double* b, double *c, std::size_t N)
{
    std::size_t vlen = get_sme_vlen(); 
    // TODO: handle tail (whilelo p0.....)
    assert(0 == N % 8*vlen);

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
        
        "mov x11, x8\n\t"

        "zero {za0.d}\n\t"
        // z15.d = {Re(b00), Re(b10), Re(b20), Re(b30), Re(b40), ...}
        // Z16.d = {Im(b00), Im(b10), Im(b20), Im(b30), Im(b40), ...}
        "ld2d {z15.d,z16.d}, p0/z, [x11]\n\t"

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
        "fmopa za0.d, p1/m,p0/m, z0.d, z15.d\n\t"

        "ld2d {z17.d,z18.d}, p0/z, [x11, #2, MUL VL]\n\t"
        "fmopa za1.d, p1/m,p0/m, z0.d, z17.d\n\t"

        "ld2d {z19.d,z20.d}, p0/z, [x11, #4, MUL VL]\n\t"
        "fmopa za2.d, p1/m,p0/m, z0.d, z19.d\n\t"

        "ld2d {z21.d,z22.d}, p0/z, [x11, #6, MUL VL]\n\t"
        "fmopa za3.d, p1/m,p0/m, z0.d, z21.d\n\t"

        "ld2d {z23.d,z24.d}, p0/z, [x11, #8, MUL VL]\n\t"
        "fmopa za4.d, p1/m,p0/m, z0.d, z23.d\n\t"

        "ld2d {z25.d,z26.d}, p0/z, [x11, #10, MUL VL]\n\t"
        "fmopa za5.d, p1/m,p0/m, z0.d, z25.d\n\t"

        "ld2d {z27.d,z28.d}, p0/z, [x11, #12, MUL VL]\n\t"
        "fmopa za6.d, p1/m,p0/m, z0.d, z27.d\n\t"

        "ld2d {z29.d,z30.d}, p0/z, [x11, #14, MUL VL]\n\t"
        "fmopa za7.d, p1/m,p0/m, z0.d, z29.d\n\t"


        "add x11, x11, x6, LSL #3\n\t" // next row

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
        "fmopa za0.d, p1/m,p0/m, z1.d, z16.d\n\t"
        "fmopa za1.d, p1/m,p0/m, z1.d, z18.d\n\t"
        "fmopa za2.d, p1/m,p0/m, z1.d, z20.d\n\t"
        "fmopa za3.d, p1/m,p0/m, z1.d, z22.d\n\t"
        "fmopa za4.d, p1/m,p0/m, z1.d, z24.d\n\t"
        "fmopa za5.d, p1/m,p0/m, z1.d, z26.d\n\t"
        "fmopa za6.d, p1/m,p0/m, z1.d, z28.d\n\t"
        "fmopa za7.d, p1/m,p0/m, z1.d, z30.d\n\t"


        // ax1 o bx1

        "ld2d {z15.d,z16.d}, p0/z, [x11]\n\t"
        "fmopa za0.d, p1/m,p0/m, z2.d, z15.d\n\t"

        "ld2d {z17.d,z18.d}, p0/z, [x11, #2, MUL VL]\n\t"
        "fmopa za1.d, p1/m,p0/m, z2.d, z17.d\n\t"

        "ld2d {z19.d,z20.d}, p0/z, [x11, #4, MUL VL]\n\t"
        "fmopa za2.d, p1/m,p0/m, z2.d, z19.d\n\t"

        "ld2d {z21.d,z22.d}, p0/z, [x11, #6, MUL VL]\n\t"
        "fmopa za3.d, p1/m,p0/m, z2.d, z21.d\n\t"

        "ld2d {z23.d,z24.d}, p0/z, [x11, #8, MUL VL]\n\t"
        "fmopa za4.d, p1/m,p0/m, z2.d, z23.d\n\t"

        "ld2d {z25.d,z26.d}, p0/z, [x11, #10, MUL VL]\n\t"
        "fmopa za5.d, p1/m,p0/m, z2.d, z25.d\n\t"

        "ld2d {z27.d,z28.d}, p0/z, [x11, #12, MUL VL]\n\t"
        "fmopa za6.d, p1/m,p0/m, z2.d, z27.d\n\t"

        "ld2d {z29.d,z30.d}, p0/z, [x11, #14, MUL VL]\n\t"
        "fmopa za7.d, p1/m,p0/m, z2.d, z29.d\n\t"

        "add x11, x11, x6, LSL #3\n\t" // next row

        "fmopa za0.d, p1/m,p0/m, z3.d, z16.d\n\t"
        "fmopa za1.d, p1/m,p0/m, z3.d, z18.d\n\t"
        "fmopa za2.d, p1/m,p0/m, z3.d, z20.d\n\t"
        "fmopa za3.d, p1/m,p0/m, z3.d, z22.d\n\t"
        "fmopa za4.d, p1/m,p0/m, z3.d, z24.d\n\t"
        "fmopa za5.d, p1/m,p0/m, z3.d, z26.d\n\t"
        "fmopa za6.d, p1/m,p0/m, z3.d, z28.d\n\t"
        "fmopa za7.d, p1/m,p0/m, z3.d, z30.d\n\t"

        // ax2 o bx2

        "ld2d {z15.d,z16.d}, p0/z, [x11]\n\t"
        "fmopa za0.d, p1/m,p0/m, z4.d, z15.d\n\t"

        "ld2d {z17.d,z18.d}, p0/z, [x11, #2, MUL VL]\n\t"
        "fmopa za1.d, p1/m,p0/m, z4.d, z17.d\n\t"

        "ld2d {z19.d,z20.d}, p0/z, [x11, #4, MUL VL]\n\t"
        "fmopa za2.d, p1/m,p0/m, z4.d, z19.d\n\t"

        "ld2d {z21.d,z22.d}, p0/z, [x11, #6, MUL VL]\n\t"
        "fmopa za3.d, p1/m,p0/m, z4.d, z21.d\n\t"

        "ld2d {z23.d,z24.d}, p0/z, [x11, #8, MUL VL]\n\t"
        "fmopa za4.d, p1/m,p0/m, z4.d, z23.d\n\t"

        "ld2d {z25.d,z26.d}, p0/z, [x11, #10, MUL VL]\n\t"
        "fmopa za5.d, p1/m,p0/m, z4.d, z25.d\n\t"

        "ld2d {z27.d,z28.d}, p0/z, [x11, #12, MUL VL]\n\t"
        "fmopa za6.d, p1/m,p0/m, z4.d, z27.d\n\t"

        "ld2d {z29.d,z30.d}, p0/z, [x11, #14, MUL VL]\n\t"
        "fmopa za7.d, p1/m,p0/m, z4.d, z29.d\n\t"

        // no next row

        // interleave last fmopas with stores
        // (We can probably reuse z registers as they are 
        //  being freed up, i.e using z16 after this fmopa,
        //  z18 after the next, etc... maximizing the 
        //  distances before a reg is reused, but it's
        //  complicated, so just use the same registers 
        //  every time for now)
        "fmopa za0.d, p1/m,p0/m, z5.d, z16.d\n\t"
        "mova z6.d, p0/m, za0h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za0h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za0h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za0h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za0h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za0h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za1.d, p1/m,p0/m, z5.d, z18.d\n\t"
        "mova z6.d, p0/m, za1h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za1h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za1h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za1h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za1h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za1h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za2.d, p1/m,p0/m, z5.d, z20.d\n\t"
        "mova z6.d, p0/m, za2h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za2h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za2h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za2h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za2h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za2h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za3.d, p1/m,p0/m, z5.d, z22.d\n\t"
        "mova z6.d, p0/m, za3h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za3h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za3h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za3h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za3h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za3h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za4.d, p1/m,p0/m, z5.d, z24.d\n\t"
        "mova z6.d, p0/m, za4h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za4h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za4h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za4h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za4h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za4h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za5.d, p1/m,p0/m, z5.d, z26.d\n\t"
        "mova z6.d, p0/m, za5h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za5h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za5h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za5h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za5h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za5h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za6.d, p1/m,p0/m, z5.d, z28.d\n\t"
        "mova z6.d, p0/m, za6h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za6h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za6h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za6h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za6h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za6h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"

        "fmopa za7.d, p1/m,p0/m, z5.d, z30.d\n\t"
        "mova z6.d, p0/m, za7h.d[w12,0]\n\t"
        "mova z7.d, p0/m, za7h.d[w12,1]\n\t"
        "mova z8.d, p0/m, za7h.d[w13,0]\n\t"
        "mova z9.d, p0/m, za7h.d[w13,1]\n\t"
        "mova z10.d, p0/m, za7h.d[w14,0]\n\t"
        "mova z11.d, p0/m, za7h.d[w14,1]\n\t"

        // Interleaved store
        "st2d {z6.d,z7.d}, p0, [x10]\n\t"
        "st2d {z8.d,z9.d}, p0, [x10, x6, LSL #3]\n\t"
        "st2d {z10.d,z11.d}, p0, [x10, x7, LSL #3]\n\t"
        "incb x10, ALL, MUL #2\n\t"


        // next 16xVLEN-size batch
        "incb x8, ALL, MUL #16\n\t"
        "cmp x8, x9\n\t"
        "b.ne .matmulloop%=\n\t"

        "smstop\n\t"
        : [dummy_c] "+m"(*(double(*)[])c)
        : [a] "r" (a), [b] "r" (b), [c] "r" (c), [N] "r" (N)
        : "z0", "z1", "z2", "z3", "z4", "z5",
          "z6", "z7", "z8", "z9", "z10", "z11",
          "z15", "z16", "z17", "z18", "z19", "z20", "z21",
          "z22", "z23", "z24", "z25", "z26", "z27", "z28",
          "z29", "z30",
          "za", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14"
    );
}



int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << argv[0] << " N\n";
        return -1;
    }
    std::uint64_t N = std::stoul(argv[1]);

    std::size_t vlen = get_sme_vlen();

    std::vector<double> A(3*3*2);
    std::vector<double> B(3*N*2);
    std::vector<double> C(3*N*2);

    std::iota(A.begin(), A.end(), 0.0);
    
    std::iota(B.begin(), B.end(), 0.0);
    std::transform(B.begin(), B.end(), B.begin(), [](double a){return a*0.3;});

    std::string bleftmargin = "                                   ";
    std::cout << bleftmargin << "B:\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            numberstr << B[(i*N+j)*2] << "+" << B[(i*N+j)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    std::cout << "A:\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        for(std::size_t j = 0; j < 3; j++)
        {
            std::stringstream numberstr;
            numberstr << A[(j*3+i)*2] << "+" << A[(j*3+i)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    std::complex<double> one{1.0,0.0};
    std::fill(C.begin(), C.end(), 0.0);
    cblas_zgemm(CBLAS_ORDER::CblasRowMajor,
                CBLAS_TRANSPOSE::CblasTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                3, N, 3, &one,
                A.data(), 3,
                B.data(), N, &one,
                C.data(), N);

    std::cout << bleftmargin << "C (reference):\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            numberstr << C[(i*N+j)*2] << "+" << C[(i*N+j)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    std::fill(C.begin(), C.end(), 0.0);
    if (vlen >= 8)
    {
        if (0 == (N % (8*vlen)))
        {
            std::cout << "Computing using 8xtile neg-a version\n";
            mm3x3t3xn_nega_8xtile(A.data(), B.data(), C.data(), N);
        }
        else
        {
            std::cout << "Computing using 1xtile neg-a version\n";
            mm3x3t3xn_nega_1xtile(A.data(), B.data(), C.data(), N);
        }
    }
    else if (vlen >= 4)
    {
        std::cout << "Computing using neg-b version\n";
        mm3x3t3xn_negb_1xtile(A.data(), B.data(), C.data(), N);
    }
    std::cout << bleftmargin << "C (ours):\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            numberstr << C[(i*N+j)*2] << "+" << C[(i*N+j)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    return 0;
}
