#include <array>
#include <iostream>
#include <iomanip>
#include <cassert>


int main(int argc, char *argv[])
{
    std::size_t vlen = 0;
    asm(
        "smstart\n\t"
        "mov %[vlen],0\n\t"
        "incd %[vlen]\n\t"
        "smstop\n\t"
        : [vlen] "=X" (vlen)
        :
        :
    );
    assert(vlen == 4);

    std::array c{
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    std::array d{
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    std::array a{1.0,2.0,3.0,4.0};
    std::array b{0.5,0.7,0.9,1.1};

    std::cout << "    B: ";
    for(std::size_t i = 0; i < vlen; i++)
    {
        std::cout << std::setw(4) << b[i] << " ";
    }
    std::cout << "\n";

    std::cout << "  A: \n";
    for(std::size_t i = 0; i < vlen; i++)
    {
        std::cout << std::setw(4) << a[i] << "\n";
    }

    asm(
        "smstart\n\t"

        "ptrue p0.d\n\t"
        "ld1d {z0.d},p0/z,[%[a]]\n\t"
        "ld1d {z1.d},p0/z,[%[b]]\n\t"

        "zero {za0.d}\n\t"

        "fmopa za0.d,p0/m,p0/m,z0.d,z1.d\n\t"

        // The correct way to store tiles is using a loop (keeps it VLA)
        "mov x12, 0\n\t"
        "mov x1, %[c]\n\t"
        ".storeloop:\n\t"
        // move za arrays to sve vectors
        "mov z0.d, p0/m, za0h.d[w12,0]\n\t"
        "mov z1.d, p0/m, za0h.d[w12,1]\n\t" // max offset for .d is 1 because SVL is a multiple of 128,
                                            // which is 2 FP64 elements, i.e also 2 rows?
        "st1d z0.d, p0,[x1, #0, MUL VL]\n\t"
        "st1d z1.d, p0,[x1, #1, MUL VL]\n\t"

        "incb x1\n\t" // lazy, better: store 2xSVLb into a reg and add it to x1
        "incb x1\n\t"
        "add x12, x12, #2\n\t"
        "cmp x12, %[vlen]\n\t"
        "b.lt .storeloop\n\t"


        // can we do something simpler if we don't care about VLA?
        "mov x12, 0\n\t"
        "str za[w12,0],[%[d], #0, MUL VL]\n\t"
        "mov x12, 7\n\t"
        "str za[w12,1],[%[d], #1, MUL VL]\n\t"
        "mov x12, 14\n\t"
        "str za[w12,2],[%[d], #2, MUL VL]\n\t"
        "mov x12, 21\n\t"
        "str za[w12,3],[%[d], #3, MUL VL]\n\t"
        // ehhh... i dunno

        "smstop\n\t"
        : [dummy_c] "+m"(*(double(*)[])c.data()),
          [dummy_d] "+m"(*(double(*)[])d.data())
        : [a] "r" (a.data()),
          [b] "r" (b.data()),
          [c] "r" (c.data()),
          [d] "r" (d.data()),
          [vlen] "r" (vlen)
        : "x1","x12","p0","z0","z1","z2","z3","za"
    );

    std::cout << "  Outer Product (proper store):\n";

    for (std::size_t i = 0; i < vlen; i++)
    {
        std::cout << "      ";
        for(std::size_t j = 0; j < vlen; j++)
        {
            std::cout << " "<< std::setw(4)  << c[i*vlen+j];
        }
        std::cout << "\n";
    }

    std::cout << "  Outer Product (hacky store):\n";

    for (std::size_t i = 0; i < vlen; i++)
    {
        std::cout << "      ";
        for(std::size_t j = 0; j < vlen; j++)
        {
            std::cout << " "<< std::setw(4)  << d[i*vlen+j];
        }
        std::cout << "\n";
    }


    return 0;
}
