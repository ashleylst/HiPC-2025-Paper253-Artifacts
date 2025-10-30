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

    std::array c22{
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    std::array c12{
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    std::array c21{
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    std::array a0{1.0,2.0,3.0,4.0};
    std::array a1{5.0,6.0,7.0,8.0};
    std::array b0{0.3,0.5,0.7,1.1};
    std::array b1{1.3,1.7,1.9,2.3};

    std::cout << "    B: ";
    for(std::size_t i = 0; i < vlen; i++)
    {
        std::cout << std::setw(4) << b0[i] << " ";
    }
    for(std::size_t i = 0; i < vlen; i++)
    {
        std::cout << std::setw(4) << b1[i] << " ";
    }
    std::cout << "\n";

    std::cout << "  A: \n";
    for(std::size_t i = 0; i < vlen; i++)
    {
        std::cout << std::setw(4) << a0[i] << "\n";
    }
    for(std::size_t i = 0; i < vlen; i++)
    {
        std::cout << std::setw(4) << a1[i] << "\n";
    }

    asm(
        "smstart\n\t"

        "ptrue p0.d\n\t"
        "ld1d {z0.d},p0/z,[%[a0]]\n\t"
        "ld1d {z1.d},p0/z,[%[a1]]\n\t"
        "ld1d {z16.d},p0/z,[%[b0]]\n\t"
        "ld1d {z17.d},p0/z,[%[b1]]\n\t"

        "zero {za0.d}\n\t"
        "zero {za1.d}\n\t"
        "zero {za2.d}\n\t"
        "zero {za3.d}\n\t"

        "fmop4a za0.d,{z0.d-z1.d},{z16.d-z17.d}\n\t"
        "fmop4a za1.d,z0.d,{z16.d-z17.d}\n\t"
        "fmop4a za2.d,{z0.d-z1.d},z16.d\n\t"
        "fmop4a za3.d,z0.d,z16.d\n\t"

        "mov x12, 0\n\t"
        "mov x1, %[c22]\n\t"
        ".storeloop0:\n\t"
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
        "b.lt .storeloop0\n\t"


        "mov x12, 0\n\t"
        "mov x1, %[c12]\n\t"
        ".storeloop1:\n\t"
        // move za arrays to sve vectors
        "mov z0.d, p0/m, za1h.d[w12,0]\n\t"
        "mov z1.d, p0/m, za1h.d[w12,1]\n\t" // max offset for .d is 1 because SVL is a multiple of 128,
                                            // which is 2 FP64 elements, i.e also 2 rows?
        "st1d z0.d, p0,[x1, #0, MUL VL]\n\t"
        "st1d z1.d, p0,[x1, #1, MUL VL]\n\t"

        "incb x1\n\t" // lazy, better: store 2xSVLb into a reg and add it to x1
        "incb x1\n\t"
        "add x12, x12, #2\n\t"
        "cmp x12, %[vlen]\n\t"
        "b.lt .storeloop1\n\t"


        "mov x12, 0\n\t"
        "mov x1, %[c21]\n\t"
        ".storeloop2:\n\t"
        // move za arrays to sve vectors
        "mov z0.d, p0/m, za2h.d[w12,0]\n\t"
        "mov z1.d, p0/m, za2h.d[w12,1]\n\t" // max offset for .d is 1 because SVL is a multiple of 128,
                                            // which is 2 FP64 elements, i.e also 2 rows?
        "st1d z0.d, p0,[x1, #0, MUL VL]\n\t"
        "st1d z1.d, p0,[x1, #1, MUL VL]\n\t"

        "incb x1\n\t" // lazy, better: store 2xSVLb into a reg and add it to x1
        "incb x1\n\t"
        "add x12, x12, #2\n\t"
        "cmp x12, %[vlen]\n\t"
        "b.lt .storeloop2\n\t"

        "smstop\n\t"
        : [dummy_c22] "+m"(*(double(*)[])c22.data()),
          [dummy_c12] "+m"(*(double(*)[])c12.data()),
          [dummy_c21] "+m"(*(double(*)[])c21.data())
        : [a0] "r" (a0.data()),
          [a1] "r" (a1.data()),
          [b0] "r" (b0.data()),
          [b1] "r" (b1.data()),
          [c22] "r" (c22.data()),
          [c12] "r" (c12.data()),
          [c21] "r" (c21.data()),
          [vlen] "r" (vlen)
        : "x1","x12","p0","z0","z1","z16","z17","za"
    );

    std::cout << "  fmop4a za, {a0-a1},{b0-b1}:\n";

    for (std::size_t i = 0; i < vlen; i++)
    {
        std::cout << "      ";
        for(std::size_t j = 0; j < vlen; j++)
        {
            std::cout << " "<< std::setw(4)  << c22[i*vlen+j];
        }
        std::cout << "\n";
    }

    std::cout << "  fmop4a za, a0,{b0-b1}:\n";

    for (std::size_t i = 0; i < vlen; i++)
    {
        std::cout << "      ";
        for(std::size_t j = 0; j < vlen; j++)
        {
            std::cout << " "<< std::setw(4)  << c12[i*vlen+j];
        }
        std::cout << "\n";
    }

    std::cout << "  fmop4a za, {a0-a1},b0:\n";

    for (std::size_t i = 0; i < vlen; i++)
    {
        std::cout << "      ";
        for(std::size_t j = 0; j < vlen; j++)
        {
            std::cout << " "<< std::setw(4)  << c21[i*vlen+j];
        }
        std::cout << "\n";
    }


    return 0;
}
