#include <cstdlib>

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
