# Modeling the performance

With the help of QEMU's 'insn' plugin we can count the occurences of certain instructions. [modelperf.py](modelperf.py) is a script that can use QEMU to estimate the performance of an implementation using a simple instruction-cost performance model.

Usage:
```
$ python3 modelperf.py --help
usage: modelperf.py [-h] --qemu-bin QEMU_BIN --qemu-insn QEMU_INSN [--sme-bits SME_BITS] [--instructions INSTRUCTIONS [INSTRUCTIONS ...]] [--dvz-instructions DVZ_INSTRUCTIONS [DVZ_INSTRUCTIONS ...]] [--cost-overrides COST_OVERRIDES [COST_OVERRIDES ...]]
                    [--roi-iterations ROI_ITERATIONS]
                    program {blas,cpp,sme_nega1,sme_nega8,sme_negb}

benchmark the SME performance using QEMU instruction counting

positional arguments:
  program               program to run
  {blas,cpp,sme_nega1,sme_nega8,sme_negb}
                        computation method to use

options:
  -h, --help            show this help message and exit
  --qemu-bin QEMU_BIN   Path to the qemu binary to launch the benchmark with
  --qemu-insn QEMU_INSN
                        Path to the qemu insn tcg plugin
  --sme-bits SME_BITS   SME vector length
  --instructions INSTRUCTIONS [INSTRUCTIONS ...]
                        instructions for which to measure counts
  --dvz-instructions DVZ_INSTRUCTIONS [DVZ_INSTRUCTIONS ...]
                        instructions that can be either fp,neon or sve for which to measure counts for each of their variants
  --cost-overrides COST_OVERRIDES [COST_OVERRIDES ...]
                        Cost overrides for example '--const-overrides fmopa=2 fmops=2', default cost is 1 for every instruction
  --roi-iterations ROI_ITERATIONS
                        Iterate ROI this many times to isolate it better from the driver application
```

The script executes the binary with N=64 and N=128, saves the differences between the instruction counts, calculates the number of instructions per N and multiplies them with the instruction costs. The sum is displayed in the end of the output.

NOTE: The insn plugin won't count instructions unless they are explicitly specified, the default instructions that the script investigated are specified in the `default_instructions` and `default_dvz_instructions` in the script source.

## Usage example

```
$ python3 modelperf.py \
  --qemu-bin $QEMU_DIRECTORY/bin/qemu-aarch64 \
  --qemu-insn $QEMU_BUILD_DIR/tests/tcg/plugins/libinsn.so \
  --cost-overrides "mov(SVE)=0" "st1d(SVE)=2" "st2d(SVE)=2" \
  -- mm3x3t3xN blas
ldr d : 61.25 per N
str d : 36.75 per N
fmla(SVE) : 11.625 per N
fmls(SVE) : 3.75 per N
fmadd(FP) : 0.375 per N
fnmsub(FP) : 0.375 per N
fmul(FP) : 0.75 per N
fmul(SVE) : 2.625 per N
mov(SVE) : 6.0 per N
ld1rd(SVE) : 8.75 per N
ld2d(SVE) : 1.875 per N
st1d(SVE) : 2.25 per N
st2d(SVE) : 2.75 per N
============================================
instruction-cost performance model: 8840.0
```
This will execute the driver using the blas implementation, override the cost of the mov instruction with 0 (register renaming) and the cost of stores with 2 (stores more expensive than loads)

## Performance estimates table


### `roi_iterations=1`

No overrides:

| implementation | cost |
| -------------- | ---- |
| `blas`         | 8904 |
| `cpp`          | 6888 |
| `sme_negb`     | 6720 |
| `sme_nega1`    | 6528 |
| `sme_nega8`    | 6528 |

`"mov(SVE)=0" "st1d(SVE)=2" "st2d(SVE)=2"`:

| implementation | cost |
| -------------- | ---- |
| `blas`         | 8840 |
| `cpp`          | 7056 |
| `sme_negb`     | 6816 |
| `sme_nega1`    | 6600 |
| `sme_nega8`    | 6600 |

`"mov(SVE)=0" "st1d(SVE)=2" "st2d(SVE)=2" "fmopa=2" "fmops=2"` (SME operations cost double of SVE operations):

| implementation | cost |
| -------------- | ---- |
| `blas`         | 8840 |
| `cpp`          | 7056 |
| `sme_negb`     | 6912 |
| `sme_nega1`    | 6648 |
| `sme_nega8`    | 6648 |

### `roi_iterations=100`

No overrides:

| implementation |  cost  |
| -------------- | ------ |
| `blas`         | 253632 |
| `cpp`          |  56832 |
| `sme_negb`     |  40032 |
| `sme_nega1`    |  20832 |
| `sme_nega8`    |  20832 |

`"mov(SVE)=0" "st1d(SVE)=2" "st2d(SVE)=2"`:

| implementation |  cost  |
| -------------- | ------ |
| `blas`         | 232976 |
| `cpp`          |  64176 |
| `sme_negb`     |  40176 |
| `sme_nega1`    |  18576 |
| `sme_nega8`    |  18576 |

`"mov(SVE)=0" "st1d(SVE)=2" "st2d(SVE)=2" "fmopa=2" "fmops=2"` (SME operations cost double of SVE operations):

| implementation |  cost  |
| -------------- | ------ |
| `blas`         | 232976 |
| `cpp`          |  64176 |
| `sme_negb`     |  49776 |
| `sme_nega1`    |  23376 |
| `sme_nega8`    |  23376 |

### `roi_iterations=100` and `N=8,16`

`"mov(SVE)=0" "st1d(SVE)=2" "st2d(SVE)=2" "fmopa=2" "fmops=2"` (SME operations cost double of SVE operations):

| implementation |  cost  |
| -------------- | ------ |
| `blas`         |  29122 |
| `cpp`          |   8022 |
| `sme_negb`     |   6222 |
| `sme_nega1`    |   2922 |

## Performance discussion

- The ROI is not isolated from the rest of the application, the instructions are also counted for the data initialization, etc... 
- This performance model cannot account for optimizations of instruction throughput, i.e. latency hiding of the `sme_nega8` version over `sme_nega1` is not reflected in the cost
- the BLAS library produces lots of overhead and doesn't seem to be worth calling for matrices this small
- the additional fneg/revd instructions significantly increase the cost of the `sme_negb` variant
- Custom SME implementations can be expected to outperform both BLAS and the C++ implementation
- The ROI is better isolated when using `--roi-iterations 100`.
  - Here, for the best SME variant, we can see a speedup of 3.45 over the C++ version and a 12.5 over the BLAS library
  - doubling the cost of SME operations still leads to a speedup of 2.75 over C++ and 9.97 over BLAS
- For the `N=8,16` case the best SME variant has a speedup of 2.75 over the C++ version and a 9.97 over the BLAS library just like the `N=64,128` case. This is expected


