import argparse

import subprocess
import re


#,match=fmla\ d,match=fmla\ v,match=fmls\ z,match=fmls\ d,match=fmls\ v,match=fmul\ z,match=fmul\ d,match=fmul\ v,match=fneg,match=revd,match=mov\ z,match=mova\ z,match=ldr\ d,match=str\ d,match=ld1rd\ {z,match=ld1d\ {z,match=ld2d\ {z,match=st1d\ {z,match=st2d\ {z,match=fmopa\ z,match=fmadd,match=fmsub \

default_dvz_instructions = ['fmla','fmls','fmadd','fmsub','fnmadd','fnmsub','fmul','fneg','revd','mov','mova','ld1rd','ld1d','ld2d','st1d','st2d']
default_instructions = ['ldr d','str d','fmopa','fmops']

def extract_instructions(output : str,
                         instructions : list[str],
                         dvz_instructions : list[str]):
    all_instructions = instructions + [inst + variant for inst in dvz_instructions for variant in [' d',' v',' z',' {d',' {v',' {z']]

    inst_counts = {}

    variants = [' d',' v',' z',' {d',' {v',' {z']
    variant_names = ['FP','NEON','SVE','FP','NEON','SVE']


    for inst in instructions:
        for line in output.split('\n'):
            inst_regex = re.compile(fr'Match: {inst}, hits (\d+)')
            inst_match = re.match(inst_regex, line)
            if inst_match:
                inst_counts[inst] = int(inst_match.group(1))
                continue

    for inst in dvz_instructions:
        for variant,vname in zip(variants,variant_names):
            for line in output.split('\n'):
                inst_regex = re.compile(fr'Match: {inst}{variant}, hits (\d+)')
                inst_match = re.match(inst_regex, line)
                if inst_match:
                    if inst+f"({vname})" in inst_counts.keys():
                        inst_counts[inst+f"({vname})"] += int(inst_match.group(1))
                    else:
                        inst_counts[inst+f"({vname})"] = int(inst_match.group(1))
                    continue

    return inst_counts


def run_program(qemu_bin : str, program : str, sme_bits : int,
                insn_path : str,
                instructions : list[str],
                dvz_instructions : list[str],
                method : str,
                N : int, roi_iterations : int):

    all_instructions = instructions + [inst + variant for inst in dvz_instructions for variant in [' d',' v',' z',' {d',' {v',' {z']]

    match_str = ',match='+",match=".join(all_instructions)



    cmd = [qemu_bin, f"-cpu", f"max,sve{sme_bits}=on,sme{sme_bits}=on", f"-plugin", f"{insn_path}{match_str}", '-d', 'plugin', program, str(N), method, str(roi_iterations)]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = result.stderr.decode()
    #print(output)
    #print(result.stdout)

    return extract_instructions(output, instructions, dvz_instructions)

def parse_arguments():

    parser = argparse.ArgumentParser(description="benchmark the SME performance using QEMU instruction counting")

    parser.add_argument("--qemu-bin", required=True, type=str, help="Path to the qemu binary to launch the benchmark with")
    parser.add_argument("--qemu-insn", required=True, type=str, help="Path to the qemu insn tcg plugin")
    parser.add_argument("--sme-bits", type=str, default="512", help="SME vector length")
    parser.add_argument("--instructions", nargs="+", type=str, default=default_instructions, help="instructions for which to measure counts")
    parser.add_argument("--dvz-instructions", nargs="+", type=str, default=default_dvz_instructions, help="instructions that can be either fp,neon or sve for which to measure counts for each of their variants")
    parser.add_argument("--cost-overrides", nargs="+", type=str, default=[], help="Cost overrides for example '--const-overrides fmopa=2 fmops=2', default cost is 1 for every instruction")
    parser.add_argument("--roi-iterations", type=int, default=1, help="Iterate ROI this many times to isolate it better from the driver application")
    parser.add_argument("program", nargs=1, type=str, help="program to run")
    parser.add_argument("method", nargs=1, type=str, choices=['blas','cpp','sme_nega1','sme_nega8','sme_negb'], help="computation method to use")

    return parser.parse_args()


def inst_count_difference(counts_more : dict[str,int], counts_less : dict[str,int]):
    diffs = {}
    for k,vm in counts_more.items():
        vl = counts_less[k]
        diff = vm-vl
        if diff != 0:
            diffs[k] = diff
    return diffs

def model_performance(counts : dict[str,float], overrides : list[str]):
    cost_sum = 0
    costs = {}
    for ov in overrides:
        split_ov = ov.split("=")
        inst = split_ov[0]
        cost = float(split_ov[1])
        costs[inst] = cost
    for inst,count in counts.items():
        if inst in costs:
            cost_sum += costs[inst]*count
        else:
            cost_sum += count

    return cost_sum

def main():

    args = parse_arguments()

    counts64 = run_program(qemu_bin=args.qemu_bin, program=args.program[0], sme_bits=args.sme_bits,
                           insn_path=args.qemu_insn, instructions=args.instructions, dvz_instructions=args.dvz_instructions,
                           method=args.method[0], N=64, roi_iterations=args.roi_iterations)

    counts128 = run_program(qemu_bin=args.qemu_bin, program=args.program[0], sme_bits=args.sme_bits,
                            insn_path=args.qemu_insn, instructions=args.instructions, dvz_instructions=args.dvz_instructions,
                            method=args.method[0], N=128, roi_iterations=args.roi_iterations)

    diffs = inst_count_difference(counts_more=counts128, counts_less=counts64)

    for inst,count in diffs.items():
        print(f"{inst} : {count/64} per N")

    cost = model_performance(counts=diffs, overrides=args.cost_overrides)

    print("============================================")
    print(f"instruction-cost performance model: {cost}")


if __name__ == "__main__":
    
    main()
