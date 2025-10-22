from cProfile import label
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

change1 = []
layout1 = []
layout2 = []

ookami_single = [0.4937, 1.0369, 2.0992, 4.3977, 8.3753]
ookami_layout1 = [0.4792, 0.8258, 1.6060, 3.0800, 5.8961]
ookami_layout2 = [2.4013, 2.5875, 2.7481, 3.5678, 6.9959]

rand_ookami_single = [0.3234, 0.5733, 1.0729, 2.0409, 4.0218]
rand_ookami_layout1 = [0.3337, 0.5209, 1.0114, 2.9115, 5.5740]
rand_ookami_layout2 = [2.3329, 2.2689, 2.0733, 3.4251, 6.4954]

o_layout1=[]
o_layout2=[]

h_layout1=[]
h_layout2=[]

cal_time = []
comm_time = []
cal_time_2 = []
comm_time_2 = []

n = 16*16*16*16
FLOP = 2574 * 4 * n / 1e9

P_peak = 24*2.7*2*8*2*2    # GFlops/s
b_s = 2*128  # GByte/s

blen = [1, 2, 4, 8, 16]
iter = 20

gmres_1 = []
gmres_2 = []
oe_1 = []
oe_2 = []

oe1 = [10.5970]
oe2 = [14.9982]

def get_sum(arr1, arr2, arr3, arr4):
    output = []
    for i in range(len(blen)):
        cnt = arr1[i]
        cnt += arr2[i]
        cnt += arr3[i]
        cnt += arr4[i]
        output.append(cnt)
    return output

def get_avg(arr, output):
    for i in range(len(blen)):
        avg = 0
        for j in range(0, iter):
            avg += arr[i*iter + j]
        output.append(avg/iter)
    return output

def draw_obs(change1, change2, filename, typ, mb, platform):
    AI_OE = [443550 * blen[i] / (16 * (96318 * blen[i] + 20742)) for i in range(0, len(change1))]
    fig, ax = plt.subplots(figsize=(6, 4))
    if typ == 'dirac':
        AI = [2574*blen[i]/(16*(168*blen[i] + 114)) for i in range(0, len(change1))]

        iter = 20

        #plt.plot(AI, [mb * AI[i] * 0.095 for i in range(0, len(single))],
        #        label='Roofline', color='grey', linestyle='dotted')
    elif typ == 'gmres':
        AI = [426918 * blen[i] / (16 * (92592 * blen[i] + 12426)) for i in range(0, len(change1))]
        FLOP = (109*2574*n*4 + 146352*n*4 + 145566) / 1e9
        iter = 1
        plt.plot(AI, [mb * 0.21 * AI[i] for i in range(0, len(change1))],
                 label='Roofline', color='grey', linestyle='dotted')
    if platform == 'juwels':
        FLOP = 2574 * 2 * n / 1e9
        plt.plot(AI, [FLOP*blen[i]/change1[i]*iter for i in range(0, len(change1))],
                 label='JUWELS, Layout 1', color='#99d8c9', marker='o')
        plt.plot(AI, [FLOP*blen[i]/change2[i]*iter for i in range(0, len(change1))],
                 label='JUWELS, Layout 2', color='#2ca25f', marker='v')
    elif platform == 'ookami':
        FLOP = 2574 * 4 * n / 1e9
        plt.plot(AI, [FLOP * blen[i] / change1[i] * iter for i in range(0, len(change1))],
                 label='Ookami, Layout 1', color='#9ebcda', marker='s')
        plt.plot(AI, [FLOP * blen[i] / change2[i] * iter for i in range(0, len(change1))],
                 label='Ookami, Layout 2', color='#8856a7', marker='*')
    plt.xticks(AI, ['AI('+ str(blen[i]) + ')' + '\n=' + str(round(AI[i], 3)) for i in range(0, len(change1))],
               fontsize=10)
    plt.xlabel('Arithmetic Intensity')
    plt.ylabel('Performance (GFLOP/s)')
    #plt.legend()
    fig.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncols=2, mode="expand", borderaxespad=0.)
    ax.set_axisbelow(True)
    plt.grid(axis='y', linestyle = '--', linewidth = 0.5)
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

def draw_gmres(change1, change2, change3, change4, filename):
    AI_OE = [(445230 * blen[i] + 9156) / (16 * (96516 * blen[i] + 21582)) for i in range(0, len(change1))]
    fig, ax = plt.subplots(figsize=(6, 4))

    AI = [398178 * blen[i] / (16 * (48636 * blen[i] + 12426)) for i in range(0, len(change1))]
    theor_juwels_perf = [AI[i] * 155 for i in range(0, len(change1))]

    FLOP_NP = (109 * 2574 * n * 4 + 146352 * n * 4 + 145566) / 1e9/2
    FLOP_P = (109 * (2574 * n * 4 + 2 * 84 * n) + 146352 * n * 2 + 145566/2) / 1e9/2

    achieved_perf_1 = [FLOP_NP * blen[i] / change1[i] for i in range(0, len(change1))]
    achieved_perf_2 = [FLOP_NP * blen[i] / change2[i] for i in range(0, len(change1))]

    print([achieved_perf_1[i] / theor_juwels_perf[i] for i in range(0, len(change1))])
    print([achieved_perf_2[i] / theor_juwels_perf[i] for i in range(0, len(change1))])

    plt.plot(AI, [FLOP_NP * blen[i] / change1[i] for i in range(0, len(change1))],
             label='JUWELS, Layout 1', color='#117733', marker='d')
    plt.plot(AI, [FLOP_NP * blen[i] / change2[i] for i in range(0, len(change1))],
             label='JUWELS, Layout 2', color='#44AA99', marker='x')
    #plt.plot(AI_OE, [FLOP_P * blen[i] / change3[i] for i in range(0, len(change1))],
    #        label='Odd-even, Layout 1', color='#f768a1', marker='p')
    #plt.plot(AI_OE, [FLOP_P * blen[i] / change4[i] for i in range(0, len(change1))],
    #         label='Odd-even, Layout 2', color='#ae017e', marker='h')

    plt.xticks(AI, ['AI(' + str(blen[i]) + ')' + '\n=' + f"{AI[i]:.3f}" for i in range(0, len(change1))])
    plt.xlabel('Arithmetic Intensity')
    plt.ylabel('Performance (GFLOP/s)')
    fig.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncols=2, mode="expand", borderaxespad=0.)    #plt.ylim(5.5, 18.5)
    ax.set_axisbelow(True)
    plt.grid(axis='y', linestyle = '--', linewidth = 0.5)
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

def read_from_file(filename, ignore=False, layout1=None, layout2=None):
    if layout2 is None:
        layout2 = []
    if layout1 is None:
        layout1 = []
    file1 = open(filename, 'r')
    cal1 = []
    cal2 = []
    cal3 = []
    comm = []
    ecal1 = []
    ecal2 = []
    ecal3 = []
    ecomm = []
    Lines = file1.readlines()

    for line in Lines:
        if line.startswith('cal time:'):
            cal1.append(float(line.split(' ')[2]))
            cal2.append(float(line.split(' ')[3]))
            cal3.append(float(line.split(' ')[4]))
        elif line.startswith('comm time:'):
            comm.append(float(re.sub("[^\d\.]", "", line)))
        elif line.startswith('layout 2 cal time'):
            ecal1.append(float(line.split(' ')[4]))
            ecal2.append(float(line.split(' ')[5]))
            ecal3.append(float(line.split(' ')[6]))
        elif line.startswith('layout 2 comm time'):
            ecomm.append(float(re.sub("[^\d\.]", "", line)[1:]))
        elif line.startswith('No blocking'):
            change1.append(float(re.sub("[^\d\.]", "", line[20:])))
        elif line.startswith('blocking'):
            if 'per entry' in line:
                layout2.append(float(re.sub("[^\d\.]", "", line[20:])))
            else:
                layout1.append(float(re.sub("[^\d\.]", "", line[20:])))

    if not ignore:
        cal_time.append(get_avg(cal1, []))
        cal_time.append(get_avg(cal2, []))
        cal_time.append(get_avg(cal3, []))
        comm_time.append(get_avg(comm, []))
        cal_time_2.append(get_avg(ecal1, []))
        cal_time_2.append(get_avg(ecal2, []))
        cal_time_2.append(get_avg(ecal3, []))
        comm_time_2.append(get_avg(ecomm, []))


def read_mem(filename):
    file1 = open(filename, 'r')
    Lines = file1.readlines()
    acc = 0

    for line in Lines:
        if not line.startswith('PAPI'):
            counters = line.split()
            cacheline = int(counters[2]) + int(counters[4])
            cacheline /= int(counters[1])
            acc += cacheline

    print(acc*2*256)

def plot_column_chart(caltime, commtime, filename):
    bsize = (
        "1",
        "2",
        "4",
        "8",
        "16",
    )
    sumarr = get_sum(caltime[0], caltime[1], caltime[2], commtime[0])
    weight_counts = {
        "clover + apply pi": np.array([caltime[0][i]/sumarr[i] for i in range(0, len(caltime[0]))]),
        "apply U^H": np.array([caltime[1][i]/sumarr[i] for i in range(0, len(caltime[1]))]),
        "apply U": np.array([caltime[2][i]/sumarr[i] for i in range(0, len(caltime[2]))]),
        "MPI Wait": np.array([commtime[0][i] / sumarr[i] for i in range(0, len(caltime[0]))]),
    }
    fig, ax = plt.subplots()
    bottom = np.zeros(5)

    for boolean, weight_count in weight_counts.items():
        p = ax.bar(bsize, weight_count, 0.5, label=boolean, bottom=bottom)
        bottom += weight_count

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=5)
    plt.xlabel('Blocking Size')
    plt.ylabel('proportion of execution time (%)')
    plt.savefig(filename)
    plt.show()

def plot_bar(change1, change2, change3, change4, change5, change6, filename):
    iter = 200
    rhs = ("1", "2", "4", "8", "16")
    perf = {
        'JUWELS, Layout 1': [FLOP*blen[i]/change1[i]*iter for i in range(0, len(change1))],
        'JUWELS, Layout 2': [FLOP*blen[i]/change2[i]*iter for i in range(0, len(change2))],
        'Ookami, Layout 1': [FLOP*blen[i]/change3[i]*iter for i in range(0, len(change3))],
        'Ookami, Layout 2': [FLOP*blen[i]/change4[i]*iter for i in range(0, len(change4))],
        'HAICGU, Layout 1': [FLOP*blen[i]/change5[i]*iter for i in range(0, len(change5))],
        'HAICGU, Layout 2': [FLOP*blen[i]/change6[i]*iter for i in range(0, len(change6))],
    }

    AI = [2574 * blen[i] / (16 * (168 * blen[i] + 114)) for i in range(0, len(change1))]
    mb_ookami = 619
    mb_juwels = 155
    mb_haicgu = 218
    mb = [155, 155, 619, 619, 218, 218]
    p = [FLOP * blen[i] / change6[i] * iter for i in range(0, len(change1))]
    print([p[i]/AI[i]/mb_haicgu for i in range(0, len(change1))])

    color = ['#117733', '#44AA99', '#332288', '#88CCEE', '#AA4499', '#CC6677']

    x = np.arange(len(rhs))  # the label locations
    width = 0.15  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(figsize=(6, 4))

    for attribute, measurement in perf.items():
        offset = width * multiplier
        ax.bar(x + offset, measurement, width, label=attribute, color=color[multiplier])
        #ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Performance (GFLOP/s)')
    ax.set_xticks(x + width, rhs)
    ax.set_xlabel('Blocking size')
    plt.legend(bbox_to_anchor=(-0.02, 1.13, 1., .102), loc='lower left',
               ncols=3, mode="expand", borderaxespad=-1.2)
    ax.set_ylim(0, 120)

    ax.set_axisbelow(True)
    plt.grid(axis='y', linestyle = '--', linewidth = 0.5)
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()


def annotate(ax, label, x, y, xytext, color):
    if x - 0.57 < 0.1:
        xytext = (0, 10)
    if 0.919 - x < 0.1:
        xytext = (10, 0)
        if label == '9%':
            xytext = (10, -10)
        if label == '37%':
            xytext = (10, -5)
    if 0.715 - x < 0.1:
        if label == '54%':
            xytext = (0, -15)
        elif label == '55%':
            xytext = (0, 10)
    if 0.819 -x < 0.1:
        if label == '51%':
            xytext = (0, 8)
    if 0.883 - x < 0.1:
        if label == '47%':
            xytext = (0, -18)
        if label == '48%':
            xytext = (10, -2)

    ax.annotate(label, xy=(x,y),
                xytext=xytext, textcoords='offset points',
                fontsize=10,
                color=color)


def draw_roofline(change1, change2, change3, change4, filename):
    AI_OE = [443550 * blen[i] / (16 * (96318 * blen[i] + 20742)) for i in range(0, len(change1))]
    fig, ax = plt.subplots(figsize=(6, 4))
    AI = [2574*blen[i]/(16*(168*blen[i] + 114)) for i in range(0, len(change1))]
    iter = 20
    FLOP = 2574 * 2 * n / 1e9

    theor_juwels_perf = [AI[i] * 155 for i in range(0, len(change1))]
    theor_ookami_perf = [AI[i] * 619 for i in range(0, len(change1))]
    achieved_juwels_perf_1 = [FLOP*blen[i]/change1[i]*iter for i in range(0, len(change1))]
    achieved_juwels_perf_2 = [FLOP*blen[i]/change2[i]*iter for i in range(0, len(change2))]
    arch_eff_1 = [f"{achieved_juwels_perf_1[i]/theor_juwels_perf[i] * 100:.0f}%" for i in range(0, len(change1))]
    arch_eff_2 = [f"{achieved_juwels_perf_2[i]/theor_juwels_perf[i] * 100:.0f}%" for i in range(0, len(change1))]
    achieved_ookami_perf_1 = [FLOP*2*blen[i]/change3[i]*iter for i in range(0, len(change3))]
    achieved_ookami_perf_2 = [FLOP*2*blen[i]/change4[i]*iter for i in range(0, len(change4))]
    arch_eff_3 = [f"{achieved_ookami_perf_1[i]/theor_ookami_perf[i] * 100:.0f}%" for i in range(0, len(change1))]
    arch_eff_4 = [f"{achieved_ookami_perf_2[i]/theor_ookami_perf[i] * 100:.0f}%" for i in range(0, len(change1))]

    plt.plot(AI, [achieved_juwels_perf_1[i] for i in range(0, len(change1))],
             label='JUWELS, Layout 1', color='#44AA99', marker='o')
    plt.plot(AI, [achieved_juwels_perf_2[i] for i in range(0, len(change1))],
             label='JUWELS, Layout 2', color='#117733', marker='v')
    plt.plot(AI, [achieved_ookami_perf_1[i] for i in range(0, len(change1))],
             label='Ookami, Layout 1', color='#88CCEE', marker='s')
    plt.plot(AI, [achieved_ookami_perf_2[i] for i in range(0, len(change1))],
             label='Ookami, Layout 2', color='#332288', marker='*')

    plt.xticks(AI, ['AI('+ str(blen[i]) + ')' + '\n=' + f"{AI[i]:.3f}" for i in range(0, len(change1))]
               )
    plt.xlabel('Arithmetic Intensity')
    plt.ylabel('Performance (GFLOP/s)')

    # conditionally position labels
    for label, x, y in zip(arch_eff_1, AI, achieved_juwels_perf_1):
        annotate(ax, label, x, y, (0, -15), '#44AA99')

    for label, x, y in zip(arch_eff_2, AI, achieved_juwels_perf_2):
        annotate(ax, label, x, y, (0, 15), '#117733')

    for label, x, y in zip(arch_eff_3, AI, achieved_ookami_perf_1):
        annotate(ax, label, x, y, (0, 10), '#88CCEE')

    for label, x, y in zip(arch_eff_4, AI, achieved_ookami_perf_2):
        annotate(ax, label, x, y, (0, 15), '#332288')
    #plt.legend()
    fig.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncols=2, mode="expand", borderaxespad=0.)
    ax.set_axisbelow(True)
    ax.set_ylim(0, 89)
    ax.set_xlim(0.56, 0.95)

    plt.grid(axis='y', linestyle = '--', linewidth = 0.5)
    plt.tight_layout()

    plt.savefig(filename)
    plt.show()


def read_data(filename, change):
    file1 = open(filename, 'r')
    Lines = file1.readlines()

    for line in Lines:
        if line.startswith('| elapsed wall clock time'):
            change.append(float(re.sub("[^\d\.]", "", line)))

    file1.close()


rc('font',**{'family':'serif','serif':['Times New Roman']})
plt.rcParams.update({'font.size': 10})

# --------- choose one ---------

# For the roofline plot
read_from_file('./comm_info/juwels256', ignore=True, layout1=layout1, layout2=layout2)
read_from_file('./comm_info/ookami128', ignore=True, layout1=o_layout1, layout2=o_layout2)
draw_roofline(layout1, layout2, o_layout1, o_layout2, filename='roofline.png')

# For bar plot
#read_from_file('./comm_info/juwels16', ignore=True, layout1=layout1, layout2=layout2)
#read_from_file('./comm_info/ookami16', ignore=True, layout1=o_layout1, layout2=o_layout2)
#read_from_file('./comm_info/haicgu16', ignore=True, layout1=h_layout1, layout2=h_layout2)
#plot_bar(layout1, layout2, o_layout1, o_layout2, h_layout1, h_layout2, 'comparison3.png')

# For GMRES plot
#read_data('./comm_info/gmres_layout1', gmres_1)
#read_data('./comm_info/gmres_layout2', gmres_2)
#read_data('./comm_info/oe_layout1', oe_1)
#read_data('./comm_info/oe_layout2', oe_2)
#draw_gmres(gmres_1, gmres_2, oe_1, oe_2, filename='bgmres_model.png')


