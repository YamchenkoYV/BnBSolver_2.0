import re
import matplotlib.pyplot as plt
import statistics as st

from argparse import ArgumentParser
from argparse import FileType


def parse(file_path):
    task_block_list = file_path.read().split("*\n\n")

    task_dict = {}
    for task_block in task_block_list[:-1]:
        t_name = re.search(r"^desciption: (.+)", task_block, re.MULTILINE).group(1)
        t_time = re.search(r"^Time: (.+) microsecond", task_block, re.MULTILINE).group(1)
        t_tps = re.search(r"^Time per subproblem: (.+) miscroseconds", task_block, re.MULTILINE).group(1)

        # check if not coverged
        try:
            t_steps = re.search(r"^Converged in (\d+) steps", task_block, re.MULTILINE).group(1)
        except:
            t_steps = -1

        t_fval = re.search(r"^BnB found = ([.0-9e+-]+)", task_block, re.MULTILINE).group(1)

        if t_name in task_dict:
            print("Warn key repeat: {}".format(t_name))

        try:
            task_dict[t_name] = {'time': float(t_time), 'tps': float(t_tps),
                             'steps': int(t_steps),
                             'fval': float(t_fval) }
        except:
            print("")
    return task_dict


def parse_str(res_str):
    task_block_list = res_str.split("*\n\n")

    task_dict = {}
    for task_block in task_block_list[:-1]:
        t_name = re.search(r"^desciption: (.+)", task_block, re.MULTILINE).group(1)
        t_time = re.search(r"^Time: (.+) microsecond", task_block, re.MULTILINE).group(1)
        t_tps = re.search(r"^Time per subproblem: (.+) miscroseconds", task_block, re.MULTILINE).group(1)

        # check if not coverged
        try:
            t_steps = re.search(r"^Converged in (\d+) steps", task_block, re.MULTILINE).group(1)
        except:
            t_steps = -1

        t_fval = re.search(r"^BnB found = ([.0-9e+-]+)", task_block, re.MULTILINE).group(1)

        if t_name in task_dict:
            print("Warn key repeat: {}".format(t_name))

        try:
            task_dict[t_name] = {'time': float(t_time), 'tps': float(t_tps),
                             'steps': int(t_steps),
                             'fval': float(t_fval) }
        except:
            print("")
    return task_dict

def analyze_time(td_1, td_2):
    print("Time analyze")
    time_diff_list = []
    for f_name, res_dict in td_1.items():
        ct_1 = float(res_dict['time'])
        ct_2 = float(td_2[f_name]['time'])
        time_diff_list.append(ct_2 / ct_1)

    pos_accel = []
    neg_accel = []
    for time in time_diff_list:
        if time > 1:
            pos_accel.append(time)
        elif time < 1:
            neg_accel.append(1.0 / time)

    print("\tPos Accel: {}, mean: {}, median: {}".format(len(pos_accel),
                                                       st.mean(pos_accel),
                                                       st.median(pos_accel)))
    print("\tNeg Accel: {}, mean: {}, median: {}".format(len(neg_accel),
                                                       st.mean(neg_accel),
                                                       st.median(neg_accel)))

    '''
    plt.plot(range(len(time_diff_list)), time_diff_list)
    plt.plot(range(len(time_diff_list)), [1]*len(time_diff_list))
    plt.show()
    '''

def analyze_tps(td_1, td_2):
    print("TPS analyze")
    time_diff_list = []
    for f_name, res_dict in td_1.items():
        ct_1 = float(res_dict['tps'])
        ct_2 = float(td_2[f_name]['tps'])
        time_diff_list.append(ct_1 / ct_2)

    pos_accel = []
    neg_accel = []
    for time in time_diff_list:
        if time > 1:
            pos_accel.append(time)
        elif time < 1:
            neg_accel.append(1.0 / time)

    print("\tPos Accel: {}, mean: {}, median: {}".format(len(pos_accel),
                                                       st.mean(pos_accel),
                                                       st.median(pos_accel)))
    print("\tNeg Accel: {}, mean: {}, median: {}".format(len(neg_accel),
                                                       st.mean(neg_accel),
                                                       st.median(neg_accel)))

    '''
    plt.plot(range(len(time_diff_list)), time_diff_list)
    plt.plot(range(len(time_diff_list)), [1]*len(time_diff_list))
    plt.show()
    '''

def analyze_fval(td_1, td_2):
    print("Func value analyze")
    pos_accel = 0
    neg_accel = 0
    for f_name, res_dict in td_1.items():
        fv_1 = float(res_dict['fval'])
        fv_2 = float(td_2[f_name]['fval'])
        if fv_1 < fv_2:
            pos_accel += 1
        elif fv_1 > fv_2:
            neg_accel += 1


    print("\tPos Accel: {}".format(pos_accel))
    print("\tNeg Accel: {}".format(neg_accel))

if __name__ == "__main__":
    parser = ArgumentParser(description="Compare bnb results")
    parser.add_argument('f1', type=FileType("r"), help="first file to compare")
    parser.add_argument('f2', type=FileType("r"), help="second file to compare")

    args = parser.parse_args()

    f1_task_dict = parse(args.f1)
    f2_task_dict = parse(args.f2)

    analyze_time(f1_task_dict, f2_task_dict)
    analyze_tps(f1_task_dict, f2_task_dict)
    analyze_fval(f1_task_dict, f2_task_dict)


