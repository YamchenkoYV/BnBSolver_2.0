import subprocess as sp
import os
import parser as prs
import shutil
import pandas as pd

OUT_DIR = "result"
EXE_DIR = ".."

RUN_COUNT = 100

benchmark_list = ['Ackley 1 function',
'Ackley 2 function',
'Brown function',
'Chung Reynolds function',
'Exponential function',
'Griewank function',
'Powell Singular 2_6s function',
'Qing function',
'Quintic function',
'Rosenbrock function',
'Schumer Steiglitz function',
'Schaffer F6 function',
'Schwefel 2.22 function',
'Streched V Sine Wave function',
'Trid 6 function',
'Trid 10 function',
'Trigonometric 1 function',
'Trigonometric 2 function',
'W Wavy function',
'Whitley function',
'Xin-She Yang 2 function',
'Xin-She Yang 3 function',
'Xin-She Yang 4function',
'Zakharov function',
'Cluster2D1 function',
'Cluster2D2 function']


def run_test_parasol():
    out_dir = OUT_DIR + os.sep + "parasol"
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    bin = EXE_DIR + os.sep + 'paralsolver' + os.sep + 'parasol.exe'
    for bench in benchmark_list:
        df = pd.DataFrame([], columns=['time', 'tps', 'steps', 'fval'])
        exec_list = [bin, bench, "' '", '0.1', '1000000', '1', '0.5', '0.5', '1000', '10']
        print("Running: {}".format(' '.join(exec_list)))
        for i in range(RUN_COUNT):
            print("Start: {}".format(i))
            p = sp.Popen(exec_list, stdout=sp.PIPE, stderr=sp.PIPE)
            output = p.communicate()
            res = prs.parse_str(output[0].decode("utf-8"))
            df_res = pd.DataFrame(res).transpose()
            df = df.append(df_res, ignore_index=True)
        res_file_path = out_dir + os.sep + bench.replace(' ', '_') + '.csv'
        df.to_csv(res_file_path)


def run_test_bnbatomic():
    out_dir = OUT_DIR + os.sep + "bnbatomic"
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    bin = EXE_DIR + os.sep + 'bnbatomic' + os.sep + 'bnbatomic.exe'
    for bench in benchmark_list:
        df = pd.DataFrame([], columns=['time', 'tps', 'steps', 'fval'])
        exec_list = [bin, bench, "' '", '0.1', '1000000', '1', '1000']
        print("Running: {}".format(' '.join(exec_list)))
        for i in range(RUN_COUNT):
            print("Start: {}".format(i))
            p = sp.Popen(exec_list, stdout=sp.PIPE, stderr=sp.PIPE)
            output = p.communicate()
            res = prs.parse_str(output[0].decode("utf-8"))
            df_res = pd.DataFrame(res).transpose()
            df = df.append(df_res, ignore_index=True)
        res_file_path = out_dir + os.sep + bench.replace(' ', '_') + '.csv'
        df.to_csv(res_file_path)


def run_test_bnbomp():
    out_dir = OUT_DIR + os.sep + "bnbomp"
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    bin = EXE_DIR + os.sep + 'bnbomp' + os.sep + 'bnbomp.exe'
    for bench in benchmark_list:
        df = pd.DataFrame([], columns=['time', 'tps', 'steps', 'fval'])
        exec_list = [bin, bench, "' '", '0.1', '1000000', '1']
        print("Running: {}".format(' '.join(exec_list)))
        for i in range(RUN_COUNT):
            print("Start: {}".format(i))
            p = sp.Popen(exec_list, stdout=sp.PIPE, stderr=sp.PIPE)
            output = p.communicate()
            res = prs.parse_str(output[0].decode("utf-8"))
            df_res = pd.DataFrame(res).transpose()
            df = df.append(df_res, ignore_index=True)
        res_file_path = out_dir + os.sep + bench.replace(' ', '_') + '.csv'
        df.to_csv(res_file_path)

if __name__=="__main__":
    run_test_bnbomp()