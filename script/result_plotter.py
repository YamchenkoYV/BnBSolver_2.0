import os
import glob
import pandas as pd
import shutil

from matplotlib import pyplot as plt

solver = 'bnbomp'

OUT_DIR = 'plots/' + solver
PROC_PATH = 'result/' + solver

if __name__=="__main__":
    if os.path.exists(OUT_DIR):
        shutil.rmtree(OUT_DIR)
    os.makedirs(OUT_DIR)
    file_list = glob.glob(PROC_PATH + os.sep + "*.csv")

    for file in file_list:
        df = pd.read_csv(file)
        func_name = file.rsplit('/')[-1].split(".")[0]
        ax = df.boxplot( column='time', return_type='axes')
        descr = df['time'].describe()
        ax.set_title(func_name)
        ax.set_ylabel('Time, ms')
        ax.set_xlabel('Box Plot')
        annotat = 'min: {:.2f}\nmax: {:.2f}\nmean: {:.2f}\nstd: {:.2f}\n25%: {:.2f}\n75%: {:.2f}\n50%: {:.2f}'.format(
            descr['min'], descr['max'], descr['mean'], descr['std'], descr['25%'], descr['75%'], descr['50%'])
        props = dict(facecolor='white', alpha=1, edgecolor='black')
        plt.text(0.8, 0.8, annotat, ha='center', va='center', transform=ax.transAxes, bbox=props)
        plt.savefig(OUT_DIR + os.sep + func_name)
        plt.close()