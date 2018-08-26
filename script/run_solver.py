import subprocess as sp
from argparse import ArgumentParser
import logging

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('bin', type=str, help="Path to exe file")
    parser.add_argument('task_name', type=str, help="Task name")
    args = parser.parse_args()

    exec_list = [args.bin, args.task_name, '0.1']
    logger.info('Running solver {} : {}'.format(args.bin, " ".join(exec_list)))
    p = sp.Popen(exec_list)
    p.wait()