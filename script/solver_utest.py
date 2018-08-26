import os
import unittest
import subprocess as sp
import parser as prs
import pandas as pd

EXE_DIR = ".."


class ParalSolverTester(unittest.TestCase):
    def setUp(self):
        self.run_count = 1000
        self.bin = EXE_DIR + os.sep + 'paralsolver' + os.sep + 'parasol.exe'
        self.bench = 'Cluster2D1 function'
        self.exec_list = [self.bin, self.bench, "'unknrec'", '0.1',
                          '1000000', '64', '0.5', '0.5', '1000', '10']

    def test_non_negative_sc(self):
        print("Running: {}".format(' '.join(self.exec_list)))
        for i in range(self.run_count):
            print("Start: {}".format(i))
            p = sp.Popen(self.exec_list, stdout=sp.PIPE, stderr=sp.PIPE)
            output = p.communicate()
            res = prs.parse_str(output[0].decode("utf-8"))
            steps = res[self.bench]['steps']
            self.assertGreater(steps, 0)

if __name__ == '__main__':
    unittest.main()