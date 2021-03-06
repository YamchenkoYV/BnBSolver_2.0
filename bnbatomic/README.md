# Multithreaded BnB interval solver

Coded by *M. Posypkin*

The program use C++ 17 threads to parallelize the branch-and-bound methods based on simple interval bounds.

**Usage:**

    bnbatomic.exe "name_of_bench" knrec|unknrec eps max_steps virtual_procs_number parallel_steps_limit

**or (to list available test functions):**
    
    bnbatomic.exe list

Generic options:

Parameter | Description
------------ | -------------
"name_of_bench" | Name of the benchmark
knrec | preset the value of the record equal to the known optimum
uknrec | don't preset the value of the record 
eps | Precision
max_steps | The maximal number of steps to perform

Algorithm specific options:

Parameter | Description
------------ | -------------
virtual_procs_number | number of virtual processors
parallel_steps_limit | the minimal number of steps to start parallelization


