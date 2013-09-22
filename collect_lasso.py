#!/usr/bin/env python
import os
import subprocess
import re
import pickle

if __name__=="__main__":
    threads = ["1","2","4","8"]

    records = \
        {k1: {'solve': [], 
              'factor': [], 
              'total': []} 
         for k1 in threads}


    regexps = {
        'solve': re.compile(r'Total solve time is ([\d\.]+)'),
        'factor':re.compile(r'Factorization time: ([\d\.]+)'),
        'total': re.compile(r'Total factorize \+ solve time ([\d\.]+)')
    }


    def record_times(output, thread):        
        for regexp, record in zip(regexps.values(), records[thread].values()):
            m = regexp.search(output)
            record.append(float(m.group(1))) 


    for thread in threads:
        os.environ["OMP_NUM_THREADS"] = thread
        problem = "matlab/lasso_data"
        result = subprocess.check_output(["./bin/pdos", problem])
        record_times(result, thread)

    pickle.dump(records, open("lasso_timings.p", "wb"))
