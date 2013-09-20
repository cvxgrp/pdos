#!/usr/bin/env python
import os
import subprocess
import re
import pickle

if __name__=="__main__":
    threads = ["1","2","4","8"]
    problem_sizes = ["10x100", "30x1000", "100x10000", "300x100000"]

    records = \
        {k1:
            {k2: {'solve': [], 
                  'factor': [], 
                  'total': []} 
             for k2 in problem_sizes}
         for k1 in threads}


    regexps = {
        'solve': re.compile(r'Total solve time is ([\d\.]+)'),
        'factor':re.compile(r'Factorization time: ([\d\.]+)'),
        'total': re.compile(r'Total factorize \+ solve time ([\d\.]+)')
    }


    def record_times(output, prob, thread):        
        for regexp, record in zip(regexps.values(), records[thread][prob].values()):
            m = regexp.search(output)
            record.append(float(m.group(1))) 


    for thread in threads:
        os.environ["OMP_NUM_THREADS"] = thread
        for prob in problem_sizes:
            for i in xrange(10):
                problem = "data/portfolio_test_%02d_%s" % (i+1, prob)
                result = subprocess.check_output(["./bin/pdos", problem])
                record_times(result, prob, thread)

    pickle.dump(records, open("timings.p", "wb"))
