import pdos_direct as p
import pdos_indirect as q
from sys import getrefcount
import cvxopt as o

#from guppy import hpy

A = o.spmatrix([1.,1.], [0,1],[0,1])
b = o.matrix([1.,1.])
c = o.matrix([1.,1.])

dims = {'f':2}
solution = p.solve(c,A,b,dims)
print solution['x']
print solution['y']
print solution['status']

print getrefcount(solution['x'])

# h = hpy()
# print h.heap()