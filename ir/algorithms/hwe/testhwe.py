import xacc

xacc.Initialize()

f = xacc.getIRGenerator('hwe').generate([2,4,'[[0,1],[1,2],[2,3]]'])
print(f.toString('q'))
print(f.nParameters())