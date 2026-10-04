#!/usr/bin/env python3
import numpy as np
# Simple reader for LAMMPS custom trajectory produced by in.md.
# Reports PTZ center-of-mass z and minimum PTZ-Pb / PTZ-I distances per frame.
def frames(path):
    with open(path) as f:
        while True:
            line=f.readline()
            if not line: return
            if line.strip()!="ITEM: TIMESTEP": continue
            step=int(f.readline())
            f.readline(); n=int(f.readline())
            f.readline(); bounds=[list(map(float,f.readline().split())) for _ in range(3)]
            f.readline()
            arr=[]
            for _ in range(n):
                p=f.readline().split()
                arr.append([int(p[0]),int(p[1]),int(p[2]),float(p[3]),float(p[4]),float(p[5]),float(p[6])])
            yield step,bounds,np.array(arr,float)
def main(path):
    for step,bounds,a in frames(path):
        ptz=a[a[:,1]==1]; pb=a[a[:,2]==2]; iod=a[a[:,2]==3]
        # COM using masses
        masses={4:12.011,5:14.007,6:32.06,7:1.008}
        w=np.array([masses[int(t)] for t in ptz[:,2]])
        com=np.average(ptz[:,4:7],axis=0,weights=w)
        dpb=np.linalg.norm(ptz[:,4:7,None,:]-pb[:,4:7][None,:,:],axis=2).min()
        di=np.linalg.norm(ptz[:,4:7,None,:]-iod[:,4:7][None,:,:],axis=2).min()
        print(step, *com, dpb, di)
if __name__=="__main__":
    import sys; main(sys.argv[1])
