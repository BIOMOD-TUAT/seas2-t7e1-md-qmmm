import numpy as np, pickle
from scipy.ndimage import map_coordinates
A="data/"   # put the elf cube files and the *_env.pdb files here
B=0.529177210903
F={"R":("R_elf_web.cube","R_env.pdb"),"TS":("TS_elf_web.cube","TS_env.pdb"),"P":("P_elf_web.cube","P_env.pdb")}
def cube(fn):
    L=open(fn).read().split('\n'); n=int(L[2].split()[0]); o=np.array(L[2].split()[1:4],float)
    N=[int(L[i].split()[0]) for i in (3,4,5)]; h=[float(L[3+i].split()[1+i]) for i in range(3)]
    d=np.array(' '.join(L[6+abs(n):]).split(),float).reshape(N)
    return o*B,np.array(h)*B,d
def pdb(fn):
    at=[]
    for l in open(fn,encoding='utf-8'):
        if l.startswith(('ATOM','HETATM')):
            at.append((l[12:16].strip(),l[76:78].strip() or l[12:16].strip()[0],np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]),float(l[54:60]),float(l[60:66])))
    return at
out={}
for k,(c,p) in F.items():
    o,h,d=cube(A+c); at=pdb(A+p)
    key={int(round(a[3])):a[2] for a in at if 0.5<a[3]<12.5}
    P,O3,OW=key[1],key[2],key[3]
    def val(x): return map_coordinates(d,((np.atleast_2d(x)-o)/h).T,order=1,mode='nearest')
    e1=(OW-O3); e1/=np.linalg.norm(e1); v=P-O3; e2=v-np.dot(v,e1)*e1; e2/=np.linalg.norm(e2); e2=-e2
    c0=(O3+OW+P)/3
    s=np.linspace(-3.5,3.5,141); X,Y=np.meshgrid(s,s)
    pts=c0+X[...,None]*e1+Y[...,None]*e2
    Z=val(pts.reshape(-1,3)).reshape(X.shape)
    n=np.cross(e1,e2)
    proj=[(a[1],np.dot(a[2]-c0,e1),np.dot(a[2]-c0,e2),np.dot(a[2]-c0,n),a[0],a[4]) for a in at if a[4]>0.5 and abs(np.dot(a[2]-c0,e1))<3.5 and abs(np.dot(a[2]-c0,e2))<3.5 and abs(np.dot(a[2]-c0,n))<0.6]
    mid=(val((P+O3)/2)[0],val((P+OW)/2)[0])
    print(k,d.shape,d.min(),d.max(),"dP-O3 %.3f dP-Ow %.3f"%(np.linalg.norm(P-O3),np.linalg.norm(P-OW)),"ELF mid P-O3 %.2f P-Ow %.2f"%mid)
    for q in proj: print("   ",q[4],q[0],"%.2f %.2f out %.2f"%q[1:4])
    out[k]=dict(s=s,Z=Z,proj=proj,mid=mid,dP3=np.linalg.norm(P-O3),dPw=np.linalg.norm(P-OW),
        key={n_:(np.dot(x-c0,e1),np.dot(x-c0,e2)) for n_,x in (("P",P),("O3",O3),("OW",OW))})
pickle.dump(out,open("elf.pkl","wb"))
