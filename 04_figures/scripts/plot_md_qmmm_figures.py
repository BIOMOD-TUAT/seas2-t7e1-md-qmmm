import os, sys, subprocess, pickle, csv, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont, ImageChops
A="img/"; O="out/"
C="../../02_md/analysis_csv/"
plt.rcParams["font.family"]="Liberation Sans"
PAL=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300"]; INK2="#52514e"; GRID="#d9d8d3"
FB="/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"; FR="/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
cm=LinearSegmentedColormap.from_list("b",["#ffffff","#cde2fb","#86b6ef","#2a78d6","#104281"])
def rd(fn): return [r for r in csv.reader(l for l in open(C+fn,encoding="utf-8-sig") if not l.startswith("#")) if r]
def sty(a):
    a.spines[["top","right"]].set_visible(False); a.yaxis.grid(True,color=GRID,lw=0.8); a.set_axisbelow(True); a.tick_params(labelsize=11,colors=INK2)
    for s in a.spines.values(): s.set_color(INK2)
CO=["0A","1A","2A","3A","4A","5A"]
# ---------- fig03
fig=plt.figure(figsize=(13,4.3),dpi=300); axs=[fig.add_axes([0.05+i*0.335,0.17,0.265,0.6]) for i in range(3)]
for c,cl in zip(CO,PAL):
    R=rd(f"Fig03_timeseries_{c}.csv"); a=np.array([[float(v) for v in r] for r in R[1:]]); t=(a[:,0]-a[0,0])/1000
    for ax,ci in zip(axs,(1,2,11)): ax.plot(t,a[:,ci],color=cl,lw=1.2,label=c)
for ax,yl,ym in zip(axs,("Protein RMSD (nm)","DNA RMSD (nm)","Dimer-interface contacts retained"),(0.5,1.0,1.0)):
    ax.set_xlim(0,50); ax.set_ylim(0,ym); ax.set_xlabel("Time (ns)",fontsize=12); ax.set_ylabel(yl,fontsize=12); sty(ax)
for ax in axs[:2]: ax.legend(frameon=False,fontsize=9.5,ncol=6,loc="upper center",bbox_to_anchor=(0.5,1.17),handlelength=1.2,columnspacing=0.8,handletextpad=0.3)
for i,l in enumerate("abc"): fig.text(0.008+i*0.335,0.9,l,fontsize=20,fontweight="bold")
fig.savefig(O+"fig03_structural_stability.png",facecolor="white"); plt.close(fig)
# ---------- fig04
RES=["Lys67","Gly68","Leu69","Glu71","Lys103"]
ANY=np.array([[np.nan,0.2,0.6,10.4,20.8,3.8],[np.nan,95.4,46.7,99.8,92.6,70.9],[np.nan,95.6,50.3,100.0,99.8,85.8],[np.nan,0.0,73.0,49.9,91.0,54.5],[np.nan,1.2,0.0,99.0,88.6,32.1]])
im=Image.open("img/contacts_500ns.png").convert("RGB")
fig=plt.figure(figsize=(13,4.2),dpi=300); ax=fig.add_axes([0.085,0.1,0.36,0.62]); cmg=cm.copy(); cmg.set_bad("#ededed")
ax.imshow(np.ma.masked_invalid(ANY),cmap=cmg,vmin=0,vmax=100,aspect="auto")
for i in range(5):
    for j in range(6):
        v=ANY[i,j]; ax.text(j,i,"–" if np.isnan(v) else f"{v:.0f}",ha="center",va="center",fontsize=13,color="white" if (not np.isnan(v) and v>60) else "#0b0b0b")
ax.set_xticks(range(6)); ax.set_xticklabels(CO,fontsize=12); ax.xaxis.tick_top(); ax.set_yticks(range(5)); ax.set_yticklabels(RES,fontsize=12); ax.tick_params(length=0)
ax.set_xticks(np.arange(-.5,6),minor=True); ax.set_yticks(np.arange(-.5,5),minor=True); ax.grid(which="minor",color="white",lw=2); ax.tick_params(which="minor",length=0)
for s in ax.spines.values(): s.set_visible(False)
fig.text(0.265,0.86,"Bulge length",fontsize=13,ha="center",color=INK2)
w=0.49; h=w*13/4.2*im.size[1]/im.size[0]; ab=fig.add_axes([0.505,0.48-h/2,w,h]); ab.imshow(im); ab.axis("off")
fig.text(0.01,0.9,"a",fontsize=20,fontweight="bold"); fig.text(0.475,0.9,"b",fontsize=20,fontweight="bold")
fig.savefig(O+"fig04_bulge_contacts.png",facecolor="white"); plt.close(fig)
# ---------- fig05 (reuse data part)
src=open("md_metal_site_and_stability.py").read(); cut=src.index("fig=plt.figure(figsize=(13,4.6),dpi=300)")
g={}; exec(src[:cut].replace('print(np.round(occ,1));','#').replace('plt.rcParams["font.family"]=["Liberation Sans",JP]','plt.rcParams["font.family"]="Liberation Sans"'),g); occ=g["occ"]; stat=g["stat"]
ROWS=[(0,"Asp55"),(1,"Glu65"),(2,"Thr66 (backbone O)"),(5,"Phosphate O"),(6,"Asp55"),(9,"Glu20"),(10,"Glu37"),(11,"Phosphate O")]
fig=plt.figure(figsize=(13,4.6),dpi=300); ax=fig.add_axes([0.06,0.17,0.38,0.66]); x=np.arange(6)
for si,(c_,lab) in enumerate(((PAL[0],"Mn①"),(PAL[1],"Mn②"))):
    q=np.array([stat[(c,si)] for c in CO]); ax.errorbar(x+(si-0.5)*0.2,q[:,1],yerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],fmt="o",ms=8,color=c_,lw=2,label=lab)
ax.set_xticks(x); ax.set_xticklabels(CO); ax.set_ylim(0,1.0); ax.set_ylabel("Mn–P distance (nm)",fontsize=12); ax.set_xlabel("Bulge length",fontsize=12); sty(ax)
plt.rcParams["font.family"]=["Liberation Sans","Noto Sans CJK JP"]
ax.legend(frameon=False,fontsize=12,loc="upper center",bbox_to_anchor=(0.5,1.18),ncol=2,handletextpad=0.2)
axb=fig.add_axes([0.68,0.1,0.30,0.72]); M=np.array([occ[i] for i,_ in ROWS]); axb.imshow(M,cmap=cm,vmin=0,vmax=100,aspect="auto")
for i in range(8):
    for j in range(6): axb.text(j,i,f"{M[i,j]:.0f}",ha="center",va="center",fontsize=11,color="white" if M[i,j]>60 else "#0b0b0b")
axb.set_xticks(range(6)); axb.set_xticklabels(CO,fontsize=11); axb.xaxis.tick_top(); axb.set_yticks(range(8)); axb.set_yticklabels([n for _,n in ROWS],fontsize=11); axb.tick_params(length=0,colors="#0b0b0b")
axb.set_xticks(np.arange(-.5,6),minor=True); axb.set_yticks(np.arange(-.5,8),minor=True); axb.grid(which="minor",color="white",lw=1.5); axb.tick_params(which="minor",length=0); axb.axhline(3.5,color="white",lw=5)
for s_ in axb.spines.values(): s_.set_visible(False)
fig.text(0.525,0.64,"Mn①",fontsize=12,va="center",color=PAL[0],fontweight="bold"); fig.text(0.525,0.28,"Mn②",fontsize=12,va="center",color=PAL[1],fontweight="bold")
fig.text(0.01,0.9,"a",fontsize=20,fontweight="bold"); fig.text(0.50,0.9,"b",fontsize=20,fontweight="bold")
fig.savefig(O+"fig05_metal_site.png",facecolor="white"); plt.close(fig); plt.rcParams["font.family"]="Liberation Sans"
# ---------- fig06
N=rd("Fig08_NAC_like.csv")[1:]; f=lambda v: None if v=="N/A" else float(v)
SER=[[f(r[k]) for r in N] for k in (3,4,8)]; T=rd("Fig08_3A_500ns_blocks_10ns.csv")[1:]; tm=[(float(r[0].split("〜")[0])+float(r[0].split("〜")[1]))/2 for r in T]
TSER=[[float(r[k]) for r in T] for k in (3,4,8)]; LAB=["Water in attack position","…and bound to Mn","…and proton relay in place"]; P3=PAL[:3]
fig=plt.figure(figsize=(13,4.3),dpi=300); ax=fig.add_axes([0.06,0.2,0.36,0.62]); x=np.arange(6); w=0.26
for i,(v,c) in enumerate(zip(SER,P3)):
    for j,val in enumerate(v):
        if val is None: ax.text(x[j]+(i-1)*w,1.5,"–",ha="center",va="bottom",fontsize=11,color=INK2)
        else: ax.bar(x[j]+(i-1)*w,val,w*0.9,color=c,label=LAB[i] if j==2 else None)
ax.set_xticks(x); ax.set_xticklabels(CO); ax.set_ylim(0,100); ax.set_ylabel("Occupancy (%)",fontsize=12); ax.set_xlabel("Bulge length",fontsize=12); sty(ax)
ax2=fig.add_axes([0.5,0.2,0.48,0.62])
for v,c in zip(TSER,P3): ax2.plot(tm,v,color=c,lw=2)
ax2.set_xlim(0,500); ax2.set_ylim(0,100); ax2.set_xlabel("Time (ns)",fontsize=12); ax2.set_ylabel("Occupancy (%)",fontsize=12); sty(ax2)
hd,lb=ax.get_legend_handles_labels(); fig.legend(hd,lb,loc="upper center",bbox_to_anchor=(0.52,1.0),ncol=3,frameon=False,fontsize=11.5,handlelength=1.2)
fig.text(0.012,0.86,"a",fontsize=20,fontweight="bold"); fig.text(0.455,0.86,"b",fontsize=20,fontweight="bold")
fig.savefig(O+"fig06_nucleophilic_water.png",facecolor="white"); plt.close(fig)
# ---------- fig07 (PIL)
wide=Image.open(A+"complex_wide.png").convert("RGB"); qm=Image.open(A+"qm_region.png").convert("RGB")
a=np.array(wide).astype(int); r,gc,b=a[...,0],a[...,1],a[...,2]
mask=((r>200)&(gc<70)&(b<70))|((r>230)&(gc>230)&(b>230))|((r>120)&(b>120)&(gc<60))
ys,xs=np.where(mask); cx,cy=np.median(xs),np.median(ys); sel=(abs(xs-cx)<220)&(abs(ys-cy)<220); x0,x1,y0,y1=xs[sel].min()-25,xs[sel].max()+25,ys[sel].min()-25,ys[sel].max()+25
RED=(255,0,0); LW=14
ah=1250; aw=round(ah*wide.size[0]/wide.size[1]); sc=ah/wide.size[1]; bh=2100; bw=round(bh*qm.size[0]/qm.size[1])
W=120+aw+420+bw+80; H=bh+420; cv=Image.new("RGB",(W,H),"white"); d=ImageDraw.Draw(cv)
ax_=120; ay_=160+(bh-ah)//2+80; bx_=120+aw+420; by_=170
cv.paste(wide.resize((aw,ah),Image.LANCZOS),(ax_,ay_)); cv.paste(qm.resize((bw,bh),Image.LANCZOS),(bx_,by_))
bx0,by0,bx1,by1=ax_+x0*sc,ay_+y0*sc,ax_+x1*sc,ay_+y1*sc
d.line([(bx1,by0),(bx_,by_)],fill=RED,width=LW); d.line([(bx1,by1),(bx_,by_+bh)],fill=RED,width=LW)
d.rectangle([bx0,by0,bx1,by1],outline=RED,width=LW); d.rectangle([bx_,by_,bx_+bw,by_+bh],outline=RED,width=LW+6)
fl=ImageFont.truetype(FB,110); fs=ImageFont.truetype(FR,84)
d.text((ax_,ay_-150),"a",font=fl,fill="black"); d.text((bx_,by_-150),"b",font=fl,fill="black")
d.text((ax_+aw/2,ay_+ah+40),"Whole enzyme–DNA complex",font=fs,fill="black",anchor="ma"); d.text((bx_+bw/2,by_+bh+40),"Region treated by quantum mechanics (QM region)",font=fs,fill="black",anchor="ma")
cv.save(O+"fig07_qmmm_model.png",dpi=(300,300))
# ---------- fig09
ims=[("a","Reactant",A+"structure_R.png"),("b","Transition state",A+"structure_TS.png"),("c","After cleavage",A+"structure_P.png")]
pw=1500; ph=round(pw*1582/2400); gap=50; m=40; top=110; W=m*2+pw*3+gap*2; H=top+ph+130; cv=Image.new("RGB",(W,H),"white"); d=ImageDraw.Draw(cv); fl=ImageFont.truetype(FB,84); fs=ImageFont.truetype(FR,64)
for i,(l,s,f_) in enumerate(ims):
    x=m+i*(pw+gap); cv.paste(Image.open(f_).convert("RGB").resize((pw,ph),Image.LANCZOS),(x,top)); d.text((x,8),l,font=fl,fill="black"); d.text((x+pw/2,top+ph+25),s,font=fs,fill="black",anchor="ma")
cv.save(O+"fig09_structures.png",dpi=(300,300))
# ---------- fig11
D=pickle.load(open("elf.pkl","rb")); L=3.0; YH=2.2; sp=Image.open("img/spin_crop.png")
def draw(ax,k):
    dd=D[k]; ax.imshow(dd["Z"],extent=[-3.5,3.5,-3.5,3.5],origin="lower",cmap="Blues",vmin=0,vmax=1,interpolation="bilinear"); ax.contour(dd["s"],dd["s"],dd["Z"],levels=[0.5],colors="#0b3a6e",linewidths=0.9)
    ax.set_xlim(-L,L); ax.set_ylim(-YH,YH); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    for n,t in (("O3","O"),("P","P"),("OW","O")):
        x,y=dd["key"][n]; ax.plot(x,y,"o",ms=5,mfc="white",mec="black",mew=1.0); ax.text(x,y-0.55,t,ha="center",va="center",fontsize=13,fontweight="bold",bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none",alpha=0.75))
FW,FH=13.0,7.6; pw_=0.285; ph_=pw_*FW/(2*L/(2*YH))/FH; fig=plt.figure(figsize=(FW,FH),dpi=300)
ax=fig.add_axes([0.5-pw_/2,0.565,pw_,ph_]); ax.imshow(sp); ax.axis("off"); fig.text(0.5-pw_/2-0.03,0.565+ph_-0.03,"a",fontsize=20,fontweight="bold"); fig.text(0.5,0.535,"Spin density (reactant)",fontsize=13,ha="center")
for i,(k,t) in enumerate((("R","ELF (reactant)"),("TS","ELF (transition state)"),("P","ELF (after cleavage)"))):
    x=0.045+i*(pw_+0.02); a_=fig.add_axes([x,0.09,pw_,ph_]); draw(a_,k); fig.text(x+pw_/2,0.055,t,fontsize=13,ha="center")
fig.text(0.012,0.09+ph_-0.03,"b",fontsize=20,fontweight="bold")
cax=fig.add_axes([0.955,0.13,0.012,0.26]); cb=fig.colorbar(plt.cm.ScalarMappable(cmap="Blues",norm=plt.Normalize(0,1)),cax=cax,ticks=[0,0.5,1]); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=11); cb.ax.set_title("ELF",fontsize=12,pad=6)
fig.savefig(O+"fig11_electronic_structure.png",facecolor="white"); plt.close(fig)
print("ok")
