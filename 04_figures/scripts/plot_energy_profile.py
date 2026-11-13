import csv, subprocess, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from PIL import Image, ImageChops
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.chart.data import XyChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_MARKER_STYLE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
A=""
SO="soffice.py"
def rows(fn):
    out=[r for r in csv.reader(l for l in open(fn,encoding="utf-8-sig") if not l.startswith("#")) if r]; return out[1:]
en={r[0]:(float(r[2]),float(r[3])) for r in rows("../../03_qmmm/results/energy_table.csv")}
ds={r[0]:[float(x) for x in r[2:]] for r in rows("../../03_qmmm/results/path_distances.csv")}
EARLY=[(1.249,7.36),(0.699,25.49),(0.044,45.55)]
keys=[k for k in en if k in ("R_opt3","rqB_TSC_km","bk040","bk050","TS","bk060","rqB_TSC_kp","bk070","sm080lg_free2") or k.startswith("ipB_")]
pts=sorted([en[k] for k in keys]+EARLY,key=lambda p:-p[0]); px=[p[0] for p in pts]; py=[p[1] for p in pts]
HL=[("Reactant",en["R_opt3"]),("Transition state",en["TS"]),("After cleavage",en["sm080lg_free2"])]
o3=sorted([k for k in ds if k in keys],key=lambda k:-ds[k][0]); x=[ds[k][0] for k in o3]
SER=[(2,"Breaking bond (P–O3′)","#2a78d6"),(1,"Forming bond (P–Ow)","#eb6834"),(6,"H to DNA oxygen (O3′)","#1baf7a"),(8,"H to relay-water oxygen","#6250d6")]
JPF="/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"; fm.fontManager.addfont(JPF); JP=fm.FontProperties(fname=JPF).get_name()
plt.rcParams["font.family"]=["Liberation Sans",JP]; INK2="#52514e"; GRID="#d9d8d3"
fig=plt.figure(figsize=(13,4.3),dpi=300)
ax=fig.add_axes([0.06,0.17,0.40,0.68]); ax.plot(px,py,color="#8a8984",lw=1.6,marker="o",ms=4.5,mfc="#8a8984",mec="white",mew=0.6,zorder=2)
for (t,(X,Y)),(dx,dy,ha) in zip(HL,((-0.02,5,"center"),(0.22,5.5,"center"),(-0.02,-8.5,"center"))):
    ax.plot(X,Y,"o",ms=9,color="#0b0b0b",mec="white",mew=1.2,zorder=4); ax.text(X+dx,Y+dy,t,ha=ha,va="center",fontsize=12)
ax.set_xlim(2.0,-1.0); ax.set_ylim(-4,75); ax.set_xlabel("Reaction coordinate (Å)",fontsize=12); ax.set_ylabel("Energy ΔE (kcal/mol)",fontsize=12)
bx=fig.add_axes([0.56,0.17,0.42,0.68])
for i,lab,c in SER: bx.plot(x,[ds[k][i] for k in o3],color=c,lw=2,marker="o",ms=4,label=lab)
bx.axvline(ds["TS"][0],color=INK2,lw=0.9,ls=(0,(3,3))); bx.set_xlim(2.0,-1.0); bx.set_ylim(0.8,3.6); bx.set_xlabel("Reaction coordinate (Å)",fontsize=12); bx.set_ylabel("Distance (Å)",fontsize=12)
bx.legend(frameon=False,fontsize=10.5,ncol=2,loc="upper center",bbox_to_anchor=(0.5,1.22),handlelength=1.6,columnspacing=1.2)
for a in (ax,bx):
    a.spines[["top","right"]].set_visible(False); a.yaxis.grid(True,color=GRID,lw=0.8); a.set_axisbelow(True); a.tick_params(labelsize=11,colors=INK2)
fig.text(0.01,0.9,"a",fontsize=20,fontweight="bold"); fig.text(0.51,0.9,"b",fontsize=20,fontweight="bold")
fig.savefig("ab.png",facecolor="white")
TAB=[("Item","Value","Meaning"),
("Activation energy ΔE‡","+57.6 kcal/mol","Height from the reactant to the top of the energy barrier (the transition state)"),
("Free-energy barrier ΔG‡ (310.15 K)","+52.7 kcal/mol","ΔE‡ plus corrections for atomic vibrations (zero-point energy, thermal energy and entropy)"),
("Number of imaginary frequencies","1 (65.5i cm⁻¹)","At the top of a barrier the energy decreases only along the reaction direction. Exactly one such direction is the condition for a transition state"),
("Largest remaining force","1.14 × 10⁻³ Ha/Bohr","Forces on the atoms vanish at the top of a barrier. This value shows how close to zero they are (slightly above the target of 1.0 × 10⁻³)"),
("Connection on both sides","Reactant ↔ after cleavage","When the structure was displaced slightly to either side and relaxed, it moved toward the reactant on one side and toward cleavage on the other"),
("Energy after cleavage","ΔE = +57.5 kcal/mol","Structure beyond the top of the barrier. It is only 0.05 kcal/mol below the top and is not a stable product")]
def table(s,x,y,w,h):
    t=s.shapes.add_table(len(TAB),3,Inches(x),Inches(y),Inches(w),Inches(h)).table; t._tbl.tblPr.set("bandRow","0")
    for j,cw in enumerate((0.27,0.19,0.54)): t.columns[j].width=Inches(w*cw)
    for i,r in enumerate(TAB):
        for j,v in enumerate(r):
            c=t.cell(i,j); c.text=v; c.vertical_anchor=MSO_ANCHOR.MIDDLE; c.fill.solid(); c.fill.fore_color.rgb=RGBColor(0xF3,0xF2,0xEF) if i==0 else RGBColor(255,255,255)
            for p in c.text_frame.paragraphs:
                for ru in p.runs: ru.font.size=Pt(11); ru.font.name="Arial"; ru.font.bold=(i==0 or j==1); ru.font.color.rgb=RGBColor(0,0,0)
def tb(s,x,y,t,sz=24):
    b=s.shapes.add_textbox(Inches(x),Inches(y),Inches(0.5),Inches(0.45)); r=b.text_frame.paragraphs[0].add_run(); r.text=t; r.font.size=Pt(sz); r.font.bold=True; r.font.name="Arial"
# table-only render
p=Presentation(); p.slide_width=Inches(13.333); p.slide_height=Inches(7.5); s=p.slides.add_slide(p.slide_layouts[6]); table(s,0.5,0.5,12.3,3.2); p.save("tab.pptx")
subprocess.run(["python3",SO,"--headless","--convert-to","pdf","tab.pptx"],capture_output=True); subprocess.run(["pdftoppm","-r","300","-png","tab.pdf","tab"])
im=Image.open("tab-1.png").convert("RGB"); bb=ImageChops.difference(im,Image.new("RGB",im.size,"white")).getbbox(); tabim=im.crop(bb)
ab=Image.open("ab.png").convert("RGB"); W=ab.size[0]; tw=int(W*0.9); tabim=tabim.resize((tw,round(tw*tabim.size[1]/tabim.size[0])),Image.LANCZOS)
cv=Image.new("RGB",(W,ab.size[1]+tabim.size[1]+160),"white"); cv.paste(ab,(0,0)); cv.paste(tabim,((W-tw)//2+40,ab.size[1]+110))
from PIL import ImageDraw, ImageFont
ImageDraw.Draw(cv).text((int(W*0.01),ab.size[1]+20),"c",font=ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",84),fill="black")
cv.save("../out/fig10_energy.png",dpi=(300,300)); cv.resize((1400,round(1400*cv.size[1]/cv.size[0]))).save("chk.png")
# editable pptx
p=Presentation(); p.slide_width=Inches(13.333); p.slide_height=Inches(7.5); s=p.slides.add_slide(p.slide_layouts[6])
def axes(ch,xt,yt,xr,yr,ymaj):
    for a,t in ((ch.category_axis,xt),(ch.value_axis,yt)):
        a.has_title=True; a.axis_title.text_frame.text=t
        for q in a.axis_title.text_frame.paragraphs:
            for r in q.runs: r.font.size=Pt(11); r.font.bold=False
        a.tick_labels.font.size=Pt(10)
    ch.category_axis.minimum_scale,ch.category_axis.maximum_scale=xr; ch.category_axis.major_unit=0.5; ch.category_axis.reverse_order=True; ch.category_axis.has_major_gridlines=False
    ch.value_axis.minimum_scale,ch.value_axis.maximum_scale=yr; ch.value_axis.major_unit=ymaj; ch.value_axis.major_gridlines.format.line.color.rgb=RGBColor(0xd9,0xd8,0xd3)
    from pptx.enum.chart import XL_AXIS_CROSSES
    ch.category_axis.crosses=XL_AXIS_CROSSES.MAXIMUM
xd=XyChartData(); se=xd.add_series("Path")
for a_,b_ in zip(px,py): se.add_data_point(a_,b_)
se=xd.add_series("Key structures")
for t,(X,Y) in HL: se.add_data_point(X,Y)
ch=s.shapes.add_chart(XL_CHART_TYPE.XY_SCATTER_LINES,Inches(0.4),Inches(0.6),Inches(6.0),Inches(3.4),xd).chart; ch.has_legend=False
s0,s1=ch.plots[0].series
s0.format.line.color.rgb=RGBColor(0x8a,0x89,0x84); s0.format.line.width=Pt(1.5); s0.smooth=False; s0.marker.style=XL_MARKER_STYLE.CIRCLE; s0.marker.size=5; s0.marker.format.fill.solid(); s0.marker.format.fill.fore_color.rgb=RGBColor(0x8a,0x89,0x84); s0.marker.format.line.fill.background()
s1.format.line.fill.background(); s1.marker.style=XL_MARKER_STYLE.CIRCLE; s1.marker.size=9; s1.marker.format.fill.solid(); s1.marker.format.fill.fore_color.rgb=RGBColor(0,0,0); s1.marker.format.line.fill.background()
axes(ch,"Reaction coordinate (Å)","Energy ΔE (kcal/mol)",(-1.0,2.0),(0,75),15)
xd=XyChartData()
for i,lab,c in SER:
    se=xd.add_series(lab)
    for a_,k in zip(x,o3): se.add_data_point(a_,ds[k][i])
c2=s.shapes.add_chart(XL_CHART_TYPE.XY_SCATTER_LINES,Inches(6.8),Inches(0.6),Inches(6.2),Inches(3.4),xd).chart
c2.has_legend=True; c2.legend.position=XL_LEGEND_POSITION.TOP; c2.legend.include_in_layout=False; c2.legend.font.size=Pt(9.5)
for se,(i,lab,c) in zip(c2.plots[0].series,SER):
    rgb=RGBColor.from_string(c[1:].upper()); se.format.line.color.rgb=rgb; se.format.line.width=Pt(2); se.smooth=False; se.marker.style=XL_MARKER_STYLE.CIRCLE; se.marker.size=5; se.marker.format.fill.solid(); se.marker.format.fill.fore_color.rgb=rgb; se.marker.format.line.fill.background()
axes(c2,"Reaction coordinate (Å)","Distance (Å)",(-1.0,2.0),(0.8,3.6),0.4)
tb(s,0.35,0.15,"a"); tb(s,6.75,0.15,"b"); tb(s,0.35,4.05,"c"); table(s,0.9,4.25,11.9,2.9)
for t,(X,Y),(dx,dy) in zip([h[0] for h in HL],[h[1] for h in HL],((0.2,-0.42),(-0.55,-0.5),(0.05,0.22))):
    fx=0.4+0.75+(2.0-X)/3.0*4.9; fy=0.6+0.15+(75-Y)/75*2.55
    b=s.shapes.add_textbox(Inches(fx+dx-0.5),Inches(fy+dy),Inches(1.0),Inches(0.3)); b.text_frame.paragraphs[0].alignment=PP_ALIGN.CENTER; r=b.text_frame.paragraphs[0].add_run(); r.text=t; r.font.size=Pt(11); r.font.name="Arial"
p.save("Fig13_QMMM_energy_editable.pptx")
