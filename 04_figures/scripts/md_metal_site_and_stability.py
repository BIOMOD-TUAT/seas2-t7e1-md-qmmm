import csv, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.chart.data import CategoryChartData, XyChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_MARKER_STYLE
from pptx.enum.text import PP_ALIGN
import openpyxl
from openpyxl.styles import Font, Alignment
C="../../02_md/analysis_csv/"
def rd(fn): return [r for r in csv.reader(l for l in open(C+fn,encoding="utf-8-sig") if not l.startswith("#")) if r]
JPF="/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"; fm.fontManager.addfont(JPF); JP=fm.FontProperties(fname=JPF).get_name()
plt.rcParams["font.family"]=["Liberation Sans",JP]
PAL=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300"]; INK2="#52514e"; GRID="#d9d8d3"
CO=["0A","1A","2A","3A","4A","5A"]
def sty(a):
    a.spines[["top","right"]].set_visible(False); a.yaxis.grid(True,color=GRID,lw=0.8); a.set_axisbelow(True); a.tick_params(labelsize=11,colors=INK2)
    for s in a.spines.values(): s.set_color(INK2)
def tbx(s,x,y,w,h,t,sz,bold=False):
    b=s.shapes.add_textbox(Inches(x),Inches(y),Inches(w),Inches(h)); tf=b.text_frame; tf.margin_left=tf.margin_right=tf.margin_top=tf.margin_bottom=0; tf.word_wrap=True
    r=tf.paragraphs[0].add_run(); r.text=t; r.font.size=Pt(sz); r.font.bold=bold; r.font.name="Arial"; r.font.color.rgb=RGBColor(0,0,0)
def axst(ch,xt,yt):
    for a,t in ((ch.category_axis,xt),(ch.value_axis,yt)):
        a.has_title=True; a.axis_title.text_frame.text=t
        for p in a.axis_title.text_frame.paragraphs:
            for r in p.runs: r.font.size=Pt(11); r.font.bold=False; r.font.name="Arial"
        a.tick_labels.font.size=Pt(10); a.tick_labels.font.name="Arial"
    ch.value_axis.major_gridlines.format.line.color.rgb=RGBColor(0xd9,0xd8,0xd3)
def fin(wb):
    for ws in wb:
        for row in ws.iter_rows():
            for c in row: c.font=Font(name="Arial",size=10); c.alignment=Alignment(wrap_text=True,vertical="top")
        for c in ws[1]: c.font=Font(name="Arial",size=10,bold=True)
# ================= Fig 7
FR={}
for c in CO:
    R=rd(f"Fig07_frames_{c}.csv"); hd=R[0]; FR[c]=(hd,np.array([[float(v) if v not in ("",) else np.nan for v in r[:len(hd)]] if True else 0 for r in [[x if x.replace('.','',1).replace('-','',1).isdigit() else "" for x in rr] for rr in R[1:]]]))
def col(c,key):
    hd,a=FR[c]; i=[k for k,h in enumerate(hd) if key(h)]; return a[:,i]
LIG=[("Asp55",["Asp55 OD1","Asp55 OD2"]),("Glu65",["Glu65 OE1","Glu65 OE2"]),("Thr66 O",["Thr66 O O"]),("Glu20",["Glu20 OE1","Glu20 OE2"]),("Glu37",["Glu37 OE1","Glu37 OE2"]),("リン酸の O",["C8 O1P","C8 O2P"])]
occ=np.zeros((12,6)); stat={}
for j,c in enumerate(CO):
    for si,site in enumerate(("活性中心1 サイト I）","活性中心1 サイト II）")):
        for li,(n,ats) in enumerate(LIG):
            m=col(c,lambda h: site in h and any(h.endswith("–"+a) for a in ats)); assert m.shape[1]==len(ats),(c,site,n,m.shape)
            occ[si*6+li,j]=100*np.mean(m.min(axis=1)<=0.26)
        p=col(c,lambda h: site in h and h.endswith("–C8 P"))[:,0]; stat[(c,si)]=np.percentile(p,[25,50,75])
    mm=col(c,lambda h: h=="活性中心 1 の Mn–Mn")[:,0]; stat[(c,"mm")]=np.percentile(mm,[25,50,75])
print(np.round(occ,1)); print({k:np.round(v,3) for k,v in stat.items() if k[0] in("3A","0A")})
fig=plt.figure(figsize=(13,4.6),dpi=300)
ax=fig.add_axes([0.055,0.17,0.25,0.66]); x=np.arange(6)
for si,(c_,lab) in enumerate(((PAL[0],"サイト I の Mn"),(PAL[1],"サイト II の Mn"))):
    q=np.array([stat[(c,si)] for c in CO]); ax.errorbar(x+(si-0.5)*0.22,q[:,1],yerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],fmt="o",ms=8,color=c_,capsize=0,lw=2,label=lab)
ax.set_xticks(x); ax.set_xticklabels(CO); ax.set_ylim(0,1.0); ax.set_ylabel("Mn とリンの距離 (nm)",fontsize=12); ax.set_xlabel("バルジの長さ",fontsize=12); sty(ax)
ax.legend(frameon=False,fontsize=11,loc="upper center",bbox_to_anchor=(0.5,1.2),ncol=2,handletextpad=0.2,columnspacing=1)
axb=fig.add_axes([0.42,0.17,0.25,0.66]); cm=LinearSegmentedColormap.from_list("b",["#ffffff","#cde2fb","#86b6ef","#2a78d6","#104281"])
axb.imshow(occ,cmap=cm,vmin=0,vmax=100,aspect="auto")
for i in range(12):
    for j in range(6): axb.text(j,i,f"{occ[i,j]:.0f}",ha="center",va="center",fontsize=9.5,color="white" if occ[i,j]>60 else "#0b0b0b")
axb.set_xticks(range(6)); axb.set_xticklabels(CO,fontsize=11); axb.set_yticks(range(12)); axb.set_yticklabels([n for n,_ in LIG]*2,fontsize=10.5)
axb.axhline(5.5,color="white",lw=4); axb.tick_params(length=0,colors=INK2)
for s in axb.spines.values(): s.set_visible(False)
fig.text(0.345,0.67,"サイト I",fontsize=11,rotation=90,va="center"); fig.text(0.345,0.335,"サイト II",fontsize=11,rotation=90,va="center")
axc=fig.add_axes([0.745,0.17,0.24,0.66]); q=np.array([stat[(c,"mm")] for c in CO])
axc.errorbar(x,q[:,1],yerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],fmt="o",ms=8,color=PAL[0],lw=2); axc.axhline(0.352,color=INK2,lw=1.2,ls="--"); axc.text(5.4,0.357,"結晶構造",fontsize=10.5,color=INK2,ha="right",va="bottom")
axc.set_xticks(x); axc.set_xticklabels(CO); axc.set_ylim(0.3,0.6); axc.set_ylabel("Mn どうしの距離 (nm)",fontsize=12); axc.set_xlabel("バルジの長さ",fontsize=12); sty(axc)
for X,l in ((0.01,"a"),(0.335,"b"),(0.695,"c")): fig.text(X,0.9,l,fontsize=20,fontweight="bold")
fig.savefig("Fig07_metal_site.png",facecolor="white"); fig.savefig("Fig07_metal_site.svg",facecolor="white"); fig.savefig("chk7.png",dpi=110,facecolor="white")
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5); s=prs.slides.add_slide(prs.slide_layouts[6])
cd=CategoryChartData(); cd.categories=CO
cd.add_series("サイト I の Mn",[float(stat[(c,0)][1]) for c in CO]); cd.add_series("サイト II の Mn",[float(stat[(c,1)][1]) for c in CO])
ch=s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS,Inches(0.3),Inches(1.6),Inches(4.0),Inches(4.4),cd).chart
ch.has_legend=True; ch.legend.position=XL_LEGEND_POSITION.TOP; ch.legend.include_in_layout=False; ch.legend.font.size=Pt(10)
for se,c_ in zip(ch.plots[0].series,PAL):
    se.format.line.fill.background(); se.marker.style=XL_MARKER_STYLE.CIRCLE; se.marker.size=9; se.marker.format.fill.solid(); se.marker.format.fill.fore_color.rgb=RGBColor.from_string(c_[1:].upper()); se.marker.format.line.fill.background(); se.smooth=False
axst(ch,"バルジの長さ","Mn とリンの距離 (nm)"); ch.value_axis.minimum_scale=0; ch.value_axis.maximum_scale=1.0; ch.value_axis.major_unit=0.2
tab=s.shapes.add_table(13,7,Inches(4.6),Inches(1.6),Inches(4.3),Inches(4.4)).table
def cell(r,c,t,fill=None,white=False,bold=False):
    ce=tab.cell(r,c); ce.text=t; p=ce.text_frame.paragraphs[0]; p.alignment=PP_ALIGN.CENTER if c else PP_ALIGN.LEFT
    for ru in p.runs: ru.font.size=Pt(9.5); ru.font.name="Arial"; ru.font.bold=bold; ru.font.color.rgb=RGBColor(255,255,255) if white else RGBColor(0,0,0)
    ce.fill.solid(); ce.fill.fore_color.rgb=RGBColor.from_string(fill or "FFFFFF"); ce.margin_top=ce.margin_bottom=Inches(0.02)
cell(0,0,"",bold=True)
for j,c in enumerate(CO): cell(0,j+1,c,bold=True)
for i in range(12):
    cell(i+1,0,("I  " if i<6 else "II  ")+LIG[i%6][0])
    for j in range(6):
        v=occ[i,j]; rgb=cm(v/100)[:3]; cell(i+1,j+1,f"{v:.0f}","%02X%02X%02X"%tuple(int(round(k*255)) for k in rgb),white=v>60)
cd=CategoryChartData(); cd.categories=CO; cd.add_series("Mn–Mn（中央値）",[float(stat[(c,"mm")][1]) for c in CO]); cd.add_series("結晶構造",[0.352]*6)
c3=s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS,Inches(9.2),Inches(1.6),Inches(3.9),Inches(4.4),cd).chart
c3.has_legend=True; c3.legend.position=XL_LEGEND_POSITION.TOP; c3.legend.include_in_layout=False; c3.legend.font.size=Pt(10)
s0,s1=c3.plots[0].series
s0.format.line.fill.background(); s0.marker.style=XL_MARKER_STYLE.CIRCLE; s0.marker.size=9; s0.marker.format.fill.solid(); s0.marker.format.fill.fore_color.rgb=RGBColor(0x2A,0x78,0xD6); s0.marker.format.line.fill.background()
s1.marker.style=XL_MARKER_STYLE.NONE; s1.format.line.color.rgb=RGBColor(0x52,0x51,0x4E); s1.format.line.width=Pt(1.25); s1.smooth=False
axst(c3,"バルジの長さ","Mn どうしの距離 (nm)"); c3.value_axis.minimum_scale=0.3; c3.value_axis.maximum_scale=0.6; c3.value_axis.major_unit=0.05
for X,l in ((0.3,"a"),(4.6,"b"),(9.2,"c")): tbx(s,X,1.1,0.5,0.45,l,24,True)

s2=prs.slides.add_slide(prs.slide_layouts[6]); s2.shapes.add_picture("Fig07_metal_site.png",Inches(0.2),Inches(1.4),Inches(12.9))
prs.save("Fig07_metal_site_editable.pptx")
wb=openpyxl.Workbook(); w=wb.active; w.title="README"
for r in [("項目","内容"),("図","Main Fig. 7（MD での金属イオンの配置：リン酸との距離、配位しているアミノ酸、Mn どうしの距離）"),("作成日","2026-11-13"),
("元データ","BIOMOD_Wiki_Fig03_07_08_20261113_v1 の 01_Analysis_CSV：Fig07_frames_0A〜5A.csv（フレームごとの距離）。各条件 0〜50 ns、501 フレーム"),
("処理","a と c：フレームごとの距離の中央値と四分位（Q1、Q3）。b：アミノ酸ごとに、その酸素（Asp は OD1・OD2、Glu は OE1・OE2、Thr66 は主鎖の O、リン酸は O1P・O2P）のどれかが Mn から 0.26 nm 以内にあったフレームの割合。値はこちらでフレームの表から計算した（集計表 Fig07_Metal_Coordination.csv は原子ごとの値）"),
("対象","切断を調べている側の活性中心（活性中心 1：Mn3 = サイト I、Mn4 = サイト II）。リンは C 鎖 C8 の P"),
("結晶構造の値","Mn–Mn 0.352 nm（図 c の点線）"),
("注意","Mn は結合を持たない単純なモデルで扱っており、配位の形は結晶構造と違うことがある。サイト II の Mn が Glu20・Glu37 に配位するのは結晶構造にはない配置。各条件は 1 本の計算")]: w.append(r)
w.column_dimensions["A"].width=16; w.column_dimensions["B"].width=115
wa=wb.create_sheet("a_Mn_P"); wa.append(["条件","サイト I：Q1 (nm)","サイト I：中央値【図】","サイト I：Q3","サイト II：Q1","サイト II：中央値【図】","サイト II：Q3"])
for c in CO: wa.append([c]+[round(float(v),3) for v in stat[(c,0)]]+[round(float(v),3) for v in stat[(c,1)]])
wb2=wb.create_sheet("b_coordination"); wb2.append(["サイト","配位子"]+CO)
for i in range(12): wb2.append(["サイト I" if i<6 else "サイト II",LIG[i%6][0]]+[round(float(v),1) for v in occ[i]])
wc=wb.create_sheet("c_Mn_Mn"); wc.append(["条件","Q1 (nm)","中央値【図】","Q3"])
for c in CO: wc.append([c]+[round(float(v),3) for v in stat[(c,"mm")]])
for ws in (wa,wb2,wc):
    for cc in "ABCDEFGH": ws.column_dimensions[cc].width=20
fin(wb); wb.save("Fig07_source_data.xlsx")
# ================= Fig 3
TS={}
for c in CO:
    R=rd(f"Fig03_timeseries_{c}.csv"); a=np.array([[float(v) for v in r] for r in R[1:]]); a[:,0]=(a[:,0]-a[0,0])/1000; TS[c]=(R[0],a)
fig=plt.figure(figsize=(13,4.3),dpi=300)
a1=fig.add_axes([0.06,0.2,0.41,0.62]); a2=fig.add_axes([0.565,0.2,0.41,0.62])
for c,cl in zip(CO,PAL):
    a=TS[c][1]; a1.plot(a[:,0],a[:,1],color=cl,lw=1.3,label=c); a2.plot(a[:,0],a[:,11],color=cl,lw=1.3)
a1.set_ylim(0,0.5); a1.set_ylabel("酵素の Cα の RMSD (nm)",fontsize=12); a2.set_ylim(0,1.0); a2.set_ylabel("二量体界面の接触の保持率",fontsize=12)
for a in (a1,a2): a.set_xlim(0,50); a.set_xlabel("時間 (ns)",fontsize=12); sty(a)
hd,lb=a1.get_legend_handles_labels(); fig.legend(hd,lb,loc="upper center",bbox_to_anchor=(0.52,1.0),ncol=6,frameon=False,fontsize=12,handlelength=1.4)
fig.text(0.012,0.86,"a",fontsize=20,fontweight="bold"); fig.text(0.515,0.86,"b",fontsize=20,fontweight="bold")
fig.savefig("Fig03_dimer_stability.png",facecolor="white"); fig.savefig("Fig03_dimer_stability.svg",facecolor="white"); fig.savefig("chk3.png",dpi=110,facecolor="white")
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5); s=prs.slides.add_slide(prs.slide_layouts[6])
for k,(ci,yt,ymax,mu,X) in enumerate(((1,"酵素の Cα の RMSD (nm)",0.5,0.1,0.4),(11,"二量体界面の接触の保持率",1.0,0.2,6.8))):
    xd=XyChartData()
    for c in CO:
        se=xd.add_series(c); a=TS[c][1]
        for t,v in zip(a[:,0],a[:,ci]): se.add_data_point(float(t),float(v))
    ch=s.shapes.add_chart(XL_CHART_TYPE.XY_SCATTER_LINES_NO_MARKERS,Inches(X),Inches(1.7),Inches(6.1),Inches(4.3),xd).chart
    ch.has_legend=(k==0)
    if k==0: ch.legend.position=XL_LEGEND_POSITION.TOP; ch.legend.include_in_layout=False; ch.legend.font.size=Pt(10)
    for se,cl in zip(ch.plots[0].series,PAL): se.format.line.color.rgb=RGBColor.from_string(cl[1:].upper()); se.format.line.width=Pt(1.25); se.smooth=False
    axst(ch,"時間 (ns)",yt); ch.value_axis.minimum_scale=0; ch.value_axis.maximum_scale=ymax; ch.value_axis.major_unit=mu; ch.category_axis.minimum_scale=0; ch.category_axis.maximum_scale=50; ch.category_axis.major_unit=10
    tbx(s,X,1.2,0.5,0.45,"ab"[k],24,True)
prs.save("Fig03_dimer_stability_editable.pptx")
wb=openpyxl.Workbook(); w=wb.active; w.title="README"
for r in [("項目","内容"),("図","Main Fig. 3（MD での二量体の安定性：酵素のずれと、二量体界面の接触の保持）"),("作成日","2026-11-13"),
("元データ","BIOMOD_Wiki_Fig03_07_08_20261113_v1 の 01_Analysis_CSV：Fig03_timeseries_0A〜5A.csv（時系列）、Fig03_Dimer_Interface.csv（平均値）。値は変えていない。時刻は各計算の解析の始点を 0 ns にそろえた（3A は 500 ns の計算の 0.4〜50.4 ns）"),
("a：Cα の RMSD","基準は解析の最初のフレーム。両方の鎖の Cα で重ね合わせて計算したずれ"),
("b：接触の保持率","最初のフレームで A 鎖と B 鎖の間に接触（重原子の距離 ≤ 0.45 nm）があった残基の組のうち、その時刻にも接触が残っている割合"),
("シート","summary：条件ごとの平均値。ts_0A〜ts_5A：時系列（図に使った列に【図】）"),
("注意","各条件は 1 本の計算で、50 ns の範囲での観察。二量体がこれより長い時間でも保たれるかは、この図からは言えない")]: w.append(r)
w.column_dimensions["A"].width=16; w.column_dimensions["B"].width=115
ws=wb.create_sheet("summary")
for r in rd("Fig03_Dimer_Interface.csv"): ws.append([ (float(v) if v.replace('.','',1).isdigit() else v) for v in r])
for c in CO:
    hd,a=TS[c]; wt=wb.create_sheet("ts_"+c); h2=list(hd); h2[0]="時間 (ns)【図の横軸】"; h2[1]+="【図 a】"; h2[11]+="【図 b】"; wt.append(h2)
    for row in a: wt.append([float(v) for v in row])
fin(wb); wb.save("Fig03_source_data.xlsx")
