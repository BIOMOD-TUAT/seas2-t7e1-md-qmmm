import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_MARKER_STYLE
import openpyxl
from openpyxl.styles import Font, Alignment
B="./"
RES=["Lys67","Gly68","Leu69","Glu71","Lys103"]; CO=["1A","2A","3A","4A","5A"]
ANY=np.array([[0.2,0.6,10.4,20.8,3.8],[95.4,46.7,99.8,92.6,70.9],[95.6,50.3,100.0,99.8,85.8],[0.0,73.0,49.9,91.0,54.5],[1.2,0.0,99.0,88.6,32.1]])
PAL=["#eb6834","#1baf7a","#eda100","#e87ba4","#008300"]; INK2="#52514e"; GRID="#d9d8d3"
JPF="/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"; fm.fontManager.addfont(JPF); JP=fm.FontProperties(fname=JPF).get_name()
plt.rcParams["font.family"]=["Liberation Sans",JP]; im=Image.open("img/contacts_500ns.png").convert("RGB")
TX={"ja":("アミノ酸残基","接触していた時間の割合 (%)"),"en":("Amino acid residue","Contact occupancy (%)")}
for lang in ("ja","en"):
    fig=plt.figure(figsize=(13,4.4),dpi=300); ax=fig.add_axes([0.06,0.17,0.36,0.63]); x=np.arange(5)
    for j,(c,cl) in enumerate(zip(CO,PAL)): ax.plot(x,ANY[:,j],color=cl,lw=2,marker="o",ms=7,label=c)
    ax.set_xticks(x); ax.set_xticklabels(RES); ax.set_xlim(-0.25,4.25); ax.set_ylim(-3,105); ax.set_yticks([0,25,50,75,100]); ax.set_xlabel(TX[lang][0],fontsize=12); ax.set_ylabel(TX[lang][1],fontsize=12)
    ax.spines[["top","right"]].set_visible(False); ax.yaxis.grid(True,color=GRID,lw=0.8); ax.set_axisbelow(True); ax.tick_params(labelsize=11,colors=INK2)
    ax.legend(frameon=False,fontsize=11,ncol=5,loc="upper center",bbox_to_anchor=(0.5,1.2),handlelength=1.3,columnspacing=1.2,handletextpad=0.4)
    w=0.49; h=w*13/4.4*im.size[1]/im.size[0]; ab=fig.add_axes([0.505,0.5-h/2,w,h]); ab.imshow(im); ab.axis("off")
    fig.text(0.01,0.9,"a",fontsize=20,fontweight="bold"); fig.text(0.475,0.9,"b",fontsize=20,fontweight="bold")
    fig.savefig(f"fig04_bulge_contacts_{lang}.png",facecolor="white"); plt.close(fig)
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5); s=prs.slides.add_slide(prs.slide_layouts[6])
cd=CategoryChartData(); cd.categories=RES
for j,c in enumerate(CO): cd.add_series(c,[float(v) for v in ANY[:,j]])
ch=s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS,Inches(0.4),Inches(1.7),Inches(5.9),Inches(4.2),cd).chart
ch.has_legend=True; ch.legend.position=XL_LEGEND_POSITION.TOP; ch.legend.include_in_layout=False; ch.legend.font.size=Pt(10)
for se,c in zip(ch.plots[0].series,PAL):
    rgb=RGBColor.from_string(c[1:].upper()); se.format.line.color.rgb=rgb; se.format.line.width=Pt(2); se.smooth=False; se.marker.style=XL_MARKER_STYLE.CIRCLE; se.marker.size=8; se.marker.format.fill.solid(); se.marker.format.fill.fore_color.rgb=rgb; se.marker.format.line.fill.background()
for a,t in ((ch.category_axis,"アミノ酸残基"),(ch.value_axis,"接触していた時間の割合 (%)")):
    a.has_title=True; a.axis_title.text_frame.text=t
    for p in a.axis_title.text_frame.paragraphs:
        for r in p.runs: r.font.size=Pt(11); r.font.bold=False
    a.tick_labels.font.size=Pt(10)
ch.value_axis.minimum_scale=0; ch.value_axis.maximum_scale=100; ch.value_axis.major_unit=25; ch.value_axis.major_gridlines.format.line.color.rgb=RGBColor(0xd9,0xd8,0xd3)
s.shapes.add_picture("img/contacts_500ns.png",Inches(6.7),Inches(2.0),Inches(6.3))
for X,l in ((0.35,"a"),(6.65,"b")):
    b=s.shapes.add_textbox(Inches(X),Inches(1.2),Inches(0.5),Inches(0.45)); r=b.text_frame.paragraphs[0].add_run(); r.text=l; r.font.size=Pt(24); r.font.bold=True; r.font.name="Arial"
prs.save("fig04_bulge_contacts_editable.pptx")
wb=openpyxl.Workbook(); w=wb.active; w.title="README"
for r in [("項目","内容"),("図","挿入した A と酵素の接触（a：アミノ酸ごとの接触の割合を、バルジの長さ別に示した。b：3A での 500 ns の時間変化）"),
("a","横軸は T7E1 の B 鎖の 5 つのアミノ酸、縦軸は挿入した A のどれかと接触していたフレームの割合（%、50 ns）。線はバルジの長さ（1A〜5A）。接触の定義：水素以外の原子間の最短距離が 0.45 nm 未満"),
("0A がない理由","0A には挿入した A がなく、「挿入した A との接触」が定義できないため"),
("a の線について","横軸のアミノ酸は連続した量ではない。線は同じ条件の点を見分けやすくするためにつないだもので、中間の値を意味しない"),
("a の数値の出どころ","元の MD 解析の図版に印字されていた数値（小数第 1 位まで）を、目で読み取って写した。元の表（CSV）とは照合していない"),
("b","3A の 500 ns の MD を 50 ns ごとに区切った、同じ接触の割合。B68 = Gly68、B69 = Leu69、B71 = Glu71、B103 = Lys103。元の MD 解析の図版から切り出した画像")]: w.append(r)
w.column_dimensions["A"].width=20; w.column_dimensions["B"].width=115
ws=wb.create_sheet("a_contacts"); ws.append(["残基（B 鎖）"]+CO)
for r,v in zip(RES,ANY): ws.append([r]+[float(x) for x in v])
for q in wb:
    for row in q.iter_rows():
        for c in row: c.font=Font(name="Arial",size=10); c.alignment=Alignment(wrap_text=True,vertical="top")
    for c in q[1]: c.font=Font(name="Arial",size=10,bold=True)
wb.save("fig04_bulge_contacts_source_data.xlsx")
