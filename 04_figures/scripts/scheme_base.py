import math, copy
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.dml import MSO_LINE
from pptx.oxml.ns import qn
from lxml import etree
S=1.1; Y0=1.75
BLACK="000000"; MAG="D4009F"; GRN="1E9E2A"; GOLD="D99A00"; PUR="7B2D9B"; GREY="6E6E6E"; RED="C8281E"
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5); s=prs.slides.add_slide(prs.slide_layouts[6])
def tb(x,y,w,h,t,sz,col=BLACK,bold=False,center=True):
    b=s.shapes.add_textbox(Inches(x),Inches(y),Inches(w),Inches(h)); tf=b.text_frame; tf.margin_left=tf.margin_right=tf.margin_top=tf.margin_bottom=0
    tf.vertical_anchor=MSO_ANCHOR.MIDDLE; tf.word_wrap=False
    if center: tf.paragraphs[0].alignment=PP_ALIGN.CENTER
    r=tf.paragraphs[0].add_run(); r.text=t; r.font.size=Pt(sz); r.font.bold=bold; r.font.name="Arial"; r.font.color.rgb=RGBColor.from_string(col); return b
def line(p,q,col=BLACK,w=2.0,dash=False,arrow=False):
    c=s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT,Inches(p[0]),Inches(p[1]),Inches(q[0]),Inches(q[1]))
    c.line.color.rgb=RGBColor.from_string(col); c.line.width=Pt(w)
    if dash: c.line.dash_style=MSO_LINE.DASH
    if arrow:
        ln=c.line._get_or_add_ln(); e=etree.SubElement(ln,qn("a:tailEnd")); e.set("type","triangle"); e.set("w","med"); e.set("len","med")
    return c
def panel(X0,atoms,bonds,arrows):
    P={}
    for n,(x,y,t,col,r) in atoms.items(): P[n]=(X0+x*S,Y0+y*S,r)
    def seg(a,b,pad=0.0):
        (x1,y1,r1),(x2,y2,r2)=P[a],P[b]; d=math.hypot(x2-x1,y2-y1); ux,uy=(x2-x1)/d,(y2-y1)/d
        return (x1+ux*(r1+pad),y1+uy*(r1+pad)),(x2-ux*(r2+pad),y2-uy*(r2+pad))
    for a,b,kind in bonds:
        p,q=seg(a,b)
        if kind=="b": line(p,q)
        elif kind=="scis": line(p,q,GOLD,5)
        elif kind=="new": line(p,q,GRN,3)
        elif kind=="hnew": line(p,q,MAG,3)
        elif kind=="coord": line(p,q,PUR,1.25,dash=True)
        elif kind=="broken": line(p,q,GOLD,3,dash=True)
    for a,b,col in arrows:
        p,q=seg(a,b,0.03); line(p,q,col,2.5,arrow=True)
    for n,(x,y,t,col,r) in atoms.items():
        cx,cy,_=P[n]
        if t=="Mn":
            sh=s.shapes.add_shape(MSO_SHAPE.OVAL,Inches(cx-r),Inches(cy-r),Inches(2*r),Inches(2*r)); sh.fill.solid(); sh.fill.fore_color.rgb=RGBColor.from_string(PUR); sh.line.fill.background()
            tf=sh.text_frame; tf.margin_left=tf.margin_right=tf.margin_top=tf.margin_bottom=0; tf.word_wrap=False; tf.paragraphs[0].alignment=PP_ALIGN.CENTER
            rr=tf.paragraphs[0].add_run(); rr.text="Mn²⁺"; rr.font.size=Pt(12); rr.font.bold=True; rr.font.name="Arial"; rr.font.color.rgb=RGBColor(255,255,255)
        else:
            w=0.9 if t=="DNA" else 0.5
            tb(cx-w/2,cy-0.2,w,0.4,t,18 if t!="DNA" else 15,col,bold=(col!=BLACK and col!=GREY) or t=="P")
H=0.13; O=0.17
base={ "DNAL":(0.6,2.4,"DNA",GREY,0.36),"O3":(1.75,2.4,"O",BLACK,O),"P":(3.0,2.4,"P",BLACK,O),"O1P":(3.0,1.3,"O",BLACK,O),
 "O2P":(3.6,3.4,"O",BLACK,O),"O5":(2.4,3.4,"O",BLACK,O),"DNAR":(2.4,4.3,"DNA",GREY,0.2),
 "Mn4":(2.2,0.45,"Mn",PUR,0.3),"Mn3":(5.05,1.4,"Mn",PUR,0.3),
 "W3O":(1.25,0.3,"O",BLACK,O),"W3H":(0.75,-0.05,"H",BLACK,H),"W2O":(0.5,1.45,"O",BLACK,O),"W2H":(-0.05,1.75,"H",BLACK,H)}
A=dict(base); A.update({"Ow":(4.5,2.4,"O",GRN,O),"OwH":(5.05,2.85,"H",BLACK,H),"Hn":(4.45,3.15,"H",MAG,H),"Hr":(0.98,1.87,"H",MAG,H),"Hm":(0.9,0.78,"H",MAG,H)})
bA=[("DNAL","O3","b"),("O3","P","scis"),("P","O1P","b"),("P","O2P","b"),("P","O5","b"),("O5","DNAR","b"),("Ow","OwH","b"),("Ow","Hn","b"),
    ("W2O","Hr","b"),("W2O","W2H","b"),("W3O","Hm","b"),("W3O","W3H","b"),("Mn4","O1P","coord"),("Mn4","W3O","coord"),("Mn3","Ow","coord")]
aA=[("Ow","P",GRN),("Hn","O2P",MAG),("Hr","O3",MAG),("Hm","W2O",MAG)]
B=dict(base); B["O3"]=(1.6,2.4,"O",BLACK,O)
B.update({"Ow":(4.2,2.4,"O",GRN,O),"OwH":(4.75,2.85,"H",BLACK,H),"Hn":(4.2,3.75,"H",MAG,H),"Hr":(1.2,1.95,"H",MAG,H),"Hm":(0.72,0.95,"H",MAG,H)})
bB=[("DNAL","O3","b"),("O3","P","broken"),("P","O1P","b"),("P","O2P","b"),("P","O5","b"),("O5","DNAR","b"),("Ow","OwH","b"),("P","Ow","new"),("O2P","Hn","hnew"),
    ("O3","Hr","hnew"),("W2O","W2H","b"),("W2O","Hm","hnew"),("W3O","W3H","b"),("Mn4","O1P","coord"),("Mn4","W3O","coord"),("Mn3","Ow","coord")]
XA=0.55; XB=7.05
panel(XA,A,bA,aA); panel(XB,B,bB,[])
for X,l,t in ((XA,"a","反応前"),(XB,"b","反応後")):
    tb(X-0.3,Y0-1.0,0.5,0.45,l,24,bold=True,center=False); tb(X-0.3,Y0+4.75*S,5.9,0.4,t,18)
prs.save("Fig_mechanism_scheme_editable.pptx")
