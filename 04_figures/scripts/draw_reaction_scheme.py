import sys, subprocess
from PIL import Image, ImageChops
B="./"
SO="soffice.py"
src=open("scheme_base.py").read(); head=src[:src.index("H=0.13; O=0.17")]
head=head.replace("S=1.1; Y0=1.75","S=0.7; Y0=2.45").replace('18 if t!="DNA" else 15','14 if t!="DNA" else 11').replace("rr.font.size=Pt(12)","rr.font.size=Pt(9)")
head=head.replace('''        elif kind=="broken": line(p,q,GOLD,3,dash=True)''','''        elif kind=="broken": line(p,q,GOLD,3,dash=True)
        elif kind=="fly": line(p,q,MAG,2.25,dash=True)''')
assert '"fly"' in head
def build(lang):
    g={}; exec(head,g); s=g["s"]; panel=g["panel"]; tb=g["tb"]; Inches=g["Inches"]; Pt=g["Pt"]; RGBColor=g["RGBColor"]; MSO_SHAPE=g["MSO_SHAPE"]
    BLACK,MAG,GRN,PUR,GREY=g["BLACK"],g["MAG"],g["GRN"],g["PUR"],g["GREY"]
    H=0.12; O=0.15
    base={"DNAL":(0.6,2.4,"DNA",GREY,0.3),"O3":(1.75,2.4,"O",BLACK,O),"P":(3.0,2.4,"P",BLACK,O),"O1P":(3.0,1.3,"O",BLACK,O),"O2P":(3.6,3.4,"O",BLACK,O),"O5":(2.4,3.4,"O",BLACK,O),"DNAR":(2.4,4.3,"DNA",GREY,0.17),
     "Mn4":(2.2,0.45,"Mn",PUR,0.24),"Mn3":(5.05,1.4,"Mn",PUR,0.24),"W3O":(1.25,0.3,"O",BLACK,O),"W3H":(0.75,-0.05,"H",BLACK,H),"W2O":(0.5,1.45,"O",BLACK,O),"W2H":(-0.05,1.75,"H",BLACK,H)}
    common=[("DNAL","O3","b"),("P","O1P","b"),("P","O2P","b"),("P","O5","b"),("O5","DNAR","b"),("Ow","OwH","b"),("W2O","W2H","b"),("W3O","W3H","b"),("Mn4","O1P","coord"),("Mn4","W3O","coord"),("Mn3","Ow","coord")]
    A=dict(base); A.update({"Ow":(4.5,2.4,"O",GRN,O),"OwH":(5.05,2.85,"H",BLACK,H),"Hn":(4.45,3.2,"H¹",MAG,H),"Hr":(1.02,1.95,"H³",MAG,H),"Hm":(0.76,0.95,"H²",MAG,H)})
    bA=common+[("O3","P","scis"),("Ow","Hn","b"),("W2O","Hr","b"),("W3O","Hm","b")]; aA=[("Ow","P",GRN),("Hn","O2P",MAG)]
    T=dict(base); T["O3"]=(1.68,2.4,"O",BLACK,O)
    T.update({"Ow":(4.2,2.4,"O",GRN,O),"OwH":(4.75,2.85,"H",BLACK,H),"Hn":(4.2,3.75,"H¹",MAG,H),"Hr":(1.12,1.93,"H³",MAG,H),"Hm":(0.86,0.9,"H²",MAG,H)})
    bT=common+[("O3","P","broken"),("P","Ow","new"),("O2P","Hn","hnew"),("W2O","Hr","fly"),("Hr","O3","fly"),("W3O","Hm","fly"),("Hm","W2O","fly")]
    Bq=dict(base); Bq["O3"]=(1.6,2.4,"O",BLACK,O)
    Bq.update({"Ow":(4.2,2.4,"O",GRN,O),"OwH":(4.75,2.85,"H",BLACK,H),"Hn":(4.2,3.75,"H¹",MAG,H),"Hr":(1.2,1.95,"H³",MAG,H),"Hm":(0.72,0.95,"H²",MAG,H)})
    bB=common+[("O3","P","broken"),("P","Ow","new"),("O2P","Hn","hnew"),("O3","Hr","hnew"),("W2O","Hm","hnew")]
    T_={"ja":("反応物","遷移状態","切断後"),"en":("Reactant","Transition state","After cleavage")}[lang]
    Y0=g["Y0"]; BW=4.1
    for k,(atoms,bonds,arrows) in enumerate(((A,bA,aA),(T,bT,[]),(Bq,bB,[]))):
        X=0.2+k*(BW+0.31)
        bx=s.shapes.add_shape(MSO_SHAPE.RECTANGLE,Inches(X),Inches(Y0-0.45),Inches(BW),Inches(3.75)); bx.fill.background(); bx.line.color.rgb=RGBColor(0,0,0); bx.line.width=Pt(1.5)
        panel(X+0.28,atoms,bonds,arrows); tb(X,Y0-0.95,BW,0.4,T_[k],18)
        if k<2:
            ar=s.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW,Inches(X+BW+0.03),Inches(Y0+1.25),Inches(0.25),Inches(0.5)); ar.fill.solid(); ar.fill.fore_color.rgb=RGBColor(0xFF,0xFF,0); ar.line.color.rgb=RGBColor(0,0,0); ar.line.width=Pt(0.75)
    f=f"fig08_reaction_mechanism_{lang}.pptx"; g["prs"].save(f)
    subprocess.run(["python3",SO,"--headless","--convert-to","pdf",f],capture_output=True); subprocess.run(["pdftoppm","-r","300","-png",f.replace(".pptx",".pdf"),f"sch_{lang}"])
    im=Image.open(f"sch_{lang}-1.png").convert("RGB"); bb=ImageChops.difference(im,Image.new("RGB",im.size,"white")).getbbox(); im.crop((bb[0]-40,bb[1]-40,bb[2]+40,bb[3]+40)).save(f"fig08_reaction_mechanism_{lang}.png",dpi=(300,300))
build("ja"); build("en")
