# -*-coding:Utf-8 -*
import Vue.Figure as figure
import Vue.Donnee as donnee
import Lecture.FonctionLectureClassique as classique
import numpy as np
import sys

########################################################################################
#-------------------------------          Input            -----------------------------
########################################################################################

if (len(sys.argv)==2):
	filename=sys.argv[1]
	print(filename+" will be printed")
else:
	print("You must give the name of the output file")
	sys.exit(1)


outputname_com="/".join(filename.split("/")[0:-1])+"/graphe_com_"+(filename.split("/")[-1]).split(".")[0]
outputname_emv="/".join(filename.split("/")[0:-1])+"/graphe_emv_"+(filename.split("/")[-1]).split(".")[0]
outputname_efr="/".join(filename.split("/")[0:-1])+"/graphe_efr_"+(filename.split("/")[-1]).split(".")[0]

########################################################################################
#-------------------------------          Figure           -----------------------------
########################################################################################

colors=["m","b","c","r","g","y","k","firebrick","purple"]
markers=["^","o",".","v"]
ncolor=0
nmarker=0


(eta,eps)   = classique.lecture(filename,0,1)
(com,emv)   = classique.lecture(filename,2,3)
(eta,efr)   = classique.lecture(filename,0,4)

Eta  = []
Eps  = []
Com  = []
Emv  = []
Efr  = []

Courbes_com = []
Courbes_emv = []
Courbes_efr = []

compt=0
ymax_com=0
ymax_emv =0
ymax_efr =0


offset = 9
for i in range(0,7):
	Eta.append(eta[compt+0:compt+offset])
	Eps.append(eps[compt+0])
	Com.append(com[compt+0:compt+offset])
	Emv.append(emv[compt+0:compt+offset])
	Efr.append(efr[compt+0:compt+offset])

	ymax_com=max(ymax_com,max(com[compt+0:compt+offset]))
	ymax_emv=max(ymax_emv,max(emv[compt+0:compt+offset]))
	ymax_efr=max(ymax_efr,max(efr[compt+0:compt+offset]))

	compt+=offset



for i in range(0,len(Eps)):
	line={"linestyle":"--","linewidth":3,"linecolor":colors[ncolor]}
	marker={"markerstyle":markers[1],"markersize":10,"fillstyle":"full"}

	Courbes_com.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Com[i],abscisse=Eta[0],line=line,marker=marker))
	Courbes_emv.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Emv[i],abscisse=Eta[0],line=line,marker=marker))
	Courbes_efr.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Efr[i],abscisse=Eta[0],line=line,marker=marker))

	ncolor+=1
	nmarker+=1

xlim=[min(Eta[0])*0.75,max(Eta[0])*1.25]
ylim_com=[0,ymax_com*1.25]
ylim_emv =[0,ymax_emv*1.25]
ylim_efr=[0,ymax_efr*1.25]

xlabel={"label":"$\eta$","fontsize":20}
ylabel_emv={"label":"MvProd Error","fontsize":20}
ylabel_com={"label":"Compression","fontsize":20}
ylabel_efr={"label":"Frobenius Error","fontsize":20}

# titre={"titre":"Test","fontsize":20,"loc":"center"}
legende={"loc":"upper right","bbox_to_anchor":None,"ncol":1,"fontsize":12}

# xscale='log'

Figure_emv=figure.Graphe1D(id=0,legende=legende,xlim=xlim,ylim=ylim_emv,xlabel=xlabel,ylabel=ylabel_emv)
Figure_com=figure.Graphe1D(id=1,legende=legende,xlim=xlim,ylim=ylim_com,xlabel=xlabel,ylabel=ylabel_com)
Figure_efr=figure.Graphe1D(id=2,legende=legende,xlim=xlim,ylim=ylim_efr,xlabel=xlabel,ylabel=ylabel_efr)

for courbe in Courbes_com:
	Figure_com.AjoutCourbe(courbe)

for courbe in Courbes_emv:
	Figure_emv.AjoutCourbe(courbe)

for courbe in Courbes_efr:
	Figure_efr.AjoutCourbe(courbe)

Figure_emv.TraceGraphe1D()
Figure_emv.EnregistreFigure(outputname_emv)

Figure_com.TraceGraphe1D()
Figure_com.EnregistreFigure(outputname_com)

Figure_efr.TraceGraphe1D()
Figure_efr.EnregistreFigure(outputname_efr)

Figure_com.FermeFigure()
Figure_emv.FermeFigure()
Figure_efr.FermeFigure()
