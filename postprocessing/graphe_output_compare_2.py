# -*-coding:Utf-8 -*
import Vue.Figure as figure
import Vue.Donnee as donnee
import Lecture.FonctionLectureClassique as classique
import numpy as np
import sys
import matplotlib.pyplot as plt

########################################################################################
#-------------------------------          Input            -----------------------------
########################################################################################

if (len(sys.argv)==3):
	filename1=sys.argv[1]
	filename2=sys.argv[2]
else:
	print("You must give the name of the output file")
	sys.exit(1)

if len(filename1.split("/"))!=1:
	outputname="/".join(filename1.split("/")[0:-1])+"/graphe_compare_2_"+(filename1.split("/")[-1]).split(".")[0]+"_"+(filename2.split("/")[-1]).split(".")[0]
else:
	outputname="graphe_compare_2_"+(filename1.split("/")[-1]).split(".")[0]+"_"+(filename2.split("/")[-1]).split(".")[0]

print(filename1+" and "+filename2+" will be printed in "+outputname)

########################################################################################
#-------------------------------          Figure           -----------------------------
########################################################################################

colors=["m","b","c","r","g","y","k","firebrick","purple"]
markers=["v","o","*","^","s","p"]
ncolor=0
nmarker=0

Eta  = []
Eps  = []
Com  = []
Efr  = []


max_com = 0
max_efr = 0
min_com = 1e30
min_efr = 1e30

for filename in [filename1,filename2]:
	(eta,eps)   = classique.lecture(filename,0,1)
	(com,efr)   = classique.lecture(filename,2,4)

	compt=0
	offset = 6
	for i in range(0,3):
		Eta.append(eta[compt+0])
		Eps.append(eps[compt+0:compt+offset])
		Com.append(com[compt+0:compt+offset])
		Efr.append(efr[compt+0:compt+offset])

		max_com=max(max_com,max(com[compt+0:compt+offset]))
		max_efr=max(max_efr,max(efr[compt+0:compt+offset]))
		min_com=min(min_com,min(com[compt+0:compt+offset]))
		min_efr=min(min_efr,min(efr[compt+0:compt+offset]))

		compt+=offset

	ylim=[min_efr*0.1,max_efr*10]
	xlim=[0,1.1]

Nuages=[]

for i in range(0,len(Eta)):
	shape=np.ones([1,len(Eps[0])])
	shape*=100
	if (i<3):
		Nuages.append(donnee.Nuage(nom=r"D3 : $\eta=$"+str(Eta[i]),abscisse=Com[i],ordonnee=Efr[i],marker=markers[i],shape=shape,color=colors[i]))
	else:
		Nuages.append(donnee.Nuage(nom=r"DN1 : $\eta=$"+str(Eta[i]),abscisse=Com[i],ordonnee=Efr[i],marker=markers[i],shape=shape,color=colors[i]))

xlabel={"label":"Compression","fontsize":20}
ylabel={"label":"Frobenius error","fontsize":20}

legende={"loc":"lower right","bbox_to_anchor":None,"ncol":1,"fontsize":12}

Figure0 = figure.ScatterPlot(id=0,xlabel=xlabel,ylabel=ylabel,legende=legende,yscale='log',xlim=xlim,ylim=ylim)
for elt in Nuages:
	Figure0.AjoutNuage(elt)
Figure0.TraceScatterPlot()

Eps_ann = Eps[0].astype('|S4')
Eps_ann=map(lambda x: '0 blocks' if x == "-1.0" else x, Eps_ann)

Figure0.Annotations(annotation_text=Eps_ann,annotation_pos=[Com[0],Efr[0]],annotation_opt=(-5,10),annotation_color=colors[0])
Figure0.Annotations(annotation_text=Eps_ann,annotation_pos=[Com[1],Efr[1]],annotation_opt=(-5,-15),annotation_color=colors[1])
Figure0.Annotations(annotation_text=Eps_ann,annotation_pos=[Com[2],Efr[2]],annotation_opt=(-5,5),annotation_color=colors[2])
Figure0.Annotations(annotation_text=Eps_ann,annotation_pos=[Com[3],Efr[3]],annotation_opt=(-5,10),annotation_color=colors[3])
Figure0.Annotations(annotation_text=Eps_ann,annotation_pos=[Com[4],Efr[4]],annotation_opt=(-5,-15),annotation_color=colors[4])
Figure0.Annotations(annotation_text=Eps_ann,annotation_pos=[Com[5],Efr[5]],annotation_opt=(-5,5),annotation_color=colors[5])

Figure0.EnregistreFigure(outputname)
Figure0.FermeFigure()
