# -*-coding:Utf-8 -*
import Vue.Figure as figure
import Vue.Donnee as donnee
import Lecture.FonctionLectureClassique as classique
import numpy as np
import sys

########################################################################################
#-------------------------------          Input            -----------------------------
########################################################################################
filename=[]
if (len(sys.argv)>=2):
	filenames=sys.argv[1:]
else:
	print("You must give the name of the output file")
	sys.exit(1)


outputname_err="/".join(filenames[1].split("/")[0:-1])+"/graphe_"+"_".join((filenames[1].split("/")[-1]).split(".")[0].split("_")[0:2])+"_comp"


########################################################################################
#-------------------------------          Figure           -----------------------------
########################################################################################

linestyle=["-","--",":"]
colors=["m","b","c","r","g","y","k","firebrick","purple"]
markers=["^","o",".","v"]

nlinestyle=0
ymax_err=0
ymin_err=1e30

Courbes_dist = []
Courbes_rank = []
Courbes_err1 = []
Courbes_err2 = []

for filename in filenames:
	test="_".join((filename.split("/")[-1]).split(".")[0].split("_")[3:])


	(dist,rank)   = classique.lecture(filename,0,1)
	(err1,err2)   = classique.lecture(filename,2,3)


	Dist  = []
	Rank  = []
	Err1  = []
	Err2  = []



	compt=0


	offset = 49
	for i in range(1,10):
		Rank.append(rank[compt+0:compt+offset])
		Dist.append(dist[compt+0])
		Err1.append(err1[compt+0:compt+offset])
		Err2.append(err2[compt+0:compt+offset])

		ymax_err=max(ymax_err,max(err1[compt+0:compt+offset]))
		ymax_err=max(ymax_err,max(err2[compt+0:compt+offset]))
		ymin_err=min(ymin_err,min(err1[compt+0:compt+offset]))
		ymin_err=min(ymin_err,min(err2[compt+0:compt+offset]))

		compt+=offset


	ncolor=0
	for i in range(0,len(Dist)):
		line1={"linestyle":linestyle[nlinestyle],"linewidth":2,"linecolor":colors[ncolor]}
		# line2={"linestyle":"--","linewidth":3,"linecolor":colors[ncolor]}
		marker={"markerstyle":"None","markersize":10,"fillstyle":"full"}

		Courbes_err1.append(donnee.Ligne(nom=test+r"$ - distance=$"+str(Dist[i]),ordonnee=Err1[i],abscisse=Rank[i],line=line1,marker=marker))

		# Courbes_err2.append(donnee.Ligne(nom=r"$SVD - distance=$"+str(Dist[i]),ordonnee=Err2[i],abscisse=Rank[i],line=line2,marker=marker))


		ncolor+=1

	nlinestyle+=1

	# xlim=[min(Rank[0])*0.75,max(Rank[0])+1]
	# ylim_erro=[1e-20,ymax_err*1.25]
	xlim=[min(Rank[0])*0.75,10]
	ylim_erro=[1e-6,1]

	xlabel={"label":"Rank","fontsize":20}
	ylabel_erro={"label":"Relative error","fontsize":20}


	# titre={"titre":"Test","fontsize":20,"loc":"center"}
	legende={"loc":"upper left","bbox_to_anchor":(1.01,1),"ncol":1,"fontsize":12}

Figure_erro=figure.Graphe1D(id=0,legende=legende,xlim=xlim,ylim=ylim_erro,xlabel=xlabel,ylabel=ylabel_erro,yscale="log",axis="off")



for courbe in Courbes_err1:
	Figure_erro.AjoutCourbe(courbe)
# for courbe in Courbes_err2:
# 	Figure_erro.AjoutCourbe(courbe)


Figure_erro.TraceGraphe1D()

Figure_erro.EnregistreFigure(outputname_err)

Figure_erro.FermeFigure()
