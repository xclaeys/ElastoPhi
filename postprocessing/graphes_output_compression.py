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
outputname_err="/".join(filename.split("/")[0:-1])+"/graphe_err_"+(filename.split("/")[-1]).split(".")[0]
outputname_tas="/".join(filename.split("/")[0:-1])+"/graphe_tas_"+(filename.split("/")[-1]).split(".")[0]
outputname_tmv="/".join(filename.split("/")[0:-1])+"/graphe_tmv_"+(filename.split("/")[-1]).split(".")[0]

########################################################################################
#-------------------------------          Figure           -----------------------------
########################################################################################

colors=["m","b","c","r"]
markers=["^","o",".","v"]
ncolor=0
nmarker=0


(eta,eps)   = classique.lecture(filename,0,1)
(com,err)   = classique.lecture(filename,2,3)
(tas,tmv)   = classique.lecture(filename,4,5)

Eta  = []
Eps  = []
Com  = []
Err  = []
Tas  = []
Tmv  = []

Courbes_com = []
Courbes_err = []
Courbes_tas = []
Courbes_tmv = []

compt=0
ymax_com=0
ymax_err =0
ymax_tas =0
ymax_tmv =0

offset = 4
for i in range(0,offset):
	Eta.append(eta[compt+0:compt+offset]) 
	Eps.append(eps[compt+0])
	Com.append(com[compt+0:compt+offset])
	Err.append(err[compt+0:compt+offset])
	Tas.append(tas[compt+0:compt+offset])
	Tmv.append(tmv[compt+0:compt+offset])
	
	ymax_com=max(ymax_com,max(com[compt+0:compt+offset]))
	ymax_err=max(ymax_err,max(err[compt+0:compt+offset]))
	ymax_tas=max(ymax_tas,max(tas[compt+0:compt+offset]))
	ymax_tmv=max(ymax_tmv,max(tmv[compt+0:compt+offset]))
	
	compt+=offset
	
	

for i in range(0,len(Eps)):
	
	line={"linestyle":"--","linewidth":3,"linecolor":colors[ncolor]}
	marker={"markerstyle":markers[nmarker],"markersize":10,"fillstyle":"full"}

	Courbes_com.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Com[i],abscisse=Eta[i],line=line,marker=marker))
	Courbes_err.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Err[i],abscisse=Eta[i],line=line,marker=marker))
	Courbes_tas.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Tas[i],abscisse=Eta[i],line=line,marker=marker))
	Courbes_tmv.append(donnee.Ligne(nom=r"$\varepsilon=$"+str(Eps[i]),ordonnee=Tmv[i],abscisse=Eta[i],line=line,marker=marker))
		
	ncolor+=1
	nmarker+=1

xlim=[min(Eta[0])*0.75,max(Eta[0])*1.25]
ylim_com=[0,ymax_com*1.25]
ylim_err =[0,ymax_err*1.25]
ylim_tas=[0,ymax_tas*1.25]
ylim_tmv =[0,ymax_tmv*1.25]

xlabel={"label":"$\eta$","fontsize":20}
ylabel_err={"label":"Error","fontsize":20}
ylabel_com={"label":"Compression","fontsize":20}
ylabel_tas={"label":"Assembling time","fontsize":20}
ylabel_tmv={"label":"Time for Matrix-vector product","fontsize":20}

# titre={"titre":"Test","fontsize":20,"loc":"center"}
legende={"loc":"upper right","bbox_to_anchor":None,"ncol":1,"fontsize":12}

Figure_err=figure.Graphe1D(id=0,legende=legende,xlim=xlim,ylim=ylim_err,xlabel=xlabel,ylabel=ylabel_err,xscale='log')
Figure_com=figure.Graphe1D(id=1,legende=legende,xlim=xlim,ylim=ylim_com,xlabel=xlabel,ylabel=ylabel_com,xscale='log')
Figure_tas=figure.Graphe1D(id=2,legende=legende,xlim=xlim,ylim=ylim_tas,xlabel=xlabel,ylabel=ylabel_tas,xscale='log')
Figure_tmv=figure.Graphe1D(id=3,legende=legende,xlim=xlim,ylim=ylim_tmv,xlabel=xlabel,ylabel=ylabel_tmv,xscale='log')


for courbe in Courbes_com:
	Figure_com.AjoutCourbe(courbe)
	
for courbe in Courbes_err:
	Figure_err.AjoutCourbe(courbe)
	
for courbe in Courbes_tas:
	Figure_tas.AjoutCourbe(courbe)
	
for courbe in Courbes_tmv:
	Figure_tmv.AjoutCourbe(courbe)
	
Figure_err.TraceGraphe1D()
Figure_err.EnregistreFigure(outputname_err)

Figure_com.TraceGraphe1D()
Figure_com.EnregistreFigure(outputname_com)

Figure_tas.TraceGraphe1D()
Figure_tas.EnregistreFigure(outputname_tas)

Figure_tmv.TraceGraphe1D()
Figure_tmv.EnregistreFigure(outputname_tmv)

Figure_com.FermeFigure()
Figure_err.FermeFigure()
Figure_tas.FermeFigure()
Figure_tmv.FermeFigure()


