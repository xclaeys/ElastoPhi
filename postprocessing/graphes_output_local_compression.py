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


outputname="/".join(filename.split("/")[0:-1])+"/graphe_mapp_"+".".join((filename.split("/")[-1]).split(".")[:-1])

########################################################################################
#-------------------------------          Figure           -----------------------------
########################################################################################

colors=["m","b","c","r"]
markers=["^","o",".","v"]
ncolor=0
nmarker=0



nc = 45
nr = 45
surface = np.ones((nr,nc))
classique.lecturesurface(filename,surface)



# 				for j in range(0,len(ordonnee)):
# 					surface[int(round(abscisse_1[j]))+(int(1./(2*lc))-1)/2,int(round(abscisse_2[j]))+(int(1./(2*lc))-1)/2]=ordonnee[j]
# 					mask[int(round(abscisse_1[j]))+(int(1./(2*lc))-1)/2,int(round(abscisse_2[j]))+(int(1./(2*lc))-1)/2]=0
					#surface[int(round(abscisse_1[j])),int(round(abscisse_2[j]))]=ordonnee[j]
					#mask[int(round(abscisse_1[j])),int(round(abscisse_2[j]))]=0
				#print(abscisse_1[j],abscisse_2[j],int(round(abscisse_1[j])),int(round(abscisse_2[j])),ordonnee[j])
				#print(mask)
				#print(surface)
				#if maillage==3:
# 				xextent=[lc*min(abscisse_1)-lc*(int(1./(2*lc))-1)/2,lc*max(abscisse_1)+lc*(int(1./(2*lc))-1)/2]
# 				yextent=[lc*min(abscisse_2)-lc*(int(1./(2*lc))-1)/2,lc*max(abscisse_2)+lc*(int(1./(2*lc))-1)/2]
				#else:
				#xextent=[lc*min(abscisse_1),lc*max(abscisse_1)]
				#yextent=[lc*min(abscisse_2),lc*max(abscisse_2)]

xlabel={"label":"x","fontsize":25}
ylabel={"label":"y","fontsize":25}
titre={"titre":"Local compression rate","fontsize":25,"loc":"center"}
				#legende={"loc":"upper left","bbox_to_anchor":None,"ncol":1,"fontsize":15}

Figure0=figure.Graphe2D(id=0,titre=titre,xlabel=xlabel,ylabel=ylabel,interpolation="none",origin="upper",axis="off")

Surface=donnee.Surface(surface=surface)

Figure0.AjoutSurface(Surface)
				#Figure0.MasqueSurface(np.transpose(mask))
				#Figure0.TraceSurfaceMasquee()
Figure0.TraceSurface()
Figure0.EnregistreFigure(outputname)
Figure0.FermeFigure()
