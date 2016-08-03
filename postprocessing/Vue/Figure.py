# -*-coding:Utf-8 -*
import Vue.Donnee as donnee
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
##########################################################################################
#-------------------------------   Classe mère 	-----------------------------------------#
##########################################################################################

class Figure:
	"""Classe définissant une figure caractérisée par :
	- titre : dict[titre,fontsize,loc] 
	- xlabel : dict[label,fontsize]
	- ylabel : dict[label,fontsize]
	- xticks : int
	- yticks : int
	- id : int 
	- format : string 
	- xlim = [int,int]
	- ylim = [int,int]"""
	
	
	def __init__(self,id=0,titre=None,xlabel=None,ylabel=None,xticks=None,yticks=None,format="eps",xlim=None,ylim=None,xscale=None,yscale=None): # Constructeur
		self.titre=titre
		self.xlabel=xlabel
		self.ylabel=ylabel
		self.xticks=xticks
		self.yticks=yticks
		self.id=id
		self.format=format
		self.xlim=xlim
		self.ylim=ylim
		self.xscale=xscale
		self.yscale=yscale
		
	def __str__(self):
		return "Figure numero {0}".format(self.id)
		
	def OptionsFigure(self):
		# Choix de la figure courante
		plt.figure(self.id)

		# Titre
		if not self.titre==None:
			plt.title(self.titre["titre"],fontsize=self.titre["fontsize"],loc=self.titre["loc"])
			
		# Labels des axes
		if self.xlabel==None:
			plt.xlabel("x",fontsize=25)
		else:
			plt.xlabel(self.xlabel["label"],fontsize=self.xlabel["fontsize"])
			
		if self.ylabel==None:
			plt.ylabel("y",fontsize=25)
		else:
			plt.ylabel(self.ylabel["label"],fontsize=self.ylabel["fontsize"])
		
		# Tickes
		if self.xticks==None:
			plt.xticks(fontsize=15)
		else:
			plt.xticks(fontsize=self.xticks)
			
		if self.yticks==None:
			plt.yticks(fontsize=15)
		else:
			plt.yticks(fontsize=self.yticks)
			
		# Limites
		if not self.xlim==None:
			plt.xlim(self.xlim[0],self.xlim[1])
		if not self.ylim==None:
			plt.ylim(self.ylim[0],self.ylim[1])	
			
		# Echelles
		if not self.xscale==None:
			plt.xscale(self.xscale)
			
		if not self.yscale==None:
			plt.yscale(self.yscale)
			
		
		
		
	def EnregistreFigure(self,outputname):
		plt.figure(self.id)
		plt.savefig(outputname+"."+self.format,bbox_inches='tight',format=self.format)
		
	def FermeFigure(self):
		plt.close(self.id)
	
	#def MontreEtFermeFigure(self):
		#plt.show()


##########################################################################################
#------------------------------- Classes filles -----------------------------------------#
##########################################################################################

class Graphe1D(Figure):
	"""Classe définissant un graphe caractérisé par :
	- titre : string 
	- xlabel : dict[label,fontsize]
	- ylabel : dict[label,fontsize]
	- des ticks	: fontsize
	- id : int 
	- format : string 
	- xlim = [int,int]
	- ylim = [int,int]
	- donnee : liste de Courbes
	- legende : dict[loc,bbox_to_anchor,ncol,fontsize]
	"""

	def __init__(self,id=0,titre=None,xlabel=None,ylabel=None,xticks=None,yticks=None,format="eps",xlim=None,ylim=None,xscale=None,yscale=None,legende=None): # Constructeur
		Figure.__init__(self,id,titre,xlabel,ylabel,xticks,yticks,format,xlim,ylim,xscale,yscale) 
		self.legende=legende
		self.donnee=[]
	
	# A appeler après avoir tracé les courbes
	def OptionsFigure(self):
		Figure.OptionsFigure(self)
		if self.legende==None:
			plt.legend()
		else:
			plt.legend(loc=self.legende["loc"],bbox_to_anchor=self.legende["bbox_to_anchor"],ncol=self.legende["ncol"], fancybox=True, shadow=False,fontsize=self.legende["fontsize"])
	
	def AjoutCourbe(self,courbe):
		self.donnee.append(courbe)
		if courbe.nom==None:
			courbe.nom="Courbe n°"+len(self.donnee)
		if courbe.line==None:
			#  linestyle '-' | '--' | '-.' | ':' | 'None' | ' ' | ''
			#  linecolor b g r c m y k w
			courbe.line={"linestyle":"s","linewidth":2,"linecolor":"b"}
		if courbe.marker==None:
			# fillstyle ('full', 'left', 'right', 'bottom', 'top', 'none')
			# markerstyle {0: 'tickleft', 1: 'tickright', 2: 'tickup', 3: 'tickdown', ',': 'pixel', '8': 'octagon', 6: 'caretup', '<': 'triangle_left', '*': 'star', 's': 'square', '1': 'tri_down', 'v': 'triangle_down', 'o': 'circle', '2': 'tri_up', '^': 'triangle_up', '.': 'point', 7: 'caretdown', '4': 'tri_right', '3': 'tri_left', ' ': 'nothing', 'h': 'hexagon1', 5: 'caretright', 'd': 'thin_diamond', 'H': 'hexagon2', '|': 'vline', 'x': 'x', '>': 'triangle_right', 'p': 'pentagon', '_': 'hline', 'D': 'diamond', None: 'nothing', '': 'nothing', '+': 'plus', 'None': 'nothing', 4: 'caretleft'}
			courbe.marker={"markerstyle":".","markersize":15,"fillstyle":"full"}
			
	def TraceGraphe1D(self):
		plt.figure(self.id)
		for elt in self.donnee:
			if isinstance(elt,donnee.Ligne):
				plt.plot(elt.abscisse,elt.ordonnee,label=elt.nom
					,linestyle=elt.line["linestyle"],color=elt.line["linecolor"],linewidth=elt.line["linewidth"],
					marker=elt.marker["markerstyle"],markersize=elt.marker["markersize"],fillstyle=elt.marker["fillstyle"])
			if isinstance(elt,donnee.LigneHorizontale):
				plt.axhline(y=elt.value,label=elt.nom
					,linestyle=elt.line["linestyle"],color=elt.line["linecolor"],linewidth=elt.line["linewidth"],
					marker=elt.marker["markerstyle"],markersize=elt.marker["markersize"],fillstyle=elt.marker["fillstyle"])
			
		self.OptionsFigure()
		
		
		
	
#----------------------------------------------------------------------------------------#	
#----------------------------------------------------------------------------------------#
class Graphe2D(Figure):
	"""Classe définissant une figure caractérisée par :
	- titre : dict[titre,fontsize,loc] 
	- xlabel : dict[label,fontsize]
	- ylabel : dict[label,fontsize]
	- xticks : int
	- yticks : int
	- id : int 
	- format : string 
	- xlim = [int,int]
	- ylim = [int,int]
	- donnee = Surface
	- surface_masquee
	- interpolation : string
	- cmap : string (colormap)
	- origin : string
	"""
	
	def __init__(self,id=0,titre=None,xlabel=None,ylabel=None,xticks=None,yticks=None,format="eps",xlim=None,ylim=None,xscale=None,yscale=None,interpolation="nearest",cmap=cm.YlGnBu,origin="lower"): # Constructeur
		Figure.__init__(self,id,titre,xlabel,ylabel,xticks,yticks,format,xlim,ylim,xscale,yscale) 
		self.donnee=None
		self.surface_masquee=None
		self.interpolation=interpolation
		self.cmap=cmap
		self.origin=origin

	def OptionsFigure(self):
		Figure.OptionsFigure(self)

	def AjoutSurface(self,donnee):
		self.donnee=donnee

	def MasqueSurface(self,mask):
		self.surface_masquee = np.ma.masked_array(self.donnee.surface[:,:],mask[:,:])
		
	#def AjoutContour(self,inline=1,fontsize=12):
		#CS=plt.contour(self.donnee.ordonnee,extent=[self.xlim[0],self.xlim[1],self.ylim[0],self.ylim[1]],colors='k')
		#plt.clabel(CS,inline=inline,fontsize=fontsize)

	def TraceSurface(self):
		plt.figure(self.id)
		#print(self.donnee)
		plt.imshow(self.donnee.surface,interpolation=self.interpolation,origin=self.origin,cmap=self.cmap,extent=[self.donnee.xextent[0],self.donnee.xextent[1],self.donnee.yextent[0],self.donnee.yextent[1]])
		plt.colorbar()
		self.OptionsFigure()
		
	def TraceSurfaceMasquee(self):
		plt.figure(self.id)

		plt.imshow(self.surface_masquee,interpolation=self.interpolation,origin=self.origin,cmap=self.cmap,extent=[self.donnee.xextent[0],self.donnee.xextent[1],self.donnee.yextent[0],self.donnee.yextent[1]])
		self.OptionsFigure()
		plt.colorbar()
		
		
		
		
		