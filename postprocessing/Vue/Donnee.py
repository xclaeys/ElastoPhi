# -*-coding:Utf-8 -*
##########################################################################################
#-------------------------------   Classe mère 	-----------------------------------------#
##########################################################################################
class Donnee:
	"""Classe définissant les données pour une graphe caractérisé par :
	- nom
	- unite
	"""
	
	
	def __init__(self,nom=None,unite=None): # Constructeur
		self.nom=nom
		self.unite=unite
		
	def __str__(self):
		return "Donnee {0}".format(self.nom)
	
##########################################################################################
#------------------------------- Classes filles -----------------------------------------#
##########################################################################################
class Courbe(Donnee):
	"""Classe définissant une courbe caractérisé par :
	- abscisse
	- ordonnee
	- nom
	- unite
	- line
	- marker"""
	
	
	def __init__(self,nom=None,unite=None,line=None,marker=None): # Constructeur
		Donnee.__init__(self,nom,unite)
		self.line=line
		self.marker=marker
		
		
		
#----------------------------------------------------------------------------------------#	
#----------------------------------------------------------------------------------------#
	
class Surface(Donnee):
	"""Classe définissant les données pour une graphe caractérisé par :
	- ordonnee
	- nom
	- unite
	- extent
	"""

	def __init__(self,surface=None,extent=None,nom=None,unite=None): # Constructeur
		Donnee.__init__(self,nom,unite)
		self.surface=surface
		self.extent=extent
		
			
		
		
		
##########################################################################################
#--------------------------- Classes petites filles -------------------------------------#
##########################################################################################
class Ligne(Courbe):
	"""Classe définissant une courbe caractérisé par :
	- abscisses
	- ordonnees
	- nom
	- unite
	- line
	- marker"""
	
	def __init__(self,abscisse=None,ordonnee=None,nom=None,unite=None,line=None,marker=None): # Constructeur
		Courbe.__init__(self,nom,unite,line,marker)
		self.abscisse=abscisse
		self.ordonnee=ordonnee
		
		
class LigneHorizontale(Courbe):
	"""Classe définissant une courbe caractérisé par :
	- value
	- nom
	- unite
	- line
	- marker"""
	
	def __init__(self,value=None,nom=None,unite=None,line=None,marker=None): # Constructeur
		Courbe.__init__(self,nom,unite,line,marker)
		self.value=value
		
		
