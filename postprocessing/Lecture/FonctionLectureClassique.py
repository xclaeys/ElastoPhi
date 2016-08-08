# -*-coding:Utf-8 -*
import numpy as np

# Lecture classique par colonne, type gnuplot. 
#Renvoie les abscisses et ordonnées
def lecture(filename, n_abs, n_ord):
	abscisses=np.array([])
	ordonnees=np.array([])
	#print(os.getcwd())
	
	with open(filename) as output:
		for ligne in output:
			ligne_s=ligne.split()
			abscisses=np.append(abscisses,[float(ligne_s[n_abs])])
			ordonnees=np.append(ordonnees,[float(ligne_s[n_ord])])
				
	return (abscisses,ordonnees)

# Lecture classique par colonne, type gnuplot mais limite le nbr de coef. 
#Renvoie les abscisses et ordonnées
def lecturetronquee(filename, n_abs, n_ord,n_term):
	abscisses=np.array([])
	ordonnees=np.array([])
	#print(os.getcwd())
	
	with open(filename) as output:
		compteur=0
		for ligne in output:
			ligne_s=ligne.split()
			if compteur<=n_term:
				abscisses=np.append(abscisses,[float(ligne_s[n_abs])])
				ordonnees=np.append(ordonnees,[float(ligne_s[n_ord])])
			compteur+=1
				
	return (abscisses,ordonnees)
	
def lecturelignespecifiee(filename, n_abs, n_ord,n_ligne):
	abscisses=np.array([])
	ordonnees=np.array([])
	#print(os.getcwd())
	compteur=0
	with open(filename) as output:
		for ligne in output:
			if compteur==0:
				ligne_s=ligne.split()
				abscisses=np.append(abscisses,[float(ligne_s[n_abs])])
				ordonnees=np.append(ordonnees,[float(ligne_s[n_ord])])
				compteur+=1
			elif compteur==n_ligne-1:
				compteur=0
			else:
				compteur+=1
	return (abscisses,ordonnees)
	
# Lecture pour un graphe 2D 
def lecturesurface(filename,surface):
	#print(os.getcwd())
	
	with open(filename) as output:
		for ligne in output:
			ligne_s=ligne.split()
			surface[ligne_s[0],ligne_s[1]] = ligne_s[2]
			print(ligne_s[0],ligne_s[1],surface[ligne_s[0],ligne_s[1]])
# 			abscisses=np.append(abscisses,[float(ligne_s[n_abs])])
# 			ordonnees=np.append(ordonnees,[float(ligne_s[n_ord])])
