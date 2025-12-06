# Programme "Traitement_Deflectometry_integration_v1.R" 
# octobre 2022


rm(list=ls())	# Clear all variables from memory

# fonctions utilisées par le programme
source("Fonctions_2.R")

####################################
#Filters0<- Filters[c("txt"),]

# Sélection des fichiers à ouvrir
if (interactive()) 
  NomFichier_a_Traiter <- choose.files(filters = c("CSV files (*.csv)","*.csv"))
Nb_Fich <-length(NomFichier_a_Traiter)
# ?read.table
if (Nb_Fich!= 0)  AmpRSignal<-ThetaDegSignal<-phaseRef<-Position<-FundamentalFreqRef<-Np<-rep(0,Nb_Fich)
op <- par()# par(op)

Fich0<-rep("",Nb_Fich)
Data<- as.list(1:Nb_Fich)

for (numFich in 1:Nb_Fich) { 
# numFich=1  
  FichierData<-NomFichier_a_Traiter[numFich]
# Nom du fichier sans le chemin
    Fich0[numFich]= substr(FichierData,nchar(getwd())+2, nchar(FichierData))
# Lecture du fichier enregistré 
    # 3 colonnes :
    # Colonne 1 : Temps en s
    # Colonne 2 : Signal mesuré sur CH1 en V
    # Colonne 3 : Signal de référence TTL sur CH2 en V
    #  Séparateur de colonnes "," (comma) 
    # Données numériques à partir de la ligne 10
    # Unités sur la ligne 9
    Signal<-Signal0 <- read.csv(Fich0[numFich],header=TRUE, skip=8)
    #N0<-length(Signal0[,1])%/%1 ; Signal<-Signal0[1:N0,]# Diminue le nombre de points
# Ouverture d'une image du graphique sauvegardé dans le répertoire de travail
  png(paste(substr(Fich0[numFich],1,nchar(Fich0[numFich])-4),".png",sep=""),width = 1024, height = 768, res=1.5*96) 

# Renommer les colonnes
  names(Signal)<-c("time","Signal","ref")
  Data[[numFich]]<-Signal # Sauvegarde dans une structure de données
  Smin<-min(Signal[,2]);Smax<-max(Signal[,2])
  tmin<-min(Signal[,1]);tmax<-max(Signal[,1])
  Refmin<-min(Signal[,3]);Refmax<-max(Signal[,3])
  #TitrePlot<- substr(strsplit(Fich0,"z_")[[numFich]][2],1,6)
  #TitrePlot<- paste(as.numeric(substr(Fich0[numFich],nchar(Fich0[numFich])-9,nchar(Fich0[numFich])-6)),"µm")
  TitrePlot<- Fich0[numFich]
 par(mar = c(4.1,4.1,3,5))# Elargissement de la marge de droite
 plot(c(tmin,tmax),c(Refmin,Refmax),cex.title=0.5,type="n",main=TitrePlot,xlab="Temps (s)",ylab="Signal TTL référence (V)")
 lines(Signal[,1],Signal[,3],col="gray")# reference # Signal de référence
 RisingPoints<-UpZeroCrossing(Data[[numFich]][c(1,3)])# Détection des fronts montants
 points(RisingPoints,rep(2.5,length(RisingPoints)),cex=2, pch="." ,col="blue") # Affichage des fronts montants

 # Graphique signal mesuré sur CH1 ; l'échelle sur l'axe vertical droit 
 par(new=TRUE) # Nouveau graphique dans la même fenêtre
 plot(Signal[,1],Signal[,2]/1e-3,type="l",col="red", axes=FALSE,xlab=NA,ylab=NA)
 abline(h=0,v=0,lty=3)
 #plot(Signal[,1],Refmin+(Refmax-Refmin)*(Signal[,2]-Smin)/(Smax-Smin),type="l",col="red", axes=FALSE,xlab=NA,ylab=NA)
 axis(side = 4,col="red") # Traçage de l'axe
 mtext(side = 4, line = 3, 'Signal mesuré (mV)',col="red")# ylab du second signal

 dev.off() 
 }
#
# Fin du chargement des données
#
###############################################

###############################
#
#  Calcul Fréquence et phase des signaux de référence 
#

ResultRef<-as.data.frame(matrix(rep(NA,Nb_Fich*2),Nb_Fich,2))

for (numFich in 1:Nb_Fich)  ResultRef[numFich,] <- FreqPhaseRef(Data[[numFich]][c(1,3)])[1:2]

names(ResultRef)<- c("FreqRef","PhaseRefDeg")
FinalResultRef<-cbind(Fichier=Fich0,ResultRef)
print(FinalResultRef)

###############################
#
#  Démodulation synchrone 
#  Calcul amplitude et phase des signaux 
#
#


ResultMesure<-data.frame(matrix(rep(NA,Nb_Fich*5),Nb_Fich,5))
names(ResultMesure)<-colnames(DSAmpPhase(Data[[1]]), do.NULL = FALSE)

for (numFich in 1:Nb_Fich)  ResultMesure[numFich,]<- DSAmpPhase(Data[[numFich]])

DSResult<-as.data.frame(cbind(Fichier=Fich0,ResultMesure))

# Ecriture fichier de résultats dans un fichier texte  "Results_Lock-in_Simulation.txt"
# Séparateur de colonne : caractère tabulation "\t"
# Séparateur décimal : ","
# Efface le fichier existant précédemment
write.table(DSResult, file = "Results_Lock-in_Simulation.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ",", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
