

####################################################
#	Fonction signal
#
#	Somme_sur_i[ A_i * Cos(2*pi*f_i*t)*exp(-alpha_i*t) ]
#
SignalFunction <- Vectorize(function(t,p) {
A <- p[[1]]			# Amplitudes
f <- Re(p[[2]])/2/pi	# Fréquences
alpha <- Im(p[[2]])	# Coefficients d'atténuation
phi <- p[[3]]		# Phases

sum(A * cos(2*pi*f*t+phi)*exp(-alpha*t))
},"t")
#
#	Fin fonction signal
#################################


###############################################
# Function : Data matrix Y
# Utilisée par la fonction "MatPencil()"

DataY <- function(x,L) {
N <- length(x)
YL=matrix(rep(NA,(N-L)*(L+1)),c(N-L,L+1))

# First column of Y
Y0=x[1:(N-L)]
YL[,1]=Y0
# Fonctions.R

# Filling the columns of Y
for (i in 2:(L+1)) {
Y0=c(Y0[-1],x[N-L+i-1])
YL[,i]=Y0
}
YL
}
# End Data matrix
################################################

# Moore-Penrose pseudo-inverse #
#
#  Utilisée par la fonction "MatPencil()"
##
MPPI=function(M) solve(Conj(t(M))%*%M)%*% Conj(t(M))

# Function moving average : vector v  
mg=function(v,n){
v0=v
N=length(v)
for (i in ((1+n):(N-n)))
v0[i]=mean(v[(i-n):(i+n)])
v0
} 
############################################

##########################################
# Matrix Pencil
#####################################
#	Filtrage du bruit avec la méthode 
#	du faisceau de matrices (matrix pencil)
#
# Ref.  Schöpfer et al., CEAS Aeronautic J. (2013)
#
####################################
# The matrix pencil number must be chosen in order to minimize
# the effect of noise in the data sequence.
# The optimum value of L is such that: (Np/3)< L <(Np/2)
# Ref. : Sarkar & Pereira, IEEE (1995)
# L=floor((N %/% 2 + N %/% 3)/2)
# L(optimum)=41 with Np=100
##############################################

# Hypothèse : Le signal à retrouver est sous la forme d'une 
# somme de  M fonctions sinusoïdales amorties:

# S(t)= Somme{m=1;M} A[m]*cos(2*pi*f[m])*t+phi[m])*exp(-alpha[m]*t)
# Ce signal est superposé à un bruit : b(t).

# Sous forme complexe, le signal s'écrit :
# 
# S(t)= Somme{m=1;M} a[m]*exp(i*omega[m])*t)+ c.c.
#    avec c.c. : complexe conjugué
#    et     a[m] = A[m]/2*exp(i*phi[m]),
#           omega[m]=2*pi*f[m]+i*alpha[m]

# Données d'entrée

# Data : tableau de données
#      colonne 1 : Vecteur Temps
#      colonne 2 : vecteur Signal
#              M : Nombre de modes à rechercher
#           Nmax : Nombre maximum de points à traiter dans le fichier
#                  afin d'éviter une trop longue durée de calcul.
#           Fmin : fréquence minimale à rechercher 

# Données de sorties 
#  Tableau de données : à M lignes et 4 colonnes
# Pour chaque ligne (mode m), on a : 
#          Freq : f[m] fréquence du mode 
#      CoeffAtt : alpha[m], le coefficient d'atténuation
#           Ampl : A[m], l'amplitude du mode
#     Phase_deg : phi[m], phase en degré
# 

MatPencil=function(Data,M=16,Nmax=2000,Fmin=0.001) {
 Time<-Data[,1] 
 Amp<-Data[,2] 
 # Echantillonnage
 Te <- (Time[length(Time)]-Time[1])/(length(Time)-1)
 (Fe <- 1/Te)	# Fréquence d'échantillonnage

 N <- length(Time)	# Nombre de points
 if (N>Nmax) N<-Nmax
 L <- floor((N %/% 2 + N %/% 3)/2)
 Y <- DataY(Amp,L)   # dim(Y)=c() (N-L) , (L+1) )
 # The Hankel matrices of the data sequence
 Y1 <- Y[,-(L+1)]	# The last column is removed
 Y2 <- Y[,-1]	# The first column is removed
 # The singular-value decomposition is used to "filter" the noise in the data
 # Singular-value decomposition (SVD) of Y1 : Y1 = U %*% D %*% V 
 svdY <- svd(Y1)
 # Diagonal matrix of singular values
 D <- diag(svdY$d)
  # Matrices of the SVD decomposition
 U <- svdY$u
 V <- svdY$v
 
 # Selection of the first M  columns of the matrices U and V
 D1 <- diag(svdY$d[1:M])
 U1 <- U[,1:M]
 V1 <- V[,1:M]
 #  Matice X2 : voir Schöpfer et al., CEAS Aeronautic J. (2013)
 X2 <- Conj(t(U1))%*% Y2 %*% V1 
  # Eigenvalue problem
 eigP <- eigen(MPPI(X2) %*% D1)
  # Eigenvalues: zi=exp(i*2*pi*fi*Te), where fi are frequencies of the data
 eigVal <- eigP$values
  # Eigenfrequencies in Hz 
 (freq <- -log(eigVal)/Te/2/pi/1i)	# OK 
  # eigenfrequencies 
 Fr <- Re(freq)
 Frequencies <- (sort(Fr[Fr>Fmin], index.return = TRUE)$x)  # 
 IX <- sort(Fr[Fr>Fmin], index.return = TRUE)$ix 
 Att <- Im(freq)*2*pi
 CoeffAtt <- abs(Att[Fr>Fmin][IX])
 ##
  #################################
 # Linear least-square problem
 #################################
  eigVect <- eigP$vectors
 A <- matrix(rep(NA,M*N),N,M)  # matrix initialization
 #
 for (i in 1:N) A[i,]=(eigVal)^(-(i-1))
 # Calcul des amplitudes complexes avec 
 Ri <- MPPI(A) %*% Amp 
 Ri0 <- Ri[Fr>0][IX]
 Amplitudes <- Mod(Ri0)# round(Mod(Ri0),1)*2 
 Phases <- Arg(Ri0)/2/pi*180
 as.data.frame(cbind(Freq=round(Frequencies,4),CoeffAtt=round(CoeffAtt,4), Ampl=Amplitudes, Phase_deg=round(Phases,4) ))

}

# MatPencil(Data[[4]][c(1,3)] )


##############################################
#
#  UpZeroCrossing(Data) 

#       Fonction trouvant les instants de passage à zéro
#       dans le sens de la montée
#       d'un signal alternatif
#       
# Données d'entrée :

# signal : tableau de données
#      colonne 1 : Vecteur Temps
#      colonne 2 : vecteur d'amplitude du signal

# Donnée de sortie : vecteur numérique des instants de "zero crossing"
# 
UpZeroCrossing=function(signal) {
  # signal<-Data[[4]][c(1,3)]
  vT<-signal[,1]# Time vector
  Refmin<-min(signal[,2]);Refmax<-max(signal[,2])
  RefMoy <-(Refmax + Refmin)/2
  ys<- signal[,2]-RefMoy# Suppression du décalage
  
  ys1<-c(ys[-1],ys[length(ys)])
  Ind<-which((ys1>0 & ys<0)==TRUE)	# points de passage à zéro : sens montée
  #Indb<-which((ys1<0 & ys>0)==TRUE)	# points de passage à zéro : sens descente
  
  # Estimation par interpolation linéaire des instants de passage par zéro,
  # sens de la montée
 return(vT[Ind]-ys[Ind]/(ys1[Ind]-ys[Ind])*(vT[Ind+1]-vT[Ind]))
  #tm<-vT[Indb]-ys[Indb]/(ys1[Indb]-ys[Indb])*(vT[Indb+1]-vT[Indb]) # Sens de la descente
  
}
# UpZeroCrossing(Data[[4]])

##############################################
#
#  FreqPhaseRef(Signal) 

#       Fonction trouvant : 
#       la fréquence, la phase initiale (t=0),
#       la variation de la période, par période
#       d'un signal rectangulaire ou sinusoïdal chirpé 
#   
#
# Données d'entrée :

# signal : tableau de données
#      colonne 1 : Vecteur Temps
#      colonne 2 : vecteur d'amplitude du signal

# Donnée de sortie : vecteur numérique c(FreqRef, PhaseRef_Deg, QuadCoef) 
#      Fréquence de référence      : FreqRef
#      Phase de référence en degré : PhaseRef_Deg (entre -180° et +180°)
#      coefficient quadratique     : QuadCoef

FreqPhaseRef<- function(Signal) {
# Signal<-Data[[4][c(1,3)]  
 tp<- UpZeroCrossing(Signal)
   Nz<-length(tp);Numz<-1:Nz#
  fm<-lm(tp~Numz+I(Numz^2))# Régression linéaire
  
  #summary(fm) # Résumé régression 
  # summary(fm)[[4]][3,4] # Pr(>|t|) pour coefficient quadratique
  coeffRegLin<-  as.numeric(coef(fm)) # Coefficients : (1) ordonnée à l'origine et (2) pente 
  Tref<-coeffRegLin[2] # pente = période du signal
  t0<-coeffRegLin[1] # instant du début de la première période
  
  if(summary(fm)[[4]][3,4] >0.05) q0<-0 else q0<-coeffRegLin[3]
  # Fréquence, phase entre -180° et +180°, et non linéarité quadratique
return(as.data.frame(cbind(FreqRef=1/Tref, PhaseRefDeg= (t0/Tref-round(t0/Tref))*360,QuadCoef=q0/Tref)) )
  
}

# FreqPhaseRef(Data[[4]][c(1,3)])



#######################################
#
#   DS_Amp_Phase(Signaux)
#
#    Traitement d'un signal 
#    Simulation d'une Détection synchrone 
#
#
# Données d'entrée : Signaux
#           - Colonne 1 : Vecteur "Temps"
#             Colonne 2 : Vecteur signal mesuré
#             Colonne 3 : Vecteur signal de référence
#
#             Ordre de l'harmonique : k (défaut k=1)
#
# Données en sortie : Vecteur c(FreqRef, PhaseRef, AmpSignal, PhaseSignal, TimeConstant)
#
#             FreqRef      : Fréquence de référence
#             PhaseRef     : Phase du signal de référence
#             AmpSignal    : Amplitude du signal à mesurer
#             PhaseSignal  : Phase du signal à mesurer
#             TimeConstant : Constante d'intégration

DSAmpPhase<- function(Signaux,k=1) {
#   Signaux<- Data[[4]]
Temps<- Signaux[,1]
Te<-Temps[2]-Temps[1]
SignalRef <-  Signaux[,c(1,3)]
SignalMes<-  Signaux[,2]
CoeffRef<-as.numeric(FreqPhaseRef(SignalRef))
FreqRef<- CoeffRef[1]
CosRef<-cos(2*pi*FreqRef*k*Temps )  
SinRef<-sin(2*pi*FreqRef*k*Temps )  

MultCos<- SignalMes*CosRef
MultSin<- SignalMes*SinRef
Tint<-diff(range(Temps))# Tint*FreqRef
X<-sum(MultCos)*Te/Tint*2
Y<-sum(MultSin)*Te/Tint*2
(R<-sqrt(X^2+Y^2))# ;AmpRSignal[numFich]<-R
(Theta=Arg(X+Y*1i)); (ThetaDegSignal=Theta/pi*180-CoeffRef[2])

return(cbind(FreqRef=FreqRef,PhaseRefDeg=CoeffRef[2],AmpMes=R,PhaseMesDeg=ThetaDegSignal,TimeConstant=Tint))  

}

# DSAmpPhase(Data[[4]])

