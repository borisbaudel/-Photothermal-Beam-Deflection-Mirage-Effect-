%pkg load optim; % pour octave

clear all;
% Model to fit
F=@(p,x) p(1)+p(2)*erf(sqrt(2)*(x-p(3))/p(4)) ;

% Lecture des données 
% Le fichier 'miroir-couteau_data.csv' est donné en exemple 
% pour tester le programme avec : a=10,b=950,y0=1,w=0.5 comme paramètres initiaux
% Le bruit a un écart-type de 5 µA rms
%dataxy = readmatrix('miroir-couteau_data.csv')
%dataxy = tabl1_I_y__1_; % "readmatrix" avec MatLab
%dataxy = dlmread('miroir-couteau_data.csv');% "dlmread" avec octave
%xdata=dataxy(:,1);% Colonne des positions du couteau, en mm
%ydata=dataxy(:,2);% Colonne des intensités mesurées, en µA

xdata=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
ydata=[-947 -944 -944 -931 -928 -897 -842 -716 -539 -288 3 304 553 749 852 922 951 961 962 966 969];

[N Nx]=size(ydata);
plot(xdata,ydata,'o');
grid;
xlabel('Position y (mm) ')
ylabel('I (µA)')
% Initialisation paramètres

pinit=[(max(ydata)+ min(ydata))/2 (ydata(N)- ydata(1))/2 (max(xdata)+ min(xdata))/2 abs((xdata(1)- xdata(N)))/4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ajustement des paramètres sur les données expérimentales
% La fonction "lsqcurvefit" du package "optim" trouve les paramètres
% qui ajustent "au mieux" le modèle sur les données expérimantales
[psol residuals]=lsqcurvefit(F, pinit,xdata,ydata);

% Ajout de la courbe ajustée sur les points expérimentaux
hold on;
xfit= linspace(0,max(xdata),101);
yfit=F(psol,xfit);
plot(xfit,yfit,'-b'); % ligne continue bleue 
a= psol(1) % Premier paramètre : valeur moyenne
b= psol(2)% Valeur de b en µA : Amplitude des variations
y0 = psol(3) % Position centrale faisceau
w= psol(4)% Valeur du rayon gaussien

text(y0+w,F(psol,y0+0.8*w),['w= ', num2str(w, 4),' mm']);
text(y0+w,F(psol,y0+w/2),['b= ', num2str(b, 4),' µA']);
legend('points expérimentaux','Meilleur ajustement','Location','Best')
% Calcul de la sensibilité (maximum de pente) 
Sensitivity2= sqrt(8/pi)*abs(psol(2))/psol(4)% µA/mm
text(y0+0.1*w,F(psol,y0),['S_2= ', num2str(Sensitivity2, 4),' µA/mm']);

% Fin
