(* ::Package:: *)

Zal=10./137.;
y=Zal;
m=1.;
NumBR=15;
BohrRadius=1/(m*Zal);
Efm=m (Sqrt[1-Zal^2]);
rinit=10^-5;

Ru1[r_,kappa_,Ee_]:=Module[{gamma=Sqrt[kappa^2-y^2]},Sqrt[2 Pi] (Sqrt[Ee] (2 Ee r)^gamma Exp[Pi y/2]*Abs[Gamma[gamma+I y]])/(2 Sqrt[(Pi Ee)] Gamma[2 gamma+1]) (Exp[-I Ee r] (-I) Sqrt[(kappa-0 I y/Ee)/(gamma+I y)] Exp[I phi] (gamma+I y) Hypergeometric1F1[gamma+1+I y,2 gamma+1,2 I Ee r]+Exp[I Ee r] (I) Sqrt[(kappa+0 I y/Ee)/(gamma-I y)] Exp[-I phi] (gamma-I y) Hypergeometric1F1[gamma+1-I y,2 gamma+1,-2 I Ee r])];
Ru2[r_,kappa_,Ee_]:=Module[{gamma=Sqrt[kappa^2-y^2]},Sqrt[2 Pi] (I Sqrt[Ee] (2 Ee r)^gamma Exp[Pi y/2] Abs[Gamma[gamma+I y]])/(2 Sqrt[(Pi Ee)] Gamma[2 gamma+1]) (Exp[-I Ee r] (-I) Sqrt[(kappa-0 I y/Ee)/(gamma+I y)] Exp[I phi] (gamma+I y) Hypergeometric1F1[gamma+1+I y,2 gamma+1,2 I Ee r]-Exp[I Ee r] (I) Sqrt[(kappa+0 I y/Ee)/(gamma-I y)] Exp[-I phi] (gamma-I y) Hypergeometric1F1[gamma+1-I y,2 gamma+1,-2 I Ee r])];

u1Electron[r_,kappa_,Ee_]:=Ru1[r,kappa,Ee]/.phi->0;
u2Electron[r_,kappa_,Ee_]:=Ru2[r,kappa,Ee]/.phi->0;
G[r_]:=(2^Sqrt[1-Zal^2] E^(-m r Zal) Zal (m r Zal)^Sqrt[1-Zal^2] Sqrt[m (1+Sqrt[1-Zal^2])])/Sqrt[Zal Gamma[1+Sqrt[4-4 Zal^2]]]/r;
F[r_]:=-((2^Sqrt[1-Zal^2] E^(-m r Zal) Zal (m r Zal)^Sqrt[1-Zal^2] Sqrt[m-m Sqrt[1-Zal^2]])/Sqrt[Zal Gamma[1+Sqrt[4-4 Zal^2]]])/r;
ge[r_,kappa_,Ee_]:=Evaluate[u1Electron[r,kappa,Ee]/r];
fe[r_,kappa_,Ee_]:=Evaluate[u2Electron[r,kappa,Ee]/r];
geI[k_,Ee_]:=geI[k,Ee]=Evaluate[Interpolation[Evaluate[Table[{r,ge[r,k,Ee]},{r,10^-6,NumBR BohrRadius,0.05}]]]];
feI[k_,Ee_]:=feI[k,Ee]=Evaluate[Interpolation[Evaluate[Table[{r,fe[r,k,Ee]},{r,10^-6,NumBR BohrRadius,0.05}]]]];

S0[k_?NumericQ,kappa_?NumericQ,J_?NumericQ,Ee_?NumericQ]:=Module[{jk=Abs[kappa]-1/2,lk},
lk=jk+1/2 kappa/Abs[kappa];1/2*NIntegrate[r^2 SphericalBesselJ[J,k r] ((1+(-1)^(lk+J)) (1+kappa) (geI[kappa,Ee][r] G[r]-feI[kappa,Ee][r] F[r])+I (1-(-1)^(lk+J)) (1-kappa) (feI[kappa,Ee][r] G[r]+geI[kappa,Ee][r] F[r])),{r,rinit,NumBR BohrRadius}]];
S1[k_?NumericQ,kappa_?NumericQ,J_?NumericQ,Ee_?NumericQ]:=Module[{jk=Abs[kappa]-1/2,lk},
lk=jk+1/2 kappa/Abs[kappa];1/2*NIntegrate[r^2 SphericalBesselJ[J+1,k r] ((1+(-1)^(lk+J)) ((J+kappa+2) geI[kappa,Ee][r] F[r]+(kappa-J) feI[kappa,Ee][r] G[r])+I (1-(-1)^(lk+J)) ((J-kappa+2) feI[kappa,Ee][r] F[r]+(J+kappa) geI[kappa,Ee][r] G[r])),{r,rinit,NumBR BohrRadius}]];
Sm1[k_?NumericQ,kappa_?NumericQ,J_?NumericQ,Ee_?NumericQ]:=Module[{jk=Abs[kappa]-1/2,lk},
lk=jk+1/2 kappa/Abs[kappa];1/2*NIntegrate[r^2 SphericalBesselJ[J-1,k r] ((1+(-1)^(lk+J)) ((1+kappa-J) geI[kappa,Ee][r] F[r]+(J+kappa+1) feI[kappa,Ee][r] G[r])+I (1-(-1)^(lk+J)) ((1-kappa-J) feI[kappa,Ee][r] F[r]-(J-kappa+1) geI[kappa,Ee][r] G[r])),{r,rinit,NumBR BohrRadius}]];
S4[k_?NumericQ,kappa_?NumericQ,J_?NumericQ,Ee_?NumericQ]:=Module[{jk=Abs[kappa]-1/2,lk},
lk=jk+1/2 kappa/Abs[kappa];1/2*NIntegrate[r^2 SphericalBesselJ[J,k r] ((1+(-1)^(lk+J)) (geI[kappa,Ee][r] G[r]+feI[kappa,Ee][r] F[r])+I (1-(-1)^(lk+J)) (feI[kappa,Ee][r] G[r]-geI[kappa,Ee][r] F[r])),{r,rinit,NumBR BohrRadius}]];
NEe[Ee_,J_,kappa_]:=Module[{jk=Abs[kappa]-1/2,lk},
lk=jk+1/2 kappa/Abs[kappa];
16/(\[Pi] m^5)*Ee^2 (2 jk+1)*NIntegrate[kk^2 (((Efm-Ee)^2-kk^2) ((Abs[S0[kk,kappa,J,Ee] ]^2)/(J (J+1))+(Abs[S1[kk,kappa,J,Ee] ]^2)/((J+1) (2 J+1))+(Abs[Sm1[kk,kappa,J,Ee]]^2 )/(J (2 J+1)))+kk^2 (Abs[S4[kk,kappa,J,Ee]]^2+1/(2 J+1)^2 (S1[kk,kappa,J,Ee]+Sm1[kk,kappa,J,Ee]) Conjugate[(S1[kk,kappa,J,Ee]+Sm1[kk,kappa,J,Ee])])+(Efm-Ee) kk (1/(2 J+1) ((Conjugate[S1[kk,kappa,J,Ee]]+Conjugate[Sm1[kk,kappa,J,Ee]]) S4[kk,kappa,J,Ee]+Conjugate[(Conjugate[S1[kk,kappa,J,Ee]]+Conjugate[Sm1[kk,kappa,J,Ee]]) S4[kk,kappa,J,Ee]]))),{kk,10^-8,Efm-Ee}]];

NEeJ0[Ee_,J_,kappa_]:=Module[{jk=Abs[kappa]-1/2,lk},
lk=jk+1/2 kappa/Abs[kappa];
16/(\[Pi] m^5)*Ee^2 (2 jk+1)*NIntegrate[kk^2*(((Efm-Ee)^2-kk^2)*((Abs[S1[kk,kappa,J,Ee] ]^2)/((J+1) (2 J+1)))+kk^2 (Abs[S4[kk,kappa,J,Ee]]^2+1/(2 J+1)^2 Abs[S1[kk,kappa,J,Ee]]^2)+(Efm-Ee) kk (1/(2 J+1) ((Conjugate[S1[kk,kappa,J,Ee]]) S4[kk,kappa,J,Ee]+Conjugate[(Conjugate[S1[kk,kappa,J,Ee]]) S4[kk,kappa,J,Ee]]))),{kk,0,Efm-Ee}]];



Emin=10^-5;
Emax=Efm;
Estep =0.005;
Jmin=1;
Jmax=31;
name="Zal.1_j";
dd=".dat";
Print["J=0"];
specJ0=ParallelTable[{Ee,1/(4 Ee^2)NEeJ0[Ee,0,-1] },{Ee,Emin,Emax,Estep}];
Export["Zal0.1_j0.dat",specJ0];
Do[Print["J=",J];
Print["kappa=",kappa];
specJ=ParallelTable[{Ee,1/(4 Ee^2)NEe[Ee,J,kappa] },{Ee,Emin,Emax,Estep}];
Export[name<>ToString[J]<>"_k"<>ToString[kappa+J]<>dd,specJ];
,{J,Jmin,Jmax,1},{kappa,-(J+1),-J,1}];


