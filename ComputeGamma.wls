#!/usr/bin/env wolframscript
(* ::Package:: *)

k12=6000
KM=0.4
{k17,k18}={10,0.0005}
{k7,k8,k19}={6,6,6}
{K5,K9}={0.1,0.0001}

Stot=33.3
Sp=5
Su=Stot-Sp

Ca=0.5
CaMtot=10
{L1,L2,L3,L4}={0.1,0.025,0.32,0.40}

p=0.1


CaM=CaMtot/(1 + L4/Ca + L3*L4/(Ca^2) + L2*L3*L4/(Ca^3) + L1*L2*L3*L4/(Ca^4))


qa[x_]:=k12*k17*K9*Sp*(1-x)
qb[x_]:=-CaM*k12*k17*K9 + CaM*k12*k17*K9*x + k17*K5*K9*x - k12*k18*K5*Sp*x + k12*k17*K9*Su*x - k12*k17*K9*Su*x^2
qc[x_]:=k12*k18*K5*x*(CaM-x*Sp)
dlt[x_]:=Sqrt[qb[x]^2 - 4*qa[x]*qc[x]]
q[x_]:=(-qb[x]+dlt[x])/(2*qa[x])

nrg[x_]:=k18*K5*x*((KM+Su+Su*(1-x)+Sp*(1-q[x]))*(CaM-x*Su-q[x]*Sp)-K5*Su*x-K9*Sp*q[x])+k12*p*(CaM-x*Su-q[x]*Sp)*(K5*x-(1-x)*(CaM-x*Su-q[x]*Sp))


Plot[q[x],{x,0,1}]


Plot[nrg[x],{x,0,0.15}]


(* ::InheritFromParent:: *)
(**)





gamma=x/.Solve[nrg[x]==0,x][[4]]
gammaSt=q[gam]
zeta=1-gamma-K5*gamma/(CaM-gamma*Su-gammaSt*Sp)
zetaSt=1-gammaSt-K9*gammaSt/(CaM-gamma*Su-gammaSt*Sp)



