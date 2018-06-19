#!/usr/bin/env wolframscript
(* ::Package:: *)

k12=6000
k17=10
k18=0.0005
K5=0.1
K=0.4
Stot=33.3


eps=k17/k18
eps18=(k18*K5)/(C-S+K5)
eps17=(k17*K5)/(C-S+K5)

pol[x_] := K+S + (eps*S + (eps-1)*(K+S) + (6*Stot-S))*x + (6*Stot-S)*(eps-1)*x^2
Solve[(1-Zeta)*pol[Zeta]-(k12*p/eps18)*(Zeta+(eps-1)*Zeta^2)==0, Zeta]



