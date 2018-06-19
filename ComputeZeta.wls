#!/usr/bin/env wolframscript
(* ::Package:: *)

eps=(k17*K9*(C-S+K5))/(k19*K5*(C-S+K9))
eps18=(k18*K5)/(C-S+K5)
eps17=(k17*K9)/(C-S+K9)

pol(_x)=K+S + (eps*S + (eps-1)*(K+S) + (6*Stot-S))*x + (6*Stot-S)*(eps-1)*x^2
Solve[(1-Zeta)*pol(Zeta)-(k12*p/eps18)*(Zeta+(eps-1)*Zeta^2)==0, Zeta]



