clc 
clear all

%% DEFINIZIONE DEL SISTEMA: dinamica laterale
%vettore dello stato laterale [beta, p, r, phi]
%vettore degli ingressi [deltaA, deltaR], ingressi sono alettone e timone.

A = [-2.3817 0 -1.0019 2.1827;
        -21.063 -16.055 0.87229 0;
        24.512 -16.651 -3.5379 0;
        0 1.0026 -0.029766 0];
    
%deltaA come ingresso, ovvero l'inclinazione dell'alettone    
B = [0 -36.263 -0.67252 0;
    -0.24719 -688.44 -67.983 0]';

%mediante la Glat decidiamo su quali paramentri agire della dinamica
%laterale
C = [1 0 0 0;0 0 0 1];

D = [0 0;0 0];

P = ss(A,B,C,D);
Px = ss(A,B,eye(4),zeros(4,2));
%clear A B C D;

%% SPECIFICHE DA RISPETTARE
wn = 2; zita = 0.5;
wb = wn*sqrt(2)/10;
eu = 0.1;
Mp = exp(-pi*zita/(sqrt(1-zita^2)));
eps = 0.01;
Ld = tf(wn^2, conv([1 eps],[1 2*zita*wn]))*eye(2); Ld = minreal(Ld);
I = eye(size(Ld));
Sd = feedback(I,Ld); Sd=minreal(Sd);
Td = I-Sd; Td = minreal(Td);
ts = 4/(zita*wn);

%% SCELTA DELLA FUNZIONE DI PESO e calcolo controllore

ref = tf(1);

W1 = minreal(1/Sd)*I; 

G1 = [W1 W1*P; -I -P];
K1 = minreal(hinfsyn(G1,2,2));
% aggiungere funzione di peso sugli ingressi per impedire gli ingressi di
% cresce troppo e di diventare non reali.
%Wu = tf(.1,[1 .5])*I;%primo aileron, secondo rudder
Mu = 1;
Wu = tf([1 wb/Mu], [eu wb])*diag([10,1]);

G2 = [W1 minreal(-W1*P); zeros(2) Wu;I -P]; G2 = minreal(G2);
K3 = minreal(h2syn(G2,2,2)); K4 = minreal(hinfsyn(G2,2,2));
Twz = minreal(lft(G2,K4,2,2)); Twzn = norm(Twz,inf)

%% L,S,T con i controllori K
[L1,S1,T1]= controlStabs(P,K1);
[L4,S4,T4]= controlStabs(P,K4);

