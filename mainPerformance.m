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
clear A B C D;

%% SPECIFICHE DA RISPETTARE
wn = 0.5; zita = 0.5;
Mp = exp(-pi*zita/(sqrt(1-zita^2)));
eps = 0.01;
Ld = tf(wn^2, conv([1 eps],[1 2*zita*wn]))*eye(2); Ld = minreal(Ld);
I = eye(size(Ld));
Sd = feedback(I,Ld); Sd=minreal(Sd);
Td = I-Sd; Td = minreal(Td);
ts = 4/(zita*wn);

%% SCELTA DELLA FUNZIONE DI PESO e calcolo controllore

ref = tf(1);

W1 = minreal(1/Sd); 

G1 = [W1 W1*P; -I -P];
K1 = minreal(hinfsyn(G1,2,2));
%K1_1 = minreal(mixsyn(P,W1,[],[]));

W3 = minreal(tf(1,[1 .01]))*I;
W3_1 = minreal(makeweight(0.5,[0.4 1],70));
W3_1 = W3_1*I;
%G2 = [W1 -W1 -W1*P; zeros(2) W3_1 W3_1*P; -W1*P -I -P];
G2 = [W1 -W1 -W1*P; zeros(2) W3 W3*P; -W1*P -I -P];
K2 = minreal(hinfsyn(G2,2,2));

%% L,S,T con i controllori K
La1 = minreal(P*K1);
Ia1 = eye(size(La1));
Sa1 = minreal(feedback(Ia1,La1)); 
Ta1 = minreal(Ia1-Sa1);

La2 = minreal(P*K2);
Ia2 = eye(size(La2));
Sa2 = minreal(feedback(Ia2,La2)); 
Ta2 = minreal(Ia2-Sa2);

