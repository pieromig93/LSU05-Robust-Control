clc 
clear all

%% DEFINIZIONE DEL SISTEMA: dinamica laterale
%vettore dello stato laterale [beta, p, r, phi]
%vettore degli ingressi [deltaA, deltaR], ingressi sono alettone e timone.

[P1,P2,P3] = createSystems();

%% SPECIFICHE DA RISPETTARE
wn = 2; zita = 0.5;
Mp = exp(-pi*zita/(sqrt(1-zita^2)));
eps = 0.01;
Ld = tf(wn^2, conv([1 eps],[1 2*zita*wn]))*eye(2); Ld = minreal(Ld);
I = eye(size(Ld));
Sd = feedback(I,Ld); Sd=minreal(Sd);
Td = I-Sd; Td = minreal(Td);
ts = 4/(zita*wn);
ref = tf(1);
%% %% DEFINIZIONE DELL'INCERTEZZA: incertezza additiva


%% DEFINIZIONE DELL'INCERTEZZA: incertezza moltiplicativa
%Secondo step

%sigma((P2-P1)/P1, (P3-P1)/P1, (P3-P2)/P2);
%max si ha con P2
Wm = minreal((P2-P1)/P1);
W1 = inv(Sd)*diag([.5,.5]);
%W1 = tf(1, [1 .01])*I;

G1 = [zeros(2) Wm*P1;-I -P1]; G1 = minreal(G1);
G2 = [W1 -W1 minreal(-W1*P1); zeros(size(I)) zeros(size(I)) Wm*P1; eye(size(I)) -eye(size(I)) -P1]; G2 = minreal(G2);

K1 = minreal(h2syn(G1,2,2)); K2 = minreal(hinfsyn(G1,2,2));
K3 = minreal(h2syn(G2,2,2)); K4 = minreal(hinfsyn(G2,2,2));

Twz1 = minreal(lft(G1,K1,2,2)); norm(Twz1,inf);
Twz2 = minreal(lft(G1,K2,2,2)); norm(Twz2,inf);
Twz3 = minreal(lft(G2,K3,2,2)); norm(Twz3,inf);
Twz4 = minreal(lft(G2,K4,2,2)); Twz4n = norm(Twz4,inf);


% Con il K4 e le funzioni di peso scelte troviamo che gli impianti risultao
% essere tutti e 3 stabili. Inoltre con questo controllore la norma Twz è
% minore di zero quindi risulta rispettata la condizione di stabilità
% robusta.

% Il prossimo step è limitare anche gli ingressi trovando un'altra
% fuinzione di peso da posizionare prima degli impianti che limita
% l'ingresso. 

% Angolo di sideslip max 2 gradi e angolo di rollio max 30 gradi.
% Modelliamo ora la funzione di peso relativa agli ingressi Wu;

Wu = tf(.1,[1 10])*diag([1, 1]);

G3 = minreal([W1 -W1 minreal(-W1*P1);zeros(size(I)) zeros(size(I)) Wu; zeros(size(I)) zeros(size(I)) minreal(Wm*P1);I -I -P1]);G3 = minreal(G3);
K5 = minreal(h2syn(G3,2,2)); K6 = minreal(hinfsyn(G3,2,2));

Twz5 = minreal(lft(G3,K5,2,2)); norm(Twz5,inf);
Twz6 = minreal(lft(G3,K6,2,2)); norm(Twz6,inf);

%% REIEZIONE DEL DISTURBO
%costruiamo ora un controllore in grado di reiettare i disturbi, definiamo
%una matrice G4

G4 = [W1 -W1 -W1 minreal(-W1*P1); zeros(2) zeros(2) zeros(2) minreal(Wm*P1); eye(2) -eye(2) -eye(2) -P1];
G4 = minreal(G4);

K7 = minreal(h2syn(G4,2,2)); K8 = minreal(hinfsyn(G4,2,2));
