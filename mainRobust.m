clc 
clear all

%% DEFINIZIONE DEL SISTEMA: dinamica laterale
%vettore dello stato laterale [beta, p, r, phi]
%vettore degli ingressi [deltaA, deltaR], ingressi sono alettone e timone.

A1 = [-2.3817 0 -1.0019 2.1827;
        -21.063 -16.055 0.87229 0;
        24.512 -16.651 -3.5379 0;
        0 1.0026 -0.029766 0];

%A2 dove U0 = U0-20%*U0    
A2 = [-2.3817 0 -0.8 2.1827;
        -21.063 -16.055 0.87229 0;
        24.512 -16.651 -3.5379 0;
        0 1.0026 -0.029766 0];

%A3 dove U0 = U0+20%*U0      
A3 = [-2.3817 0 -1.20 2.1827;
        -21.063 -16.055 0.87229 0;
        24.512 -16.651 -3.5379 0;
        0 1.0026 -0.029766 0];    
      
B = [0 -36.263 -0.67252 0;
    -0.24719 -688.44 -67.983 0]';

C = [1 0 0 0;0 0 0 1];

D = [0 0;0 0];

P1 = minreal(ss(A1,B,C,D));
P2 = minreal(ss(A2,B,C,D));
P3 = minreal(ss(A3,B,C,D));
clear A1 A2 A3 B C D;

P1.StateName = {'sideslip angle';'roll rate';'yaw rate';'roll angle'};
P1.InputName = {'aileron';'rudder'};
P1.OutputName = {'sideslip angle';'roll angle'};

P2.StateName = {'sideslip angle';'roll rate';'yaw rate';'roll angle'};
P2.InputName = {'aileron';'rudder'};
P2.OutputName = {'sideslip angle';'roll angle'};

P3.StateName = {'sideslip angle';'roll rate';'yaw rate';'roll angle'};
P3.InputName = {'aileron';'rudder'};
P3.OutputName = {'sideslip angle';'roll angle'};

%% SPECIFICHE DA RISPETTARE
wn = 2; zita = 0.5;
wb = wn*sqrt(2);
eu = 0.1;
Mp = exp(-pi*zita/(sqrt(1-zita^2)));
eps = 0.01;
Ld = tf(wn^2, conv([1 eps],[1 2*zita*wn]))*eye(2); Ld = minreal(Ld);
I = eye(size(Ld));
Sd = feedback(I,Ld); Sd=minreal(Sd);
Td = I-Sd; Td = minreal(Td);
ts = 4/(zita*wn);

%% DEFINIZIONE DELL'INCERTEZZA
%Secondo step

%sigma((P2-P1)/P1, (P3-P1)/P1, (P3-P2)/P2);
%max si ha con P2
Wm = minreal((P3-P1)/P1);
W1 = inv(Sd)*0.5*eye(2);
%W1 = tf(1, [1 .01])*I;

G1 = [zeros(2) Wm*P1;-I -P1]; G1 = minreal(G1);
G2 = [W1 -W1 minreal(-W1*P1); zeros(size(I)) zeros(size(I)) Wm*P1; eye(size(I)) -eye(size(I)) -P1]; G2 = minreal(G2);

K1 = minreal(h2syn(G1,2,2)); K2 = minreal(hinfsyn(G1,2,2));
K3 = minreal(h2syn(G2,2,2)); K4 = minreal(hinfsyn(G2,2,2));


Twz1 = minreal(lft(G1,K1,2,2)); norm(Twz1,inf);
Twz2 = minreal(lft(G1,K2,2,2)); norm(Twz2,inf);
Twz3 = minreal(lft(G2,K3,2,2)); norm(Twz3,inf);
Twz4 = minreal(lft(G2,K4,2,2)); Twz4n = norm(Twz4,inf);


% if Twz4n<1
%     Twz4n
%     L1 = minreal(P1*K4);
%     S1 = minreal(I/(I+L1));
%     pole(S1)
%     T1 = minreal(I-S1);
%     step(T1); 
% end

% Con il K4 e le funzioni di peso scelte troviamo che gli impianti risultao
% essere tutti e 3 stabili. Inoltre con questo controllore la norma Twz è
% minore di zero quindi risulta rispettata la condizione di stabilità
% robusta.

% Il prossimo step è limitare anche gli ingressi trovando un'altra
% fuinzione di peso da posizionare prima degli impianti che limita
% l'ingresso. 

% Angolo di sideslip max 2 gradi e angolo di rollio max 30 gradi.
% Modelliamo ora la funzione di peso relativa agli ingressi Wu;

%Wu = tf([1 wb/12],[eu wb]); Wu = Wu*[1 0;0 1];
%Wu = tf(.5,[1 20])*I;
%Wu = tf(.5,[1 20])*[10 0; 0 1];
Wu = tf(.5,[1 20])*[0.005 0; 0 1];%Passa-basso

G3 = minreal([W1 -W1 minreal(-W1*P1);zeros(size(I)) zeros(size(I)) Wu; zeros(size(I)) zeros(size(I)) minreal(Wm*P1);I -I -P1]);G3 = minreal(G3);
K5 = minreal(h2syn(G3,2,2)); K6 = minreal(hinfsyn(G3,2,2));

Twz5 = minreal(lft(G3,K5,2,2)); norm(Twz5,inf)
Twz6 = minreal(lft(G3,K6,2,2)); norm(Twz6,inf)
