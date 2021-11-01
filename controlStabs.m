function [L1,S1,T1] = controlStabs(P,K)
    I = eye(2);
    L1 = P*K;
    S1 = minreal(I/(I+L1));
    pole(S1)
    T1 = I-S1;
end
