function [P1, P2, P3] = createSystems()

    A1 = [-2.3817 0 -1.0019 2.1827;
            -21.063 -16.055 0.87229 0;
            24.512 -16.651 -3.5379 0;
            0 1.0026 0.029766 0];
    
    B = [0 -36.263 -0.67252 0;
        -0.24719 -688.44 -67.983 0]';
    
    teta0 = atan(A1(4,3));
    
    U0n = 9.8*cos(teta0)/A1(1,4);
    U02 = U0n-(0.15*U0n);
    U03 = U0n+(0.15*U0n);
    Ysr_n = B(1,2)*U0n;
    Ysr_2 = Ysr_n/U02;
    Ysr_3 = Ysr_n/U03;
        
    %A2 dove U0 = U0-15%*U0    
    A2 = [-2.3817 0 -1.0019 9.8*cos(teta0)/U02;
            -21.063 -16.055 0.87229 0;
            24.512 -16.651 -3.5379 0;
            0 1.0026 0.029766 0];
    
    B2 = [0 -36.263 -0.67252 0;
          Ysr_2 -688.44 -67.983 0]';
    %A3 dove U0 = U0+15%*U0      
    A3 = [-2.3817 0 -1.0019 9.8*cos(teta0)/U03;
            -21.063 -16.055 0.87229 0;
            24.512 -16.651 -3.5379 0;
            0 1.0026 0.029766 0];    

    B3 = [0 -36.263 -0.67252 0;
          Ysr_3 -688.44 -67.983 0]';

    C = [1 0 0 0;0 0 0 1];

    D = [0 0;0 0];

    P1 = minreal(ss(A1,B,C,D));
    P2 = minreal(ss(A2,B2,C,D));
    P3 = minreal(ss(A3,B3,C,D));

    P1.StateName = {'sideslip angle';'roll rate';'yaw rate';'roll angle'};
    P1.InputName = {'aileron';'rudder'};
    P1.OutputName = {'sideslip angle';'roll angle'};

    P2.StateName = {'sideslip angle';'roll rate';'yaw rate';'roll angle'};
    P2.InputName = {'aileron';'rudder'};
    P2.OutputName = {'sideslip angle';'roll angle'};

    P3.StateName = {'sideslip angle';'roll rate';'yaw rate';'roll angle'};
    P3.InputName = {'aileron';'rudder'};
    P3.OutputName = {'sideslip angle';'roll angle'};
end