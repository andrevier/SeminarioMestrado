% Inicialização do modelo de maquina de indução no regime permanente:
% Cálculo dos estados iniciais para os Flux Linkages (concatenação do
% fluxo).
clear all
clc
% Parâmetros da máquina:
Rs = 1.77;
Rr = 1.34;
Xls = 5.25;
Xlr = 4.57;
Xm = 139.0;
J = 0.025;
p = 4;

% Regime permanente
f = 60;
VLLrms = 460;
s = 0.0172;
Wsyn = 2*pi*f;
WdA = s*Wsyn;

% Pico do fasor Va em regime permanente
Va = VLLrms*sqrt(2/3);
Vs_0 = 3/2*Va;
theta_Vs_0 = angle(Vs_0);
theta_da_0 = 0; % supor que o ângulo inicial coincide com o eixo a.

Vsd_0 = sqrt(2/3)*abs(Vs_0)*cos(theta_Vs_0 - theta_da_0);
Vsq_0 = sqrt(2/3)*abs(Vs_0)*sin(theta_Vs_0 - theta_da_0);

% indutâncias
Ls = (Xls + Xm)/Wsyn;
Lm = Xm/Wsyn;
Lr = (Xlr + Xm)/Wsyn;

M = [Ls, 0,  Lm, 0; ...
     0,  Ls, 0,  Lm; ...
     Lm, 0,  Lr, 0; ...
     0,  Lm, 0,  Lr];

Minv = inv(M);

A = [Rs,         -Wsyn*Ls,   0,         -Wsyn*Lm; ...
     Wsyn*Ls,     Rs,        Wsyn*Lm,    0      ; ...
     0,          -s*Wsyn*Lm, Rr,        -s*Wsyn*Lr; ...
     s*Wsyn*Lm,   0,         s*Wsyn*Lr,  Rr];

% Cálculo das correntes iniciais em regime permanente.
Idq_0 = A\[Vsd_0;Vsq_0;0;0];
Isd_0 = Idq_0(1)
Isq_0 = Idq_0(2)
Ird_0 = Idq_0(3)
Irq_0 = Idq_0(4)
% Cálculo da concatenação do fluxo inicial, estados iniciais, para a
% o sistema.
FL_dq = M*[Isd_0;Isq_0;Ird_0;Irq_0];
FLsd_0 = FL_dq(1)
FLsq_0 = FL_dq(2)
FLrd_0 = FL_dq(3)
FLrq_0 = FL_dq(4)

Tem_0 = (p/2)*Lm*(Isq_0*Ird_0-Isd_0*Irq_0)
TL_0 = Tem_0
Wmech_0 = (2/p)*(1-s)*Wsyn
