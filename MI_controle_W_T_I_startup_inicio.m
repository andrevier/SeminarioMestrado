% Inicialização do modelo de maquina de indução com controle de velocidade e de
% corrente para uma máquina partindo do zero.

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
tau_r = Lr/Rr;

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

Tem_n = (p/2)*Lm*(Isq_0*Ird_0-Isd_0*Irq_0) % torque elétrico nominal.
TL_0 = 0.0                                 % torque mecânico inicial.
Wmech_n = (2/p)*(1-s)*Wsyn                 % velocidade mecânica nominal.
Wmech_0 = 0.0                              % velocidade mecânica inicial.

% Cálculo da concatenação do fluxo inicial, estados iniciais, para a
% o sistema.
FluxLinkage_dq = M*[Isd_0;Isq_0;Ird_0;Irq_0];
FLsd_0 = FluxLinkage_dq(1);
FLsq_0 = FluxLinkage_dq(2);
FLrd_0 = FluxLinkage_dq(3);
FLrq_0 = FluxLinkage_dq(4);

[thetas, FLs_dq_0] = cart2pol(FLsd_0,FLsq_0);
[thetar, FLr_dq_0] = cart2pol(FLrd_0,FLrq_0);
[thetaIs, Is_0] = cart2pol(Isd_0,Isq_0);
[theta_Vs_dq, Vs_dq_0]=cart2pol(Vsd_0, Vsq_0)
% O eixo d é alinhado com o fluxo enlaçado pelo rotor.
FLrd_0 = FLr_dq_0
[FLsd_0, FLsq_0] = pol2cart(thetas-thetar,FLs_dq_0);
[Isd_0, Isq_0] = pol2cart(thetaIs-thetar,Is_0);
[Vsd_0, Vsq_0] = pol2cart(theta_Vs_dq - thetar, Vs_dq_0);

disp('controle de velocidade')
% Cálculo dos parâmetros do PI para a malha de velocidade.
Wc = 25                       % rad/s frequência de crossover.
k = (p/2)*(Lm^2/Lr)*Isd_0     % constante do torque.
MF = 60*pi/180;               % margem de fase em rad.

Ki = (Wc^2*J)/(k*sqrt(1+tan(MF)^2))
Kp = Ki*tan(MF)/Wc

% Conferir os parâmetros na função de transf. de malha aberta.
GL = (Kp + Ki/(j*Wc))*k/(J*j*Wc)
GLabs = abs(GL)
GLang = angle(GL)*180/pi

disp('controle de torque')
% Cálculo dos parâmetros do PI para a malha de torque.
sigma = 1 - Lm^2/(Lr*Ls);
Wci = 10*Wc
MFi = 60*pi/180;
tanMFi = tan(MFi - pi);      % Fórmula da margem de fase

Kii = Wci*(Wci*sigma*Ls + Rs*tanMFi)/sqrt(tanMFi^2 + 1)
Kpi = (Kii/Wci)*(Wci*sigma*Ls*tanMFi - Rs)/(Wci*sigma*Ls + Rs*tanMFi)

% Conferir a função de transferência:
GL_i = (Kpi + Kii/(j*Wci))/(Rs + (j*Wci)*sigma*Ls)
GLabs_i = abs(GL_i)
GLang_i = angle(GL_i)*180/pi

disp('controle de tensão')
% Cálculo dos parâmetros do PI para a malha de tensão.
% máximo comprimento de banda considerando a constante de tempo do conversor 4ms.
Wcii = 2*pi/(4*10^-3)          
MFii = 60*pi/180;
tanMFii = tan(MFii - pi);      % Fórmula da margem de fase

Kiii = Wcii*(Wcii*sigma*Ls + Rs*tanMFii)/sqrt(tanMFii^2 + 1)
Kpii = (Kiii/Wcii)*(Wcii*sigma*Ls*tanMFii - Rs)/(Wcii*sigma*Ls + Rs*tanMFii)

% Conferir a função de transferência:
GL_ii = (Kpii + Kiii/(j*Wcii))/(Rs + (j*Wcii)*sigma*Ls)
GLabs_ii = abs(GL_ii)
GLang_ii = angle(GL_ii)*180/pi

disp('controle do fluxo enlaçado pelo rotor')
% Cálculo dos parâmetros do PI para a malha do fluxo enlaçado pelo rotor.
Wcf = Wci                     % rad/s frequência de crossover.
MFf = 60*pi/180;               % margem de fase em rad.
Kif = Kii
Kpf = Kpi

% Executar a simulação
tempo = 4.5;
sim('MI_controle_W_T_I_startup',tempo)

%gráficos
figure(1)
subplot(311)
plot(torque)
title('Torque da carga x Torque eletrico PU.')

subplot(312)
plot(Velocidade_mec)
title('Velocidade mecanica de referencia x medida PU.')

subplot(313)
plot(Fluxo_r)
title('Fluxo enlaçado pelo rotor PU.')