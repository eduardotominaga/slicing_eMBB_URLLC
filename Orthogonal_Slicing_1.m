%% Description
% Simulation code for the orthogonal slicing between eMBB and URLLC for the
% case where only one active URLLC devices transmits in each minislot.
% This script plots the eMBB sum rate (in bits/s/Hz) versus the URLLC sum
% rate (in bits/s/Hz) given the reliability requirements for eMBB and URLLC
% and the average received SNRs for the eMBB and URLLC devices.

%% Reset
clearvars
close all
clc

%% Parameters
N=1e7;                          % Number of Monte Carlo realizations
Gamma_u_dB=10;                  % Average received SNR of the URLLC devices [dB]
Gamma_u=10^(Gamma_u_dB/10);     % Average received SNR of the URLLC devices (linear scale)
Gamma_b_dB=20;                  % Average received SNR of the eMBB devices [dB]
Gamma_b=10^(Gamma_b_dB/10);     % Average received SNR of the eMBB devices (linear scale)
Eu=1e-5;                        % Reliability requirement for URLLC
Eb=1e-3;                        % Reliability requirement for eMBB
F=10;                           % Total number of frequency channels
Fu=0:F;                         % Number of frequency channels allocated for URLLC traffic

%% eMBB Parameters
Gb_min=Gamma_b*log(1/(1-Eb));           % Threshold SNR
Gb_tar=Gamma_b/expint(Gb_min/Gamma_b);  % Target SNR
rb_orth=log2(1+Gb_tar);                 % eMBB orthogonal data rate [bits/s/Hz]

%% Channel Realizations
Hu1_max=raylrnd(1/sqrt(2),F,N);    % Channel coefficients - URLLC device 1
Gu1_max=Gamma_u*Hu1_max.^2;        % Channel gains - URLLC device 1

%% Computation of ru for one URLLC device transmitting in a minislot
parfor x=1:length(Fu)
    Gu1=Gu1_max(1:Fu(x),:);
    for ru=0:0.01:5
        Du=zeros(1,N);
        for j=1:N
           Du(j)=(1/Fu(x))*sum(log2(1+Gu1(:,j)))>=ru;
        end
        Pr_Eu=1-mean(Du);
        if Pr_Eu<=Eu
           ru_max(x)=ru; 
        else
            break;
        end
    end
end

%% Computation of the eMBB sum rate and the URLLC sum rate
rb_sum_1=(F-Fu)*rb_orth;
ru_sum_1=ru_max;

%% Plotting the curves
figure(1)
    plot(rb_sum_1(length(rb_sum_1):-1:1),ru_sum_1(length(ru_max):-1:1),'-o','LineWidth',1.5)
    grid on
    xlabel('$r_{B,sum}$','Interpreter','latex','FontSize',12)
    ylabel('$r_{U,NOMA}$','Interpreter','latex','FontSize',12)
    title('$\Gamma_B = 10$ dB, $\Gamma_U = 20$ dB, $\epsilon_B = 10^{-3}$, $\epsilon_U = 10^{-5}$',...
        'Interpreter','latex','FontSize',14)
    leg=legend('1 URLLC device/minislot','2 URLLC devices/minislot');
    set(leg,'Interpreter','latex','FontSize',12)

%% Saving the results
save('Results_1_eMBB_1_URLLC.mat','rb_sum_1','ru_sum_1')
