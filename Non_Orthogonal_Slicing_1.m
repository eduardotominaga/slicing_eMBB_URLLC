%% Description
% Simulation code for the non-orthogonal slicing between eMBB and URLLC for
% the case where only one active URLLC devices transmits in each minislot.
% This script plots the eMBB sum rate (in bits/s/Hz) versus the URLLC sum
% rate (in bits/s/Hz) given the reliability requirements for eMBB and URLLC
% and the average received SNRs for the eMBB and URLLC devices.

%% Reset
clearvars
close all
clc

%% Parameters
N=1e6;                          % Number of Monte Carlo realizations
Gamma_u_dB=10;                  % Average channel gain of the URLLC devices [dB]
Gamma_u=10^(Gamma_u_dB/10);     % Average channel gain of the URLLC devices (linear scale)
Gamma_b_dB=20;                  % Average channel gain of the eMBB devices [dB]
Gamma_b=10^(Gamma_b_dB/10);     % Average channel gain of the eMBB devices (linear scale)
Eb=1e-3;                        % Reliability requirement for eMBB
Eu=1e-5;                        % Reliability requirement for URLLC
F=10;                           % Total number of frequency channels
Fu=F;                           % Number of channels allocated for URLLC traffic   

%% eMBB Parameters
Gb_min=Gamma_b*log(1/(1-Eb));               % Threshold SNR
Gb_tar_max=Gamma_b/expint(Gb_min/Gamma_b);  % Maximum Target SNR
rb_max=log2(1+Gb_tar_max);                  % Maximum eMBB data rate [bits/s/Hz]

%% Channel Realizations
Hu=raylrnd(1/sqrt(2),Fu,N);           % Channel coefficients - URLLC device
Gu=Gamma_u*Hu.^2;                     % Channel gains - URLLC device

%% H-NOMA with SIC Decoding - 1 eMBB device and 1 URLLC device
rb=0:0.1:round(rb_max+0.05,1);
ru_max=zeros(1,length(rb));
parfor x=1:length(rb)
   for ru=0:0.01:4
      flag=0;
      for Gb_tar=Gb_min:0.1:Gb_tar_max
         Sigma_u=Gu./(1+Gb_tar);
         Db=zeros(F,N);
         Du=zeros(1,N);
         Pr_Eb=zeros(F,1);
         for j=1:N
            Du(j)=(1/Fu)*sum(log2(1+Sigma_u(:,j)))>=ru;
            if Du(j)==1
                Db(:,j)=log2(1+Gb_tar)>=rb(x);
            else
                for i=1:F
                    Db(i,j)=log2(1+Gb_tar/(1+Gu(i,j)))>=rb(x);
                end
            end
         end
         Pr_Eu=1-mean(Du);
         for i=1:F
            Pr_Eb(i)=1-mean(Db(i,:));
         end
         if Pr_Eu<=Eu && sum(Pr_Eb(:)<=Eb)==F
            ru_max(x)=ru;
            Gb_tar_required(x)=Gb_tar;
            flag=1;
            break;
         end
      end
      if flag==0
           break; 
      end
   end
end

%% Computation of ru_sum and rb_sum
ru_sum_1_Non=ru_max;
rb_sum_1_Non=F*rb;

%% Plotting the curve
figure(1)    
    plot(rb_sum_1,ru_sum_1,'LineWidth',1.5)
    grid on
    xlabel('$r_{B,sum}$ [bits/s/Hz]','Interpreter','latex','fontsize',12)
    ylabel('$r_{U,sum}$ [bits/s/Hz]','Interpreter','latex','fontsize',12)
    leg=legend('1 URLLC device/minislot','2 URLLC devices/minislot','3 URLLC devices/minislot');
    set(leg,'color','none','Interpreter','latex','fontsize',10)
    title('$\Gamma_B = 20$ dB, $\Gamma_U = 10$ dB, $\epsilon_B = 10^{-3}$, $\epsilon_U = 10^{-5}$','Interpreter','latex','fontsize',14)
    
%% Saving the results
save('Results_1_Non.mat','rb_sum_1_Non','ru_sum_1_Non')
