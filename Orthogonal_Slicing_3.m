%% Description
% Simulation code for the orthogonal slicing between eMBB and URLLC for the
% case where three active URLLC devices transmit in each minislot.
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
Hu2_max=raylrnd(1/sqrt(2),F,N);    % Channel coefficients - URLLC device 2
Gu2_max=Gamma_u*Hu2_max.^2;        % Channel gains - URLLC device 2
Hu3_max=raylrnd(1/sqrt(2),F,N);    % Channel coefficients - URLLC device 3
Gu3_max=Gamma_u*Hu3_max.^2;        % Channel gains - URLLC device 3

%% Computation of ru_orth
parfor x=1:length(Fu)
    Gu1=Gu1_max(1:Fu(x),:);
    Gu2=Gu2_max(1:Fu(x),:);
    Gu3=Gu3_max(1:Fu(x),:);
    for ru=0:0.01:5
        Du1=zeros(1,N);
        Du2=zeros(1,N);
        Du3=zeros(1,N);
        Sigma_u1=zeros(Fu(x),N);
        Sigma_u2=zeros(Fu(x),N);
        Sigma_u3=zeros(Fu(x),N);
        for j=1:N
            Sigma_u1(1:Fu(x),j)=Gu1(:,j)./(1+Gu2(:,j)+Gu3(:,j));          % SINR - URLLC device 1
            Sigma_u2(1:Fu(x),j)=Gu2(:,j)./(1+Gu1(:,j)+Gu3(:,j));          % SINR - URLLC device 2
            Sigma_u3(1:Fu(x),j)=Gu3(:,j)./(1+Gu1(:,j)+Gu2(:,j));          % SINR - URLLC device 3
            sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
            sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
            sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
            if sum_MutualInfo_u1>sum_MutualInfo_u2 && sum_MutualInfo_u1>sum_MutualInfo_u3
                Du1(j)=(1/Fu(x))*sum(log2(1+Sigma_u1(:,j)))>=ru;
                if Du1(j)==1
                    Sigma_u2(1:Fu(x),j)=Gu2(:,j)./(1+Gu3(:,j));               % SINR Update - URLLC device 2
                    Sigma_u3(1:Fu(x),j)=Gu3(:,j)./(1+Gu2(:,j));               % SINR Update - URLLC device 3
                    sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                    sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                    if sum_MutualInfo_u2>sum_MutualInfo_u3
                        Du2(j)=(1/Fu(x))*sum(log2(1+Sigma_u2(:,j)))>=ru;
                        if Du2(j)==1
                            Du3(j)=(1/Fu(x))*sum(log2(1+Gu3(:,j)))>=ru; 
                        end
                    else
                        Du3(j)=(1/Fu(x))*sum(log2(1+Sigma_u3(:,j)))>=ru;
                        if Du3(j)==1
                            Du2(j)=(1/Fu(x))*sum(log2(1+Gu2(:,j)))>=ru; 
                        end
                    end
                end
            elseif sum_MutualInfo_u2>sum_MutualInfo_u1 && sum_MutualInfo_u2>sum_MutualInfo_u3
                Du2(j)=(1/Fu(x))*sum(log2(1+Sigma_u2(:,j)))>=ru;
                if Du2(j)==1
                    Sigma_u1(1:Fu(x),j)=Gu1(:,j)./(1+Gu3(:,j));               % SINR Update - URLLC device 2
                    Sigma_u3(1:Fu(x),j)=Gu3(:,j)./(1+Gu1(:,j));               % SINR Update - URLLC device 3
                    sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                    sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                    if sum_MutualInfo_u1>sum_MutualInfo_u3
                        Du1(j)=(1/Fu(x))*sum(log2(1+Sigma_u1(:,j)))>=ru;
                        if Du1(j)==1
                            Du3(j)=(1/Fu(x))*sum(log2(1+Gu3(:,j)))>=ru; 
                        end
                    else
                        Du3(j)=(1/Fu(x))*sum(log2(1+Sigma_u3(:,j)))>=ru;
                        if Du3(j)==1
                            Du1(j)=(1/Fu(x))*sum(log2(1+Gu1(:,j)))>=ru; 
                        end
                    end
                end
            else
                Du3(j)=(1/Fu(x))*sum(log2(1+Sigma_u3(:,j)))>=ru;
                if Du3(j)==1
                    Sigma_u1(1:Fu(x),j)=Gu1(:,j)./(1+Gu2(:,j));               % SINR Update - URLLC device 2
                    Sigma_u2(1:Fu(x),j)=Gu2(:,j)./(1+Gu1(:,j));               % SINR Update - URLLC device 3
                    sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                    sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                    if sum_MutualInfo_u1>sum_MutualInfo_u2
                        Du1(j)=(1/Fu(x))*sum(log2(1+Sigma_u1(:,j)))>=ru;
                        if Du1(j)==1
                            Du2(j)=(1/Fu(x))*sum(log2(1+Gu2(:,j)))>=ru; 
                        end
                    else
                        Du2(j)=(1/Fu(x))*sum(log2(1+Sigma_u2(:,j)))>=ru;
                        if Du2(j)==1
                            Du1(j)=(1/Fu(x))*sum(log2(1+Gu1(:,j)))>=ru; 
                        end
                    end
                end
            end
            
        end            
        Pr_Eu1=1-mean(Du1);
        Pr_Eu2=1-mean(Du2);
        Pr_Eu3=1-mean(Du3);
        if Pr_Eu1<=Eu && Pr_Eu2<=Eu && Pr_Eu3<=Eu
            ru_max_3(x)=ru;
        else
           break; 
        end
    end
end

%% Computation of rb_sum and ru_sum
rb_sum_3=(F-Fu)*rb_orth;
ru_sum_3=3*ru_max_3;

%% Plotting the curves
figure(1)
    plot(rb_sum(length(rb_sum):-1:1),ru_max_3(length(ru_max_3):-1:1),'-o','LineWidth',1.5)
    grid on    
    xlabel('$r_{B,sum}$','Interpreter','latex','FontSize',12)
    ylabel('$r_{U,NOMA}$','Interpreter','latex','FontSize',12)
    title('$\Gamma_B = 20$ dB, $\Gamma_U = 10$ dB, $\epsilon_B = 10^{-3}$, $\epsilon_U = 10^{-5}$',...
        'Interpreter','latex','FontSize',14)

%% Saving the results
save('Results_1_eMBB_3_URLLC.mat','rb_sum_3','ru_sum_3')
