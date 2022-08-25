%% Description
% Simulation code for the non-orthogonal slicing between eMBB and URLLC for
% the case where four active URLLC devices transmit in each minislot.
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
Hu1=raylrnd(1/sqrt(2),F,N);     % Channel coefficients - URLLC device 1
Gu1=Gamma_u*Hu1.^2;             % Channel gains - URLLC device 1
Hu2=raylrnd(1/sqrt(2),F,N);     % Channel coefficients - URLLC device 2
Gu2=Gamma_u*Hu2.^2;             % Channel gains - URLLC device 2
Hu3=raylrnd(1/sqrt(2),F,N);     % Channel coefficients - URLLC device 3
Gu3=Gamma_u*Hu3.^2;             % Channel gains - URLLC device 3
Hu4=raylrnd(1/sqrt(2),F,N);     % Channel coefficients - URLLC device 4
Gu4=Gamma_u*Hu4.^2;             % Channel gains - URLLC device 4

%% Computation of ru_orth
rb=0:0.1:round(rb_max+0.05,1);     % Espectral efficiency - eMBB devices
ru_max=zeros(1,length(rb));
parfor x=1:length(rb)
    for ru=0:0.01:5
        flag=0;
        for Gb_tar=Gb_min:0.1:Gb_tar_max
            Sigma_u1=Gu1./(1+Gu2+Gu3+Gu4+Gb_tar);          % SINR - URLLC device 1
            Sigma_u2=Gu2./(1+Gu1+Gu3+Gu4+Gb_tar);          % SINR - URLLC device 2
            Sigma_u3=Gu3./(1+Gu1+Gu2+Gu4+Gb_tar);          % SINR - URLLC device 3
            Sigma_u4=Gu4./(1+Gu1+Gu2+Gu3+Gb_tar);          % SINR - URLLC device 4
            Db=zeros(F,N);
            Du1=zeros(1,N);
            Du2=zeros(1,N);
            Du3=zeros(1,N);
            Du4=zeros(1,N);
            Pr_Eb=zeros(F,1);
            for j=1:N
                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                if sum_MutualInfo_u1>sum_MutualInfo_u2 && sum_MutualInfo_u1>sum_MutualInfo_u3 && sum_MutualInfo_u1>sum_MutualInfo_u4
                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                    if Du1(j)==1
                        
                        Sigma_u2(:,j)=Gu2(:,j)./(1+Gu3(:,j)+Gu4(:,j)+Gb_tar);
                        Sigma_u3(:,j)=Gu3(:,j)./(1+Gu2(:,j)+Gu4(:,j)+Gb_tar);
                        Sigma_u4(:,j)=Gu4(:,j)./(1+Gu2(:,j)+Gu3(:,j)+Gb_tar);
                        sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                        sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                        sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                        if sum_MutualInfo_u2>sum_MutualInfo_u3 && sum_MutualInfo_u2>sum_MutualInfo_u4
                            Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                            if Du2(j)==1
                                Sigma_u3(:,j)=Gu3(:,j)./(1+Gu4(:,j)+Gb_tar);              
                                Sigma_u4(:,j)=Gu4(:,j)./(1+Gu3(:,j)+Gb_tar); 
                                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                                if sum_MutualInfo_u3>sum_MutualInfo_u4
                                    Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                                    if Du3(j)==1
                                       Du4(j)=(1/Fu)*sum(log2(1+Gu4(:,j)./(1+Gb_tar)))>=ru;
                                       if Du4(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                                    if Du4(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                    
                        elseif sum_MutualInfo_u3>sum_MutualInfo_u2 && sum_MutualInfo_u3>sum_MutualInfo_u4
                            Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                            if Du3(j)==1
                                Sigma_u2(:,j)=Gu2(:,j)./(1+Gu4(:,j)+Gb_tar);              
                                Sigma_u4(:,j)=Gu4(:,j)./(1+Gu2(:,j)+Gb_tar);
                                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                                if sum_MutualInfo_u2>sum_MutualInfo_u4
                                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                                    if Du2(j)==1
                                       Du4(j)=(1/Fu)*sum(log2(1+Gu4(:,j)./(1+Gb_tar)))>=ru;
                                       if Du4(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                                    if Du4(j)==1
                                       Du2(j)=(1/Fu)*sum(log2(1+Gu2(:,j)./(1+Gb_tar)))>=ru;
                                       if Du2(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                      
                        else
                            Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                            if Du4(j)==1
                                Sigma_u2(:,j)=Gu2(:,j)./(1+Gu3(:,j)+Gb_tar);              
                                Sigma_u3(:,j)=Gu3(:,j)./(1+Gu2(:,j)+Gb_tar);
                                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                                if sum_MutualInfo_u2>sum_MutualInfo_u3
                                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                                    if Du2(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                                    if Du3(j)==1
                                       Du2(j)=(1/Fu)*sum(log2(1+Gu2(:,j)./(1+Gb_tar)))>=ru;
                                       if Du2(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end
                        end
                        
                    end                    
                elseif sum_MutualInfo_u2>sum_MutualInfo_u1 && sum_MutualInfo_u2>sum_MutualInfo_u3 && sum_MutualInfo_u2>sum_MutualInfo_u4
                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                    if Du2(j)==1
                        
                        Sigma_u1(:,j)=Gu1(:,j)./(1+Gu3(:,j)+Gu4(:,j)+Gb_tar);
                        Sigma_u3(:,j)=Gu3(:,j)./(1+Gu1(:,j)+Gu4(:,j)+Gb_tar);
                        Sigma_u4(:,j)=Gu4(:,j)./(1+Gu1(:,j)+Gu3(:,j)+Gb_tar);
                        sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                        sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                        sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                        if sum_MutualInfo_u1>sum_MutualInfo_u3 && sum_MutualInfo_u1>sum_MutualInfo_u4
                            Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                            if Du1(j)==1
                                Sigma_u3(:,j)=Gu3(:,j)./(1+Gu4(:,j)+Gb_tar);              
                                Sigma_u4(:,j)=Gu4(:,j)./(1+Gu3(:,j)+Gb_tar); 
                                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                                if sum_MutualInfo_u3>sum_MutualInfo_u4
                                    Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                                    if Du3(j)==1
                                       Du4(j)=(1/Fu)*sum(log2(1+Gu4(:,j)./(1+Gb_tar)))>=ru;
                                       if Du4(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                                    if Du4(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                    
                        elseif sum_MutualInfo_u3>sum_MutualInfo_u1 && sum_MutualInfo_u3>sum_MutualInfo_u4
                            Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                            if Du3(j)==1
                                Sigma_u1(:,j)=Gu1(:,j)./(1+Gu4(:,j)+Gb_tar);              
                                Sigma_u4(:,j)=Gu4(:,j)./(1+Gu1(:,j)+Gb_tar);
                                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                                if sum_MutualInfo_u1>sum_MutualInfo_u4
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du4(j)=(1/Fu)*sum(log2(1+Gu4(:,j)./(1+Gb_tar)))>=ru;
                                       if Du4(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                                    if Du4(j)==1
                                       Du1(j)=(1/Fu)*sum(log2(1+Gu1(:,j)./(1+Gb_tar)))>=ru;
                                       if Du1(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                      
                        else
                            Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                            if Du4(j)==1
                                Sigma_u1(:,j)=Gu1(:,j)./(1+Gu3(:,j)+Gb_tar);              
                                Sigma_u3(:,j)=Gu3(:,j)./(1+Gu1(:,j)+Gb_tar);
                                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                                if sum_MutualInfo_u1>sum_MutualInfo_u3
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                                    if Du3(j)==1
                                       Du1(j)=(1/Fu)*sum(log2(1+Gu1(:,j)./(1+Gb_tar)))>=ru;
                                       if Du1(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end
                        end
                        
                    end
                elseif sum_MutualInfo_u3>sum_MutualInfo_u1 && sum_MutualInfo_u3>sum_MutualInfo_u2 && sum_MutualInfo_u3>sum_MutualInfo_u4
                    Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                    if Du3(j)==1
                        Sigma_u1(:,j)=Gu1(:,j)./(1+Gu2(:,j)+Gu4(:,j)+Gb_tar);
                        Sigma_u2(:,j)=Gu2(:,j)./(1+Gu1(:,j)+Gu4(:,j)+Gb_tar);
                        Sigma_u4(:,j)=Gu4(:,j)./(1+Gu1(:,j)+Gu2(:,j)+Gb_tar);
                        sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                        sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                        sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                        if sum_MutualInfo_u1>sum_MutualInfo_u2 && sum_MutualInfo_u1>sum_MutualInfo_u4
                            Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                            if Du1(j)==1
                                Sigma_u2(:,j)=Gu2(:,j)./(1+Gu4(:,j)+Gb_tar);              
                                Sigma_u4(:,j)=Gu4(:,j)./(1+Gu2(:,j)+Gb_tar); 
                                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                                if sum_MutualInfo_u2>sum_MutualInfo_u4
                                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                                    if Du2(j)==1
                                       Du4(j)=(1/Fu)*sum(log2(1+Gu4(:,j)./(1+Gb_tar)))>=ru;
                                       if Du4(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                                    if Du4(j)==1
                                       Du2(j)=(1/Fu)*sum(log2(1+Gu2(:,j)./(1+Gb_tar)))>=ru;
                                       if Du2(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                    
                        elseif sum_MutualInfo_u2>sum_MutualInfo_u1 && sum_MutualInfo_u2>sum_MutualInfo_u4
                            Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                            if Du2(j)==1
                                Sigma_u1(:,j)=Gu1(:,j)./(1+Gu4(:,j)+Gb_tar);              
                                Sigma_u4(:,j)=Gu4(:,j)./(1+Gu1(:,j)+Gb_tar);
                                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                                sum_MutualInfo_u4=sum(log2(1+Sigma_u4(:,j)));
                                if sum_MutualInfo_u1>sum_MutualInfo_u4
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du4(j)=(1/Fu)*sum(log2(1+Gu4(:,j)./(1+Gb_tar)))>=ru;
                                       if Du4(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                                    if Du4(j)==1
                                       Du1(j)=(1/Fu)*sum(log2(1+Gu1(:,j)./(1+Gb_tar)))>=ru;
                                       if Du1(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                      
                        else
                            Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                            if Du4(j)==1
                                Sigma_u1(:,j)=Gu1(:,j)./(1+Gu2(:,j)+Gb_tar);              
                                Sigma_u2(:,j)=Gu2(:,j)./(1+Gu1(:,j)+Gb_tar);
                                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                                if sum_MutualInfo_u1>sum_MutualInfo_u2
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du2(j)=(1/Fu)*sum(log2(1+Gu2(:,j)./(1+Gb_tar)))>=ru;
                                       if Du2(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                                    if Du2(j)==1
                                       Du1(j)=(1/Fu)*sum(log2(1+Gu1(:,j)./(1+Gb_tar)))>=ru;
                                       if Du1(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end
                        end   
                    end
                else
                    Du4(j)=(1/Fu)*sum(log2(1+Sigma_u4(:,j)))>=ru;
                    if Du4(j)==1
                        Sigma_u1(:,j)=Gu1(:,j)./(1+Gu2(:,j)+Gu3(:,j)+Gb_tar);
                        Sigma_u2(:,j)=Gu2(:,j)./(1+Gu1(:,j)+Gu3(:,j)+Gb_tar);
                        Sigma_u3(:,j)=Gu3(:,j)./(1+Gu1(:,j)+Gu2(:,j)+Gb_tar);
                        sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                        sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                        sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                        if sum_MutualInfo_u1>sum_MutualInfo_u2 && sum_MutualInfo_u1>sum_MutualInfo_u3
                            Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                            if Du1(j)==1
                                Sigma_u2(:,j)=Gu2(:,j)./(1+Gu3(:,j)+Gb_tar);              
                                Sigma_u3(:,j)=Gu3(:,j)./(1+Gu2(:,j)+Gb_tar); 
                                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                                if sum_MutualInfo_u2>sum_MutualInfo_u3
                                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                                    if Du2(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                                    if Du3(j)==1
                                       Du2(j)=(1/Fu)*sum(log2(1+Gu2(:,j)./(1+Gb_tar)))>=ru;
                                       if Du2(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                    
                        elseif sum_MutualInfo_u2>sum_MutualInfo_u1 && sum_MutualInfo_u2>sum_MutualInfo_u3
                            Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                            if Du2(j)==1
                                Sigma_u1(:,j)=Gu1(:,j)./(1+Gu3(:,j)+Gb_tar);              
                                Sigma_u3(:,j)=Gu3(:,j)./(1+Gu1(:,j)+Gb_tar);
                                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                                sum_MutualInfo_u3=sum(log2(1+Sigma_u3(:,j)));
                                if sum_MutualInfo_u1>sum_MutualInfo_u3
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du3(j)=(1/Fu)*sum(log2(1+Gu3(:,j)./(1+Gb_tar)))>=ru;
                                       if Du3(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end                      
                        else
                            Du3(j)=(1/Fu)*sum(log2(1+Sigma_u3(:,j)))>=ru;
                            if Du3(j)==1
                                Sigma_u1(:,j)=Gu1(:,j)./(1+Gu2(:,j)+Gb_tar);              
                                Sigma_u2(:,j)=Gu2(:,j)./(1+Gu1(:,j)+Gb_tar);
                                sum_MutualInfo_u1=sum(log2(1+Sigma_u1(:,j)));
                                sum_MutualInfo_u2=sum(log2(1+Sigma_u2(:,j)));
                                if sum_MutualInfo_u1>sum_MutualInfo_u2
                                    Du1(j)=(1/Fu)*sum(log2(1+Sigma_u1(:,j)))>=ru;
                                    if Du1(j)==1
                                       Du2(j)=(1/Fu)*sum(log2(1+Gu2(:,j)./(1+Gb_tar)))>=ru;
                                       if Du2(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                else
                                    Du2(j)=(1/Fu)*sum(log2(1+Sigma_u2(:,j)))>=ru;
                                    if Du2(j)==1
                                       Du1(j)=(1/Fu)*sum(log2(1+Gu1(:,j)./(1+Gb_tar)))>=ru;
                                       if Du1(j)==1
                                           Db(:,j)=log2(1+Gb_tar)>=rb(x);
                                       end
                                    end
                                end
                            end
                        end                                                 
                    end
                end
                                                                                                               
            end
            Pr_Eu1=1-mean(Du1);
            Pr_Eu2=1-mean(Du2);
            Pr_Eu3=1-mean(Du3);
            Pr_Eu4=1-mean(Du4);
            for i=1:F
                Pr_Eb(i)=1-mean(Db(i,:));
            end
            if Pr_Eu1<=Eu && Pr_Eu2<=Eu && Pr_Eu3<=Eu && Pr_Eu4<=Eu && sum(Pr_Eb(:)<=Eb)==F
                ru_max(x)=ru;
                flag=1;
                break;
            end
        end
        if flag==0
           break; 
        end
    end
end

%% Computation of rb_sum and ru_sum
rb_sum_4_Non=F*rb;
ru_sum_4_Non=4*ru_max;

%% Plotting the curves
figure(1)
    plot(rb_sum_4,ru_sum_4,'LineWidth',1.5)
    grid on    
    xlabel('$r_{B,sum}$','Interpreter','latex','FontSize',12)
    ylabel('$r_{U,sum}$','Interpreter','latex','FontSize',12)
    title('$\Gamma_B = 20$ dB, $\Gamma_U = 10$ dB, $\epsilon_B = 10^{-3}$, $\epsilon_U = 10^{-5}$',...
        'Interpreter','latex','FontSize',14)

%% Saving the results
save('Result_4_Non.mat','rb_sum_4_Non','ru_sum_4_Non')
