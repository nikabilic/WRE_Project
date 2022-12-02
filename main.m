load P.txt;
load Q_obs.txt;
load temperature.txt;
load kc.txt;
load area_rating_curve.txt;

%PARAMETERS
%hydrological model
sw=0.25;     %wilting point
s1=0.4;      %soil moisture stress threshold
n=0.3;       %porosity
Qb=7;        %m^3/s, baseflow
tsup=22;     %h, mean superficial residence time
A=4000;      %km^2, catchment area
phi=38;       %degree, latitude% 
n_years=6;  

%ARBITRARY VALUES
Ksat=5*10^-6;         %saturated hydraulic conductivity, K(s)=Ksat*s^c [m/s]
c=7;          
tsub=10*24;        %h, mean sub-superficial residence time
z=1;          %m, root zone thickness

dt=1;
doTest=1;

%RUN HYDROLOGICAL MODEL
[Q, R, I,s,L,ET] = f_hydromodel(phi,temperature,n_years,P,Ksat,sw,s1,kc,n,Qb,tsup,tsub,A,c,dt,z, doTest)

%CALIBRATION MODEL

%NS index
Q=Q';
NS_function=1-sum((Q_obs-Q).^2)/sum((Q_obs-mean(Q_obs)).^2); %[-]
%@(Q) was in the equation in front of 1????
NS_function=NS_old  %first NS 

%preallocation of vectors
Ksat=zeros(N_it,1);                 %[m/s]
c=zeros(N_it,1);                    %[-]
tsub=zeros(N_it,1);                 %[h]
z=zeros(N_it,1);                    %[mm]
NS=zeros(N_it,1);                   %[-]
%T_SA=zeros(N_it,1);                 %[-]

%initial values
Ksat(1)=1e-5;                       %[m/s]  Hydraulic conductivity for saturated soil
c(1)=5;                            %[-]    Exponent for power-law relation L(s)
tsub(1)=400;                        %[h]    Mean sub-superficial residence time
z(1)=2000;                          %[mm]   Root zone depth

%st dev
sigma_Ksat=0.05*(10^(-5)-10^(-7));  %[m/s]
sigma_c=0.05*(20-1);                %[-]
sigma_tsub=0.05*(400-1);            %[h]
sigma_z=0.05*(2000-1);              %[mm]

%run the model
%Nyears=n_years
%[Q1, R, I,soil_saturation,L,ET]= f_hydromodel(phi,temperature,Nyears,P,Ksat(1),sw,s1,kc,n,Qb,tsup,tsub(1),A,c(1),dt,z(1),doTest);

%NS old
%Q1=Q1';
%NS(1)=NS_function(Q1); %[-]
%NS(1)=1-sum((Q_obs-Q1).^2)/sum((Q_obs-mean(Q_obs)).^2);
%NS_old=NS(1);                       %[-]

%stimulated annealing temp-first T
T_SA(1)=exp(-cr*1);                 %[-]

%iteration until convergence NS>0.87 - can it be as a conditions in a loop?
N_it=100;                         %[-]    Number of iterations for the simulated annealing algorithm %should be 10 000 but because of long running time i put it to 100
cr=1/1200;       %cooling rate
Nb_NS_stored=1;                     %[-], accepted values od NS

for i=1:N_it
      T_SA(i) = exp(-cr*i);            %[-]
       Ksat_new=TruncNormRnd(Ksat(Nb_NS_stored),sigma_Ksat,10^(-7),10^(-5));   %[m/s]
       c_new=TruncNormRnd(c(Nb_NS_stored),sigma_c,1,20);                       %[-]
       tsub_new=TruncNormRnd(tsub(Nb_NS_stored),sigma_tsub,1,400);             %[h]
       z_new=TruncNormRnd(z(Nb_NS_stored),sigma_z,1,2000);                     %[mm]
       
       %run model
       [Q_new, R, I,soil_saturation,L,ET]= f_hydromodel(phi,temperature,Nyears,P,Ksat_new,sw,s1,kc,n,Qb,tsup,tsub_new,A,c_new,dt,z_new,doTest);
       
       %NS new
       Q_new=Q_new';
       %NS_new=NS_function(Q_new);      %[-]
       NS_new=1-sum((Q_obs-Q_new).^2)/sum((Q_obs-mean(Q_obs)).^2); 

       if NS_new>NS_old
            %Update index counting accepted values in preallocated vector
            Nb_NS_stored=Nb_NS_stored+1;
            %Save the current parameter set
            Ksat(Nb_NS_stored)=Ksat_new;
            c(Nb_NS_stored)=c_new;
            tsub(Nb_NS_stored)=tsub_new;
            z(Nb_NS_stored)=z_new;
            NS(Nb_NS_stored)=NS_new;
        
            NS_old=NS_new;
            
            %T_SA(Nb_NS_stored)=T_SA_1;     
       else
        %Accept with probability based on simulated annealing
  %not sure how to structure the probability condition
       if rand<exp((NS_new-NS_old)/T_SA(i))
            %Update index counting accepted values in preallocated vector
            Nb_NS_stored=Nb_NS_stored+1;
            %Save the current parameter set
            Ksat(Nb_NS_stored)=Ksat_new;
            c(Nb_NS_stored)=c_new;
            tsub(Nb_NS_stored)=tsub_new;
            z(Nb_NS_stored)=z_new;
            NS(Nb_NS_stored)=NS_new;
            
            NS_old=NS_new;
            
            %T_SA(Nb_NS_stored)=T_SA_1;
        end
       end
end

%clean empty vectors
Ksat(Nb_NS_stored+1:end)=[];
c(Nb_NS_stored+1:end)=[];
z(Nb_NS_stored+1:end)=[];
tsub(Nb_NS_stored+1:end)=[];
%T_SA(Nb_NS_stored+1:end)=[];
NS(Nb_NS_stored+1:end)=[];

save("calibration_parameter.mat","Ksat","c","tsub","z","T_SA","NS")

load("calibration_parameter.mat");

%DISCHARGE GENERATION
Nyears_gen = 100;               %[years]
max_Kc = Ksat(NS==max(NS));     %[m/s]  hydraulic conductivity (saturated soil) for the maximum Nash-Sutcliffe index NS
max_z = z(NS==max(NS));         %[mm]   root zone depth for the maximum Nash-Sutcliffe index NS
max_c = c(NS==max(NS));         %[-]    exponent for power-law relation for the maximum Nash-Sutcliffe index NS
max_tsub = tsub(NS==max(NS));   %[h]    mean sub-superficial residence time for the maximum Nash-Sutcliffe index NS

[time_gen,Q_gen,R_gen,I_gen,soil_saturation_gen,L_gen,ET_gen,P_gen_mean_monthly] = f_discharge(Nyears_gen, P, temperature, kc, dt, n, max_z,max_Kc,sw,s1,max_c,tsup,max_tsub,A,Qb,phi,doTest);

%plots
figure                                             
plot(time_gen,Q_gen);                 
ylabel('Generated hourly discharge  [m^{3}/s]','fontsize',14) 
xlabel('Time [hours]','fontsize',14)
box off

figure
plot(1:365,Q_gen_daily)
ylabel('Mean daily discharge [m^3/s]','fontsize',14)
xlabel('Days','fontsize',14)
box off

figure
plot(1:12,Q_gen_mean_monthly)
ylabel('Mean monthly discharge [m^3/s]','fontsize',14)
xlabel('Months','fontsize',14)
box off
