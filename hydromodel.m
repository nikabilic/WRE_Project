function [Q, R, I,s,L,ET] = hydromodel(phi,temperature,n_years,P,Ksat,sw,s1,kc,n,Qb,tsup,tsub,A,c,dt,z, doTest)

% Day/month parameters
day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; 
month_end=cumsum(day_month);                     
month_start=month_end-day_month+1;

% Potential evaportransporation is the monthly potential evapotransiration 
% and is computed using the Thornwaith equation: 
% using the mean daylight hours of month m, phi  being the phi and , D being the Julian day

lat_rad = deg2rad(phi);                                        %[rad] phi in radians
delta = 0.409*sin(2*pi*(1:365)/365 -1.39);                          %delta parameter, necessary to compute the Ws
Ws = acos(-tan(lat_rad)*tan(delta));                                %parameter necessary to compute ND

ND = 24*Ws/pi;                                                      %[h] number of daylight hours of day D
NM = zeros(12,1);                                                   %preallocation of Nm vector
for i=1:12             
       NM(i,1)  = mean(ND(month_start(i):month_end(i)));            %[h] mean daylight hours for each month
end
% Heat index: 
% Experimental exponent: 
% Computation of I and a

heat_index_I = sum((temperature/5).^1.514);
a=6.75e-7*heat_index_I^3-7.71e-5*heat_index_I^2+1.79e-2*heat_index_I+0.49;

% Computation of potential evapotranspiration
ET0 =  16*NM./12 .* (10*temperature./heat_index_I).^a; %[mm/month]
ET0_daily=ET0./day_month';                                          %[mm/day]
ET0_hourly = ET0_daily./24;%[mm/h]

% Runoff generation and evapotranspirtation computation
% Initialisation of values and space preallocation for vectors Precipitation/Runoff/Infiltration

t=0;

I=zeros(size(P));       %[L/T]  Infiltration
R=zeros(size(P));       %[L/T]  Runoff
s=zeros(size(P));       %[-]    Soil moisture
ET=zeros(size(P));      %[L/T]  Evapotranspiration
L=zeros(size(P));       %[L/]   Leaching
qsub=zeros(size(P));    %[L/T]  Sub-surface specific discharge
Vsub=zeros(size(P));    %[L]    Water stored in sub-surface layer
Qsub=zeros(size(P));    %[L3/T] Sub-surface discharge
qsup=zeros(size(P));    %[L/T]  Superficial specifc discharge
Vsup=zeros(size(P));    %[L]    Water stored in superficial layer
Qsup=zeros(size(P));    %[L3/T] Discharge in superficial layer
Q=zeros(size(P));       %[L3/T] Total discharge arriving to the reservoir

%initial moisture
s(1)=0.5;               %[-]

%initial volumes
Vsub(1)=0;              %[mm]
Vsup(1)=0;              %[mm]

%error variable to check if there is an error in evapotranspiration ET calculations
error_in_ET_calculation=0;

% Iterations for hourly calculations
% for loops for years, months and hours 
for year=1:n_years
    for month=1:12 
        for hour=1:day_month(month)*24
            t=t+1;
            
% Runoff generation calculations
            %Runoff R and infiltration I
            I(t) =  min(P(t),Ksat*1000*3600);   %[mm/h] because P is in [mm/h] et Ksat must be converted from [m/s] to [mm/h] 
            R(t) = P(t)-I(t);                   %[mm/h]
            
% Evapotranspiration ET computation in [mm/h]
 
            %Evapotranspiration computation ET [mm/h]
            if s(t)>=0 && s(t)<=sw
                ET(t) = 0;                                                                                          %[mm/h]
            elseif s(t)>sw && s(t)<s1
                ET(t) = kc(month)*ET0_hourly(month) *(s(t)-sw)/(s1-sw);   %[mm/h]
            elseif s(t)>s1 && s(t)<=1
                ET(t) = kc(month)*ET0_hourly(month);                                                %[mm/h]
            else
                error_in_ET_calculation=1;
            end
% Leaching L calculations
            %Leaching
            L(t) = Ksat*1000*3600*s(t)^c;   %[mm/h] (Ksat must be converted from [m/s] to [mm/h])
            
%Soil moisture calculation from soily moisture dynamics
            %Soil moisture calculations
            s(t+1)=s(t)+dt*(I(t)-ET(t)-L(t))/(n*z); %[-]
            
            %Boundary conditions (defition of saturation)
            if s(t+1)>1
                s(t+1)=1;
            elseif s(t+1)<0
                s(t+1)=0;
            end
            
% For a linear reservoie scheme, calculations of total discharge arriving in reservoir
            %Superficial specific discharge, storage and discharge calculations
            qsup(t)=Vsup(t)/tsup;         %[mm/h]
            Vsup(t+1)=Vsup(t)+dt*(R(t)-qsup(t));        %[mm]
            Qsup(t)=A*1e6*qsup(t)*1e-3/3600;  %[m3/s] (bassin_eare must be converted from [km2] to [m2] and qsup from [mm/h] to [m/s])
            
            %Sub-surface specific discharge, storage and discharge calculations
            qsub(t)=Vsub(t)/tsub;         %[mm/h]
            Vsub(t+1)=Vsub(t)+dt*(L(t)-qsub(t));        %[mm]
            Qsub(t)=A*1e6*qsub(t)*1e-3/3600;  %[m3/s] (bassin_eare must be converted from [km2] to [m2] and qsup from [mm/h] to [m/s])
            
            %Total discharge calculations
            Q(t)=Qsup(t)+Qsub(t)+Qb;                    %[m3/s]
            
        end
    end
end

%just transposing the Q vector to avoid a 52k x 52k matrix
Q=Q';

% Model Testing

%compute total discharge Q
if doTest % if doTest=1, check mass balance
P_tot=sum(P)*dt;
R_tot=sum(R)*dt;
L_tot=sum(L)*dt;
ET_tot=sum(ET)*dt;
%"testS" balance for the root zone (input/output). testS = 1: necessary (but not sufficient)
%condition for the implementation to be correct
testS=P_tot/(ET_tot+R_tot+L_tot+n*z*(s(end)-s(1)))
%"testQ" balance for the whole system (input/output). testQ = 1: necessary (but not sufficient)
%condition for the implementation to be correct
testQ=sum(P-ET)/(sum(qsup+qsub)+n*z*(s(end)-s(1))+Vsup(end)+Vsub(end))
end

end

