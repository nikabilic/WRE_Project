function [time_gen,Q_gen,R_gen,I_gen,soil_saturation_gen,L_gen,ET_gen,Pgen_hourly,Q_gen_mean_monthly] = discharge(Nyears_gen, P, temperature, kc, dt, n, max_z,max_Kc,sw,s1,max_c,tsup,max_tsub,A,Qb,phi,doTest)

N_days=length(P)/24;                                     %[-]    Number of days of the given rainfall series
P_daily=zeros(N_days,1);

%upscaling - is this right?
for i=1:N_days
    P_daily(i)=sum(P(24*(i-1)+1:24*i));                  %[mm/d] Average daily rainfall
end

day_month=[31 28 31 30 31 30 31 31 30 31 30 31];         %       Number of days for each month
month_end=cumsum(day_month);                             %       Last day of each month
month_start=month_end-day_month+1;                       %       First day of each month

NyearsP=length(P)/365/24;                                %[-]    Number of years of precipitation

%monthly statistics
Pmean_monthly=zeros(12,1);                               %[mm]   Monthly mean precipitation (counting also non rainy days)
Pstd_monthly=zeros(12,1);                                %[mm]   Monthly standard deviation of precipitation (counting also non rainy days)
Lambda_monthly=zeros(12,1);                              %[1/d]  Monthly rainfall events arrival rate
Alpha_monthly=zeros(12,1);                               %[mm]   Monthly mean daily rainfall depth (counting only the rainy days)

P_daily_matrix=reshape(P_daily,365,NyearsP);              %      Matrix with P_daily values organized in day(row) and years(column)
for m=1:12                                                %      Loop on months 
    P_temp=P_daily_matrix(month_start(m):month_end(m),:); %      Temporary matrix. Stores the average daily rainfall of a specific month, for each year  
    Pmean_monthly(m)=mean(P_temp,'all');
    Pstd_monthly(m)=std(P_temp,0,"all");
    Lambda_monthly(m)=numel(P_temp(P_temp>0))/numel(P_temp);
    Alpha_monthly(m)=sum(P_temp,'all')/numel(P_temp(P_temp>0));
end

save("Pobs_statistics_parameter.mat","Pmean_monthly","Pstd_monthly","Alpha_monthly","Lambda_monthly")

%generate 100 year daily rainfall
Pgen_daily_matrix=zeros(365,Nyears_gen);                  %      Matrix with P_daily values organized in day(row) and years(column)

for i=1:Nyears_gen                                        %      Loop on the years 
    day_counter=0;                                        %      Counter for the day of the year, initialized every year
    for m=1:12                                            %      Loop on the months
        for d=1:day_month(m)                              %      Loop on the days of the given month
            day_counter=day_counter+1;                    %      Update days counter
       %does it have to be lambda/alpha daily here or is lambda/aplha monthly fine? 
       %also there is a matlab code for Poisson process which is a bit different
            if rand<Lambda_monthly(m)*1
                Pgen_daily_matrix(day_counter,i)=-Alpha_monthly(m)*log(-rand+1);
            end
        end
    end
end

%recompute monthly stats
P_gen_mean_monthly=zeros(12,1);
Pgenstd_monthly=zeros(12,1);
Lambdagen_monthly=zeros(12,1);
Alphagen_monthly=zeros(12,1);

for m=1:12                                                %      Loop on the months 
    P_temp=Pgen_daily_matrix(month_start(m):month_end(m),:);
    P_gen_mean_monthly(m)=mean(P_temp,'all');
    Pgenstd_monthly(m)=std(P_temp,0,"all");
    Lambdagen_monthly(m)=numel(P_temp(P_temp>0))/numel(P_temp);
    Alphagen_monthly(m)=sum(P_temp,'all')/numel(P_temp(P_temp>0));
end

save("Pgen_statistics_parameter.mat","P_gen_mean_monthly","Pgenstd_monthly","Alphagen_monthly","Lambdagen_monthly")

%downscaling P HOURLY IS EMPTY-I cannot solve this problem
Pgen_daily_matrix=Pgen_daily_matrix(:);     %[mm/d] Transforms the matrix into a column vector
Pgen_hourly=downscaling(Pgen_daily_matrix)                     %[mm/h] Hourly generated rainfall

%run model with new generated rainfall
[Q_gen, R_gen, I_gen,soil_saturation_gen,L_gen,ET_gen]=f_hydromodel(phi,temperature,Nyears_gen,Pgen_hourly,max_Kc,sw,s1,kc,n,Qb,tsup,max_tsub,A,max_c,dt,max_z,doTest);

time_gen=2016+(0:length(Pgen_hourly)-1)/(24*365);

%mean monthly discharge
Q_gen_mean_daily = zeros(length(Q_gen)/24,1); %mean discharge per day in [m^3/s] along the generation time
Q_gen_mean_monthly = zeros(12,1); %mean discharge per month in [m^3/s] for a strict year

for i=1:length(Q_gen)/24
    Q_gen_mean_daily(i) = mean(Q_gen(1+24*(i-1) : i*24));
end

Q_gen_daily_matrix = reshape(Q_gen_mean_daily,365,Nyears_gen);
Q_gen_daily=mean(Q_gen_daily_matrix,2)

for m=1:12
    Q_temp=Q_gen_daily_matrix(month_start(m):month_end(m),:);
    Q_gen_mean_monthly(m)=mean(Q_temp,'all');
end

end

