%This script computes and saves the adjacency matrices A_sched(t) and A_real(t),
%relative to the scheduled and realized flights for one day. 

dt=20; %20 minutes time snapshot
N_AIR=14; %number of airlines
N_A=322; %number of airports

load('day.mat')
%The table 'day' contains the following information for all the flight of one day: 
%N = unique flight ID, Month and Day when the flight was operated, Airline (number refers to list in 'airline.mat'
%file), flight number, origin and destination airports (number refers to
%the order in the list in 'airports.mat' file), departure and arrivals time
%(scheduled and actual) in minutes from midnight, minutes of delay, and
%flags for diverted and cancelled flights. 

%If the flight arrived early, we set its arrival time to its scheduled
%one, otherwise we could allow connections that were not possible according to
%schedule (see methods and section S6)
day.ARRIVAL_TIME=max(day.SCHEDULED_ARRIVAL, day.ARRIVAL_TIME);
%if the flight departs early, we set its departure time to the scheduled
%one. This is needed so that the correction described in section S6 works
%properly 
day.DEPARTURE_TIME=max(day.SCHEDULED_DEPARTURE, day.DEPARTURE_TIME);

%Here, if we want we can impose a minimum connecting time of Dt, by simply
%adding Dt to all scheduled and arrival times:
%Dt=15;
%day.ARRIVAL_TIME=day.ARRIVAL_TIME+Dt;
%day.SCHEDULED_ARRIVAL=day.SCHEDULED_ARRIVAL+Dt;

t0=0; %beginning of the day, in this case it is 0 because schedules are measured in minutes from midnight

N=height(day); %number of flights
day.N=(1:N)';

N_step=ceil((max(day.ARRIVAL_TIME)-t0)/60)*(60/dt);  %number of hours to analyze for each day, to count all landings of flight departing that day
A_SCHED=cell(1,N_step);
A_REAL=cell(1,N_step);

%SCHEDULED
for j=1:N_step
    % The array index will contain a line for each link that we want to have in the adjacency matrix
    index=[];
    %flights departing in this snapshot
    I_SCHED=day(day.SCHEDULED_DEPARTURE>=t0+20*(j-1) & day.SCHEDULED_DEPARTURE<t0+20*(j),:);
    for k=1:length(I_SCHED.DAY)
        index=cat(1,index,[(I_SCHED.AIRLINE(k)-1)*N_A+I_SCHED.ORIGIN_AIRPORT(k),N_A*N_AIR+I_SCHED.N(k),1]); %for the airport i and airline j i use node (j-1)*N_A+i
    end
    %flights landing in this snapshot
    I_SCHED=day(day.SCHEDULED_ARRIVAL>=t0+20*(j-1) & day.SCHEDULED_ARRIVAL<=t0+20*(j),:);
    for k=1:length(I_SCHED.DAY)
        index=cat(1,index,[N_A*N_AIR+I_SCHED.N(k),(I_SCHED.AIRLINE(k)-1)*N_A+I_SCHED.DESTINATION_AIRPORT(k),1]);
    end
    index=[N_A*N_AIR+N,N_A*N_AIR+N,0]; %to assure that all matrix have the same size N_A*N_AIR+N, with one entry per airport copy plus one per flight
   
     A_SCHED{1,j}=sparse(index(:,1),index(:,2),index(:,3));
end

%REALIZED

for j=1:N_step
    % The array index will contain a line for each link that we want to have in the adjacency matrix
    index=[];
    %flights departing in this snapshot
    I_REAL=day(day.DEPARTURE_TIME>=t0+20*(j-1) & day.DEPARTURE_TIME<t0+20*(j) & day.CANCELLED==0 & day.DIVERTED==0,:);
    for k=1:length(I_REAL.DAY)
        index=cat(1,index,[(I_REAL.AIRLINE(k)-1)*N_A+I_REAL.ORIGIN_AIRPORT(k),N_A*N_AIR+I_REAL.N(k),1]); %for the airport i and airline j i use node (j-1)*N_A+i
    end
    %flights landing in this snapshot
    I_REAL=day(day.ARRIVAL_TIME>=t0+20*(j-1) & day.ARRIVAL_TIME<=t0+20*(j)& day.CANCELLED==0,:);
    for k=1:length(I_REAL.DAY)
        index=cat(1,index,[N_A*N_AIR+I_REAL.N(k),(I_REAL.AIRLINE(k)-1)*N_A+I_REAL.DESTINATION_AIRPORT(k),1]);
    end
    index=[N_A*N_AIR+N,N_A*N_AIR+N,0]; %to assure that all matrix have the same size
   
    A_REAL{1,j}=sparse(index(:,1),index(:,2),index(:,3));
end

save('A_SCHED.mat','A_SCHED')
save('A_REAL.mat','A_REAL')





