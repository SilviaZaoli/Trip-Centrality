%This script computes the incoming and outgoing Trip centralities for each
%airport (both the layer-specific ones and the aggregated onse), as well as
%the flight centralities (i.e. centrality of secondary nodes)

dt=20; %20 minutes time snapshot
N_AIR=14; %number of airlines
N_A=322; %number of airports

load('A_SCHED')
load('A_REAL')

N_SCHED=length(A_SCHED{1,1})-N_A*N_AIR; %number of flights
t0=0; %beginning of the day, in this case it is 0 because schedules are measured in minutes from midnight

N_step=length(A_SCHED); %number of hours to analyze for each day, to count all landings of flight departing that day

%compute K matrix
e=0.1; %value of parameter epsilon

index=zeros(N_A*N_AIR^2+N,3); %will contain links of matrix K
count=1;
for i=1:N_A
    for j=1:N_AIR
        for k=1:N_AIR
            if k~=j
                index(count,:)=[(j-1)*N_A+i,(k-1)*N_A+i,e];   %connect copies of the same airport on different layers with weight e
                count=count+1;
            end
        end
    end
end
for i=1:N_A*N_AIR+N_SCHED
    index(count+i-1,:)=[i,i,1]; %ones on the diagonal
end
K=sparse(index(:,1),index(:,2),index(:,3));
G=sparse(index(:,1),index(:,2),ones(length(index(:,1)),1)); %this matrix is used later when we want to sum the rows or columns of Q (see later) over all copies of one airport in different layers

clearvars index count

a=sqrt(0.2);   %value of parameter \tilda(alpha)=sqrt(alpha_
D=N_A*N_AIR+N; %size of matrices

%Computation of Q matrices 
Q_SCHED=speye(D);
Q_REAL=speye(D);
for i=1:N_step
    Q_SCHED=Q_SCHED*(speye(D)+a*A_SCHED{1,i}*K); 
    Q_REAL=Q_REAL*(speye(D)+a*A_REAL{1,i}*K);
    Q_REAL=min(Q_REAL,Q_SCHED);  %correction for possible new walk opened by delay (see section S6)
end
if e~=1 
    Q_SCHED=(Q_SCHED-speye(D))/K;
    Q_REAL_cloned=(Q_REAL_cloned-speye(D))/K;
else  %if e=1, K is not invertible, as explained in the paper, so we do divide by K in this case. The aggregated centralities can still be used
    Q_SCHED=(Q_SCHED-speye(D));
    Q_REAL=(Q_REAL-speye(D));
end

%Computation of layer-specific centralities
c_SCHED_in=full(sum(Q_SCHED(1:N_A*N_AIR,1:N_A*N_AIR),1));  %incoming centrality vector for scheduled network
c_SCHED_out=full(sum(Q_SCHED(1:N_A*N_AIR,1:N_A*N_AIR),2));  %outgoing centrality vector for scheduled network
c_REAL_in=full(sum(Q_REAL(1:N_A*N_AIR,1:N_A*N_AIR),1));    %incoming centrality vector for realized network
c_REAL_out=full(sum(Q_REAL(1:N_A*N_AIR,1:N_A*N_AIR),2));    %outgoing centrality vector for realized network

save(strcat('c_SCHED_in_lay-spec_e_',num2str(e),'.mat'),'c_SCHED_in')
save(strcat('c_SCHED_out_lay-spec_e_',num2str(e),'.mat'),'c_SCHED_out')
save(strcat('c_REAL_in_lay-spec_e_',num2str(e),'.mat'),'c_REAL_in')
save(strcat('c_REAL_out_lay-spec_e_',num2str(e),'.mat'),'c_REAL_out')

%centrality of flights
c_F_SCHED_in=full(sum(Q_SCHED,1));
c_F_SCHED_in=c_F_SCHED_in(N_A*N_AIR+1:end);
c_F_SCHED_out=full(sum(Q_SCHED,2));
c_F_SCHED_out=c_F_SCHED_out(N_A*N_AIR+1:end);
c_F_REAL_in=full(sum(Q_REAL,1));
c_F_REAL_in=c_F_REAL_in(N_A*N_AIR+1:end);
c_F_REAL_out=full(sum(Q_REAL,2));
c_F_REAL_out=c_F_REAL_out(N_A*N_AIR+1:end);
save(strcat('c_F_SCHED_in_e_',num2str(e),'.mat'),'c_F_SCHED_in')
save(strcat('c_F_SCHED_out_e_',num2str(e),'.mat'),'c_F_SCHED_out')
save(strcat('c_F_REAL_in_e_',num2str(e),'.mat'),'c_F_REAL_in')
save(strcat('c_F_REAL_out_e_',num2str(e),'.mat'),'c_F_REAL_out')

%Computation of aggregate centralities
Q_SCHED=G*Q_SCHED*G; %This double multiplication by G has the effect of summing each element over all rows and columns corresponding to copies of the same airport in other layers
Q_REAL=G*Q_REAL*G;
Q_SCHED=Q_SCHED(1:N_A,1:N_A); %The matrix Q_SCHED has repeted block of size N_AxN_A that are all the same after the above multiplication by G, so I keep just one
Q_REAL=Q_REAL(1:N_A,1:N_A);
c_SCHED_in=full(sum(Q_SCHED,1));  %incoming centrality vector for scheduled network
c_SCHED_out=full(sum(Q_SCHED,2));  %outgoing centrality vector for scheduled network
c_REAL_in=full(sum(Q_REAL,1));    %incoming centrality vector for realized network
c_REAL_out=full(sum(Q_REAL,2));    %outgoing centrality vector for realized network

save(strcat('c_SCHED_in_aggr_e_',num2str(h),'.mat'),'c_SCHED_in')
save(strcat('c_SCHED_out_aggr_e_',num2str(h),'.mat'),'c_SCHED_out')
save(strcat('c_REAL_in_aggr_e_',num2str(h),'.mat'),'c_REAL_in')
save(strcat('c_REAL_out_aggr_e_',num2str(h),'.mat'),'c_REAL_out')




