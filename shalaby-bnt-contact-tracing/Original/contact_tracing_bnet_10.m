function contact_tracing_bnet_10(file, infected_person_1,infected_person_2,infected_person_3,infected_person_4,infected_person_5,infected_person_6,infected_person_7,infected_person_8,infected_person_9,infected_person_10, T_infection)
tic;
addpath(genpath(pwd));
fid = fopen(file);
C = textscan(fid, '%s %s %s', 'HeaderLines', 0);
fclose(fid);
time_index_original=str2double(C{1});
time_index=unique(time_index_original);
first_contact=str2double(C{2});
second_contact=str2double(C{3});
first_contact_modified=first_contact;
second_contact_modified=second_contact;
time_index_modified=time_index;
C= cat(1, first_contact, second_contact);
C= unique(C);
C=cat(2,C,C);
for i=1:length(C)
    C(i,2)=i;
end

for j=1:length(first_contact)
    first_contact(j)=C(find(ismember(C(:,1),first_contact(j))),2);
    second_contact(j)=C(find(ismember(C(:,1),second_contact(j))),2);
end

%%Ts=unique(time_index);
Ts=time_index(1):600:time_index(end);
Ts=Ts';
Ts=cat(2,Ts, Ts);
for i=1:length(Ts)
    Ts(i,2)=i;
end

for t=1:length(time_index_original)
    time_index_original(t)=Ts(find(ismember(Ts(:,1),time_index_original(t))),2);
end




data=cat(2, first_contact, second_contact, time_index_original);
data=unique(data,'rows');
data=sortrows(data,3);
%load listcontacts_2009_05_16.txt_pruned_1_90.mat

%% DBN definition
%% 1=Not Infected
%  2=Infected
a=size(C);
N = a(1);
a=size(Ts);
T = a(1)+1;


node_sizes=2*ones(1,N*T);
dag=false(N*T, N*T);
idx = @(n,t) N*(t-1)+n;

for t= 1:T-1
    for n= 1:N
        dag(idx(n,t), idx(n,t+1))=true;
    end
end

time_steps=[min(data(:,3)):max(data(:,3))]';
contacts=cell(1,length(time_steps));

for t=1:length(time_steps)
    a=[];
    for i=1:length(data)
       if data(i,3)==time_steps(t)
           a=[a data(i,1)];
           a=[a data(i,2)];
       end
    end
    contacts{t}=unique(a);
end



for t=1:length(time_steps)
    if not(isempty(contacts{t}))
        contact_pair=combntns(contacts{t},2);
        k=size(contact_pair);
        k=k(1);
        for c=1:k
            i=contact_pair(c,1);
            j=contact_pair(c,2); 
            dag(idx(i,t),idx(j,t+1)) = true;
            dag(idx(j,t),idx(i,t+1)) = true;
        end
    end
end

bnet = mk_bnet(dag, node_sizes);


 
%%Defining the CPDS 
CPT_single=[1.0 0 0 1.0]';
CPT_root=[ 0.95 0.05]';

for n=1:N
    bnet.CPD{idx(n,1)}=tabular_CPD(bnet,idx(n,1),'CPT',CPT_root);
end

for n=1:N
    for t=2:T
        ps=[];
        node_parents=parents(dag,idx(n,t));
        w=length(node_parents);
        if w==1
            bnet.CPD{idx(n,t)}=tabular_CPD(bnet,idx(n,t),'CPT',CPT_single);
        else
            for np=1:length(node_parents)
                if node_parents(np)==idx(n,t-1)
                    ps=[ps 0];
                else
                    ps=[ps 0.9615];
                end
            end
            bnet.CPD{idx(n,t)}=noisyor_CPD(bnet,idx(n,t), 1.0, ps);
        end
    end
end





evidence=cell(1,N*T);

%engine=jtree_inf_engine(bnet);
%engine=var_elim_inf_engine(bnet);
engine=pearl_inf_engine(bnet);
%save workspace.mat
evidence{idx(infected_person_1,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_2,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_3,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_4,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_5,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_6,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_7,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_8,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_9,T_infection)}=2; %Person 1 infected at last time step
evidence{idx(infected_person_10,T_infection)}=2; %Person 1 infected at last time step

[engine]=enter_evidence(engine, evidence);

results_p=cell(N,T);
for t=1:T
    t;
    for n=1:N
        n;
        marg=marginal_nodes(engine, [idx(n,t)]);
        results_p{n,t}=marg.T';
    end
end

save(strcat(file, '_', num2str(T_infection),'_' , num2str(infected_person_1), '_10.mat'), 'results_p');
toc



%Most are noisy-OR CPDS except the first time step which are tabular priors

        
end





