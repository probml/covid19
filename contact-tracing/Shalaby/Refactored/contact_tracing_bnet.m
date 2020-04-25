function results_p = contact_tracing_bnet(data, N, T, infected_person, T_infection)
 % data(:, first-contact, second-contact, time-index)


%% DBN definition
%% 1=Not Infected
%  2=Infected

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
evidence{idx(infected_person,T_infection)}=2; %Person 1 infected at last time step


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


        
end




