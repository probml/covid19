function myBNTTest
tic
N = 150;
T = N+2;

node_sizes = 2*ones(1,N*T);
dag = false(N*T,N*T);

idx = @(n,t) N*(t-1)+n;

for t = 1:T-1
  for n = 1:N
    dag(idx(n,t),idx(n,t+1)) = true;
  end
end

twoparents = zeros(1,N*T);

for i = 1:N-1
  interact_dag(i,i+1,i+1); %Person i and i+1 interact at time i+1
end

bnet = mk_bnet(dag, node_sizes);

CPT_interact1 = [ 1.0 0.8 1.0 0.0 0.0 0.2 0.0 1.0 ]'; %If two parents, direct ancestor comes first, match value
CPT_interact2 = [ 1.0 1.0 0.8 0.0 0.0 0.0 0.2 1.0 ]'; %If two parents, direct ancestor comes second, match value
CPT_root = [ 0.05 0.95 ]'; %Prior probability of infection is 0.05
CPT_single = [ 1.0 0.0 0.0 1.0 ]'; %If single parent, match value

for i = 1:N-1
  interact_cpds(i,i+1,i+1); %Person i and i+1 interact at time i+1
end

for n = 1:N
  bnet.CPD{idx(n,1)} = tabular_CPD(bnet,idx(n,1),'CPT',CPT_root);
  for t = 2:T
    if(twoparents(idx(n,t))), continue; end;
    bnet.CPD{idx(n,t)} = tabular_CPD(bnet,idx(n,t),'CPT',CPT_single);
  end
end

evidence = cell(1,N*T);

%engine = jtree_inf_engine(bnet);
%engine = var_elim_inf_engine(bnet);
engine = pearl_inf_engine(bnet);
evidence{idx(1,T)} = 1; %Person 1 infected at last time step

[engine] = enter_evidence(engine,evidence);
marg = marginal_nodes(engine, [idx(N,T)]); %Person N infected at last time step?

marg.T

function interact_dag(i,j,t)
  dag(idx(i,t),idx(j,t+1)) = true;
  dag(idx(j,t),idx(i,t+1)) = true;
  twoparents(idx(j,t+1)) = true;
  twoparents(idx(i,t+1)) = true;
end

function interact_cpds(i,j,t)  
  if(i < j)
    bnet.CPD{idx(i,t+1)} = tabular_CPD(bnet,idx(i,t+1),'CPT',CPT_interact1);
    bnet.CPD{idx(j,t+1)} = tabular_CPD(bnet,idx(j,t+1),'CPT',CPT_interact2);
  else
    bnet.CPD{idx(i,t+1)} = tabular_CPD(bnet,idx(i,t+1),'CPT',CPT_interact2);
    bnet.CPD{idx(j,t+1)} = tabular_CPD(bnet,idx(j,t+1),'CPT',CPT_interact1);
  end
end

toc
end
