%% Processing the Data

%/usr/local/bin/matlab -nodesktop << EOF  
addpath(genpath(pwd))

fid = fopen('listcontacts_2009_06_27_pruned.txt');
C = textscan(fid, '%s %s %s', 'HeaderLines', 0);
fid = fopen('listcontacts_2009_06_27_pruned.txt');
C = textscan(fid, '%s %s %s', 'HeaderLines', 0);
fclose(fid);
time_index=str2double(C{1});
first_contact=str2double(C{2});
second_contact=str2double(C{3});
first_contact_modified=first_contact;
second_contact_modified=second_contact;
time_index_modified=time_index;
F= cat(1, first_contact, second_contact);
F= unique(F);
F=cat(2,F,F);
for i=1:length(F)
    F(i,2)=i;
end

for j=1:length(first_contact)
    first_contact(j)=F(find(ismember(F(:,1),first_contact(j))),2);
    second_contact(j)=F(find(ismember(F(:,1),second_contact(j))),2);
end

T=unique(time_index);
T=cat(2,T, T);
for t=1:length(T)
    T(t,2)=t;
end

for t=1:length(time_index)
    time_index(t)=T(find(ismember(T(:,1),time_index(t))),2);
end


data=cat(2, first_contact, second_contact, time_index);


%% DBN definition

%%P=inputdlg('Input the  size of the population:');
P=(length(F));
%%P=str2num(cell2num(P));
N=combination(P,2);
names=names_generation(N,P);
     
     


ss=length(names);

    



intrac=intra_generation(P);

[intra, names] = mk_adj_mat(intrac, names, 1);
interc=inter_generation(P);



inter = mk_adj_mat(interc, names, 0);
obs=observed_nodes_generation(N, P);

for i=1:length(obs)
  onodes(i) = strmatch(obs{i}, names ,'exact'); %stringmatch(obs{i}, names);
end
onodes = sort(onodes);

dnodes=1:ss;
ns =zeros(1,ss);
for i= 1:P
    ns(strmatch(sprintf('Y%d',i),names, 'exact')) = 4;
    ns(strmatch(sprintf('F%d',i),names , 'exact')) = 4;
    for j=1:P
        if i ~=j
            ns(strmatch(sprintf('S%d,%d',i,j),names ,'exact'))=2;
            ns(strmatch(sprintf('T%d,%d',i,j),names,'exact'))=4;
        end
    end
end
            
    
    

eclass1=1:ss;
 
eclass2=(1:ss)+ss;
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1',eclass1, 'eclass2', eclass2, 'observed', onodes);

%bnet.CPD{strmatch('F2',names)}=tabular_CPD(bnet, strmatch('F2',names));
%% STEP2: DEFINING the CPTs.

% CPD generation parameters
IRLS_iter = 100;

Y=CPTY();

F1=CPTF();
Fin= [0.25;0.25;0.25;0.25];
Sin=[0.5;0.5];
% F2=CPTF2();
% F3=CPTF3();
T=CPTT();
% for i=1: 2*ss
%     bnet.CPD{i}=tabular_CPD(bnet, i);
% end
load('softmax_weights_3nodes','o');
W=[ o(1:4,1:4);repmat(o(5:8,1:4),P-1,1)]; %softmax weights matrix
clamped = 0;

for i=1:P
bnet.CPD{strmatch(sprintf('Y%d',i),names,'exact')}=tabular_CPD(bnet, strmatch(sprintf('Y%d',i),names,'exact'), Y); %% The Ys of time slice 1
bnet.CPD{strmatch(sprintf('Y%d',i),names,'exact')+ss}=tabular_CPD(bnet, strmatch(sprintf('Y%d',i),names,'exact')+ss, Y); %% The Ys of time slice 2 

bnet.CPD{strmatch(sprintf('F%d',i),names,'exact')}=tabular_CPD(bnet, strmatch(sprintf('F%d',i),names,'exact'),Fin); %% The Fs of  time slice 1
bnet.CPD{strmatch(sprintf('F%d',i),names,'exact')+ss}=softmax_CPD(bnet, strmatch(sprintf('F%d',i),names,'exact')+ss, 'discrete' , parents(inter, strmatch(sprintf('F%d',i),names,'exact')), 'weights', W,  'offset', zeros(1,4), 'clamped', clamped, 'max_iter',IRLS_iter );

for j=1:P
    if i~=j
        bnet.CPD{strmatch(sprintf('T%d,%d',i,j),names,'exact')+ss}=tabular_CPD(bnet, strmatch(sprintf('T%d,%d',i,j),names,'exact')+ss, T); %% The Ts of time slice 2
        bnet.CPD{strmatch(sprintf('T%d,%d',i,j),names,'exact')}=tabular_CPD(bnet, strmatch(sprintf('T%d,%d',i,j),names,'exact'), T); %% the Ts of time slice 1
    end
    if i<j
        bnet.CPD{strmatch(sprintf('S%d,%d',i,j),names,'exact')+ss}=tabular_CPD(bnet, strmatch(sprintf('S%d,%d',i,j),names,'exact')+ss, Sin); %% the Ss of time slice 2
        bnet.CPD{strmatch(sprintf('S%d,%d',i,j),names,'exact')}=tabular_CPD(bnet, strmatch(sprintf('S%d,%d',i,j),names,'exact'),Sin); %% the Ss of time slice 1
    end
end

end




% make rnd params
%for i=1:2*ss
 % bnet.CPD{i} = tabular_CPD(bnet, i);
%end


%% STEP3:  INFERENCE ENGINE DEFINITION
% T = 3;

 %engine = {};
 %engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
  %engine{end+1} = jtree_dbn_inf_engine(bnet);
 %engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
 %engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet)); % observed nodes have children
 
 %inf_time = cmp_inference_dbn(bnet, engine, T)
 %learning_time = cmp_learning_dbn(bnet, engine, T)
 %%T=inputdlg('Input the  time duration:');
 T=length(T);
%%T=str2num(cell2num(T));
%  T= 10;%fixed length sequences

%engine = {};
%engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);


  engine = jtree_unrolled_dbn_inf_engine(bnet, T);
%engine = smoother_engine(jtree_2TBN_inf_engine(bnet));


 %engine =smoother_engine(hmm_2TBN_inf_engine(bnet));

%  engine=jtree_dbn_inf_engine(bnet);
% engine=bk_inf_engine(bnet, 'clusters', 'ff');





%% STEP4: ENTERING EVIDENCE
%   TESTING
%%
% for i =1:4
%     for j=1:4
%         for k=1:4
%             
% ev = sample_dbn(bnet, T);
evidence = cell(ss,T);


%observed_ID_health=str2num(cell2num(inputdlg('Enter the ID of the patients who were tested')));
observed_ID_health=num2str(8);
for u=1:size(observed_ID_health')
prompt = {['Enter the time steps at which the health state of Patient '  num2str(observed_ID_health(u)) ' being observed:'],['Enter the health state of Patient ' num2str(observed_ID_health(u)) ' at these time steps:']};
dlg_title = ['Health data of Patient ' num2str(observed_ID_health(u))'];
num_lines = 1;
def = {'','', '', ''};
%answer = inputdlg(prompt,dlg_title,num_lines,def);
answer=cell(2,1);
answer{1}=num2str(3);
answer{2}=num2str(3);
evidence(  strmatch(sprintf('Y%d',(observed_ID_health(u))),names,'exact'),str2num((answer{1})))=num2cell(str2num(answer{2}));

end

% t=1;
% prompt = {['Enter the IDs of individuals who have been in contact at time step] ' num2str(t) ' :']};  
% clear t 
% def={''};
% for t=2:T
% x={['Enter the IDs of individuals who have been in contact at time step] ' num2str(t) ' :']};
% x=char(x);
% prompt{t}=x;
% def{t}=def{1};
% end
% dlg_title = ['Contact Data of The Population'];
% num_lines = 1;

%answer2 = inputdlg(prompt,dlg_title,num_lines,def);
evidence(strmatch('S',names),:)={1};
answer2=cell(T,1);
for t=1:T
    b=[];
    for i=1:length(data)
        if data(i,1)==t 
            d=data(i,3);
            b=cat(1,b,d); 
            d=data(i,2);
            b=cat(1,b,d);     
            answer2{t}=num2str((unique(b))')
        %    a(t)=b';
        end
    end
end
for t=1:T
    E=str2num(answer2{t});
    for l=1:length(E)
        for k=1:length(E)
            if l<k
                evidence(strmatch(sprintf('S%d,%d',E(l),E(k)),names,'exact'),t)={2};
            end
        end    
    end
end 
% save('free_memory','engine')
% save('free_memory','evidence')
% save('free_memory','names')
% clear
% load('free_memory','engine')
% load('free_memory','evidence')
% load('free_memory','names')
% 
% evidence={[],[],[],[],[],[],[],[],[],[];4,[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];4,[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];4,[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];1,[],[],[],[],[],[],[],[],[];1,1,1,1,1,1,1,2,2,2;[],[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];1,1,1,1,1,1,1,1,1,1;[],[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];1,1,1,1,2,2,2,1,1,1;[],[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];1,1,1,1,1,1,1,1,1,1;[],[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];1,1,1,1,1,1,1,1,1,1;[],[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];2,2,2,2,1,1,1,1,1,1;[],[],[],[],[],[],[],[],[],[];[],[],[],[],[],[],[],[],[],[];};

[engine, ll] = enter_evidence(engine, evidence);
save(output_variables);
for i=1:T
    marg=marginal_nodes(engine, strmatch('F1',names, 'exact'), i);
%         marg=marginal_nodes(engine,5, 2);
    marg.T'
end

%         end 
%     end
% end

%     i
%     
    

% end

% for i=1:10
%     marg=marginal_nodes(engine, strmatch('F3',names), i);
%     marg.T
%     i
%     
% 
% end
% % 
% % G = bnet.dag;
% % draw_graph(G);
%         end
%     end
% end
% 
exit
 
