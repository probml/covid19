
folder = '../Data'; 
file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
fname = sprintf('%s/%s', folder, file);
 [data, C, Ts] = read_data(file); % data: 402x3, C: 80x2, Ts: 38x2
 
 a=size(C);
N = a(1); % 80 people
a=size(Ts);
T = a(1)+1; % 39 steps

 infected_person=51; %most social
 %T_infection=38; % end
 T_infection=1; % start 
 results =  contact_tracing_bnet(data, N, T,infected_person,T_infection); % {80,39}

bel  = zeros(N, T);
for n=1:N
    for t=1:T
        %fprintf('n=%d, t=%d\n',n,t);
        tmp = results{n,t};
        if (n==infected_person) && (t==T_infection)
            tmp = [0 1];
        end
        bel(n,t) = tmp(2);
    end
end

heatmap(bel')
colormap('jet')

%{
Tinfect=1
> bel(:,end)'
ans =
  Columns 1 through 17
    0.4459    0.2611    0.8417    0.6021    0.6021    0.6021    0.1309    0.6243    0.2901    0.7066    0.8823    0.8136    0.6021    0.6405    0.1309    0.7312    0.2901
  Columns 18 through 34
    0.4721    0.1858    0.7654    0.1809    0.6021    0.1495    0.6243    0.5401    0.2332    0.4459    0.0645    0.8840    0.8823    0.2332    0.1213    0.6405    0.1495
  Columns 35 through 51
    0.8840    0.1213    0.8792    0.6021    0.0645    0.4357    0.1858    0.5688    0.8136    0.0868    0.0990    0.4459    0.1284    0.1495    0.5975    0.0645    1.0000
  Columns 52 through 68
    0.0641    0.1265    0.9017    0.5401    0.7287    0.1005    0.4459    0.5688    0.0855    0.1106    0.4459    0.1284    0.7312    0.0633    0.0748    0.0990    0.4459
  Columns 69 through 80
%}

%{
Tinfect=38
bel(:,end)'
ans =
  Columns 1 through 17
         0    0.2611    0.2139         0         0         0    0.1309    0.1161    0.2901    0.4196    0.0355    0.0355         0         0    0.1309    0.4465    0.2901
  Columns 18 through 34
    0.2476    0.1858         0    0.1809         0         0    0.1161    0.2994    0.2332         0         0    0.2139    0.0355    0.2332         0         0         0
  Columns 35 through 51
    0.2139         0    0.0355         0         0    0.1957    0.1858         0    0.1634    0.0868    0.0990         0    0.1284         0    0.3215         0    1.0000
  Columns 52 through 68
         0         0    0.1161    0.2994         0    0.1005         0         0    0.0855    0.1106         0    0.1284    0.4465    0.0633    0.0748    0.0990         0
  Columns 69 through 80
         0         0         0         0         0         0         0         0         0         0    0.1809         0

%}


% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=51; %most social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=51; %most social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=73; %medium social
% T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=73; %medium social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=73; %medium social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=28; %least social
% % T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=28; %least social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
% infected_person=28; %least social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);


%Another File

% 
% file='listcontacts_2009_07_04.txt_pruned_last'; %127 individuals and 36 time_steps
% infected_person=10; %most social
% T_infection=38;
% contact_tracing_bnet(file,infected_person,T_infection);

% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=10; %most social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=10; %most social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=86; %medium social
% T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=86; %medium social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=86; %medium social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=125; %least social
% T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);

% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=125; %least social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_07_04.txt_pruned'; %127 individuals and 36 time_steps
% infected_person=125; %least social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% 
% 
% %%Another file %%
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=6; %most social
% T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);

% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=6; %most social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=6; %most social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=108; %medium social
% T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=108; %medium social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=108; %medium social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=90; %least social
% T_infection=1;
% contact_tracing_bnet(file,infected_person,T_infection);

% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=90; %least social
% T_infection=18;
% contact_tracing_bnet(file,infected_person,T_infection);
% 
% 
% file='listcontacts_2009_05_16.txt_pruned'; %241 individuals and 36 time_steps
% infected_person=90; %least social
% T_infection=36;
% contact_tracing_bnet(file,infected_person,T_infection);


