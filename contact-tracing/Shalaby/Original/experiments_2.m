file2='listcontacts_2009_06_27.txt_pruned'; %80 individuals and 36 time_steps
N=80;
for n=1:N
    infected_person=n; %most social
    T_infection=1;
    contact_tracing_bnet(file2,infected_person,T_infection);
end