function [data, C, Ts] = read_data(file)

fid = fopen(file);
C = textscan(fid, '%s %s %s', 'HeaderLines', 0);
fclose(fid);
time_index_original=str2double(C{1});
time_index=unique(time_index_original);
first_contact=str2double(C{2});
second_contact=str2double(C{3});
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

end
