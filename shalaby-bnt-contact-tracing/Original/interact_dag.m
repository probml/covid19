function interact_dag(contacts_t,t)
    contact_pair=combntns(contacts_t,2);
    for c=1:length(contact_pair)
       i=contact_pair(1,1);
       j=contact_pair(1,2) ;
       dag(idx(i,t),idx(j,t+1)) = true;
       dag(idx(j,t),idx(i,t+1)) = true;
    end
end