function exp_data = import_dd(fname)
T = readtable(fname);

last = strcmp(T.last_ask, 'true');
del_rew = T.del_rew;
last_10k = T(last & del_rew == 10000,:);
last_500 = T(last & del_rew == 500,:);
exp_data = [];
dell = [0; sort(last_500.del_m)];

exp_data.delays = dell;
value_hi = sortrows([last_10k.del_m last_10k.adj_rew],1);
exp_data.value_hi = [ 10000; value_hi(:,2) ];
value_low = sortrows([last_500.del_m last_500.adj_rew],1);
exp_data.value_low = [ 500; value_low(:,2) ];
exp_data.response_time_hi = T(del_rew == 10000,:).rt;
exp_data.response_time_low = T(del_rew == 500,:).rt;
