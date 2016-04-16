load Output_Table

theta = Output_Table(:,1);
Output_Table = Output_Table(:,2:end);
temp(:,1:12) = Output_Table(:,1:12);
temp(:,13:24) = Output_Table(:,19:30);
temp(:,25:26) = Output_Table(:,13:14);
temp(:,27:28) = Output_Table(:,15:16);
temp(:,33:34) = Output_Table(:,17:18);
temp(:,29:30) = Output_Table(:,31:32);
temp(:,31:32) = Output_Table(:,33:34);
temp(:,35:36) = Output_Table(:,35:36);

s7val =  sort((e1-e8)*temp',2);
s8val =  sort((e2-e8)*temp',2);
s9val  =  sort((e3-e8)*temp',2);
s10val  =  sort((e4-e8)*temp',2);
s11val  =  sort((e5-e8)*temp',2);
s12val  =  sort((e6-e8)*temp',2);
s13val  =  sort((e7-e8)*temp',2);

figure()
stem(theta,s7val);
hold on
stem(theta,s8val);
hold on
stem(theta,s9val);
hold on
stem(theta,s10val);
hold on
stem(theta,s11val);
hold on
stem(theta,s12val);
hold on
stem(theta,s13val);