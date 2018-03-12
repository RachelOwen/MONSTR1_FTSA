%% average pump-probe data and place in format for Tianhao's m-file
clear all; clc;
%% parameters
% file folder parameters
folder = './';                              % folder
FFIndex = '03';
PProot1 = strcat('PuPr-withA',FFIndex,'.dat');                 % pump-probe filename root
PProot2 = strcat('PuPr-withB',FFIndex,'.dat');                 % pump-probe filename root
PProot3 = strcat('PuPr-withC',FFIndex,'.dat');                 % pump-probe filename root
RRrefroot = strcat('PuPr-without',FFIndex,'.dat');             % probe filename root
PPbkgd = strcat('PuPr-bkgd',FFIndex,'.dat');
% bkgd = strcat('bkgd0',FFIndex,'.dat');
bkgd = PPbkgd;
% data parameters
Nbegin0 = 1;                          % points to cut at beinging
Nend0 = 1023;                            % points to cut at end
a = 75; b = 500;                         % filter variables
conCat = 1;                             % write out ? "0" = no, "1" = yes
lowendcut = 00;                      % lower end of zero truncation
highendcut = 1000;                     % higher end of zero truncation
N = 0.25;
offset1 = 0.0;           % add / subtract a constant offset
offset2 = 0.0;           % add / subtract a constant offset
offset3 = 0.0;           % add / subtract a constant offset
slope1 = 0.000;        % add / subtract a linear slope from data
slope2 = 0.000;        % add / subtract a linear slope from data
slope3 = 0.000;        % add / subtract a linear slope from data
%% load all spectra into matrix
filename11 = strcat(folder,PProot1);
filename12 = strcat(folder,PProot2);
filename13 = strcat(folder,PProot3);
filename2 = strcat(folder,RRrefroot);
filename3 = strcat(folder,PPbkgd);
% filename4 = strcat(folder,bkgd);
A1 = Aread(filename11);               % read files
A2 = Aread(filename12);               % read files
A3 = Aread(filename13);               % read files
B = Aread(filename2);
C = Aread(filename3);
% D = Aread(filename4);
PPSpec1 = round(A1(2,:)');       % make arrays
PPSpec2 = round(A2(2,:)');       % make arrays
PPSpec3 = round(A3(2,:)');       % make arrays
PPrefSpec = round(B(2,:)');
bkgdSpec = round(C(2,:)');
PPbkgdSpec = bkgdSpec;
clear A; clear B; clear C; clear D;
disp('Reading files complete!')
%% calculate DT and remove noise
% DTSpec = ( (PPSpec - PPbkgdSpec.*1.0) - (PPrefSpec - bkgdSpec.*1.0) ) ...
%     ./ (PPrefSpec - bkgdSpec.*1.0);
DTSpec1 = ( (PPSpec1 - PPbkgdSpec.*1.0) - (PPrefSpec - bkgdSpec.*1.0));
DTSpec2 = ( (PPSpec2 - PPbkgdSpec.*1.0) - (PPrefSpec - bkgdSpec.*1.0));
DTSpec3 = ( (PPSpec3 - PPbkgdSpec.*1.0) - (PPrefSpec - bkgdSpec.*1.0));
for n = (1:1024)
    DTSpec1(n) = (DTSpec1(n) + slope1.*(n-512));
    DTSpec2(n) = (DTSpec2(n) + slope2.*(n-512));
    DTSpec3(n) = (DTSpec3(n) + slope3.*(n-512));
end
DTSpec1 = DTSpec1 + offset1;
DTSpec2 = DTSpec2 + offset2;
DTSpec3 = DTSpec3 + offset3;
DTSpec1(1 : lowendcut) = 0; 
DTSpec2(1 : lowendcut) = 0; 
DTSpec3(1 : lowendcut) = 0; 
DTSpec1(highendcut : 1024) = 0;
DTSpec2(highendcut : 1024) = 0;
DTSpec3(highendcut : 1024) = 0;
DTSpec1 = DTSpec1/max(abs(DTSpec1));           % normalize
DTSpec2 = DTSpec2/max(abs(DTSpec2));           % normalize
DTSpec3 = DTSpec3/max(abs(DTSpec3));           % normalize
DTSpecFiltered1 = shave(DTSpec1,conCat,a,b);    % examine transform and filter
DTSpecFiltered2 = shave(DTSpec2,conCat,a,b);    % examine transform and filter
DTSpecFiltered3 = shave(DTSpec3,conCat,a,b);    % examine transform and filter
%DTSpec1(1:Nbegin0) = 0;                 % removing ends
%DTSpec1(end-Nend0:end) = 0;
figure(20)
if conCat == 0
    subplot(311)
    plot(real(DTSpecFiltered1))  
    hold on
    plot( real(DTSpec1), '-k') % plot
    hold off
    subplot(312)
    plot(real(DTSpecFiltered2))  
    hold on
    plot( real(DTSpec2), '-k') % plot
    hold off
    subplot(313)
    plot(real(DTSpecFiltered3))  
    hold on
    plot( real(DTSpec3), '-k') % plot
    hold off
end
subplot(311)
plot(DTSpec1)
hold on
plot(DTSpecFiltered1, '-k')
hold off
subplot(312)
plot(DTSpec2)
hold on
plot(DTSpecFiltered2, '-k')
hold off
subplot(313)
plot(DTSpec3)
hold on
plot(DTSpecFiltered3, '-k')
hold off
saveas(gcf, strcat('./2D',num2str(str2num(FFIndex)), '/DT_Data'), 'pdf');
%% concatonate and save matrix to file
if conCat == 0
    disp('no file created!');
elseif conCat == 1
    data = [DTSpecFiltered1 round(PPrefSpec(1:1024)) round(PPSpec1(1:1024)) DTSpec1(1:1024)];
    dlmwrite(strcat('./2D',num2str(str2num(FFIndex)),'/DT_data_A.dat'), data, '\t');
    data = [DTSpecFiltered2 round(PPrefSpec(1:1024)) round(PPSpec2(1:1024)) DTSpec2(1:1024)];
    dlmwrite(strcat('./2D',num2str(str2num(FFIndex)),'/DT_data_B.dat'), data, '\t');
    data = [DTSpecFiltered3 round(PPrefSpec(1:1024)) round(PPSpec3(1:1024)) DTSpec3(1:1024)];
    dlmwrite(strcat('./2D',num2str(str2num(FFIndex)),'/DT_data_C.dat'), data, '\t');
    disp('File written!');
else
    disp('wrong choice for "conCat"!');
end
    


