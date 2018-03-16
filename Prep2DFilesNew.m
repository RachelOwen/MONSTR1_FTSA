%% Set directory for files and filenames
clear all; clc;
FFIndex = '11';
% Number of scans saved simultaneously
NScans = 1;
% calibration slope is negative (0) or positive (1)
slope = 0;
RootDir = strcat('.\2D',num2str(str2num(FFIndex)));
SpecDir = strcat('.\',RootDir,'\spectra\spec');
FileName = {'FWM';'REF';'SI';'tra';'refa';'traref'};
%% Convert files recorded at starting position. Please notice the program 
% will skip the missing files.
for i=1:6
fid = fopen(strcat(FileName{i},FFIndex,'.dat'),'r');
if (fid==-1), continue, end
[A,count] = fscanf(fid,'%f %f',[2 inf]);
fclose(fid);
% flip the matrix left to right if the calibration slope is positive.
% A = fliplr(A);
dlmwrite(strcat(RootDir,'\',FileName{i},'00.dat'), round(A(2,:)'));
end

for j = 1 : NScans
    %% load all spectra into one matrix
    specfile = strcat(SpecDir,num2str(j),'.dat');
    A = dlmread(specfile,'\t');
%     if j == 1
%         N = size(A,1);
%         Res = mod(N,1024);
%     end
    
    N = size(A,1);
    tausteps = 0;
    if tausteps
        Ntau = tausteps * 1024;
        ACut = A(1:Ntau);
    else
        Res = mod(N,1024);
        Ntau = N - Res;
        ACut = A(1:Ntau);
    end
    
    % write to matrix
    M2D = reshape(ACut,1024,[]);
    
%     X = size(M2D);
%     Y = floor(X(2)/2);
%     
%     M2D2 = zeros(1024,Y);
%     for n = 1 : Y
%         M2D2(:,n) = M2D(:,2*n-1);
%     end
%     dlmwrite(strcat('.\',RootDir,'\SI0total',num2str(j),'.dat'), M2D2, '\t');
   
    % flip the matrix up side down if the calibration slope is positive.
    if slope
     M2D=flipud(M2D);
    end
    % save matrix to file
    dlmwrite(strcat('.\',RootDir,'\SI0total',num2str(j),'.dat'), M2D, '\t');
%         dlmwrite(strcat('.\',RootDir,'\SI0total','.dat'), M2D, '\t');

end
% dlmwrite(strcat('.\',RootDir,'\SI0total.dat'), fliplr(M2D), '\t'); %Backwards scan
clear A; clear M2D; 
msgbox('All files for 2D program have been generated', 'Mission Completed');
