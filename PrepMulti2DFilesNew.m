%% Build multiple 2D files for 3D scans, added by Takeshi

%% Set directory for files and filenames
tic;
clear all; clc;
m = 0;
FFIndex = '13';
PorT = 'T'; %changing prepulse power (P) or delay (T)
RootDir = strcat('./3D',num2str(str2num(FFIndex)));

FileName = {'FWM';'REF';'SI';'tra';'refa';'traref'};

% Determine whether this scan is rephasing (1) or non-rephasing (0)
isRephasing = 1;

if isRephasing
    MatrixDir = '/2dmatrixS1';
    SpecDir = strcat('./',RootDir,'/spectra/S1spec');
else
    MatrixDir = '/2dmatrixS2';
    SpecDir = strcat('./',RootDir,'/spectra/S2spec');
end

%% Convert files recorded at starting position. Please notice the program 
% will skip the missing files.
for i=1:6
fid = fopen(strcat(FileName{i},FFIndex,'.dat'),'r');
if (fid==-1), continue, end
[A,count] = fscanf(fid,'%f %f',[2 inf]);
fclose(fid);
dlmwrite(strcat(RootDir,'/',FileName{i},'00.dat'), round(A(2,:)'));
end

% making directory
if ~isdir(strcat('./',RootDir,MatrixDir)) mkdir(strcat('./',RootDir,MatrixDir)); end

%% 1. writing data to a new directry
%% 2. swap the data file if its number is odd (scanning reverse way)
while 1
    % number to string, mth scan of P or T
    mStr = num2str(m);
    
    fid=fopen(strcat(SpecDir,PorT, mStr, '.dat'), 'r');
    if (fid==-1), break, end
    fclose(fid);
    
    specfile = strcat(SpecDir,PorT,mStr,'.dat');
    A = dlmread(specfile,'\t');
    N = size(A,1);
    Res = mod(N,1024);
    Ntau = N - Res;
    ACut = A(1:Ntau);

    % write to matrix
    M2D = reshape(ACut,1024,[]);
    % If m is odd number, flip in the left-right direction
    if isRephasing
        if mod(m,2)
            M2D = fliplr(M2D);
        end
    else
        if mod(m+1,2)
            M2D = fliplr(M2D);
        end
    end
    dlmwrite(strcat('./',RootDir, MatrixDir, '/SItotal',PorT,mStr,'.dat'), M2D, '\t');
    clear A; clear M2D; 
    msgbox(strcat('The 2D file for', PorT,' ', mStr,' has been generated'),'replace');
    m=m+1;
end
msgbox('All files for 3D program have been generated', 'Mission Completed');
toc;