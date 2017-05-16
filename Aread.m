function Amatrix = Aread(Astring)

% Simple function to open two column data produced by MCP.vi suit
% from which loads a matrix and closes file
% ADB Oct 2007

clear Amatrix;   
fid=fopen(Astring,'r');
[Amatrix,count] = fscanf(fid,'%f %f',[2 inf]);
fclose(fid);
clear count;
clear Astring;