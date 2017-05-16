function newdata=shave(data,fig,a,b)

% function to remove unwanted fourier component
% data = data set (vector)
% fig = switch between figure or predefined values
% a, b are the start and stop values

ttild=fft(data);
Pt=ttild.*conj(ttild);

if fig==0
     % plot select start and stop points
     plot(Pt(1:length(Pt)/2));title('orig power spectrum')
     a=ginput(2)
     start=floor(a(1,1))
     stop=floor(a(2,1))
 else
     start = a;
     stop = b;
 end

newdata=remo(start,stop,ttild.');
redo=ifft(newdata);

length(redo);
newdata=real(redo);
newdata = newdata.';

function newdata=remo(start,stop,data)
fmax=length(data);
data(1:start-1);
qq=fmax-stop-1;
length(data);
data(stop:(fmax-stop-1));
data((fmax-start):length(data));

newdata = [data(1:start-1),zeros(1,stop-start),data(stop:(fmax-stop-1)),zeros(1,stop-start),data((fmax-start):length(data))];