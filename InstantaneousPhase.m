%Hilbert transform program for kymograph written by Kumiko Yoshioka-Kobayashi
%Lab. of Growth Regulation System, Infront, Kyoto University
%181213
%Requirements: Matlab2016a or later, Signal Processing Toolbox, Image
%processing Toolbox

clear all; close all; clc;
%Parameters
mov_w=61;  %Window size for calculating the average, must be odd
h_mov_w=(mov_w-1)/2;
sg_l=15; %Window size for SG filtering, must be odd
sg_d=4; %n-Degree polynominal fillting for SG filtering


%Options
sg_o=1; %1:Performing filtering 0:Using non-filtered image
saveimage_o=0; %1:Saving images 0:Not saving


%Main process
%Import kymograph image (time dimension must be horizontal
[inputname,inputpath] = uigetfile('*.tif');   %Open standard dialog box for retrieving a tif file
infullpath=strcat(inputpath,inputname); %Create full path of the file
raw0=double(imread(infullpath));   %Read the file data and put it into a matrix, convert it into double to prepera for following processes
raw=rot90(raw0,1);
[r,c]=size(raw);    %Get the image size

%Moving_average for de-trending (along the time axis)
mov=zeros(r,c);
karnel=repmat(1/mov_w,mov_w,1);

for i=1:r
    mov(i,:)=conv(padarray(raw(i,:),[0 h_mov_w],'symmetric'),karnel,'valid'); %Take moving average, while applying mirror padding
end

%De-trending
det=raw-mov;

%%
%Savitzky-Golay filtering
if sg_o==1
    det_sg=sgolayfilt(det,sg_d,sg_l,[],2);
else
    det_sg=det';
end


%Hilbert transform (along the time axis)
hil= hilbert(det_sg'); 
instphase=angle(hil);   %P = angle(Z) returns the phase angles, in radians, for each element of complex array Z.
instphase2=rot90(instphase',-1);

%Display
scrsz = get(groot,'ScreenSize');
%figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)]);
figure,
subplot(5,1,1)
imshow(raw,[min(min(raw)) max(max(raw))]) ;
title('Raw image','FontSize',12);
subplot(5,1,2)
imshow(mov,[min(min(mov)) max(max(mov))]) ;
title('Moving average','FontSize',12);
subplot(5,1,3)
imshow(det,[min(min(det)) max(max(det))]) ;
title('De-trended','FontSize',12);
subplot(5,1,4)
imshow(det_sg,[min(min(det_sg)) max(max(det_sg))]) ;
title('SG-filtered','FontSize',12);
subplot(5,1,5)
imshow(instphase',[-3.14 3.14],'Colormap',jet) ;
title('Phase','FontSize',12);
grid on

figure,
imshow(raw0, [0 10000], 'Colormap',jet) ; %colorbar
if saveimage_o==1
    newname = strrep(inputname,'.tif','_raw.pdf');
    saveas(gcf,newname);
end;
figure, 
imshow(instphase2,[-3.14 3.14],'Colormap',jet) ; %colorbar
if saveimage_o==1
    newname = strrep(inputname,'.tif','_phase.pdf');
    newname2 = strrep(inputname,'.tif','.mat');
    saveas(gcf,newname);
    save(newname2,'instphase2');
end;
%%
%Save the image matrix into a text file
if saveimage_o==1
    newname = strrep(inputname,'.tif','.txt')
    %[outoutname,outputpath] = uiputfile(newname,'Save the instantaneous phase image as'); %Open dialog to obtain file name from imput
    outfullpath=strcat(inputpath, newname); 
    dlmwrite(outfullpath,instphase2); %Write matrix to ASCII-delimited file

end

%%
figure
time=(0:1:253)';
data1=mean(instphase2(:,[175:185]),2);
plot(time, instphase2(:,[180,230,280]));legend;