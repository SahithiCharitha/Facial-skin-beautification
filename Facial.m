clc;
close all;
clear all;
warning off;
 
[filename,pathname]=uigetfile('*.jpg','*.png');
I=imread([pathname,filename]);
I=imresize(I,[256,256]);
imshow(I);title('INPUT IMAGE');
a=rgb2gray(I);
 
%%LUMINANCE(Intensity)
cform = makecform('srgb2lab');
lab = applycform(I,cform); 
figure,imshow(lab);title('LAB IMAGE');
 v = lab(:,:,1);
figure,imshow(v);title('L IMAGE');
 g = lab(:,:,2);
figure,imshow(g);title('a* IMAGE');
 
d = lab(:,:,3);
figure,imshow(d);title('b* IMAGE');
 
%%EDGE PRESERVING SMOOTHING OPERATOR (Lighting layer)
z=wlsfil(v);
z=cast(z,'uint8');
figure,imshow(z);title('Lighting layer IMAGE');
 
%% DETAIL LAYER
s= imsubtract(v,z);
figure,imshow(s,[]);title('DETAIL LAYER IMAGE');
 
%%COLOUR LAYER
g1=g(:,1:(256/2));
d1=d(:,(256-(256/2))+1:256);
gc=[g1,d1];
figure,imshow(gc);
title('COLOUR LAYER IMAGE 1');
 
I1=final();

I1=imresize(I1,[256,256]);
figure,imshow(I1);
[R1,z1,gc12]=bret(I1);
 
cd ../

bw=edge(z,'canny',[],2); 
bww=~bw;
figure,imshow(bww);
title('bw');
bw1 = bwareaopen(bw,30); %remove small objects from image - extra edges
dup=I;
se = strel('disk',5); %creates a morphological structuring disk element with radius of 2 pixels
bw2 = imclose(bw1,se); % Morpholigically closes the image. The morphological close operation is a dilation followed by an erosion, using the same structuring element for both operations.
bw3 = imdilate(bw2,se);
bw3 = imfill(bw3,'holes');
bw4=~bw3;
epps=0.0001;
ALPHAL=1.2;
LAMBDAL=0.4;
RL2=im2bw(z,0.5);
RL1=bwareaopen(RL2,3500);
RL=uint8(~RL1);
figure,imshow(RL,[]);
 
ML=uint8(RL);
 
G=log(double(z));
[Mxl ,Myl]=imgradient(double(ML));

[Gxl ,Gyl]=imgradient((G));
 
 [AL,BL]=size(ML);
count =1;
for il=1:AL
for jl=1:BL
gradL(il,jl)=LAMBDAL.*(((Mxl(il,jl))/(((Gxl(il,jl)).^(ALPHAL))+epps)) + ((Myl(il,jl))/(((Gyl(il,jl)).^(ALPHAL))+epps)));
WL(il,jl)=(((ML(il,jl))-(RL(il,jl))).^2);
ML1(il,jl) =(WL(il,jl).*(((ML(il,jl))-(RL(il,jl))).^2));
ML2(il,jl) = double(ML1(il,jl)) + gradL(il,jl);
count =count+1
end
end
MIL=min(ML2);
%%
RS=(uint8(bw4)).*(RL);
figure,imshow(RS);
MS=RS;
[l,k]=size(RS);
for i=1:l
    for j=1:k
if ((RS(i,j))==1)
    MS(i,j)=a(i,j);
end
    end
end
figure,imshow(MS);
[x ,y]=size(v);
for i=1:x
for j=1:y
GS1(i,j)=norm(double(v(i,j)));
GS2(i,j)=norm(double(g(i,j)));
GS3(i,j)=norm(double(d(i,j)));
end
end

[AS,BS]=size(MS);
for is=1:AS
for js=1:BS
if ((v(is,js))==0)
ALPHAS=1;
LAMBDAS=0.1;
GS=cat(3,GS1,GS2,GS3);
else
ALPHAS=0.01;
LAMBDAS=0.1;
GS=v;
end
end
end
 
[Mxs ,Mys]=imgradient((MS));
[Gxs ,Gys]=imgradient((GS));
count1 =1;
 
for il=1:AS
for jl=1:BS

gradS(is,js)=LAMBDAS.*(((Mxs(is,js))/(((Gxs(is,js)).^(ALPHAS))+epps)) + ((Mys(is,js))/(((Gys(is,js)).^(ALPHAS))+epps)));
WS(is,js)=(((MS(is,js))-(RS(is,js))).^2);
MS1(is,js) =(WS(is,js).*(((MS(is,js))-(RS(is,js))).^2));
MS2(is,js) = (MS1(is,js)) + real(gradS(is,js));
count1 =count1+1
end
end
MIS=min(MS2);
 %%
cw=edge(z,'canny',0.25,3); % find edges in the grey scale image
cw1 = bwareaopen(cw,10); %remove small objects from image - extra edges
sec = strel('disk',8); %creates a morphological structuring disk element with radius of 2 pixels
cw2 = imclose(cw1,sec); % Morpholigically closes the image. The morphological close operation is a dilation followed by an erosion, using the same structuring element for both operations.
dupe=dup;
cw3 = imdilate(cw2,sec);
cw3 = imfill(cw3,'holes');
cw4=~cw3;
 
RC=uint8(cw4);
figure,imshow(RC);
 
[l1,k1]=size(RC);
MC=RC;
for i=1:l1
    for j=1:k1
if ((RC(i,j))==1)

MC(i,j)=a(i,j);
end
    end
end
figure,imshow(MC);
 
[x1 ,y1]=size(v);
for i=1:x1
for j=1:y1
GC1(i,j)=norm(double(v(i,j)));
GC2(i,j)=norm(double(g(i,j)));
GC3(i,j)=norm(double(d(i,j)));
end
end
[AC,BC]=size(MC);
for ic=1:AC
for jc=1:BC
if ((v(ic,jc))==0)
ALPHAC=1;
LAMBDAC=0.8;
GC=cat(3,GC1,GC2,GC3);
else
ALPHAC=1;
LAMBDAC=0.8;
GC=v;
end
end
end
[Mxc ,Myc]=imgradient((MC));
	[Gxc ,Gyc]=imgradient((GC));
count2 =1; 	
for ic=1:AC
for jc=1:BC
gradC(ic,jc)=LAMBDAC.*(((Mxc(ic,jc))/(((Gxc(ic,jc)).^(ALPHAC))+epps)) + ((Myc(ic,jc))/(((Gyc(ic,jc)).^(ALPHAC))+epps)));
WC(ic,jc)=(((MC(ic,jc))-(RC(ic,jc))).^2);
MC1(ic,jc) =(WC(ic,jc).*(((MC(ic,jc))-(RC(ic,jc))).^2));
MC2(ic,jc) = MC1(ic,jc) + gradC(ic,jc);
count2 =count2+1
end
end
MIC=min(MC2);
%%
[im]=interactive(s,z,gc,MS,ML,MC,dupe);
figure,imshow(im);
