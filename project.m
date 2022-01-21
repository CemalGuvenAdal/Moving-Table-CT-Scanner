%loading images for the chestvolume
load chestVolume;
load mri.mat
D1=V;
B = squeeze(D);
%displaying a 3d image of the chest

%scanning table length
tablelength=0.76;
%tablemoving velocity
tablevelocity=0.005;
%how many slices
sizearray=size(D1);
slicenumber=sizearray(3);
%lengthofeachslice
slicelength=tablelength/slicenumber;
%time for scanning each slice
scanduration=slicelength/tablevelocity;
scannerframe=200;
numberofangles=scannerframe*(scanduration);
rps=90*scanduration;

%downrounding number of angles for each slice
numberofangles=floor(numberofangles)
x=tomo(D1,numberofangles,rps);
%output
D1=V(:,:,:);
volumeViewer(D1)
figure
subplot(2,2,1)
imshow(D1(:,:,100),[])
title('Original for Slice 100')
subplot(2,2,2)
imshow(D1(:,:,150),[])
title('Original for Slice 150')
subplot(2,2,3)
imshow(D1(:,:,200),[])
title('Original for Slice 200')
subplot(2,2,4)
imshow(D1(:,:,250),[])
title('Original for Slice 250')
figure;

% Reconstructed Lung
subplot(2,2,1)
imshow(x(:,:,100),[])
title('Reconstructed Lung CT for Slice 100')
subplot(2,2,2)
imshow(x(:,:,150),[])
title('Reconstructed Lung CT for Slice 150')
subplot(2,2,3)
imshow(x(:,:,200),[])
title('Reconstructed Lung CT for Slice 200')
subplot(2,2,4)
imshow(x(:,:,250),[])
title('Reconstructed Lung CT for Slice 250')
figure;
volumeViewer(x);
thtrial=linspace(0,45,95);
trial=radon(D1(:,:,200),thtrial);
trialreconstructed=iradon(trial,thtrial);
imshow(trialreconstructed,[]);
title("Reconstruction of slice 200 with 95 degrees and 0:45 interval")
figure;
%brain dataset
volumeViewer(B)

subplot(2,2,1)
imshow(B(:,:,5),[])
title('Original for Slice 5')
subplot(2,2,2)
imshow(B(:,:,10),[])
title('Original for Slice 10')
subplot(2,2,3)
imshow(B(:,:,15),[])
title('Original for Slice 15')
subplot(2,2,4)
imshow(B(:,:,20),[])
title('Original for Slice 20')
figure;

%new parameters assigned for brain dataset
%scanning table length
tablelength=0.35;
%tablemoving velocity
tablevelocity=0.005;
%how many slices
sizearray=size(B);
slicenumber=sizearray(3);
%lengthofeachslice
slicelength=tablelength/(slicenumber*2);
%time for scanning each slice
scanduration=slicelength/tablevelocity;
scannerframe=200;
numberofangles=scannerframe*(scanduration);
rps=90*scanduration;
%------------------------------
y=tomo(B,numberofangles,rps);
subplot(2,2,1)
imshow(y(:,:,5),[])
title('Reconstructed Brain slice  for Slice 5')
subplot(2,2,2)
imshow(y(:,:,10),[])
title('Reconstructed Brain slice for Slice 10')
subplot(2,2,3)
imshow(y(:,:,15),[])
title('Reconstructed Brain slice  for Slice 15')
subplot(2,2,4)
imshow(y(:,:,20),[])
title('Reconstructed Brain slice for Slice 20')
figure;
thtrial=linspace(0,116,259);
trial=radon(B(:,:,5),thtrial);
trialreconstructed=iradon(trial,thtrial);
imshow(trialreconstructed,[]);
title("Reconstruction of slice 5 with 259 degrees and 0:116 interval")
figure
volumeViewer(y);


function fintable=tomo(matrix,numberofangles,rps)
sizearray=size(matrix);
slicenumber=sizearray(3);
thetatr=linspace(0,rps,numberofangles);
Rb=radon(matrix(:,:,1),thetatr)
radonsize=size(Rb);
rowsize=radonsize(1);
colsize=numberofangles*slicenumber;
colsize=int16(colsize)
dptable=zeros(rowsize,colsize);
size(dptable)
numberofangles
%Reconstruction of the dptable
startdegree=0
for i=1:slicenumber
    theta=linspace(startdegree,startdegree+rps,numberofangles);
    Rad=radon(matrix(:,:,i),theta);
    dptable(:,(i-1)*numberofangles+1:i*numberofangles)=Rad;
    startdegree=floor(startdegree+rps);
 
end
groupnumber=360/rps;
projectionnumber=groupnumber*numberofangles;
projectionnumber=floor(projectionnumber);
thetareconstruction=linspace(0,360,projectionnumber);
trial=iradon(dptable(:,1:projectionnumber),thetareconstruction);
trialsize=size(trial);

reconstructiontable=zeros(radonsize(1),radonsize(2),sizearray(3));
for j=1:sizearray(3)
    j
    ref=j*numberofangles;
    search=(j-groupnumber)*numberofangles
    if j<=groupnumber
        search=1;
    end

    for x=search:search+groupnumber*numberofangles
        x=floor(x);
        if or(x>ref-projectionnumber/2, x<=ref+projectionnumber/2)
            y=rem(x,projectionnumber);
            y=floor(y);
            if y==0
                y=1;
            end
            reconstructiontable(:,y,j)=dptable(:,x);
        end
        
    
    end

end

%-----------
d=size(reconstructiontable);
reconstheta=linspace(0,360,d(2)) ;
denem=iradon(reconstructiontable(:,:,5),reconstheta);
fintablesize=size(denem);
fintable=zeros(fintablesize(1),fintablesize(2),sizearray(3));
for h=1:sizearray(3)
    fintable(:,:,h)=iradon(reconstructiontable(:,:,h),reconstheta);
end
end



