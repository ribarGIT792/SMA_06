
%  Single layer ANN

function SL_ANN

clc;close all

Ninputs=2       % iejimu skaicius
Nneurons=10      % vienintelio sluoksnio neuronu skaicius  
Noutputs=1      % isejimu skaicius
NLearning=10    % apmokymo aibe

W0=ones(Nneurons,Ninputs)     % daugikliu matrica funkciju i neuronus is iejimu
W1=ones(Noutputs,Nneurons)     % daugikliu matrica funkciju i isejimus is neuronu
B0=zeros(Nneurons,1)           % pastovios dalys i neuronus

for i=1:Ninputs, aL(i,:)=[0:0.45:0.45*(NLearning-1)];end
for i=1:Noutputs, sL(i,1:NLearning)=Learning(aL(i,:)); end
aL
sL
figure(1);hold on;plot(aL,sL,'b*')

s=OutputANN(aL,W0,B0,W1,@phiNeuron)
psi0=TargetFunction(s,sL)
for i=1:100
    D=GradVector(aL,sL,s,W0,B0,W1,@phiNeuron,@DphiNeuron)
    [GradW0,GradB0,GradW1]=GradVectorToMatrices(D,Ninputs,Nneurons,Noutputs) 
    ds=0.1;
    W0=W0+GradW0(:,:,1)*ds;
    B0=B0+GradB0(:,:,1)*ds;
    W1=W1+GradW1(:,:,1)*ds;
    s=OutputANN(aL,W0,B0,W1,@phiNeuron)
    psi=TargetFunction(s,sL) 
    W0,B0,W1
    if psi<psi0, psi0=psi; else, break, end
        stest=OutputANN(aL,W0,B0,W1,@phiNeuron)
        figure(1);hold on;plot(aL,stest,'r*')
    pause
end
end
    
function phi=phiNeuron(a)     %     Neurono funkcija
% a -  i neurona ateinanti informacija is iejimo ir pastovios dalies
%   
phi=a;
return
end

function dphi=DphiNeuron(a)   %     Neurono funkcijos isvestine
% a -  i neurona ateinanti informacija is iejimo ir pastovios dalies
%   
dphi=1;
return, end

function s=OutputANN(a,W0,B0,W1,phiNeuron)   %   ANN isejimo funkcija
% a  -  iejimas
% W0(Nneurons,Ninputs) -  perdavimo koeficientai is iejimu i neuronus
% W1(Noutputs,Nneurons) -  perdavimo koeficientai is neuronus i isejimus
% B0(Nneurons) -  neuronu iejimu pastovios dalys
% phiANN - neurono funkcijos vardas
% 
    s=W1*(phiNeuron(W0*a+B0));
return, end

function sL=Learning(a)    %        ANN apmokymo funkcija
%   a  -  iejimu vektorius
% 
sL=a.^2;
return, end

function D=GradVector(aL,sL,s,W0,B0,W1,phiNeuron,DphiNeuron)    % Tikslo funkcijos gradientas
% aL,sL - apmokymo iejimau ir isejimau matricos, stulpeliai yra apmokymo atvejai
% s - esamo ANN atsakas i apmokymo iejimus
% W0(Nneurons,Ninputs) -  perdavimo koeficientai is iejimu i neuronus
% W1(Noutputs,Nneurons) -  perdavimo koeficientai is neuronu i isejimus
% B0(Nneurons) -  neuronu iejimu pastovios dalys
% phiANN - neurono funkcijos vardas
% ***********  skaiciuojama vienam isejimui D=[dF/dW1, dF/dW0, dF/dB0]
% ***********  kai isejimu daugiau D, turi kelias eilutes

disp('---------------------GradVector   pradeda')
Ninputs=size(W0,2)       % iejimu skaicius
Nneurons=size(W0,1)      % vienintelio sluoksnio neuronu skaicius  
Noutputs=size(W1,1)      % isejimu skaicius  
NLearning=size(aL,2)     % apmokymo aibes dydis  

if size(aL,2)~=size(sL,2), disp('nesuderinti aL ir sL stulpeliu skaiciai ');end

D=zeros(Nneurons+Ninputs*Nneurons+Nneurons,Noutputs);  % paruosiamas gradiento vektorius

for iii=1:Noutputs  %---------------------

    for k=1:NLearning 

        D(1:Nneurons,iii)=D(1:Nneurons,iii)+reshape(phiNeuron(W0*aL(:,k)+B0),[],1);
        
        for j=1:Ninputs
            D(Nneurons+(j-1)*Nneurons+1:Nneurons+j*Nneurons,iii)=...
            D(Nneurons+(j-1)*Nneurons+1:Nneurons+j*Nneurons,iii)+...
            reshape(DphiNeuron(W0*aL(:,k)+B0),[],1)*aL(j,k).*W1(iii,:)';
        end
        D(Nneurons+Ninputs*Nneurons+1:Nneurons+Ninputs*Nneurons+Nneurons,iii)=...
        D(Nneurons+Ninputs*Nneurons+1:Nneurons+Ninputs*Nneurons+Nneurons,iii)+reshape(DphiNeuron(W0*aL(:,k)+B0),[],1);
        D(:,iii)=D(:,iii)*(s(iii)-sL(iii,k));
    end
    % gradiento normalizavimas:
    D(:,iii)= D(:,iii)/norm(D(:,iii));
end  %---------------------

return,end

function [GradW0,GradB0,GradW1]=GradVectorToMatrices(D,Ninputs,Nneurons,Noutputs)
% gradiento vektoriu isdesto ANN matricu pozicijose
for iii=1:Noutputs
    GradW1(iii,:)=D(1:Nneurons,iii);
    GradW0(:,:,iii)=reshape(D(Nneurons+1:Nneurons+Ninputs*Nneurons,iii),Nneurons,Ninputs);
    GradB0(:,iii)=D(Nneurons+Ninputs*Nneurons+1:Nneurons+Ninputs*Nneurons+Nneurons,iii);
end
return,end

function psi=TargetFunction(s,sL) 
% sL - apmokymo iejimau matrica, stulpeliai yra apmokymo atvejai
% s - esamo ANN atsakas i apmokymo iejimus
for iii=1:size(sL,1)
    psi(iii)=sum((s(iii,:)-sL(iii,:)).^2)/2; 
end
return,end
