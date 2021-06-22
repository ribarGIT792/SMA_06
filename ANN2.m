
%  Single layer ANN

function SL_ANN

clc;close all

Ninputs=2       % iejimu skaicius
Nneurons=5      % vienintelio sluoksnio neuronu skaicius  
Noutputs=1      % isejimu skaicius
NLearning=10    % apmokymo aibe

W0=ones(Nneurons,Ninputs)*0.1;     % daugikliu matrica funkciju i neuronus is iejimu
W1=ones(Noutputs,Nneurons)*0.1;     % daugikliu matrica funkciju i isejimus is neuronu
B0=zeros(Nneurons,1)*0.1;           % pastovios dalys i neuronus

for i=1:Ninputs, aL(i,:)=[0:0.58:0.58*(NLearning-1)];end
for i=1:Noutputs, sL(i,1:NLearning)=Learning(aL(i,:)); end
aL
sL
figure(1);hold on;plot(aL,sL,'b*')

s=OutputANN(aL,W0,B0,W1,@PhiNeuron);
psi0=TargetFunction(s,sL);

for i=1:100
    [DW1,DW0,DB0]=GradVector(aL,sL,s,W0,B0,W1,@PhiNeuron)
    ds=0.01;
    W0=W0-DW0(:,:)*ds;
    B0=B0-DB0*ds;
    W1=W1-DW1*ds;
    s=OutputANN(aL,W0,B0,W1,@PhiNeuron);
    psi=TargetFunction(s,sL); 
    W0,B0,W1
    if psi<psi0, psi0=psi; else, break, end
        stest=OutputANN(aL,W0,B0,W1,@PhiNeuron)
        figure(1);hold on;plot(aL,stest,'r*')
    pause
end
end
    
function phi=PhiNeuron(a)     %     Neurono funkcija
% a -  i neurona ateinanti informacija is iejimo ir pastovios dalies
%  
phi=sin(a);
return
end

function s=OutputANN(a,W0,B0,W1,PhiNeuron)   %   ANN isejimo funkcija
% a  -  iejimas
% W0(Nneurons,Ninputs) -  perdavimo koeficientai is iejimu i neuronus
% W1(Noutputs,Nneurons) -  perdavimo koeficientai is neuronus i isejimus
% B0(Nneurons) -  neuronu iejimu pastovios dalys
% phiANN - neurono funkcijos vardas
% 
    s=W1*PhiNeuron(W0*a+B0);
return, end 

function sL=Learning(a)    %        ANN apmokymo funkcija
%   a  -  iejimu matrica [Ninputs,Nlearning]
% 
sL=atan(a);
% sL=a;
return, end

function [DW1,DW0,DB0]=GradVector(aL,sL,s,W0,B0,W1,PhiNeuron,DphiNeuron)    % Tikslo funkcijos gradientas
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

DB0=zeros( size(B0));DW0=zeros( size(W0));DW1=zeros( size(W1));

% for iii=1:Noutputs  %---------------------

    for k=1:NLearning 
dds=0.01;
s1=OutputANN(aL(:,k),W0,B0,W1,PhiNeuron)
for i=1:Nneurons
    W11=W1;W11(i)=W11(i)+dds;
    DW1(i)=DW1(i)+(s1-sL(k))*(OutputANN(aL(:,k),W0,B0,W11,PhiNeuron)-s1)/dds;
end

for j=1:Ninputs
    for i=1:Nneurons
        W01=W0;W01(i,j)=W01(i,j)+dds;
        DW0(i,j)=DW0(i,j)+(s1-sL(k))*(OutputANN(aL(:,k),W01,B0,W1,PhiNeuron)-s1)/dds;
    end
end

for i=1:Nneurons
    B01=B0;B01(i)=B01(i)+dds;
    DB0(i)=DB0(i)+(s1-sL(k))*(OutputANN(aL(:,k),W0,B01,W1,PhiNeuron)-s1)/dds;
end
    
    end
    
    % gradiento normalizavimas:
    max(abs(DW1))
    max(max(abs(DW0)))
    max(abs(DB0))
    amax=max([max(abs(DW1)),max(max(abs(DW0))),max(abs(DB0))]);
    DW1=DW1/amax;  DW0=DW0/amax; DB0=DB0/amax;  
    
    
    
% end  %---------------------

return,end


function psi=TargetFunction(s,sL) 
% sL - apmokymo iejimau matrica, stulpeliai yra apmokymo atvejai
% s - esamo ANN atsakas i apmokymo iejimus
for iii=1:size(sL,1)
    psi(iii)=sum((s(iii,:)-sL(iii,:)).^2)/2; 
end
return,end
