
%  Single layer ANN

function SL_ANN

clc;close all

Ninputs=2       % iejimu skaicius
Nneurons=5      % vienintelio sluoksnio neuronu skaicius  
Noutputs=1      % isejimu skaicius
NLearning=10    % apmokymo aibe

W0=ones(Nneurons,Ninputs)*0.1     % daugikliu matrica funkciju i neuronus is iejimu
W1=ones(Noutputs,Nneurons)*0.1     % daugikliu matrica funkciju i isejimus is neuronu
B0=zeros(Nneurons,1)*0.1           % pastovios dalys i neuronus

for i=1:Ninputs, aL(i,:)=[0:0.58:0.58*(NLearning-1)];end
for i=1:Noutputs, sL(i,1:NLearning)=Learning(aL(i,:)); end
aL
sL
figure(1);hold on;plot(aL,sL,'b*')

NF0=NeuronFunction(aL,W0,B0,@PhiNeuron);
s=OutputANN(NF0,W1)
psi0=TargetFunction(s,sL)

for i=1:100
    [DW1,DW0,DB0]=GradVector(aL,sL,s,W0,B0,W1,@PhiNeuron,@DphiNeuron)
    ds=0.05 ;
    W0=W0-DW0(:,:)*ds;
    B0=B0-DB0(:,:)*ds;
    W1=W1-DW1(:,:)*ds;
    NF0=NeuronFunction(aL,W0,B0,@PhiNeuron);
    s=OutputANN(NF0,W1)
    psi=TargetFunction(s,sL) 
    W0,B0,W1
    if psi<psi0, psi0=psi; else, break, end
        NF0=NeuronFunction(aL,W0,B0,@PhiNeuron);
        stest=OutputANN(NF0,W1)
        figure(1);hold on;plot(aL,stest,'r*')
    pause
end
end
    
function phi=PhiNeuron(a)     %     Neurono funkcija
% a -  i neurona ateinanti informacija is iejimo ir pastovios dalies
%  
phi=a.^3;
return
end

function dphi=DphiNeuron(a)   %     Neurono funkcijos isvestine
% a -  i neurona ateinanti informacija is iejimo ir pastovios dalies
%   
dphi=2*a.^2;
return, end

function s=OutputANN(NF0,W1)   %   ANN isejimo funkcija
% a  -  iejimas
% W0(Nneurons,Ninputs) -  perdavimo koeficientai is iejimu i neuronus
% W1(Noutputs,Nneurons) -  perdavimo koeficientai is neuronus i isejimus
% B0(Nneurons) -  neuronu iejimu pastovios dalys
% phiANN - neurono funkcijos vardas
% 
    s=W1*NF0;
return, end 

function sL=Learning(a)    %        ANN apmokymo funkcija
%   a  -  iejimu vektorius
% 
sL=5*a.^3;
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

for i=1:Nneurons
    zzz=zeros(1,Nneurons);zzz(i)=1;
    NF0=NeuronFunction(aL(:,k),W0,B0,PhiNeuron);
    DW1(i)=DW1(i)+OutputANN(NF0,zzz)*(s(k)-sL(k)); 
end

for j=1:Ninputs
    for i=1:Nneurons
        zzz=zeros(size(W0));
        deriv=DphiNeuron(W0(i,:)*aL(:,k)+B0(i))
        zzz(i,:)=deriv*aL(:,k)';
        DW0(i,:)=DW0(i,:)+OutputANN(zzz,W1)*(s(k)-sL(k)); 
    end
end


for i=1:Nneurons
        zzz=zeros(Nneurons,1);
        zzz(i)=DphiNeuron(aL(j,k)*W0(i,j)+B0(i));
        DB0(i)=DB0(i)+OutputANN(zzz,W1)*(s(k)-sL(k));
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


function NF=NeuronFunction(a,W0,B0,PhiNeuron)

'NeuronFunction pradeda'
W0
B0
a
W0*a+B0
    NF=PhiNeuron(W0*a+B0);
return,end


function psi=TargetFunction(s,sL) 
% sL - apmokymo iejimau matrica, stulpeliai yra apmokymo atvejai
% s - esamo ANN atsakas i apmokymo iejimus
for iii=1:size(sL,1)
    psi(iii)=sum((s(iii,:)-sL(iii,:)).^2)/2; 
end
return,end
