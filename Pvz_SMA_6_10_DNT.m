
%  Single layer ANN

function SL_ANN

clc;close all
hndl=[];
Ninputs=1       % iejimu skaicius
Nneurons=5     % vienintelio sluoksnio neuronu skaicius
Noutputs=1      % isejimu skaicius
NLearning=40    % apmokymo aibe
% 
% W0=ones(Nneurons,Ninputs)*0.1;     % daugikliu matrica funkciju i neuronus is iejimu
% W1=ones(Noutputs,Nneurons)*0.1;    % daugikliu matrica funkciju i isejimus is neuronu
% B0=ones(Nneurons,1)*0.1;           % pastovios dalys i neuronus

W0=randn(Nneurons,Ninputs)*0.1;     % daugikliu matrica funkciju i neuronus is iejimu
W1=randn(Noutputs,Nneurons)*0.1;    % daugikliu matrica funkciju i isejimus is neuronu
B0=randn(Nneurons,1)*0.1;           % pastovios dalys i neuronu

xmin=-1;xmax=1;dx=(xmax-xmin)/(NLearning-1);
for i=1:Ninputs, aL(i,:)=[xmin:dx:xmax];end
for i=1:Noutputs, sL(i,1:NLearning)=Learning(aL(i,:)); end
phi(1,1:NLearning)=PhiNeuron(aL(1,:));
aL   % apmokymo iejimai
sL   % reikalaujama apmokymo priklausomybe
figure(2);plot(aL,phi,'g*');  % neurono funkcija
figure(1);hold on;plot(aL,sL,'b*'); % apmokymo priklausomybe
yspan=(max(sL(1,:))-min(sL(1,:)))
axis([xmin xmax -yspan*0.15 max(sL(1,:))+yspan*0.15])
s=OutputANN(aL,W0,B0,W1,@PhiNeuron); % ANN isejimai esant apmokymo iejimams
psi0=TargetFunction(s,sL);

ds=1;  % pradinis zingsnis pagal gradienta

for i=1:200000   % svoriu optimizavimo pagal gradienta ciklas
    [DW1,DW0,DB0]=GradVector(aL,sL,W0,B0,W1,@PhiNeuron);
    W0=W0-DW0(:,:)*ds; B0=B0-DB0*ds; W1=W1-DW1*ds;

    
    s=OutputANN(aL,W0,B0,W1,@PhiNeuron);
    psi=TargetFunction(s,sL);
    if psi<psi0, psi0=psi;
    else, W0=W0+DW0(:,:)*ds; B0=B0+DB0*ds;  W1=W1+DW1*ds; ds=ds/2;  % mazinamas zingsnis
    end
    stest=OutputANN(aL,W0,B0,W1,@PhiNeuron);
    figure(1);hold on;
    if round(i/100)==i/100
%         W1,W0,B0,
        if ~isempty(hndl);delete(hndl);end
        hndl=plot(aL(1,:),stest,'r*');
        fprintf(1,'\niter = %10g, psi=%10.6f, step = %10.3e ',i,psi,ds);pause(0.5)
    end
%             pause
end

end


function phi=PhiNeuron(a)     %     Neurono funkcija
% a -  i neurona ateinanti informacija is iejimu ir pastovios dalies
%
% phi=a;
% phi=a.^2;
% phi=2.^a;
%  phi=atan(6*a);
phi=1./(1+exp(-10*a));
% phi=exp(-15*a.^2);
return
end


function sL=Learning(a)    %        ANN apmokymo funkcija
%   a  -  iejimu matrica [Ninputs,Nlearning]
%
% sL=-5*a.^3+10;
% sL=atan(3*a);
% sL=abs(a-0.2);
% sL=sin(3*a)+3;
% sL=sin(10*a)+3;   % phi=exp(-15*a.^2);
sL=sign(3*(a+0.5))-sign(2*(a-0.7));
return, end


function s=OutputANN(a,W0,B0,W1,PhiNeuron)   %   ANN isejimo funkcija
% a  -  iejimas
% W0(Nneurons,Ninputs) -  perdavimo koeficientai is iejimu i neuronus
% W1(Noutputs,Nneurons) -  perdavimo koeficientai is neuronus i isejimus
% B0(Nneurons) -  neuronu iejimu pastovios dalys
% phiANN - neurono funkcijos vardas
%
% disp('OutputANN pradeda')
%     s=W1*(PhiNeuron(W0*a)+B0);
s=W1*(PhiNeuron(W0*a+B0));
% disp('OutputANN baige')
return, end


function [DW1,DW0,DB0]=GradVector(aL,sL,W0,B0,W1,PhiNeuron)    % Tikslo funkcijos gradientas
% aL,sL - apmokymo iejimai ir isejimai matricos, stulpeliai yra apmokymo atvejai
% s - esamo ANN atsakas i apmokymo iejimus
% W0(Nneurons,Ninputs) -  perdavimo koeficientai is iejimu i neuronus
% W1(Noutputs,Nneurons) -  perdavimo koeficientai is neuronu i isejimus
% B0(Nneurons) -  neuronu iejimu pastovios dalys
% phiANN - neurono funkcijos vardas
% ***********  skaiciuojama vienam isejimui D=[dF/dW1, dF/dW0, dF/dB0]
% ***********  kai isejimu daugiau D, turi kelias eilutes

% disp('---------------------GradVector   pradeda')
Ninputs=size(W0,2);       % iejimu skaicius
Nneurons=size(W0,1);      % vienintelio sluoksnio neuronu skaicius
Noutputs=size(W1,1);      % isejimu skaicius
NLearning=size(aL,2);     % apmokymo aibes dydis

if size(aL,2)~=size(sL,2), disp('nesuderinti aL ir sL stulpeliu skaiciai ');end

DB0=zeros( size(B0));DW0=zeros( size(W0));DW1=zeros( size(W1));

% for iii=1:Noutputs  %---------------------

for k=1:NLearning
    dds=0.000001;
    s1=OutputANN(aL(:,k),W0,B0,W1,PhiNeuron);
    
    for i=1:Nneurons
        W11=W1;W11(i)=W11(i)+dds;
%         OutputANN(aL(:,k),W0,B0,W11,PhiNeuron);
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
amax=max([max(abs(DW1)),max(max(abs(DW0))),max(abs(DB0))]);
DW1=DW1/amax;  DW0=DW0/amax; DB0=DB0/amax;

%     DW1=DW1/max(abs(DW1));DW0=DW0/max(max(abs(DW0)));DB0=DB0/max(abs(DB0));



% end  %---------------------
% disp('---------------------GradVector   baige')
return,end


function psi=TargetFunction(s,sL)
% sL - apmokymo iejimau matrica, stulpeliai yra apmokymo atvejai
% s - esamo ANN atsakas i apmokymo iejimus
% for iii=1:size(sL,1)
%     psi(iii)=sum((s(iii,:)-sL(iii,:)).^2)/2;
% end
psi = sum((s-sL).^2)/2; 
return,end
