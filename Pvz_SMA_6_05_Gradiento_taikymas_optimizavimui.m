
% Greiciausio nusileidimo metodas, kai gradientas apskaiciuojamas
%   diferencijuojant analitiskai (option=1)
%   diferencijuojant skaitiskai (option=2)
%   Broideno artiniu, pagal reikalà atnaujinant gradiento reikðmæ skaitiniu
%                                                diferencijavimu (option=3)
function kelio_optimizavimas
clc,close all
scrsz = get(0,'ScreenSize')
        height=[2    3    2   2]
        alfax=[ 0.5 0.5 0.25 0.5];
        ax=   [  0   3   -2   -2];
        alfay=[ 0.5 0.5 0.5 0.25];
        ay=   [  0   -2   3  -3 ];
        
%         height=[ 3     2]
%         alfax=[ 0.5  0.5];
%         ax=   [  3   -2];
%         alfay=[ 0.5  0.25];
%         ay=   [  -2   -3 ];
%         
xx=[-5:0.1:5];yy=[-6:0.1:6]; % grid
[X,Y]=meshgrid(xx,yy); Z=kalnas(X,Y); % mesh and mountains function
fig1=figure(1);set(fig1,'Position',[50 scrsz(4)/1.8 scrsz(3)/3 scrsz(4)/3]); 
hold on,grid on,axis equal,axis([min(xx)-1 max(xx)+1 min(yy)-1 max(yy)+1 -0.1 max(1,max(height)*1.2)]);view([1 1 1]);xlabel('x'),ylabel('y');
surf(X,Y,Z,'FaceAlpha',0.8,'EdgeColor',0.8*[1 1 1],'FaceColor',0.9*[0.6 1 0],'FaceLighting','gouraud');light('Position',[5,5,5]);


if 0  % du pradinio kelio varijantai 
    npoints1=30;xp=5;yp=-5;xg=0;yg=2;xgf=-5;ygf=5;
    npoints2=30;
    x01=[xp:(xg-xp)/npoints1:xg];x02=[xg:(xgf-xg)/npoints2:xgf];x0=[x01(1:end-1),x02]
    y01=[yp:(yg-yp)/npoints1:yg];y02=[yg:(ygf-yg)/npoints2:ygf];y0=[y01(1:end-1),y02]
    lng=min(length(x0),length(y0));
    x0=x0(1:lng);y0=y0(1:lng);
else
    npoints=60;xp=5;yp=-5;xg=-5;yg=5.5; 
    x0=[xp:(xg-xp)/npoints:xg];y0=[yp:(yg-yp)/npoints:yg];
    lng=min(length(x0),length(y0));
    x0=x0(1:lng);y0=y0(1:lng);
end

option=3;       % gradiento skaiciavimo budas
switch option    
    case 1,  eilute='analitinis diferencijavimas'
    case 2,  eilute='skaitinis diferencijavimas'
    case 3,  eilute='Broyden+skaitinis diferencijavimas'
    otherwise,  'nezinomas gradiento skaiciavimo budas',return
end          

z0=kalnas(x0,y0);  % pradinis artinys
n=length(x0);  % number of points including the first and the last
hndl=plot3(x0,y0,z0,'r-*');
ABF=5;    % penalty coefficient due to altitude change 
pradinis_kelias=kelio_ilgis(x0,y0);
pradine_tikslo_funkcija=target(x0,y0,ABF);
title(sprintf('%s: kelio ilgis= %g\ntikslo funkcija= %g,  ABF=%d',eilute,pradinis_kelias,pradine_tikslo_funkcija,ABF))
pause

eps=1e-5; % tikslumas optimaliam sprendiniui
EPS=1e-1; % tikslumas, prie kurio galima padidinti zingsni
step0=0.2;step=step0; % zingsnio gradiento kryptimi pradinis dydis
itmax=1000  % kiek kartu galima keisti krypti
x=x0;y=y0;  % pradinis artinys
hndl1=[];
          
fff=target(x,y,ABF);% pradine funkcijos reiksme

if option == 3 % pradine gradiento reiksme (skaitinis diferencijavimas):
    hhh=1e-6;
    for i=2:n-1,xxx=x;xxx(i)=xxx(i)+hhh;fffx=target(xxx,y,ABF);gradx(i-1)=(fffx-fff)/hhh;end
    for i=2:n-1,yyy=y;yyy(i)=yyy(i)+hhh;fffy=target(x,yyy,ABF);grady(i-1)=(fffy-fff)/hhh;end    
    grad=[gradx,grady];
    iSD=0;
end

for iii=1:itmax   % ciklas per nusileidimo kryptis         
    xx=x;yy=y;fff0=target(xx,yy,ABF);  % argumentu ir funkcijos reiksmes pries 
        % pradedant nusileidima, naudojama Broideno formuleje po to, kai 
        % bus baigtas nusileidimas viena kryptimi
    
switch option    % gradiento skaiciavimo budas
    case 1   % jeigu gradientas skaiciuojamas analitiskai 
          [gradx,grady]=gradient(x,y,ABF); grad=[gradx,grady]; 
    case 2   % jeigu gradientas skaiciuojamas skaitiniu diferencijavimu 
        hhh=1e-4;
        for i=2:n-1,xxx=x;xxx(i)=xxx(i)+hhh;fffx=target(xxx,y,ABF);gradx(i-1)=(fffx-fff0)/hhh;end
        for i=2:n-1,yyy=y;yyy(i)=yyy(i)+hhh;fffy=target(x,yyy,ABF);grady(i-1)=(fffy-fff0)/hhh;end  
        grad=[gradx,grady];
    case 3   % jeigu taikoma Broideno formule
          ; % galima skaiciuoti tik atlikus zingsni, taikoma ciklo gale   
    otherwise
        'nezinomas gradiento skaiciavimo budas'
end

    gradn=grad/norm(grad); % gradiento norma
            % gradiento projekcija i statmena pradiniam keliui krypti:
    nrm= [y0(3:end)-y0(1:end-2);x0(1:end-2)-x0(3:end)];  % normal to the path right 
    grd=[gradn(1:n-2);gradn(n-1:2*n-4)]; 
    for ii=1:n-2, grd(:,ii)=dot(nrm(:,ii),grd(:,ii))*nrm(:,ii)/norm(nrm(:,ii)); end
    gradn=[grd(1,:),grd(2,:)]; 
       
    fff=fff0;
    step=step0; % kiekvieno nusileidimo pradzioje imama pradine zingsnio reiksme
    nsteps=15;if option == 3, nsteps=4;end 
    for iiii=1:nsteps  % ---------------nusileidimas funkcijos mazejimo kryptimi
        x(2:n-1)=x(2:n-1)-gradn(1:n-2)*step; 
        y(2:n-1)=y(2:n-1)-gradn(n-1:2*n-4)*step;
        fff1=target(x,y,ABF);     
        tikslumas=abs(fff-fff1)/(abs(fff)+abs(fff1));
        fprintf(1,'\n %s : kryptis %d  T.f.reiksme %g tikslumas %g  zingsnis % g',eilute,iii,fff1,tikslumas,step);
        if tikslumas < eps,break,end
        if fff1 > fff, x(2:end-1)=x(2:end-1)+gradn(1:n-2)*step; y(2:end-1)=y(2:end-1)+gradn(n-1:2*n-4)*step; step=step/2; 
        else, fff=fff1; if tikslumas < EPS, step=step*2;end
        end
    end % ---------------- nusileidimo pagal krypti pabaiga
    
                % vizualizavimas:
            if ~isempty(hndl1),delete(hndl1);end 
            zz=kalnas(x,y);
            hndl1=plot3(x,y,zz,'k*');
            kelias=kelio_ilgis(x,y);
            title(sprintf('%s :kelio ilgis= %g  --> %g\n tikslo funkcija= %g --> %g,  ABF=%g',...
            eilute,pradinis_kelias,kelias,pradine_tikslo_funkcija,fff1,ABF));
        pause(0.01);
        
        if option == 3 % jeigu taikoma Broideno gradiento aproksimacija
            s=[x(2:end-1)-xx(2:end-1),y(2:end-1)-yy(2:end-1)];
            if abs(fff1-fff0)/(abs(fff1)+abs(fff0))<eps | norm(s) < eps, % Gradiento atnaujinimas skaitiniu diferencijavimu
                iSD=iSD+1;fprintf(1,'\n  %d  Gradiento atnaujinimas skaitiniu diferencijavimu',iSD);
                hhh=1e-6;
                for i=2:n-1,xxx=x;xxx(i)=xxx(i)+hhh;fffx=target(xxx,y,ABF);gradx(i-1)=(fffx-fff)/hhh;end
                for i=2:n-1,yyy=y;yyy(i)=yyy(i)+hhh;fffy=target(x,yyy,ABF);grady(i-1)=(fffy-fff)/hhh;end  
                grad=[gradx,grady];
            else  % Gradiento atnaujinimas pagal Broideno formule 
                grad=grad+(fff1-fff0-grad*s')*s/(s*s');
            end
        end

end
    return

%   Mountains function 
    function fff=kalnas(x,y)
        fff=0;for jj=1:length(alfax), fff=fff+height(jj)./(1+alfax(jj)*(x-ax(jj)).^2+alfay(jj).*(y-ay(jj)).^2); end
    return
    end
    
%   Target function
    function rez=target(x,y,AB)
        z=AB*kalnas(x,y);
%         rez=sum(sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2));
                rez=sum((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);
    return
    end

%   Kelio ilgio funkcija
    function rez=kelio_ilgis(x,y)
        z=kalnas(x,y);
        rez=sum(sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2));          
    return
    end

 
%   Derivatives of mountain function
    function [dfffx,dfffy]=Dkalnas(x,y)
        dfffx=0;dfffy=0;
        for j=1:length(alfax), 
            dfffx=dfffx-height(j)*(1+alfax(j)*(x-ax(j)).^2+alfay(j).*(y-ay(j)).^2).^(-2).*2*alfax(j).*(x-ax(j));
            dfffy=dfffy-height(j)*(1+alfax(j)*(x-ax(j)).^2+alfay(j).*(y-ay(j)).^2).^(-2).*2*alfay(j).*(y-ay(j));
        end
    return
    end
 
    % Gradient
    function [rezx,rezy]=gradient(x,y,AB)
        % returns 2 vectors of length n-2
        z=kalnas(x,y);
        [dfffx1,dfffy1]=Dkalnas(x,y);
        dfffx=dfffx1(2:end-1);dfffy=dfffy1(2:end-1);
        AB=AB^2;
        rezx=2*(x(2:end-1)-x(1:end-2))-2*(x(3:end)-x(2:end-1))+2*AB*(z(2:end-1)-z(1:end-2)).*dfffx-2*AB*(z(3:end)-z(2:end-1)).*dfffx;
        rezy=2*(y(2:end-1)-y(1:end-2))-2*(y(3:end)-y(2:end-1))+2*AB*(z(2:end-1)-z(1:end-2)).*dfffy-2*AB*(z(3:end)-z(2:end-1)).*dfffy; 
    return
    end
end