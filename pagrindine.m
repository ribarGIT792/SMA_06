
% Greiciausio nusileidimo metodas, kai gradientas apskaiciuojamas
% kaip Broideno artinys

function pagrindine
clc,
close all
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

xx=[-5:0.1:5];yy=[-6:0.1:6]; % grid
[X,Y]=meshgrid(xx,yy); Z=kalnas(X,Y); % mesh and mountains function
fig1=figure(2);
set(fig1,'Position',[50 scrsz(4)/1.8 scrsz(3)/3 scrsz(4)/3]);
hold on,grid on,axis equal,axis([min(xx) max(xx) min(yy) max(yy) -0.1 max(1,max(height)*1.2)]);view([1 1 1]);xlabel('x'),ylabel('y');
surf(X,Y,Z,'FaceAlpha',0.8,'EdgeColor',0.8*[1 1 1],'FaceColor',0.9*[0.6 1 0],'FaceLighting','gouraud');light('Position',[5,5,5]);

npoints1=30;xp=5;yp=-5;xg=0;yg=2;xgf=-5;ygf=5;
% xp=3;yp=-5;xg=3;yg=2;xgf=3-1e-6;ygf=5;
npoints2=30;
x01=[xp:(xg-xp)/npoints1:xg];x02=[xg:(xgf-xg)/npoints2:xgf];x0=[x01(1:end-1),x02]
y01=[yp:(yg-yp)/npoints1:yg];y02=[yg:(ygf-yg)/npoints2:ygf];y0=[y01(1:end-1),y02]
lng=min(length(x0),length(y0));
x0=x0(1:lng);y0=y0(1:lng);
z0=kalnas(x0,y0);  % intial trial solution
n=length(x0);  % number of points including the first and the last
hndl=plot3(x0,y0,z0,'r-*');
ABF=5;%30;    % penalty coefficient due to altidude change
pradinis_kelio_ilgis=sqrt(target(x0,y0,1));
pradine_tikslo_funkcija=target(x0,y0,ABF);
title(sprintf('kelio ilgis= %g\ntikslo funkcija= %g,  ABF=%d',pradinis_kelio_ilgis,pradine_tikslo_funkcija,ABF))
pause

eps=1e-5; % precision for optimal solution
EPS=1e-1;
step0=0.2;step=step0;
itmax=4000  % kiek kartu galima keisti krypti
x=x0;y=y0;  % pradinis artinys
hndl1=[];

fff=target(x,y,ABF);% pradine funkcijos reiksme
% vieno gradiento elemento (nr) patikrinimas su skaitiniu iverciu
% nr = 10;
% for i = -5:0.1:5
%     figure(5); hold on; 
%     xx = x; xx(nr)=i;
%     fx = target(xx,y,ABF);
%     plot(i, fx,'b.');
%     
%     [gradx,grady]=gradient(xx,y,ABF);
%     figure(6); hold on; 
%     plot(i,gradx(nr-1),'bo');
% %     plot(i,grady(nr-1),'r.');
%     dx = 1e-9;
%     xx(nr)=i+dx;
%      fx1 = target(xx,y,ABF);
%      plot(i,(fx1-fx)/dx,'m*');
% end
% pradine gradiento reiksme (skaitinis diferencijavimas):
hhh=0.01;
for i=2:n-1,xx=x;xx(i)=xx(i)+hhh;fffx=target(xx,y,ABF);gradx(i-1)=(fffx-fff)/hhh;end
for i=2:n-1,yy=y;yy(i)=yy(i)+hhh;fffy=target(x,yy,ABF);grady(i-1)=(fffy-fff)/hhh;end
grad=[gradx,grady];

grad = ones(size(grad)); % initial values for Broyden gradient

%     step=step0;
for iii=1:itmax
    xxx=x;yyy=y;
    fff0=target(x,y,ABF);  % funkcijos reiksme pries pradedant nusileidima
    
    switch 4   % gradiento skaiciavimo budas
        case 1   % jeigu gradientas skaiciuojamas analitiskai
            [gradx,grady]=gradient(x,y,ABF); 
            grad=[gradx,grady];
        case 2   % jeigu gradientas skaiciuojamas skaitiniu diferencijavimu
%             hhh=0.01;
             hhh = 1e-3;
            for i=2:n-1,xx=x;xx(i)=xx(i)+hhh;fffx=target(xx,y,ABF);gradx(i-1)=(fffx-fff0)/hhh;end
            for i=2:n-1,yy=y;yy(i)=yy(i)+hhh;fffy=target(x,yy,ABF);grady(i-1)=(fffy-fff0)/hhh;end
            grad=[gradx,grady];

        case 3
            % aukstesnes tikslumo eiles isvestines ivertis
            hhh = 1e-6;
            for i=2:n-1,
                xxb2 = x; xxb1 = x; xxa1 = x; xxa2 = x; 
                xxb2(i) = xxb2(i)-2*hhh; xxb1(i) = xxb1(i)-hhh; 
                xxa2(i) = xxa2(i)+2*hhh; xxa1(i) = xxa1(i)+hhh; 
                fffxb2=target(xxb2,y,ABF); fffxb1=target(xxb1,y,ABF);
                fffxa2=target(xxa2,y,ABF); fffxa1=target(xxa1,y,ABF);
                gradx(i-1) = (fffxb2-8*fffxb1+8*fffxa1-fffxa2)/(12*hhh);
            end
            for i=2:n-1,
                yyb2 = y; yyb1 = y; yya1 = y; yya2 = y; 
                yyb2(i) = yyb2(i)-2*hhh; yyb1(i) = yyb1(i)-hhh; 
                yya2(i) = yya2(i)+2*hhh; yya1(i) = yya1(i)+hhh; 
                fffyb2=target(x, yyb2,ABF); fffyb1=target(x, yyb1,ABF);
                fffya2=target(x, yya2,ABF); fffya1=target(x, yya1,ABF);
                grady(i-1) = (fffyb2-8*fffyb1+8*fffya1-fffya2)/(12*hhh);
            end
            grad=[gradx,grady];
        otherwise   % jeigu taikoma Broideno formule
            ;
    end
    gradn=grad/norm(grad); % gradiento norma
    % gradiento projekcija i statmena pradiniam keliui krypti:
    nrm= [y0(3:end)-y0(1:end-2);x0(1:end-2)-x0(3:end)];  % normal to the path right
    % CX: normal to current path
%     nrm= [y(3:end)-y(1:end-2);x(1:end-2)-x(3:end)];  % normal to the path right 
    grd=[gradn(1:n-2);gradn(n-1:2*n-4)];
    for ii=1:n-2, grd(:,ii)=dot(nrm(:,ii),grd(:,ii))*nrm(:,ii)/norm(nrm(:,ii)); end
    gradn=[grd(1,:),grd(2,:)];
    fff=fff0;
    step=step0;
    for iiii=1:20  % ejimas pagal mazejimo krypti
        x(2:n-1)=x(2:n-1)-gradn(1:n-2)*step;
        y(2:n-1)=y(2:n-1)-gradn(n-1:2*n-4)*step;
        zz=kalnas(x,y);
        
        % vizualizavimas:
%                     if ~isempty(hndl1),delete(hndl1);end
%                     hndl1=plot3(x,y,zz,'k*');
        fff1=target(x,y,ABF);
        %             kelio_ilgis=sqrt(target(x,y,1));
%                 title(sprintf('kelio ilgis= %g  --> %g\n tikslo funkcija= %g --> %g,  ABF=%g',...
%                     pradinis_kelio_ilgis,kelio_ilgis,pradine_tikslo_funkcija,fff1,ABF));
%                 pause(0.01);
        
        tikslumas=norm(fff-fff1)/(norm(fff)+norm(fff1));
        fprintf(1,'\n kryptis %d  T.f.reiksme %g tikslumas %g  zingsnis % g',iii,fff1,tikslumas,step);
%         if tikslumas < eps,break,end
        if fff1 > fff, 
            x(2:end-1)=x(2:end-1)+gradn(1:n-2)*step; 
            y(2:end-1)=y(2:end-1)+gradn(n-1:2*n-4)*step;
            step=step/2;
        else, 
            fff=fff1;
            if tikslumas < EPS,
                step=step*2;
            end
        end
    end 
    % vizualizavimas:
    if ~isempty(hndl1),delete(hndl1);end
    zz=kalnas(x,y);
    hndl1=plot3(x,y,zz,'k*-');
    fff1=target(x,y,ABF);
    kelio_ilgis=sqrt(target(x,y,1));
    title(sprintf('kelio ilgis= %g  --> %g\n tikslo funkcija= %g --> %g,  ABF=%g',...
        pradinis_kelio_ilgis,kelio_ilgis,pradine_tikslo_funkcija,fff1,ABF));
    pause(0.01);
    
    sss=[x(2:end-1)-xxx(2:end-1),y(2:end-1)-yyy(2:end-1)];
    if norm(sss) < 1e-12
        break;
    else
        grad=grad+(fff1-fff0-grad*sss')*sss/(sss*sss');    % Broideno formule
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
%         CX: target function no sqrt
        rez=sum((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);
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
%         rezx=2*(x(2:end-1)-x(1:end-2))-2*(x(3:end)-x(2:end-1))+2*AB*(z(2:end-1)-z(1:end-2)).*dfffx-2*AB*(z(3:end)-z(2:end-1)).*dfffx;
%         rezy=2*(y(2:end-1)-y(1:end-2))-2*(y(3:end)-y(2:end-1))+2*AB*(z(2:end-1)-z(1:end-2)).*dfffy-2*AB*(z(3:end)-z(2:end-1)).*dfffy;
        % CX: AB^2
        rezx=2*(x(2:end-1)-x(1:end-2))-2*(x(3:end)-x(2:end-1))+2*AB^2*(z(2:end-1)-z(1:end-2)).*dfffx-2*AB^2*(z(3:end)-z(2:end-1)).*dfffx;
        rezy=2*(y(2:end-1)-y(1:end-2))-2*(y(3:end)-y(2:end-1))+2*AB^2*(z(2:end-1)-z(1:end-2)).*dfffy-2*AB^2*(z(3:end)-z(2:end-1)).*dfffy;
        return
    end
end
%---------------------------------------------***********************-----------------------------
% % Broideno metodas
% function pagrindine
% clc,close all
%
% eps=1e-10
% itmax=100000
% x=10*[1;1;1;1];
% % x=[1;1;1;0];
% n=length(x);
%
% % Pradines Jakobio matricos reiksmes apskaiciavimas:
% dx=sum(abs(x))*1e-5;
% f0=f(x);
% for i=1:n, x1=x; x1(i)=x1(i)+dx; f1=f(x1); A(:,i)=(f1-f0)/dx; end
% % A=-eye(n)*10   %*10  *(-10)
%
% % Broideno metodo iteracijos:
% fi=f(x);  % pradine funkcijos reiksme
% for iii=1:itmax
%     deltax=-A\f(x); x=x+deltax; fi1=f(x); A=A+(fi1-fi-A*deltax)*deltax'/(deltax'*deltax);
%
%     tikslumas=norm(deltax)/(norm(x)+norm(deltax));
%     fprintf(1,'\n iteracija %d  tikslumas %g',iii,tikslumas);
%     if tikslumas < eps
%         fprintf(1,'\n sprendinys x ='); fprintf(1,'  %g',x);
%         fprintf(1,'\n funkcijos reiksme f ='); fprintf(1,'  %g',f(x));
%         break
%     elseif iii == itmax
%         fprintf(1,'\n ****tikslumas nepasiektas. Paskutinis artinys x ='); fprintf(1,'  %g',x);
%         fprintf(1,'\n funkcijos reiksme f ='); fprintf(1,'  %g',f(x));
%         break
%     end
% fi=fi1;
% end
%
%     return
% end
