
% Broideno metodas
function pagrindine
clc,close all
scrsz = get(0,'ScreenSize')

x=[-8:0.6:8];y=[-8:0.6:8];
Z=pavirsius(@f,x,y);
% [left, bottom, width, height]
fig1=figure(1);set(fig1,'Position',[50 scrsz(4)/1.8 scrsz(3)/3 scrsz(4)/3],'Color','w');
hold on,grid on,axis equal,axis([min(x) max(x) min(y) max(y) -5 5]);view([0 0 1]);xlabel('x'),ylabel('y');
mesh(x,y,Z(:,:,1)','FaceAlpha',0.1,'FaceColor','r','EdgeColor','r');contour(x,y,Z(:,:,1)',[0 0],'LineWidth',1.5,'LineColor','r');
mesh(x,y,Z(:,:,2)','FaceAlpha',0.1,'FaceColor','b','EdgeColor','b');contour(x,y,Z(:,:,2)',[0 0],'LineWidth',1.5,'LineColor','b');
xx=axis; fill([xx(1),xx(1),xx(2),xx(2)],[xx(3),xx(4),xx(4),xx(3)],'b','FaceAlpha',0.1);

eps=1e-5;itmax=200;
x=[-3.5;-6];  % pradinis artinys
% x=[-3;-1.2];
ff=f(x);
figure(1);plot3(x(1),x(2),0,'b*');line([x(1),x(1)],[x(2),x(2)],[0,ff(1)],'Color','black');line([x(1),x(1)],[x(2),x(2)],[0,ff(2)],'Color','black');

n=length(x);

 
% Pradinio Jakobio matricos artinio apskaiciavimas: **********
if 1  % pradinis A artinys skaiciuojamas pagal skaitinio diferencijavimo formule;
    dx=sum(abs(x))*1e-5; f0=f(x);
    for i=1:n
        x1=x;x1(i)=x1(i)+dx; f1=f(x1); A(:,i)=(f1-f0)/dx;
    end
else,  % pradinis A artinys yra parinkta istrizainine matrica 
    A=eye(n)*10  %*10  *15 *8) 
end  
    % pradines gradientu linijos
    dff=A;
    line([x(1),x(1)+dff(1,1)],[x(2),x(2)+dff(1,2)],[0 0],'Color','magenta','LineWidth',2.5)
    line([x(1),x(1)+dff(2,1)],[x(2),x(2)+dff(2,2)],[0 0],'Color','cyan','LineWidth',2.5)

% Broideno metodo iteracijos:
ff=f(x);  % pradine funkcijos reiksme
for iii=1:itmax
    deltax=-A\ff; x1=x+deltax, ff1=f(x1); A=A+(ff1-ff-A*deltax)*deltax'/(deltax'*deltax),
    figure(1);plot3(x1(1),x1(2),0,'k*');line([x(1),x1(1)],[x(2),x1(2)],[0,0],'Color','k','LineWidth',1.5);
    line([x(1),x1(1)],[x(2),x1(2)],[ff(1),0*ff1(1)],'Color','k','LineWidth',2.5);
    line([x1(1),x1(1)],[x1(2),x1(2)],[0,ff1(1)],'Color','k');
    line([x(1),x1(1)],[x(2),x1(2)],[ff(2),0*ff1(2)],'Color','k','LineWidth',2.5);
    line([x1(1),x1(1)],[x1(2),x1(2)],[0,ff1(2)],'Color','k');

    % gradientu linijos
    dff=A; 
    line([x1(1),x1(1)+dff(1,1)],[x1(2),x1(2)+dff(1,2)],[0 0],'Color','magenta','LineWidth',2.5)
    line([x1(1),x1(1)+dff(2,1)],[x1(2),x1(2)+dff(2,2)],[0 0],'Color','cyan','LineWidth',2.5)
    figure(1);plot3(x1(1),x1(2),0,'r*');line([x(1),x1(1)],[x(2),x1(2)],[0,0],'Color','red');
    pause;
 
    tikslumas=norm(deltax)/(norm(x)+norm(deltax));
    
    fprintf(1,'\n iteracija %d  tikslumas %g',iii,tikslumas);
    if tikslumas < eps,         fprintf(1,'\n sprendinys x =');fprintf(1,'  %g',x); plot3(x(1),x(2),0,'rp'); break;
    elseif iii == itmax,fprintf(1,'\n ****tikslumas nepasiektas. Paskutinis artinys x =');fprintf(1,'  %g',x); plot3(x(1),x(2),0,'gp'); break;
    end
    
    ff=ff1;x=x1;
    
end

return
end


%   Lygciu sistemos funkcija 
    function fff=f(x)
    fff=0.1*[x(1)^2+x(2)^2-2;
         -(x(1)-3)^2-x(2)^2+3];
    return
    end
    

    function Z=pavirsius(funk,x,y)
        for i=1:length(x), for j=1:length(y)
              Z(i,j,1:2)=funk([x(i),y(j)]);
        end, end
    return
    end