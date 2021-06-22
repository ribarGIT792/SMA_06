
% Greiciausio nusileidimo metodas

function pagrindine
clc,close all
scrsz = get(0,'ScreenSize')

x=[-5:0.15:5];y=[-6:0.15:6];
Z=pavirsius(@target,x,y);
fig1=figure(1);set(fig1,'Position',[50 scrsz(4)/1.8 scrsz(3)/3 scrsz(4)/3]);
hold on,grid on,axis equal,axis([min(x) max(x) min(y) max(y) 0 5]);view([1 1 1]);xlabel('x'),ylabel('y');
mesh(x,y,Z(:,:,1)','FaceAlpha',0.2);contour(x,y,Z(:,:,1)',[0 0],'LineWidth',1.5);

Zf=pavirsius(@f,x,y);
contour(x,y,Zf(:,:,1)',[0 0],'LineWidth',1.5,'LineColor','r');
contour(x,y,Zf(:,:,2)',[0 0],'LineWidth',1.5,'LineColor','b');

xx=axis;
fill([xx(1),xx(1),xx(2),xx(2)],[xx(3),xx(4),xx(4),xx(3)],'c','FaceAlpha',0.2);

eps=1e-5
step0=0.5;step=step0;
itmax=30  % kiek kartu galima keisti krypti
% x=[6;4];  % patenka i lokalu minimuma
x=[-1;4];  
% x=[-2;4];
% x=[-1;-4];

plot3(x(1),x(2),0,'g*'); 
pause
for iii=1:itmax
    
    grad=gradient(x);   
    fff=target(x);
    for j=1:30  % ejimas pagal mazejimo krypti
        deltax=grad/norm(grad)*step; 
        x=x-deltax';
            % vizualizavimas:
            plot3(x(1),x(2),0,'ro');
            h = findobj(gca,'Type','line');
            pause(0.01);delete(h(1));    
            plot3(x(1),x(2),0,'ko');
        fff1=target(x);
        if fff1 > fff, x=x+deltax'; step=step/10; else, fff=fff1;end   
    end
    step=step0;


    tikslumas=norm(fff);
    fprintf(1,'\n kryptis %d  tikslumas %g',iii,tikslumas);
    if tikslumas < eps,         fprintf(1,'\n sprendinys x =');fprintf(1,'  %g',x); plot3(x(1),x(2),0,'rp'); break;
    elseif iii == itmax,fprintf(1,'\n ****tikslumas nepasiektas. Paskutinis artinys x =');fprintf(1,'  %g',x); plot3(x(1),x(2),0,'gp'); break;
    end
    
end


    return
end

%   Lygciu sistemos funkcija 
    function fff=f(x)
          fff=[sin(4*x(1))+x(2)^2-2;
           x(1)^2+x(2)^2-(0.5*sin(6*atan(x(2)/x(1)))+1.5)^2];
    return
    end
    
%  Jakobio matrica
    function dfff=df(x)
        dfff=[4*cos(4*x(1)), 2*x(2);
              2*x(1)-2*(0.5*sin(6*atan(x(2)/x(1)))+1.5)*0.5*cos(6*atan(x(2)/x(1)))*6*1/(1+(x(2)/x(1))^2)*(-x(2)/x(1)^2),...
              2*x(2)-2*(0.5*sin(6*atan(x(2)/x(1)))+1.5)*0.5*cos(6*atan(x(2)/x(1)))*6*1/(1+(x(2)/x(1))^2)*(-1/x(1))];
    return
    end
    
%     Tikslo funkcija
    function rez=target(x)
    rez=f(x)'*f(x)/2;
    return
    end
    
    % Gradientas
    function rez=gradient(x)
    rez=f(x)'*df(x);
    return
    end
    
    function ZZ=pavirsius(funk,x,y)
    for i=1:length(x)
        for j=1:length(y)
            aa=funk([x(i),y(j)]); n=length(aa);
            for iii=1:n,  ZZ(i,j,iii)=aa(iii); end
        end
    end
        
    return
    end