
% Greiciausio nusileidimo metodas su Broideno aproksimacija Jakobio
% matricai

function pagrindine
clc,close all
scrsz = get(0,'ScreenSize')

x=[-5:0.1:5];y=[-6:0.1:6];
Z=pavirsius(@target,x,y);
fig1=figure(1);set(fig1,'Position',[50 scrsz(4)/1.8 scrsz(3)/3 scrsz(4)/3],'Color','w');
hold on,grid on,axis equal,axis([min(x) max(x) min(y) max(y) 0 5]);view([0 0 1]);xlabel('x'),ylabel('y');
mesh(x,y,Z(:,:,1)','FaceAlpha',0.2);
xx=axis; fill([xx(1),xx(1),xx(2),xx(2)],[xx(3),xx(4),xx(4),xx(3)],'c','FaceAlpha',0.2);

Zf=pavirsius(@f,x,y);
contour(x,y,Zf(:,:,1)',[0 0],'LineWidth',1.5,'LineColor','r');
contour(x,y,Zf(:,:,2)',[0 0],'LineWidth',1.5,'LineColor','b');

eps=1e-5
step0=0.5;step=step0;
itmax=200
 
% x=[3;3]; 
% x=[-3;-3];
% x=[3;-3];
% x=[-3;-3];
x=[0;25];
% x=[-1;4];
% x=[12;4];
% x=[6;4]; 
plot3(x(1),x(2),0,'g*');

% Pradines Jakobio matricos reiksmes apskaiciavimas
% Jos reikia Broideno aproksimacijos pirmam zingsniui
dx=sum(abs(x))*1e-5; ff=f(x); n=length(ff);
for i=1:n, x1=x;x1(i)=x1(i)+dx; ff1=f(x1); A(:,i)=(ff1-ff)/dx; end

for iii=1:itmax  % ciklas nusileidimo paieskos krypciu keitimui 
    grad=gradient(x,A);    fff=target(x);
    xp=x; % prisimename dabartini x, iki greiciausio nusileidimo 
    for j=1:30  % ejimas funkcijos mazejimo kryptimi
        deltax=grad/norm(grad)*step;  x=x-deltax';
            % vizualizavimas:
            plot3(x(1),x(2),0,'ro'); h = findobj(gca,'Type','line'); pause(0.01);
            delete(h(1)); plot3(x(1),x(2),0,'ko');
        fff1=target(x);
        if fff1 > fff, x=x+deltax'; step=step/10; else, fff=fff1;end   
    end
    step=step0;

    tikslumas=norm(fff);
    fprintf(1,'\n iteracija %d  tikslumas %g',iii,tikslumas);
    if tikslumas < eps,         fprintf(1,'\n sprendinys x =');fprintf(1,'  %g',x); plot3(x(1),x(2),0,'rp'); break;
    elseif iii == itmax,fprintf(1,'\n ****tikslumas nepasiektas. Paskutinis artinys x =');fprintf(1,'  %g',x); plot3(x(1),x(2),0,'gp'); break;
    end

   if xp == x,  % Kai Broideno aproksimacijos pritaikyti negalima, 
                %Jakobio matricos reiksme isimties tvarka tenka apskaiciuoti, skaitiskai diferencijuojant: 
        fprintf(1,'\n nepavyksta pakeisti krypties  \n paskutinis artinys x ='); fprintf(1,'  %g',x); 
        dx=sum(abs(x))*1e-5; f0=f(x);for i=1:n, x1=x; x1(i)=x1(i)+dx; f1=f(x1); A(:,i)=(f1-f0)/dx; end
   else
        A=df_Broiden(xp,x,A); 
   end
end

    return
end

%   Lygciu sistemos funkcija 
    function fff=f(x)
    fff=[x(1)^2+x(2)^2-2;
         x(1)^2-x(2)^2];
    return
    end
    
%   Jakobio matrica
    function AA=df_Broiden(x,x1,A)
    ff=f(x); ff1=f(x1); s=x1-x;
    if s == 0, AA=A; else, AA=A+(ff1-ff-A*s)*s'/(s'*s); end
    return
    end
    
%   Tikslo funkcija
    function rez=target(x)
    rez=f(x)'*f(x)/2;
    return
    end
    
    % Gradientas
    function rez=gradient(x,A)
    rez=f(x)'*A;
    return
    end
    
    function ZZ=pavirsius(funk,x,y)
    for i=1:length(x), for j=1:length(y)
            aa=funk([x(i),y(j)]); n=length(aa);
            for iii=1:n,  ZZ(i,j,iii)=aa(iii); end
    end, end      
    return
    end