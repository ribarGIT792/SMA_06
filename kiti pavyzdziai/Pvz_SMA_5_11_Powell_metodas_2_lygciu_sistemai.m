
% Powell metodas

function pagrindine
clc,close all

x=[-5:0.1:5];y=[-6:0.1:6];
Z=pavirsius(@target,x,y);
figure(1),hold on,grid on,axis equal   
contour(x,y,Z(:,:,1)');
xlabel('x'),ylabel('y')
colorbar;

eps=1e-5
step=3.
itmax=30
x=[-12;-4]
n=length(x)
plot3(x(1),x(2),0,'g*');
gradientai=diag(ones(n,1))

x0=x; % prisimename pradini artini
for iii=1:itmax
    for i=1:n
        grad=gradientai(i,:);
        fff=target(x);
        xxx0=x;fff0=fff; % reiksme pries minimizavima
        for j=1:100  % ejimas pagal j krypti
            deltax=grad/norm(grad)*step;
            x=x-deltax';
            fff1=target(x);
            if fff1>fff && j==1, x=x+deltax';step=-step;continue,end
                % vizualizavimas:
                plot3(x(1),x(2),0,'ro');
                h = findobj(gca,'Type','line');
                pause(1);delete(h(1));    
                plot3(x(1),x(2),0,'ko');
            if fff1 > fff, x=x+deltax';deltaf(i)=fff-fff0;break,end
            fff=fff1;
        end
    end

    tikslumas=norm(fff);
    fprintf(1,'\n iteracija %d  tikslumas %g',iii,tikslumas);
    if tikslumas < eps
        fprintf(1,'\n sprendinys x =');
        fprintf(1,'  %g',x);
        break
    elseif iii == itmax
        fprintf(1,'\n ****tikslumas nepasiektas. Paskutinis artinys x =');
        fprintf(1,'  %g',x);
        break
    end
    [a,ind]=min(deltaf);
    step=step/2;
    if a < 0, 
        gradientai(ind,:)=(x-x0)/norm((x-x0)); 
    else, x0=x; gradientai=diag(ones(n,1));        
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

%     Tikslo funkcija
    function rez=target(x)
    rez=f(x)'*f(x)/2;
    return
    end
    
    
    function Z=pavirsius(funk,x,y)
    for i=1:length(x)
        for j=1:length(y)
            Z(i,j)=funk([x(i),y(j)]);
        end
    end
        
    return
    end