
% Powell metodas

function pagrindine
clc,close all


step=0.01
eps=1e-6
itmax=200
x=[-1;-1;-1;-1];
% x=[-0.828147  6.95332  -4.33552  2.98925]'

n=length(x)
gradientai=diag(ones(n,1))

x0=x; % prisimename pradini artini
fff=target(x)

for iii=1:itmax
    for i=1:n
        grad=gradientai(i,:);
        fff=target(x);
        xxx0=x;fff0=fff; % reiksme pries minimizavima
        for j=1:100  % ejimas pagal j krypti
            deltax=grad/norm(grad)*step;
            x=x+deltax';
            fff1=target(x);
            if fff1>fff && j==1, x=x-deltax';step=-step;continue,end

            if fff1 > fff, x=x-deltax';deltaf(i)=fff-fff0;break,end
            fff=fff1;
        end
    end

    tikslumas=norm(fff);
    fprintf(1,'\n iteracija %d  tikslumas %g zingsnis %g',iii,tikslumas,step);
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
    step=step/1.02;
    if a < 0, 
        gradientai(ind,:)=(x-x0)/norm((x-x0)); 
    else,  x0=x; gradientai=diag(ones(n,1));        
    end
    
end


    return
end

%   Lygciu sistemos funkcija 
function F=f(X) 
 F(1)=X(1)+2*X(2)+X(3)+4*X(4)-20.7;
 F(2)=X(1)^2+2*X(1)*X(2)+X(4)^3-15.88;
 F(3)=X(1)^3+X(3)^2+X(4)-21.218;
 F(4)=3*X(2)+X(3)*X(4)-7.9;
 F=F(:);
 return
end 

%     Tikslo funkcija
    function rez=target(x)
    rez=f(x)'*f(x)/2;
    return
    end
    

    