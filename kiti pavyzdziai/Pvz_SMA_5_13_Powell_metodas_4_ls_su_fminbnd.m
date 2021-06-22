
% Powell metodas

function pagrindine
clc,close all
global grad
global x
eps=1e-6
itmax=300
x=[1000;100;2000;0];
% x=[-0.828147  6.953632  -4.33552  2.98925]'
% x=[1.88415  -0.849997  -3.52124  -2.1538]'

n=length(x)
gradientai=diag(ones(n,1))

x0=x; % prisimename pradini artini
fff=target(x)

for iii=1:itmax
    for i=1:n
        grad=gradientai(i,:)';
        fff0=target(x); % reiksme pries minimizavima
        [alpha,fff]=fminbnd(@fun,-500,500);
        x=x+alpha*grad;
        deltaf(i)=fff-fff0;
    end
    tikslumas=norm(fff);
    fprintf(1,'\n iteracija %d  tikslumas %g ',iii,tikslumas);
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
    gradientai(ind,:)=(x-x0)/norm((x-x0));
    x0=x;
  
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
    
% Funkcija minimizavimui pagal krypti

function rez=fun(alpha)
global grad
global x
rez=target(x+alpha*grad);
return
end

    