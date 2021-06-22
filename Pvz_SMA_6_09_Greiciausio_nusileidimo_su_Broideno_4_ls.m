
% Greiciausio nusileidimo metodas su Broideno aproksimacija Jakobio
% matricai   4 ls

function pagrindine
clc,close all

eps=1e-10;
step0=0.1;step=step0;
itmax=1000;

x=[-1;1;-2;1];
%   x=[-0.828147  6.95332  -4.33552  2.98925]'

n=length(x);

% Pradines Jakobio matricos reiksmes apytikslis apskaiciavimas, skaitiskai diferencijuojant :
dx=sum(abs(x))*1e-5; f0=f(x); for i=1:n, x1=x; x1(i)=x1(i)+dx; f1=f(x1); A(:,i)=(f1-f0)/dx; end

for iii=1:itmax  % ciklas nusileidimo paieskos krypciu keitimui 
    grad=gradient(x,A);    fff=target(x);
    xp=x; % prisimename x sios krypties pradzioje 
    for j=1:10  % funkcijos minimizavimas pagal parinkta krypti
        deltax=grad/norm(grad)*step;  x=x-deltax'; fff1=target(x);
        
        if fff1 > fff, x=x+deltax'; step=step/1.05; break; else, fff=fff1;end  %*****var1****
%         if fff1 > fff, x=x+deltax'; step=step/10; break; else, fff=fff1;end  %*****var2****
    end
%     step=step0;  %***var2******

    tikslumas=norm(fff);
    fprintf(1,'\n kryptis %d  tikslumas %g  zingsnis %g',iii,tikslumas,step);
    if tikslumas < eps, fprintf(1,'\n sprendinys x =');fprintf(1,'  %g',x);break;
    elseif iii == itmax,fprintf(1,'\n ****tikslumas nepasiektas. Paskutinis artinys x =');fprintf(1,'  %g',x); break;
    end

   if xp == x,  % Kai Broideno aproksimacijos pritaikyti negalima, 
                %Jakobio matricos reiksme isimties tvarka tenka apskaiciuoti, skaitiskai diferencijuojant: 
        fprintf(1,'\n nepavyksta pakeisti krypties  \n paskutinis artinys x ='); fprintf(1,'  %g',x); 
        dx=sum(abs(x))*1e-5; f0=f(x);for i=1:n, x1=x; x1(i)=x1(i)+dx; f1=f(x1); A(:,i)=(f1-f0)/dx; end
   else
        A=df_Broiden(xp,x,A); 
   end
    
   
end
fprintf(1,'\n funkcijos reiksme f =');  fprintf(1,'  %g',f(x));fprintf(1,'\n');

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

%  Jakobio matrica
    function AA=df_Broiden(x,x1,A)
    ff=f(x); ff1=f(x1); s=x1-x;
    if s == 0, AA=A; 'A matrica nepakeista' 
    else, AA=A+(ff1-ff-A*s)*s'/(s'*s);
    end
    return
    end

%     Tikslo funkcija
    function rez=target(x)
    rez=f(x)'*f(x)/2;
    return
    end
    
    % Gradientas
    function rez=gradient(x,A)
    rez=f(x)'*A;
    return
    end