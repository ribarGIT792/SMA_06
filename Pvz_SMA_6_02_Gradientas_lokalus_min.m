
% Gradientas  lokalus minimumas

function pagrindine
clc,close all

global X F
syms x1 x2 X F 
    X=[x1, x2];
    F=-exp(-0.1*x1^2-0.1*x2^2)-2*exp(-0.1*(x1-4)^2-0.1*(x2-6)^2)+4; 
    DF=jacobian(F,X)  % gradiento vektorius


opt=[]
opt='pagal gradienta'
% opt='greiciausio nusileidimo'

if ~isempty(opt)
% funkcija minimizuojama:
x=[-4:0.5:9];y=[-3:0.5:8];
Z=pavirsius(@f,x,y);
figure(2),hold on,grid on,axis equal; view([1 1 1]); 
title(opt);xlabel('x');ylabel('y');zlabel('f');
size(Z) 
mesh(x,y,Z','FaceAlpha',0.6); 
xlabel('x'),ylabel('y')
xx=axis;
fill([xx(1),xx(1),xx(2),xx(2)],[xx(3),xx(4),xx(4),xx(3)],'c','FaceAlpha',0.2)

step=0.6;
P=[9,4];  % taskas
P=[6,-1];
% P=[6,-0.65];
P=[6,-0.5];

% P=[6,0]; 
fnk=eval(subs(F,X,P));
plot3([P(1)],[P(2)],[fnk],'ko')
plot3([P(1)],[P(2)],0,'ro')
plot3([P(1),P(1)],[P(2),P(2)],[0 fnk],'k-')

fnk1=-1e10;
fnk=subs(F,X,P)
for i=1:50
    if fnk > fnk1 | strcmp(opt,'pagal gradienta') 
        grad=eval(subs(DF,X,P));
        grad=-step*grad/norm(grad);
    end
    quiver(P(1),P(2),grad(1),grad(2),'Color','r','LineWidth',2,'Autoscale','off')
    quiver(P(1),P(2),grad(1),0,'Color','r','LineWidth',0.5,'Autoscale','off')
    quiver(P(1),P(2),0,grad(2),'Color','r','LineWidth',0.5,'Autoscale','off')
    fnk1=fnk;
    fnk=eval(subs(F,X,P+grad));
    plot3([P(1)+grad(1)],[P(2)+grad(2)],[fnk],'k*')
    plot3([P(1),P(1)+grad(1)],[P(2),P(2)+grad(2)],[fnk1,fnk],'k-')
    plot3([P(1)+grad(1),P(1)+grad(1)],[P(2)+grad(2),P(2)+grad(2)],[0,fnk],'k--')
    P=P+grad;
    pause;
end
end
    return
end

% ---------------------------------------------------------------------
    function fff=f(x)
    %   Si funkcija reikalinga, kad butu galima jos varda perduoti 
    %   kitai funkcijai faktiniu parametru sarase 
    global X F
    fff=eval(subs(F,X,x));
    return
    end
    
 % ---------------------------------------------------------------------
    function Z=pavirsius(funk,x,y)
    % fukcija suformuoja dvieju kintamuju funkcijos masyva vaizdavimui
        for i=1:length(x), for j=1:length(y), Z(i,j)=funk([x(i),y(j)]);end,end
    return,end
% ---------------------------------------------------------------------