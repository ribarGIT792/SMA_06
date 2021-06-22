
% Broideno metodas: ortogonalaus poerdvio paaiskinimas
clc;close all
s=[1;1;1];
Z=eye(length(s))-s*s'/(s'*s)
figure(1),hold on,grid on,axis equal
for i=1:3
    plot3([0 Z(1,i)],[0 Z(2,i)],[0 Z(3,i)]);
    text(Z(1,i)/2,Z(2,i)/2,Z(3,i)/2,sprintf('Z(:,%1d)',i),'Color','b','Fontsize',12)
end

    plot3([0 s(1)],[0 s(2)],[0 s(3)],'r');
    text(s(1)/2,s(2)/2,s(3)/2,'s','Color','r','Fontsize',12)