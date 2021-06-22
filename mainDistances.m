% reikia rasti x ir y, kad taskus jungianciu liniju (kiekvienas su
% kiekvienu) ilgiu suma butu kuo vienodesne ir artima duotam stygu ilgiui totalLength
function main
clear all;
close all;
clc;
n = 10;
totalLength = 80;
rng(1); % kad atsikartotu taskai (butu galima palyginti)
positions = rand(2*n,1);

figure(1); hold on; grid on;
axis([-3 3 -3 3]);

plot(positions(1:2:end),positions(2:2:end),'ro','MarkerSize',6);
h = [];
step0 = 1;
eps = 1e-9;
maxIt = 500;
x0 = positions;
gradIn = quasiGradientInitial(x0, totalLength, n, 1e-6);


f0 =  target(x0, totalLength, n);
xI = positions +f0*gradIn'/norm(gradIn); % apskaiciuojamas taskas, is kurio
% pateko á positions
step = step0;
grad = gradIn;
prec = 1;
stepmin = step0/1e5;
countChange = 0;
for i = 1:maxIt
    
    
    switch 2
        case 1
            grad = f0*quasiGradientInitial(x0, totalLength, n, step/1e6);
        case 2
            
            if step < stepmin
                % jei reikejo "grubiai" perskaiciuoti krypti daug kartu
                if countChange > 50
                    fprintf(1,'gradientas (grubus) perskaiciuotas %d kartu\n',countChange);
                    break;
                end
                % nepavyksta pakeisti krypties
                grad = f0*quasiGradientInitial(x0, totalLength, n, step/1e6);
                gradIn = grad;
                fprintf(1,'perskaiciuojamas gradientas\n');
                step = step0;
                countChange = countChange + 1;
            else
                grad = f0*quasiGradient(x0, xI, gradIn, totalLength, n);
            end
        otherwise;
    end
    if norm(grad) < eps
        
        fprintf(1,'gradiento norma %6.5f\n',norm(grad));
        break;
    end
    gradPlus = grad/norm(grad);
    x1 = x0-step*gradPlus';
    f1 =  target(x1, totalLength, n);
    
    if(f1>f0)
        step = step / 10; % sumazinamas zingsnis
    else
        prec = norm(x1-x0)/(norm(x1)+norm(x0));
        f0 = f1;
        xI = x0;
        x0 = x1;
        gradIn = grad; % pakeiciama gradiento reiksme
    end
    fprintf(1,'iteration %d, funkcija %6.5f, step %6.8f, totalLength %6.5f, prec %15.10f\n',...
        i,f1, step,sum(calculateDistances(x1,n)), prec);
    if ~isempty(h)
        delete(h);
    end
    h = visualization(x0,n);
    pause(0.0001);
    
    if prec < eps
        fprintf(1,'pasiektas tikslumas\n');
        break;
    end
    
end

end

function f = target(positions, totalLength, n)
distances = calculateDistances(positions,n);
% meanDistance = mean(distances);
% f = sum(distances);
% f = (f-totalLength)^2 + sum((distances-meanDistance).^2);
f = sum((distances-totalLength/(size(distances,1))).^2);
end

function distances = calculateDistances(positions,n)
x = positions(1:2:end);
y = positions(2:2:end);
distances = [];
for i = 1:n-1
    dist = sqrt((x(i)-x(i+1:end)).^2+(y(i)-y(i+1:end)).^2);
    distances = [distances; dist];
end
end

function A = quasiGradient(positions, positionsI, A, totalLength, n)
s = positions - positionsI;
y = target(positions, totalLength, n)-target(positionsI, totalLength, n);
if norm(s)>1e-9
    A = A + (y-A*s)*s'/(s'*s);
end

end

function df = quasiGradientInitial(positions, totalLength, n, dx)
df = zeros(size(positions));
df = df';
f0 = target(positions, totalLength, n);
for i = 1:numel(positions)
    XNew = positions;
    XNew(i) = XNew(i)+dx;
    f1 = target(XNew, totalLength, n);
    df(i) = (f1-f0)/dx;
end
end

function h = visualization(positions,n)
x = positions(1:2:end);
y = positions(2:2:end);

h =  plot(x,y,'bo','MarkerFaceColor',[0 0 1],'MarkerSize',6);
for i = 1:n-1
    for j = i+1:n
        h = [h; plot([x(i) x(j)],[y(i) y(j)],'Color',[0.5 0.5 0.5])];
    end
end
end