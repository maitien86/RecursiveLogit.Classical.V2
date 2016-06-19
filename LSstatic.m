globalVar;
t = 1;
A = zeros(100,1);
B = zeros(1832,1);
for n = 1: nbobs
    LS = LSatt(n).value;
    path = Obs(n,:);
    lpath = size(find(path),2);
    for i = 2:lpath - 1
       A(t) = Uturn(path(i),path(i+1));% LS(path(i),path(i+1));
       t = t+1;
    end
end

for n = 1: nbobs
    LS = LSatt(n).value;
    path = Obs(n,:);
    lpath = size(find(path),2);
    for i = 2:lpath - 1
       B(n) = B(n) + LS(path(i),path(i+1));% LS(path(i),path(i+1));
    end
end
