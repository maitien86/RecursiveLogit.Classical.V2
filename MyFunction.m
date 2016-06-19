function value = MyFunction(x)
    global Msave;
    n = size(Msave,1);
    b = zeros(n,1);
    b(n)=1;
    b_sp = sparse(b);    
    value = x - Msave * x - b_sp;
end