function [f,g] = LL(x)
    global Op;
    Op.x = x;
    [f,g] = getLL_v2();
    Op.nFev  = Op.nFev + 1;
end