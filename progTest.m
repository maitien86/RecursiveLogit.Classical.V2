global Op;
[val1, Op.grad ] = getLL_test();
Op.grad
h = 0.000000001 * [0 1 0 0 0]';
Op.x = Op.x +h;
[val2 gr] = getLL_test();
(val2 - val1)/0.000000001