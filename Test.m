clear all clc
A = [1 2 3 4;
     5 6 7 8;
     9 10 11 12;
     13 14 15 16];
B = [10; 20; 30; 40];
C = [0 1 0 1;
    1 0 1 0;
    0 1 0 1;
    1 0 1 0];
D = max(max(abs(C-A)));
% Ax = b
x = B\A;
max(max(A))