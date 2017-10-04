load ex5-2-1;

[m,l] = size(A)
B = [3 4 6]
 N = setdiff(1:l,B) 
 [L,U] = lu(A(:,B)) 
 x_B = U\(L\b)
u = L'\(U'\p(B))
c=p(N)-A(:,N)'*u
s = 2
d=U\(L\A(:,N(s)))

clear;
load ex5-2-6
B = [1 2]
N = [-3 4 -5]
[m,l] = size(A)
x = zeros(l,1)
x(3:5) = [0 5 0]'
[L,U] = lu(A(:,B))
x(B) = U\(L\(b-A(:,abs(N))*x(abs(N))))
 u = L'\(U'\p(B))
 c = p(abs(N))-A(:,abs(N))'*u
 s = 1
 d = U\(L\(A(:,abs(N(s)))))
 
 load ex5-2-6
[m,l] = size(A)
A = [A zeros(m,l); eye(l) eye(l)]%%%??????
b = [b; ub]
p = [p; zeros(l,1)]
B = [1 2 4 6 7 8 10]
[x_B,B] = rsm(A,b,p,B)
x = zeros(2*l,1) 
x(B) = x_B

%5-2-2
load ex5-2-2;
p = [-p;0;0;0;0]
A = horzcat(A,eye(4))
[m l]=size(A)
B= [4 5 6 7]
N=setdiff(1:l,B)
[x_B,B] = rsm(A,b,p,B)
z = p(B)' * x_B

%5-2-8
clear;
load ex5-2-8;
B = [2, 5];
N = [-1 -3, 4, -6, -7, -8];
[m l] = size(A)
%verify the basis corresponds to basic feasible solution
x = zeros(m,1)
x(abs(N)) = [0 0 15 0 0 0]'
[L U] = lu(A(:,B))
x(B) = U\(L\(b-A(:,abs(N))*x(abs(N))))
xstart = [0, 5.6, 0, 15, 1.6, 0, 0, 0]'
sum(abs(xstart-x))

%calc c values for N set
u = L'\(U'\p(B))
c = p(abs(N))-A(:,abs(N))'*u

%only x6 and x7 are eligible to enter the basis
%pick x7
s = 5
d = U\(L\(A(:,abs(N(s)))))

%Let x7=lambda>0
%x2 = 5.6 - 0.2*lambda
%x5 = 1.6 + 0.8*lambda
%As lambda increases, x2 approaches lb, which is -inf
%Thus, x5 is the blocking variable.
%At its ub of 5, lambda=1/2. Update values.
lambda=1/2
x(7)=x(7)+lambda
x(B)=x(B)-lambda*d
B(2)=7
N(5)=5

%Calculate c
[L U] = lu(A(:,B));
u = L'\(U'\p(B));
c = p(abs(N))-A(:,abs(N))'*u,

%pick x6 enter basis, only option
s = 4;
d = U\(L\(A(:,abs(N(s))))),

%Let x6=lambda>0
%x2 = 5.5 + 0.5*lambda
%x7 = 0.5 + 0.5*lambda
%As lambda increases, x7 approaches ub of 10, lambda=19
%As lambda increases, x2 approaches ub of 6, lambda=1
%Thus, x2 is the blocking variable
lambda=1;

x(6)=x(6)+lambda;
x(B)=x(B)-lambda*d,
B(1)=6, N(4)=2,


%Calculate c
[L U] = lu(A(:,B));
u = L'\(U'\p(B));
c = p(abs(N))-A(:,abs(N))'*u,

%optimal solution x, all signs of c and N differ
%solution x=[0, 6, 0, 15, 2, 1, 1, 0]
%Optional code to check below
%B = [2, 5];N = [-1 -3, 4, -6, -7, -8];
%[x2 B N] = rsmbdd(A,b,p,lb,ub,B,N);
%max(abs(x2-x))



%5-2-10
clear;
load ex5-2-10
v = lb
N = [-2 -3 -4]
v(1) = b(2) - A(2,2:4)*v(2:4)
d = sign(A(1,:)*v - b(1))
A = [A; [-d; 0]]
lb(5) = 0; ub(5) = inf
B=[1 5]
w = [zeros(4,1); 1]
[x,B,N] = rsmbdd(A,b,w,lb,ub,B,N)

%5-1-11
clear;
load ex3-4-1
p = [p; zeros(5,1)]
x = simplex([A -eye(5)],b,p)
p'*x
%5-1-111
clear;
load ex3-6-5
%add slacking varibles
A = [A [-eye(2); zeros(1,2)]]
p = [p; zeros(2,1)]
lb = [0; 0; -inf; 0; 0]
ub = inf*ones(5,1)
x = simplex(A,b,p,lb,ub)
%5-2-12
clear;
load ex5-2-12;
x = simplex(A,b,p,lb,ub)

%9-1-2
A = [0 -1 2 0 ; 1 1 3 0; 3 -1 -3 0]
b = [-3; 2; 0]
c1 = [7 -1 -3]
c2= [-7 -1 -3]
p = [0 0 0 1 zeros(1,5)]'
H = [-c1 1; -c2 1; A]
h = [0; 0; b]
H=[H, -eye(5)]
lb = [-inf 0 0 -inf 0 0 0 0 0]'
ub = [inf inf inf inf inf inf inf inf inf]'
simplex(H,h,p,lb,ub)

%9-2-3
clear all;
load ex9-2-3;
[m,n] = size(A)
p = [zeros(n,1);1]
augA=[A ones(m,1); -A ones(m,1)]
T=totbl(augA,[b;-b],p)
T=relabel(T,'x3','eps')
T=ljx(T,5,1);
T=ljx(T,6,2);
T=permrows(T,[1 2 3 4 5 6 6]) 
T=ljx(T,4,3)

%9-2-5
clear;
load ex9-2-3
H= [A eye(3); -A eye(3)]
h= [b;-b] 
p= [zeros(2,1);ones(3,1)]
T=totbl(H,h,p)
T= relabel(T,'x3','y1','x4','y2','x5','y3')
T=ljx(T,5,1)
ljx(T,6,2)
T=permrows(T,[1 2 3 4 7 5 6]); 
T=ljx(T,4,4);
T=ljx(T,3,5);

%9-2-1
clear;
load ex9-2-1
[m,n] = size(A)
p = [zeros(n,1);1]
augA = [A ones(m,1);-A ones(m,1)]
T = totbl(augA,[b;-b],p)
T = relabel(T,'x3','eps')
T=ljx(T,1,1)
T=ljx(T,2,2)
T = permrows(T,[3:7 1 2])
T=ljx(T,4,3)

%6-2-4
load ex6-2-4 
B=[2 4]
N=[1 3]
nb=[1 1]'
Ab = A(:,B)
An= A(:,N)
Pb=p(B)
Pn=p(N)
Ab\b
Pn'- Pb'*(Ab\An)
z= Pb' * (Ab\nb)

%6-3-3
clear;
load ex6-3-3
T=totbl(A,b,p)
T=ljx(T,1,1)
T=ljx(T,1,2)



