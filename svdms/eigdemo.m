B = [6 5 -5
     2 6 -2
     2 5 -1];

V = [1 1 1; 0 1 1; 1 0 1]';
S = diag([6 4 1]);

%should be the same
B*V
V*S

x = [1 2 3]';

c = V\x

[x V*c]

for i = 1:15
  x = B*x;
  x = x/norm(x);
end

x
x'*B*x /(x'*x)

for i = 1:100
  x = B*x;
  x = x/norm(x);
end

x
x'*B*x /(x'*x)



X = rand(3,3);
for i = 1:15
  X = B*X;
  %X = X ./ repmat(sqrt(sum(X.^2,1)),3,1);
end

x = X(:,1);
x'*B*x /(x'*x)

x = X(:,2);
x'*B*x /(x'*x)

x = X(:,3);
x'*B*x /(x'*x)
