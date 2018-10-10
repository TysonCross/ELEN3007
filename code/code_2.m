t(1)=0; % causal time
X(1)=2; % initial population
b=2;    % birth rate
d=1.5;  % death rate
u=rand; % uniform distribution

t(2)=t(1)-log(u)/(b+d)*X(1);

clc;
fprintf('t: %i -> %i \n',t) % interevent T_i change from random number u
