t(1)=0; % causal time
X(1)=2; % initial population
b=2;    % birth rate
d=1.5;  % death rate
u=rand; % uniform distribution

if u<b/(b+d)
    X(2)=X(1)+1; % birth
else
    X(2)=X(1)-1; %death
end

clc;
fprintf('X: %i -> %i \n',X) % X either increases to 3, or decreases to 1;
