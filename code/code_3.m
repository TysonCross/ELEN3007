clc; clear all;  delete(get(0,'Children'));

%%% Display setting and output setup
scr = get(groot,'ScreenSize');                              % screen resolution
fig2 =  figure('Position',...                               % draw figure
        [scr(3)/9 scr(4)-(scr(3)/5) scr(3)/5 scr(3)/5]);
set(fig2,'numbertitle','off',...                            % Give figure useful title
        'name','Figure 2',...
        'Color','white',...
        'Units','normalized');
fontName='Helvetica';
fontsize=16;
set(0,'defaultAxesFontName', fontName);                     % Make fonts pretty
set(0,'defaultTextFontName', fontName);
set(0,'DefaultAxesFontSize', fontsize)
set(groot,'FixedWidthFontName', 'ElroNet Monospace')  


% FIGURE 2: Exponential Growth & Birth and Death
X0=2; b=2; d=1.5; % Initial value and parameters
x=[0:.1:8];
y=X0*exp((b-d) .*x);% Deterministic solution
plot(x,y,'k--', 'Linewidth' ,2);
axis( [0,8,0,100]);
title({'Exponential Growth & Birth and Death'},'Fontsize',14);
hold on

for k=1:5 % Five Sample Paths
    clear t n
    t(1)=0; n(1)=X0;
    j=1;
    while n(j)>0 & n(j)<100 % Stop when hits zero or reaches size=100
        u1=rand; u2=rand; % Two uniform random numbers
        t(j+1)=-log(u1)/(b*n(j)+d*n(j))+t(j);% Interevent time
        if u2<b/ (b+d)
            n(j+1)=n(j)+1; %Birth
        else
            n(j+1)=n(j)-1; %Death
        end
        j=j+1;
    end
    stairs(t,n,'-','linewidth',2,'color',rand(1,3));
end

hold off
xlabel('Time'); ylabel('Population size');
axis( [0,8,-0.5,100]);

%%% Display setting and output setup
scr = get(groot,'ScreenSize');                              % screen resolution
fig4 =  figure('Position',...                               % draw figure
        [scr(3)/3 scr(4)-(scr(4)/2.8) scr(3)/1.5 scr(4)/2.8]);
set(fig4,'numbertitle','off',...                            % Give figure useful title
        'name','Figure 4',...
        'Color','white',...
        'Units','normalized');
set(0,'defaultAxesFontName', fontName);                     % Make fonts pretty
set(0,'defaultTextFontName', fontName);
set(0,'DefaultAxesFontSize', fontsize)
set(groot,'FixedWidthFontName', 'ElroNet Monospace')  

% FIGURE 4: Example 1 Logistic Growth
b1=1.0; d1=0.5; r=b1-d1; K=250; X0=2; time=40; % Parameters and Initial values
x=[0:0.1:time];
y=X0*K./(X0+(K-X0)*exp(-r*x));
for k1=1:2
    subplot(1,2,k1)
    plot(x,y, 'k--','linewidth',2);
    axis([0,40/k1^2,0,300/k1^(log(6)/log(2))]);
    xlabel('Time'); ylabel('Population size');
    hold on
end
for k2=1:3 % 3 sample paths
    clear t n
    t(1)=0; n(1)=X0;
    j=1;
    while n(j)>0.0 & t(j)<time; % Stop when size hits zero or time reaches 40
        u1=rand;u2=rand; % Two uniform random numbers
        lam=b1*n(j);
        mu=d1*n(j)+r*n(j)^2/K;
        tot=lam+mu;
        t(j+1)=t(j)-log(u1)/(tot); % Interevent time
        if u2<lam/tot
            n(j+1) =n(j) +1;
        else
            n (j+1)=n(j)-1 ;
        end
        j=j+1;
    end
    for k3=1:2
        sp(k3)=subplot(1,2,k3)
        stairs(t,n, '-','linewidth',2,'color',rand(1,3));
    end
end
title(sp(1),{'Example 2 Logistic Growth'},'Fontsize',14);
title(sp(2),{'3 Sample Paths'},'Fontsize',14);


% Estimate Probability of Extinction in Logistic Growth
count=0;
for k4=1:10000 % Number of sample paths
    clear n
    n=X0;
    while n>0 & n<50; % Stop when size either hits zero or 50
        u=rand;
        lam=b1*n;
        mu=d1*n+r*n^2/K;
        tot=lam+mu;
        if u<lam/tot
            n=n+1; % Birth
        else
            n=n-1; % Death
        end
    end
    if n==0
        count=count+1;
    end
end
probext=count/10000
estext=(d1/b1)^(X0)

%%% Display setting and output setup
scr = get(groot,'ScreenSize');                              % screen resolution
fig5 =  figure('Position',...                               % draw figure
        [scr(3)/3 scr(4)/2-(scr(4)/2.8) scr(3)/1.5 scr(4)/2.8]);
set(fig5,'numbertitle','off',...                            % Give figure useful title
        'name','Figure 5',...
        'Color','white',...
        'Units','normalized');
set(0,'defaultAxesFontName', fontName);                     % Make fonts pretty
set(0,'defaultTextFontName', fontName);
set(0,'DefaultAxesFontSize', fontsize)
set(groot,'FixedWidthFontName', 'ElroNet Monospace')  

% FIGURE 5: Example 2 SIR Epidemic Model
gam=0.2; beta=0.4; N=800; i0=2 %Parameters and initial values
dt=0.05; tim=100; time=tim/dt; i(1)=i0; s(1)=N-i0; sumc=0;
for tt=1:time %Euler's method for solving ODE
    i(tt+1)=i(tt)+dt*((beta/N)*i(tt)*s(tt)-gam*i(tt));
    s(tt+1)=s(tt)+dt*(-(beta/N)*i(tt)*s(tt));
    sumc=sumc+dt*(beta/N)*i(tt)*s(tt);
end
TotalCases=round(sumc)+i0 % Count number of cases in ODE
for k1=1:2
    subplot(1,2,k1)
    plot([0:dt:tim],i,'k--','linewidth',2);
    axis([0,100/k1^(log(5)/log(2)),0,150/k1^(log(7)/log(2))]);
    xlabel('Time'); ylabel('Infectives');
    hold on
end

for k2=1:3 % Sample paths for SIR Markov chain
    clear t i s
    t(1)=0; i(1)=i0; s(1)=N-i0;
    j=1;
    while i(j)>0 & t(j)<time % Stop when infectives hit zero or time=lOO
        ul=rand;u2=rand;
        t(j+1)=-log(u1)/((beta/N)*i(j)*s(j)+gam*i(j))+t(j);
        if (u2<=(beta/N)*s(j)/(beta/N*s(j)+gam))
            i(j+1)=i(j)+1;
            s(j+1)=s(j)-1;
        else
            i(j+1)=i(j)-1;
            s(j+1)=s(j);
        end
        j=j+1;
    end
    for k3=1:2
        sp(k3)=subplot(1,2,k3)
        stairs(t,i,'-','linewidth',2,'color',rand(1,3));
    end
end
title(sp(1),{'Example 2 SIR Epidemic Model'},'Fontsize',14);
title(sp(2),{'Sample paths for SIR Markov chain'},'Fontsize',14);

% Estimate Probability of Epidemic Extinction in SIR model
count=0;
for k4=1:10000 % Number of sample paths
    clear t s i
    i=i0; s=N-i0;
    j=1;
    while i>0 & i<25 % Stop when size either hits zero or reaches 25
        ul=rand; u2=rand;
        if (u2<=(beta/N)*s/(beta/N*s+gam))
            i=i+1;
            s=s-1;
        else
            i=i-1;
            s=s;
        end
    j=j+1;
    end
    if i==0
        count=count+1;
    end
end
probext=count/10000
estext=(gam/beta)^i0