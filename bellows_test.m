% init_density = 0:1:300;
% yVec1 = (init_density)./(1+(0.01*init_density).^5);
% yVec2 = (init_density)./(1+(0.02*init_density).^5);
% plot(init_density, yVec1);
% grid on
% hold on
% plot(init_density, yVec2,'r');
% xlabel('Initial density, $N_0$','interpreter','latex');
% ylabel('No. surviving @ $t = 1$','interpreter','latex');
% legend('II','III');
N0Vec = 1:150;
ePop = zeros(6,length(N0Vec));
tVec = 0:0.001:5;
ode_options = odeset('RelTol',1e-6);
% mu expression

% N0 = 100; % initial density
for j = 1:length(N0Vec)
    disp(j);
    IC = N0Vec(j);
    % [~, xout1] = ode45(@bellows_1,[0:0.01:1],IC,ode_options, 1);
    % [~, xout2] = ode45(@bellows_2,[0:0.01:1],IC,ode_options, 0.01);
    % [~, xout3] = ode45(@bellows_3,[0:0.01:1],IC,ode_options, 0.01);
    % [~, xout4] = ode45(@bellows_4,[0:0.01:1],IC,ode_options, 0.01, 1);
    [~, xout5a] = ode45(@bellows_5,[0:0.01:1],IC,ode_options, 0.01, 1);
    [~, xout5b] = ode45(@bellows_5,[0:0.01:1],IC,ode_options, 0.01, 5);
    [~, xout5c] = ode45(@bellows_5,[0:0.01:1],IC,ode_options, 0.02, 5);
    % [~, xout6] = ode45(@bellows_6,[0:0.01:1],IC,ode_options, 0.01, 5);
    % [~, xout7] = ode45(@bellows_7,[0:0.01:1],IC,ode_options, 0.01, 1);
    
    % ePop(1,j) = xout1(end); 
    % ePop(2,j) = xout2(end); 
    % ePop(3,j) = xout3(end); 
    % ePop(4,j) = xout4(end); 
    ePop(1,j) = xout5a(end); 
    ePop(2,j) = xout5b(end); 
    ePop(3,j) = xout5c(end); 
    ePop(4,j) = IC*exp(-log(1+(0.01*IC)^1)); % xout5d(end); 
    ePop(5,j) = IC*exp(-log(1+(0.01*IC)^5)); 
    ePop(6,j) = IC*exp(-log(1+(0.02*IC)^5)); 
    % ePop(6,j) = xout6(end); 
    % ePop(7,j) = xout7(end); 
end


subplot(1,2,2)
plot(N0Vec,ePop(1,:),'k--','Linewidth',1.5);
ylim([0,80])
title('(b)');
xlabel('initial density','Interpreter','latex'); 
ylabel('number surviving','Interpreter','latex'); 
set(gca,'Fontsize',[15]);
hold on
plot(N0Vec,ePop(2,:),'k-.','Linewidth',1.5);
plot(N0Vec,ePop(3,:),'k-','Linewidth',1.5);

subplot(1,2,1)
plot(N0Vec,ePop(4,:),'k--','Linewidth',1.5);
ylim([0,80])
title('(a)');
xlabel('initial density','Interpreter','latex'); 
ylabel('number surviving','Interpreter','latex'); 
set(gca,'Fontsize',[15]);
hold on
plot(N0Vec,ePop(5,:),'k-.','Linewidth',1.5);
plot(N0Vec,ePop(6,:),'k-','Linewidth',1.5);





function dx = bellows_1(t,x,b)
% define mu expression

dx = -x*b*log(x);
end

function dx = bellows_2(t,x,a)
% define mu expression

dx = -x*a*x;
end

function dx = bellows_3(t,x,a)
% define mu expression

dx = -x*log(1+(a*x));
end

function dx = bellows_4(t,x,a,b)
% define mu expression

dx = -x*a*(x^b);
end

function dx = bellows_5(t,x,a,b)
% define mu expression

dx = -x*log(1+(a*x)^b);
end

function dx = bellows_6(t,x,a,b)
% define mu expression

dx = -x*b*log(1+(a*x));
end

function dx = bellows_7(t,x,a,b)
% define mu expression

dx = -x*log(1+exp(b*x-a));
end

