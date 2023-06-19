fr = 10.^linspace(-1,2,100);
u = @(t,fr) sin(2*pi*t*fr);
wn = 10*2;

aRange = linspace(0,1,10);
figureNB;
for j = 1:length(aRange)
    a = aRange(j);
    for i = 1:length(fr)
        dx = @(t,x) [x(2); ...
                    wn^2*u(t,fr(i)) - 2*a*wn*x(2) - wn^2*x(1)];
        [T,Y] = ode15s(dx,[0,10/fr(i)],[0;0]);
        % plot(T,Y(:,1));
        % hold on;
        % drawnow;
        A(i) = range(Y(:,1));
    end
    plot(fr,A);
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    hold on;
    drawnow;
end



a = 0;
dt = 1e-5;
T = 0:dt:10;
wn = 10;

dx = @(t,x) [x(2); ...
            wn^2*u(t,10) - 2*a*wn*x(2) - wn^2*x(1)];
dx = @(t,x) [x(2); ...
            wn^2*randn - 2*a*wn*x(2) - wn^2*x(1)];
Y = zeros(2,length(T));
N = zeros(length(T),1);
rho = 0.1;

for i = 1:length(T)-1
    Y(1,i+1) = Y(1,i) + dt*Y(2,i);
    Y(2,i+1) = Y(2,i) + dt*(wn^2*N(i) - 2*a*wn*Y(2,i) - wn^2*Y(1,i));
    N(i+1) = N(i) - rho*dt*N(i) + sqrt(2*dt*rho)*randn;
    % N(i+1) = sin(2*pi*T(i)*10);
    % N(i+1) = randn;
end

% figureNB;
% plot(T,Y(1,:));
pspectrum(Y(1,:),1/dt,'FrequencyResolution',1,'FrequencyLimits',[0,100])
set(gca,'xscale','log')