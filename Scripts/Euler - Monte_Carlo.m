%% Name: Abhishek Joel George
%  Matriculation Number: 250900

%% Winnowing 4.6

clear all
close all 
clc


%% Input data

% Initial conditions
IC.Xp_0 = [0.5, 0.5];   % [m]   (Particle position)
IC.u_0 = 0.2;           % [m/s] (Air velocity)
IC.up_0 = [0, 0];       % [m/s] (Particle velocity)

%% Boundary conditions

h = 0.1;               
IC.Xp_c = [0.55, -0.5]; % Bin separation
Bin1x = [0 0 IC.Xp_c(1) IC.Xp_c(1)];
Bin1y = [IC.Xp_c(2) -1 -1 IC.Xp_c(2)];
Bin2x = [IC.Xp_c(1) IC.Xp_c(1) 1 1];
Bin2y = [IC.Xp_c(2) -1 -1 IC.Xp_c(2)];

% BC Plots
subplot(1, 2, 1)
yline(h / 2)
text(0.48, 0, 'Jet')
hold on
yline(-h / 2)
hold on
area(Bin1x, Bin1y, IC.Xp_c(2))
text(0.48, -0.58, 'Bin 1')
text(0.6, -0.58, 'Bin 2')
hold on
area(Bin2x, Bin2y, IC.Xp_c(2))
hold on

subplot(1, 2, 2)
yline(h / 2)
text(0.48, 0, 'Jet')
hold on
yline(-h / 2)
hold on
area(Bin1x, Bin1y, IC.Xp_c(2))
text(0.48, -0.58, 'Bin 1')
text(0.6, -0.58, 'Bin 2')
hold on
area(Bin2x, Bin2y, IC.Xp_c(2))
hold on

%% Air properties
rho.f = 1.2;            % Air density [kg/m^3]
mu_f = 1.8e-5;          % Dynamic viscosity [Pa·s]

%% Particle properties
rho.grain = 750;        % Grain density [kg/m^3]
d.grain = 2.5e-3;       % Grain diameter [m]

%% Chaff properties
rho.chaff = 50;         % Chaff density [kg/m^3]
d.chaff = 3.25e-3;      % Chaff diameter [m]

%% Time interval

tstart = 0;
tend = 2;


%% Define step sizes

M = 100;                % Number of random samples to be generated
N = [100000, 10000, 1000, 100]; % Number of timesteps

%% Forces
F.g = [0; -9.81];

%% Functions

f_Vp = @(diameter) pi * (diameter)^3 / 6; % Volume
f_Fb = @(density) (-1) * rho.f / density * F.g; % Buoyancy
f_Re = @(uf, up, diameter) (rho.f * norm(uf - up) * diameter) / mu_f; % Reynolds number
f_Fd = @(diameter, density, volume, Cd_coeff, uf, up) (pi * (diameter^2) * rho.f * Cd_coeff) / (2 * density * volume) *norm(uf - up) * (uf - up); % Drag

%% ODE solution functions
f_d_up = @(Fb, Fd) F.g + Fb + Fd;

%% Errors - Solution arrays

%   Euler - error
error.euler.grain_local=zeros(1,length(N));
error.euler.grain_global=zeros(1,length(N));
error.euler.chaff_local=zeros(1,length(N));
error.euler.chaff_global=zeros(1,length(N));

% Monte Carlo
error.MC.grain_local = zeros(1, length(N));
error.MC.grain_global = zeros(1, length(N));
error.MC.chaff_local = zeros(1, length(N));
error.MC.chaff_global = zeros(1, length(N));

%% Solution for varying timesteps

for j = 1:length(N)
    dt(j) = (tend - tstart) / N(j);
    t = tstart:dt(j):tend;

    %% Initialization

    Xp = zeros(2, N(j));
    up = zeros(2, N(j));
    uf = zeros(2, N(j));
    
    % Initial conditions
    Xp(1, 1) = IC.Xp_0(1);
    Xp(2, 1) = IC.Xp_0(2);
    up(1, 1) = IC.up_0(1);
    up(1, 1) = IC.up_0(2);
    uf(1, 1) = IC.u_0;
    MC_uf = [IC.u_0;0];
    
    %% Euler Solution 
    
    for m = 1:2 % Select particle properties
        if m == 1
            density = rho.grain;
            dia = d.grain;
        elseif m == 2
            density = rho.chaff;
            dia = d.chaff;
        end

        Vp = f_Vp(dia);
        Fb = f_Fb(density);
    
        for i = 1:N(j)-1
      

            Re = f_Re(uf(:, i), up(:, i), dia);
            Cd = Calc_Cd(Re);
            Fd = f_Fd(dia, density, Vp, Cd, uf(:, i), up(:, i));
    
            Xp(:, i + 1) = Xp(:, i) + dt(j) * up(:, i); % Euler Solution
            up(:, i + 1) = up(:, i) + dt(j) * f_d_up(Fb, Fd);
            uf(1, i + 1) = Calc_uf(Xp(1, i + 1), Xp(2, i + 1), h, IC.u_0);
        end 

        % Euler Error calculation

        if m ==1
            if N(j) == max(N)
                Xp_best.grain_euler = Xp;
                error.euler.grain_local(1,j) = 0;
                error.euler.grain_global(1,j) = 0;
            else
                w=0;
                i=1;
                while w==0 || i>(length(t)-1)
                    diff=abs(Xp(1,i)-IC.Xp_0(1));
                    if diff>=0.001
                        w=i;
                    else
                        i=i+1;
                    end
                end
                o= round(t(w)/dt(1),0);
                error.euler.grain_local(1,j) = abs(Xp_best.grain_euler(1,o)-Xp(1,w))/Xp_best.grain_euler(1,o);
                
                error.euler.grain_global(1,j) = abs(Xp_best.grain_euler(1,end)-Xp(1,end))/Xp_best.grain_euler(1,end);
                
            end

        elseif m ==2
            if N(j) == max(N)
                Xp_best.chaff_euler = Xp;
                error.euler.chaff_local(1,j) = 0;
                error.euler.chaff_global(1,j) = 0;
            else
                w=0;
                i=1;
                while w==0 || i>(length(t)-1)
                    diff=abs(Xp(1,i)-IC.Xp_0(1));
                    if diff>=0.001
                        w=i;
                    else
                        i=i+1;
                    end
                end
                o= round(t(w)/dt(1),0);
                error.euler.chaff_local(1,j) = abs(Xp_best.chaff_euler(1,o)-Xp(1,w))/Xp_best.chaff_euler(1,o);
                
                error.euler.chaff_global(1,j) = abs(Xp_best.chaff_euler(1,end)-Xp(1,end))/Xp_best.chaff_euler(1,end);
                
            end
        end
        %% Plot Euler trajaectories

        subplot(1, 2, 1)
        hold on
        plot(Xp(1, :), Xp(2, :));
        xlim([0.45 0.65]);
        ylim([-0.55 0.55]);
        title("Euler");
    end

  


    %% Monte Carlo Solution
    
    for m = 1:2
        if m == 1
            dia = d.grain;
            density = rho.grain;
        else
            dia = d.chaff;
            density = rho.chaff;
        end
    
        for i = 1:N(j)-1
            MC_uf = [IC.u_0;0];
            MC_Fd = [0; 0];
            Vp = f_Vp(dia);
            Fb = f_Fb(density);
            Re = f_Re(uf(:, i), up(:, i), dia);
            Cd = Calc_Cd(Re);
            Fd = f_Fd(dia, density, Vp, Cd, uf(:, i), up(:, i));
    
            for k = 1:M
                t_n = t(i) + dt(j) * rand();
                MC_dt = t_n - t(i);
                MC_up = up(:, i) + MC_dt .* f_d_up(Fb, Fd);
                MC_Xp = Xp(:, i) + MC_dt .* up(:, i);
                MC_uf = Calc_uf(MC_Xp(1), MC_Xp(2), h, IC.u_0);
    
                if MC_uf > MC_up 
                    MC_Re = f_Re(MC_uf, MC_up, dia);
                    MC_Cd = Calc_Cd(MC_Re);
                    MC_Fd = MC_Fd + f_Fd(dia, density, Vp, MC_Cd, MC_uf, MC_up);
                else
                    k = k - 1;
                end
            end
    
            up(:, i + 1) = up(:, i) + dt(j) * f_d_up(Fb, MC_Fd ./ M);
            Xp(:, i + 1) = Xp(:, i) + dt(j) * up(:, i);
            uf(:, i + 1) = Calc_uf(Xp(1, i + 1), Xp(2, i + 1), h, IC.u_0);
        end

        %% Error Calculation - Monte Carlo
        
        if m ==1
            if N(j) == max(N)
                Xp_best.grain_MC = Xp;
                error.MC.grain_local(1,j) = 0;
                error.MC.grain_global(1,j) = 0;
            else

                w=0;
                i=1;
                while w==0 || i>(length(t)-1)
                    diff=abs(Xp(1,i)-IC.Xp_0(1));
                    if diff>=0.001
                        w=i;
                    else
                        i=i+1;
                    end
                end
                o= round(t(w)/dt(1),0);
                    
                
                error.MC.grain_local(1,j) = abs(Xp_best.grain_MC(1,o)-Xp(1,w))/Xp_best.grain_MC(1,o);
                
                error.MC.grain_global(1,j) = abs(Xp_best.grain_MC(1,end)-Xp(1,end))/Xp_best.grain_MC(1,end);
                
            end

        elseif m ==2
            if N(j) == max(N)
                Xp_best.chaff_MC = Xp;
                error.MC.chaff_local(1,j) = 0;
                error.MC.chaff_global(1,j) = 0;
            else

                w=0;
                i=1;
                while w==0 || i>(length(t)-1)
                    diff=abs(Xp(1,i)-IC.Xp_0(1));
                    if diff>=0.001
                        w=i;
                    else
                        i=i+1;
                    end
                end
                o= round(t(w)/dt(1),0);
                error.MC.chaff_local(1,j) = abs(Xp_best.chaff_MC(1,o)-Xp(1,w))/Xp_best.chaff_MC(1,o);
                
                error.MC.chaff_global(1,j) = abs(Xp_best.chaff_MC(1,end)-Xp(1,end))/Xp_best.chaff_MC(1,end);
                
            end
        end
        
        % Trajectories - Monte Carlo

        subplot(1, 2, 2)
        plot(Xp(1, :), Xp(2, :))
        xlim([0.45 0.8]);
        ylim([-0.55 0.55]);
        title("Monte Carlo")
    end
end

%% Format trajectories plots

subplot(1,2,1)
title('Euler Solution')
legend({'','','','','Grain - N=1e5','Chaff - N=1e5','Grain - N=1e4','Chaff - N=1e4','Grain - N=1e3','Chaff - N=1e3','Grain - N=1e2','Chaff - N=1e2'},'NumColumns',2);
subplot(1,2,2)
title('Monte Carlo Solution')
legend({'','','','','Grain - N=1e5','Chaff - N=1e5','Grain - N=1e4','Chaff - N=1e4','Grain - N=1e3','Chaff - N=1e3','Grain - N=1e2','Chaff - N=1e2'},'NumColumns',2);

%% Error vs. Time Step Plot
figure

% Grain x-component global error plot
subplot(2, 2, 1)
loglog(N, error.euler.grain_global(1, :));
hold on

loglog(N, error.MC.grain_global(1, :));


xlabel('Number of Timesteps (N)')
ylabel('Error')
title('Grain x-component Global Error')
legend('Euler','Monte-Carlo')
grid on

% Chaff x-component global error plot
subplot(2, 2, 2)
loglog(N, error.euler.chaff_global(1, :));
hold on

loglog(N, error.MC.chaff_global(1, :));


xlabel('Number of Timesteps (N)')
ylabel('Error')
title('Chaff x-component Global Error')
legend('Euler','Monte-Carlo')
grid on

% Grain x-component local error plot
subplot(2, 2, 3)
loglog(N, error.euler.grain_local(1, :));
hold on

loglog(N, error.MC.grain_local(1, :));


xlabel('Number of Timesteps (N)')
ylabel('Error')
title('Grain x-component Local Error')
legend('Euler','Monte-Carlo')
grid on

% Chaff x-component local error plot
subplot(2, 2, 4)
loglog(N, error.euler.chaff_local(1, :));
hold on

loglog(N, error.MC.chaff_local(1, :));
xlabel('Number of Timesteps (N)')
ylabel('Error')
title('Chaff x-component local Error')
legend('Euler','Monte-Carlo')
grid on


%% Terminal velocity 

f_ut= @(Cd,density,volume,dia) sqrt((2*(density-rho.f)*volume*abs(F.g(2)))/(rho.f*Cd*(pi*(dia)^2)));
dt = 0.01;
t = tstart:dt:tend;
thres = 0.01;

%% Solution arrays

Xp = zeros(2,length(t));
up = zeros(2,length(t));
uf = zeros(2,length(t));


%% Initial conditions

Xp(:,1) = IC.Xp_0;
uf(1,1) = IC.u_0;
up(:,1) = IC.up_0;
IC.angle = [-95,-85];

%% Solution - Euler

for m = 1:2
    if m==1 
        dia = d.grain;
        density = rho.grain;
    else
        dia = d.chaff;
        density = rho.chaff;
    end

    ut = zeros(1,length(t));

    Vp = f_Vp(dia);
    Fb = f_Fb(density);

    y=1;
    i=1;
    
    while y>thres
        Re = f_Re(uf(:,i),up(:,i),dia);
        Cd = Calc_Cd(Re);
        Fd = f_Fd(dia,density,Vp,Cd,uf(:,i),up(:,i));

        up(:,i+1) = up(:,i) + dt.*f_d_up(Fb,Fd);
        Xp(:,i+1) = Xp(:,i) + dt.*up(:,i);
        uf(:,i+1) = Calc_uf(Xp(1,i+1),Xp(1,i+1),h,IC.u_0);
        ut(i) = f_ut(Cd,density,Vp,dia);
        %i = i+1;

        if i==1
            y =1;
        else
            y = abs(ut(i)-ut(i-1));
        end
        i = i+1;
    end
    if diff<thres && m==1
        u_t.max_grain = max(ut);
    elseif diff<thres && m==2
        u_t.max_chaff = max(ut);
    end
end



%% Simulate trajectories of 1000 heavy and light particles

R = 1000;
dt = 0.001;
tend =3;
t = tstart:dt:tend;

%% Boundary plots

figure
% BC Plots
subplot(1,2,1)
yline(h/2)
text(0.35,0,'Jet')
hold on
yline(-h/2)
hold on
area(Bin1x,Bin1y,IC.Xp_c(2))
text(0.48,-0.58,'Bin 1')
hold on
area(Bin2x, Bin2y, IC.Xp_c(2))
text(0.6,-0.58,'Bin 2')

subplot(1,2,2)
yline(h/2)
text(0.35,0,'Jet')
hold on
yline(-h/2)
hold on
area(Bin1x,Bin1y,IC.Xp_c(2))
text(0.48,-0.58,'Bin 1')
hold on
area(Bin2x, Bin2y, IC.Xp_c(2))
text(0.6,-0.58,'Bin 2')

IC.Xpc = 0.45:0.02:0.65;
Bin1 = zeros(1,length(IC.Xpc)+1);
Bin2 = zeros(1,length(IC.Xpc)+1);



%% Distribution parameters

d.grain_s = 1e-3;


rho.chaff_s = 20;
d.chaff_min = 2e-3;
d.chaff_max = 5e-3;





for m =1:2
    Xp = zeros(2,length(t));
    up = zeros(2,length(t));
    uf = zeros(2,length(t));
    ut = zeros(1,length(t));
    
    
    % Initial conditions
    
    Xp(1,1) = IC.Xp_0(1);
    Xp(2,1) = IC.Xp_0(2);
    uf(1,1) = IC.u_0;

    % Particle property selection
    c=1;
    if m==1 % Grain
        rho.samples=rho.grain*ones(1,R); % (Constant grain density)
        while c<(R+1)
            d.samples(c)=d.grain+d.grain_s*randn(); % (Normally distributed grain diameter)
            if d.samples(c)>(d.grain-d.grain_s*2) && d.samples(c)<(d.grain+d.grain_s*2)
                c=c+1;
            end
        end
        rho.samples_g=rho.samples;
        d.samples_g=d.samples;
    elseif m==2 % Chaff
        while c<(R+1)
            rho.samples(c)=rho.chaff+rho.chaff_s*randn(); % (Normally distributed chaff density)
            if rho.samples(c)>(rho.chaff-rho.chaff_s*2) && rho.samples(c)<(rho.chaff+rho.chaff_s*2)
                c=c+1;
            end
        end
        d.samples=d.chaff_min+(abs(d.chaff_max-d.chaff_min))*rand(1,R); % (Uniformlly distributed chaff diameter)
        rho.samples_c=rho.samples;
        d.samples_c=d.samples;
    end
    Vp=f_Vp(dia);
    Fb=f_Fb(density);
    
    w=1;

    for r = 1:R
        if m==1
            up(1,1) = u_t.max_grain*cos((IC.angle(1)+(w-1))*pi/180);
            up(2,1) = u_t.max_chaff*sin((IC.angle(1)+(w-1))*pi/180);
        elseif m==2
            up(1,1) = u_t.max_grain*cos((IC.angle(1)+(w-1))*pi/180);
            up(2,1) = u_t.max_chaff*sin((IC.angle(1)+(w-1))*pi/180);
        end

        density = rho.samples(r);
        dia = d.samples(r);

        %Euler solution

        for i = 1:length(t)-1
            Re = f_Re(uf(:,i),up(:,i),dia);
            Cd = Calc_Cd(Re);
            Fd = f_Fd(dia,density,Vp,Cd,uf(:,i),up(:,i));

            up(:,i+1) = up(:,i) + dt.*f_d_up(Fb,Fd);
            Xp(:,i+1) = Xp(:,i) + dt.*up(:,i);
            uf(:,i+1) = Calc_uf(Xp(1,i+1),Xp(2,i+1),h,IC.u_0);
            if Xp(2,i+1)>1
                error=1;
            end

            

        end
        if w==10
            w=1;
        else
            w=w+1;
        end

        if m==1
           subplot(1,2,1)
           hold on
           plot(Xp(1,:),Xp(2,:));
           y=0;
           i=1;
           while y==0 || i>length(t)-1
               if Xp(2,i)<=IC.Xp_c(2)
                  y = i-1;
               else
                   i=i+1;
               end
           end
           if Xp(1,y)>=IC.Xp_c(1)
               Bin2(1) = Bin2(1)+1;
           end

           for k=1:length(IC.Xpc)
               if Xp(1,y)>=IC.Xpc(k)
                   Bin2(k+1) = Bin2(k+1)+1;
               end
           end
        elseif m==2
                subplot(1,2,2)
                hold on
                plot(Xp(1,:),Xp(2,:));
               
               if Xp(1,end)<=IC.Xp_c(1)
                   Bin1(1) = Bin1(1)+1;
               end
    
               for k=1:length(IC.Xpc)
                   if Xp(1,end)<=IC.Xpc(k)
                       Bin1(k+1) = Bin1(k+1)+1;
                   end
               end
        end
    end
end


%% Proportions 

Bin_ratio(1,:) = Bin2.*(100/R);
Bin_ratio(2,:) = Bin1.*(100/R);
Bin_opt = sum(Bin_ratio);
k = find(Bin_opt == min(min(Bin_opt)));
Xc_opt = IC.Xpc(k-1);

subplot(1,2,1)
hold on
xline(Xc_opt,'LineStyle','- . ')
Xc_opt1=['Optimal x_c =' num2str(Xc_opt)];
text(Xc_opt-0.05,0.50,Xc_opt1)
Xopt=[num2str(Bin_ratio(1,k)) ' % of grain fell in Bin 2'];
text(0.6,0.34,'Optimal Bin placement')
text(0.6,0.3,Xopt)
txt1=[num2str(Bin_ratio(1,1)) ' % of grain fell in Bin 2'];
text(0.6,0.4,txt1)
xlim([0.3 0.8])
ylim([-0.75 0.55])
title('Grain')

subplot(1,2,2)
hold on
xline(Xc_opt,'LineStyle','- . ')
Xc_opt1=['Optimal x_c =' num2str(Xc_opt)];
text(Xc_opt-0.05,0.50,Xc_opt1)
Xopt=[num2str(Bin_ratio(2,k)) ' % of chaff fell in Bin 1'];
text(0.6,0.34,'Optimal Bin placement')
text(0.6,0.3,Xopt)
txt2=[num2str(Bin_ratio(2,1)) ' % of chaff fell in Bin 1'];
text(0.6,0.4,txt2)
xlim([0.3 0.8])
ylim([-0.75 0.55])
title('Chaff')

figure
subplot(2,2,1)
histogram(rho.samples_g,'Normalization','pdf')
title('Grain density distribution [Constant]')
subplot(2,2,3)
histogram(d.samples_g,'Normalization','pdf')
title('Grain diameter normal distribution')
subplot(2,2,2)
histogram(rho.samples_c,'Normalization','pdf')
title('Chaff density normal distribution')
subplot(2,2,4)
histogram(d.samples_c,'Normalization','pdf')
title('Chaff diameter uniform distribution')



                   



        

         







        









       


%% Functions
function Cd = Calc_Cd(Re)
    if Re < 800
        Cd = (24 / Re) * (1 + 0.15 * Re^(0.687));
    else
        Cd = 0.44;
    end
end

function uf = Calc_uf(x, y, h, u_0)
    if x >= 5 * h
        uf = 6.2 * u_0 * sqrt(h / x) * exp(-50 * (y^2 / x^2));
    else
        uf = 0;
    end
end





         
