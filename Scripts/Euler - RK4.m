%% Name: Abhishek Joel George
%  Matriculation Number: 250900


%% Winnowing - 3.5

clear all

%% Input data

% Initial conditions
IC.Xp_0=[0.5,0.5];    % [m]   (Particle position)
IC.u_0=0.2;          % [m/s] (Air velocity)
IC.up_0=[0,0];       % [m/s] (Particle velocity)

%% Boundary conditions

h=0.1;               
IC.Xp_c=[0.55,-0.5]; % Bin separation
Bin1x=[0 0 IC.Xp_c(1) IC.Xp_c(1)];
Bin1y=[IC.Xp_c(2) -1 -1 IC.Xp_c(2)];
Bin2x=[IC.Xp_c(1) IC.Xp_c(1) 1 1];
Bin2y=[IC.Xp_c(2) -1 -1 IC.Xp_c(2)];

% BC Plots
subplot(1,2,1)
yline(h/2)
text(0.48,0,'Jet')
hold on
yline(-h/2)
hold on
area(Bin1x,Bin1y,IC.Xp_c(2))
text(0.48,-0.58,'Bin 1')
text(0.6,-0.58,'Bin 2')
hold on
area(Bin2x,Bin2y,IC.Xp_c(2))
text(0.48,-0.58,'Bin 1')
text(0.6,-0.58,'Bin 2')
subplot(1,2,2)
yline(h/2)
text(0.48,0,'Jet')
hold on
yline(-h/2)
hold on
area(Bin1x,Bin1y,IC.Xp_c(2))
text(0.48,-0.58,'Bin 1')
text(0.6,-0.58,'Bin 2')
hold on
area(Bin2x,Bin2y,IC.Xp_c(2))
text(0.48,-0.58,'Bin 1')
text(0.6,-0.58,'Bin 2')
hold on

%% Air properties
rho.f = 1.2;
mu_f = 1.8e-5;

%% Particle properties
rho.grain = 750;
d.grain = 2.5e-3;

%% Chaff properties
rho.chaff = 50;
d.chaff = 3.25e-3;

%% Time interval

tstart =0;
tend =2;


%% Initialization 
n = 5;
dt =zeros(1,n);


%   Euler - error
error.euler.grain_local=zeros(2,n);
error.euler.grain_global=zeros(2,n);
error.euler.chaff_local=zeros(2,n);
error.euler.chaff_global=zeros(2,n);

%   Runge-kutta
k.u_p=zeros(2,4);
k.X_p=zeros(2,4);
error.rkutta.grain_local=zeros(2,n);
error.rkutta.grain_global=zeros(2,n);
error.rkutta.chaff_local=zeros(2,n);
error.rkutta.chaff_global=zeros(2,n);

%% Forces
F.g = [0;-9.81];

%% Functions

f_Vp= @(diameter) pi*(diameter)^3/6; % Volume
f_Fb= @(density) (-1)*rho.f/density*F.g; %Buoyancy
f_Re= @(uf,up,diameter) (rho.f*norm(uf-up)*diameter)/mu_f; %Reynolds number
f_Fd= @(diameter,density,volume,Cd_coeff,uf,up) (pi*(diameter^2)*rho.f*Cd_coeff)/(2*density*volume)*norm(uf-up)*(uf-up);%Drag

%% ODE solution functions
f_d_up= @(Fb,Fd) F.g+Fb+Fd;


%% Solutions 

for j = 1:n
    dt(j) = 10^(-(n+1)+j);
end

[dt_best, idx_best] = min(dt);

for j = 1:n
    t = tstart:dt(j):tend;
    
        % Initialization
        Xp = zeros(2,length(t));
        up = zeros(2,length(t));
        uf = zeros(2,length(t));
    
        % Initial conditions
        Xp(:,1) = IC.Xp_0;
        up(:,1) = IC.up_0;
        uf(1,1) = IC.u_0;
        
        %% Choose chaff or particles 
        
        for m = 1:2
            if m ==1
                dia = d.grain;
                density = rho.grain;
            elseif m==2
                dia = d.chaff;
                density = rho.chaff;
            end
            Vp = f_Vp(dia);
            Fb = f_Fb(density);
            %% Euler solution
            for i = 1:length(t)-1
                Re = f_Re(uf(:,i),up(:,i),dia);
                Cd = Calc_Cd(Re);
                Fd = f_Fd(dia,density,Vp,Cd,uf(:,i),up(:,i));
                up(:,i+1) = up(:,i) + dt(j).*f_d_up(Fb,Fd); % velocity
                Xp(:,i+1) = Xp(:,i) + dt(j).*up(:,i);
                uf(1,i+1) = Calc_uf(Xp(1,i+1),Xp(2,i+1),h,IC.u_0);
            end
            
            %% Euler error calculations
            if m == 1 % grain
                if j == idx_best
                    Xp_grain_best = Xp;
                    error.euler.grain_local(:,j)=[0 0];
                    error.euler.grain_global(:,j)= [0 0];
                    
                else
                    error.euler.grain_local(1,j) = abs((Xp(1,3)-Xp_grain_best(1,3))/Xp_grain_best(1,3)); %x
                    error.euler.grain_local(2,j) = abs((Xp(2,3)-Xp_grain_best(2,3))/Xp_grain_best(2,3));%y
                    error.euler.grain_global(1,j)=abs((Xp(1,end)-Xp_grain_best(1,end))/Xp_grain_best(1,end)); %x
                    error.euler.grain_global(2,j)=abs((Xp(2,end)-Xp_grain_best(2,end))/Xp_grain_best(2,end)); %y
                end
            elseif m == 2 %chaff
                if j == idx_best
                    Xp_chaff_best=Xp;
                    error.euler.chaff_local(:,j)=[0,0];
                    error.euler.chaff_global(:,j)=[0,0];
                else
                    error.euler.chaff_local(1,j)=abs((Xp(1,3)-Xp_chaff_best(1,3))/Xp_chaff_best(1,3)); % Error in x
                    error.euler.chaff_local(2,j)=abs((Xp(2,3)-Xp_chaff_best(2,3))/Xp_chaff_best(2,3)); % Error in y
                    error.euler.chaff_global(1,j)=abs((Xp(1,end)-Xp_chaff_best(1,end))/Xp_chaff_best(1,end)); % Error in x
                    error.euler.chaff_global(2,j)=abs((Xp(2,end)-Xp_chaff_best(2,end))/Xp_chaff_best(2,end)); % Error in y
                end
            end
    
            %% Plots
            subplot(1,2,1)
            hold on 
            if dt(j) < 0.1
                plot(Xp(1,:), Xp(2,:));
            end
    
          
            xlim([0.45 0.65])
            ylim([-0.5 0.55])
    
           
            
            %% Rk4 solution
            for i = 1:length(t)-1
                %       k1
                Re=f_Re(uf(:,i),up(:,i),dia);
                Cd=Calc_Cd(Re);
                Fd=f_Fd(dia,density,Vp,Cd,uf(:,i),up(:,i));
                k.up(:,1)=f_d_up(Fb,Fd); % Velocity derivative at f(tn,yn)
                k.Xp(:,1)=up(:,i); % Position derivative at f(tn,yn)
                %       k2
                k.Xp(:,2)=up(:,i)+k.up(:,1)*dt(j)/2; % Position derivative at f(tn+dt/2,yn+dt/2*k1)
                Re=f_Re(uf(:,i),k.Xp(:,2),dia);
                Cd=Calc_Cd(Re);
                Fd=f_Fd(dia,density,Vp,Cd,uf(:,i),k.Xp(:,2));
                k.up(:,2)=f_d_up(Fb,Fd); % Velocity derivative at f(tn+dt/2,yn+dt/2*k1)
                %       k3
                k.Xp(:,3)=up(:,i)+k.up(:,2)*dt(j)/2; % Position derivative at f(tn+dt/2,yn+dt/2*k2)
                Re=f_Re(uf(:,i),k.Xp(:,3),dia);
                Cd=Calc_Cd(Re);
                Fd=f_Fd(dia,density,Vp,Cd,uf(:,i),k.Xp(:,3));
                k.up(:,3)=f_d_up(Fb,Fd); % Velocity derivative at f(tn+dt/2,yn+dt/2*k2)
                %       k4
                k.Xp(:,4)=up(:,i)+k.up(:,3)*dt(j); % Position derivative at f(tn+dt,yn+dt*k3)
                Re=f_Re(uf(:,i),k.Xp(:,4),dia);
                Cd=Calc_Cd(Re);
                Fd=f_Fd(dia,density,Vp,Cd,uf(:,i),k.Xp(:,4));
                k.up(:,4)=f_d_up(Fb,Fd); % Velocity derivative at f(tn+dt,yn+dt*k3)
                %       Update
                up(:,i+1)=up(:,i)+dt(j)/6.*(k.up(:,1)+2*k.up(:,2)+2*k.up(:,3)+k.up(:,4));
                Xp(:,i+1)=Xp(:,i)+dt(j)/6.*(k.Xp(:,1)+2*k.Xp(:,2)+2*k.Xp(:,3)+k.Xp(:,4));
                uf(1,i+1)=Calc_uf(Xp(1,i+1),Xp(2,i+1),h,IC.u_0);
            end
            %% RK4 Error Calculations
            if m == 1 % grain
                if j == idx_best
                    Xp_grain_best_rk = Xp;
                    error.rkutta.grain_local(:,j) = [0, 0];
                    error.rkutta.grain_global(:,j) = [0, 0];
                else
                    error.rkutta.grain_local(1,j) = abs((Xp(1,3) - Xp_grain_best_rk(1,3)) / Xp_grain_best_rk(1,3)); % x
                    error.rkutta.grain_local(2,j) = abs((Xp(2,3) - Xp_grain_best_rk(2,3)) / Xp_grain_best_rk(2,3)); % y
                    error.rkutta.grain_global(1,j) = abs((Xp(1,end) - Xp_grain_best_rk(1,end)) / Xp_grain_best_rk(1,end)); % x
                    error.rkutta.grain_global(2,j) = abs((Xp(2,end) - Xp_grain_best_rk(2,end)) / Xp_grain_best_rk(2,end)); % y
                end
            elseif m == 2 % chaff
                if j == idx_best
                    Xp_chaff_best_rk = Xp;
                    error.rkutta.chaff_local(:,j) = [0, 0];
                    error.rkutta.chaff_global(:,j) = [0, 0];
                else
                    error.rkutta.chaff_local(1,j) = abs((Xp(1,3) - Xp_chaff_best_rk(1,3)) / Xp_chaff_best_rk(1,3)); % x
                    error.rkutta.chaff_local(2,j) = abs((Xp(2,3) - Xp_chaff_best_rk(2,3)) / Xp_chaff_best_rk(2,3)); % y
                    error.rkutta.chaff_global(1,j) = abs((Xp(1,end) - Xp_chaff_best_rk(1,end)) / Xp_chaff_best_rk(1,end)); % x
                    error.rkutta.chaff_global(2,j) = abs((Xp(2,end) - Xp_chaff_best_rk(2,end)) / Xp_chaff_best_rk(2,end)); % y
                end
            end
            
            %% Plots for RK4
            subplot(1,2,2)
            hold on
            if dt(j) < 0.1
                plot(Xp(1,:), Xp(2,:));
            end
            
            xlim([0.45 0.65])
            ylim([-0.5 0.55])
        end
    end

    %%  Format plots

% Solutionplots
subplot(1,2,1)
title('Euler Solution')
legend({'','','','','Grain - dt=1e-5','Chaff - dt=1e-5','Grain - dt=1e-4','Chaff - dt=1e-4','Grain - dt=1e-3','Chaff - dt=1e-3','Grain - dt=1e-2','Chaff - dt=1e-2'},'NumColumns',2);
subplot(1,2,2)
title('RK4 Solution')
legend({'','','','','Grain - dt=1e-5','Chaff - dt=1e-5','Grain - dt=1e-4','Chaff - dt=1e-4','Grain - dt=1e-3','Chaff - dt=1e-3','Grain - dt=1e-2','Chaff - dt=1e-2'},'NumColumns',2);

%% Error vs. Time Step Plot
figure
subplot(2,2,1) % Grain x-component error plot
loglog(dt, error.euler.grain_global(1,:), '-o');
hold on
loglog(dt,error.euler.grain_local(1,:),'-x');
hold on
loglog(dt, error.rkutta.grain_global(1,:), '-s');
hold on
loglog(dt,error.rkutta.grain_local(1,:),'-d');
hold on
loglog(dt,dt,dt,dt.^2,dt,dt.^3,dt,dt.^5);
legend('Euler Grain (Global Error - x)','Euler Grain (Local Error - x)','RK4 Grain (Global Error - x)','RK4 Grain (Local Error - x)','dt','dt^2','dt^3','dt^4')
hold on
xlabel('dt')
ylabel(' Error')
legend()
title(' grain x')

subplot(2,2,2)
loglog(dt, error.euler.chaff_global(1,:), '-^');
hold on
loglog(dt,error.euler.chaff_local(1,:),'-s');
hold on
loglog(dt, error.rkutta.chaff_global(1,:), '-d');
hold on
loglog(dt,error.rkutta.chaff_local(1,:),'-x');
hold on
loglog(dt,dt,dt,dt.^2,dt,dt.^3,dt,dt.^5);
legend('Euler Chaff (Global Error - x)','Euler Chaff (Local Error - x)','RK4 Chaff (Global Error - x)','RK4 Chaff (Local Error - x)','dt','dt^2','dt^3','dt^4')
xlabel('dt')
ylabel(' Error')
legend
title(' chaff x')

subplot(2,2,3)
loglog(dt, error.euler.grain_global(2,:), '-o');
hold on
loglog(dt,error.euler.grain_local(2,:),'-d');
hold on
loglog(dt, error.rkutta.grain_global(2,:), '-s');
hold on
loglog(dt,error.rkutta.grain_local(2,:),'-x');
hold on
loglog(dt,dt,dt,dt.^2,dt,dt.^3,dt,dt.^4);
legend('Euler Grain (Global Error - y)','Euler Grain (Local Error - y)','RK4 Grain (Global Error - y)','RK4 Grain (Local Error - y)','dt','dt^2','dt^3','dt^4')
xlabel('dt')
ylabel(' Error')
legend
title(' grain y')

subplot(2,2,4)
loglog(dt, error.euler.chaff_global(2,:), '-^');
hold on
loglog(dt,error.euler.chaff_local(2,:),'-o');
hold on
loglog(dt, error.rkutta.chaff_global(2,:), '-d');
hold on
loglog(dt,error.rkutta.chaff_local(2,:),'-x');
hold on
loglog(dt,dt,dt,dt.^2,dt,dt.^3,dt,dt.^4);
legend('Euler Chaff (Global Error - y)','Euler Chaff (Local Error - y)','RK4 Chaff (Global Error - y)','RK4 Chaff (Local Error - y)','dt','dt^2','dt^3','dt^4')
%grid on
xlabel('dt')
ylabel(' Error')
legend
title(' chaff y')


%% Min velocity for separation of grain and chaff - (Euler)
%% Definition of time vector
u_0 = 0.05:0.001:0.2; % Jet velocities to iterate over
tstart = 0;
tend = 2;
t = linspace(tstart, tend, 151); % Define time vector
dt = t(2) - t(1); % Calculate time step

%% Definition of solution arrays
uf = zeros(2, length(t)); % Air velocity
up = zeros(2, length(t)); % Particle velocity
Xp = zeros(2, length(t)); % Particle position

%% Initial conditions
Xp(:, 1) = IC.Xp_0; % Initial position
up(:, 1) = IC.up_0; % Initial velocity

%% Figure setup
figure;
hold on;
xlabel('x-position [m]');
ylabel('y-position [m]');
xlim([0.45 0.65]);
ylim([-0.5 0.55]);

%% Iterate over different jet velocities (u_0)

min_u0_chaff = NaN; % Minimum u_0 for chaff separation

for l = 1:length(u_0)
    % Update initial air velocity for the current iteration
    uf(1, 1) = u_0(l);
    
    % Solve for each type of particle (grain or chaff)
    for m = 1:2
        if m == 1 % Grain
            dia = d.grain;
            density = rho.grain;
            particle_name = 'Grain';
            color = 'b'; % Blue for grain
        elseif m == 2 % Chaff
            dia = d.chaff;
            density = rho.chaff;
            particle_name = 'Chaff';
            color = 'r'; % Red for chaff
        end

        % Calculate particle properties
        Vp = f_Vp(dia); % Volume of particle
        Fb = f_Fb(density); % Buoyancy force

        % Euler solution for this jet velocity and particle type
        for i = 1:(length(t) - 1)
            
            Re = f_Re(uf(:, i), up(:, i), dia); % Reynolds number
            Cd = Calc_Cd(Re); % Drag coefficient
            Fd = f_Fd(dia, density, Vp, Cd, uf(:, i), up(:, i)); 

            % Update velocity and position using Euler method
            up(:, i + 1) = up(:, i) + dt * f_d_up(Fb, Fd); % Update velocity
            Xp(:, i + 1) = Xp(:, i) + dt * up(:, i); % Update position

            % Update air velocity
            uf(1, i + 1) = Calc_uf(Xp(1, i + 1), Xp(2, i + 1), h, u_0(l)); 
        end

        % Check separation conditions
        if m == 1  && Xp(1, end) <= IC.Xp_c(1)
            min_u0_grain = u_0(l); % Record minimum u_0 for grain
        elseif m == 2 && isnan(min_u0_chaff) && Xp(1, end) >= IC.Xp_c(1)
            min_u0_chaff = u_0(l); % Record minimum u_0 for chaff
        end

        % Plot trajectory
        plot(Xp(1, :), Xp(2, :), 'Color', color, 'DisplayName', [particle_name ', u_0=' num2str(u_0(l))]);
    end

    % Break loop if separation is achieved for both particle types
    if ~isnan(min_u0_grain) && ~isnan(min_u0_chaff)
        break;
    end
end

%% Display Results
legend('show'); % Add legend
title('Particle Separation Trajectories');
if ~isnan(min_u0_grain) && ~isnan(min_u0_chaff)
    
    fprintf('Minimum u_0 : %.3f m/s\n', min_u0_chaff);

end



      














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
        uf = 6.2 * u_0 * sqrt(h / x) * exp(-50 * (y^2 / x^2)); % [m/s] (Air velocity)
    else
        uf = 0;
    end
end


            


    
















