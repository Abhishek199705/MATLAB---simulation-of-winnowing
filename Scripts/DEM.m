%% Problem 5.7: Collisions of particles, comparison of hard-sphere and DEM
% Name: Abhishek Joel George
% Matriculation Number: 250900


%% Refresh

clc
clear all
close all

%% Collision time and hard sphere algorithm

x1 = [0, 0, 0];
x2 = [1.1, 1.3, 0];
v1 = [0, 0, 0];
v2 = [-1, -1, 0];
d1 = 1;
d2 = 1;
m1 = 0.05;
m2 = 0.05;
e = [1, 0.75];
dissipation = zeros(length(e),1);



% Compute collision time
tcol = collisiontime(x1,x2,v1,v2,d1,d2);

% Advance particles to collison time

x1 = x1 + tcol*v1;
x2 = x2 + tcol*v2;
x1
x2

% Contact point before collision

contact_hard_bcol = 0.5*(x1 + x2); 

% Compute post-collision state

for i = 1:length(e)
    [v1f(i,:), v2f(i,:), dissipation(i,:)] = performcollision(m1,m2,e(i),x1,x2,v1,v2);
end

%% Refresh - Soft Sphere
%{
clc
clear all
close all
%}


%% Initial conditions of particles

x1 = [0, 0, 0];
x2 = [1.1, 1.3, 0];
v1 = [0, 0, 0];
v2 = [-1, -1, 0];
d1 = 1;
d2 = 1;
m = 0.05;
k = 1;


%% Timestep

dt = 1e-4;
n = 1/dt;

%% Leap frog 


[delta,normal] = overlap(x1,x2,d1,d2);
F = force(delta,normal);

for i = 1:n
    
    
    %F = force(delta,normal,k);
    % particle 2
    vpred2 = v2 + 0.5*dt*(F/m);
    x2 = x2 + vpred2*dt;

    % particle 1
    vpred1 = v1 + 0.5*dt*(-F/m);
    x1 = x1 + vpred1*dt;

    [delta,normal] = overlap(x1,x2,d1,d2);
    
    F = force(delta,normal,k);

    % velocity corrector
    v2 = vpred2 + 0.5*dt*(F/m);
    v1 = vpred1 + 0.5*dt*(-F/m);
    v12 = v1 - v2;


end


%% 
%{
clc
clear all
%}

%% Account for co-efficient of restitution

e1 = 1;
e2 = 0.75;
k_load = 750;
x1 = [0, 0, 0];
x2 = [1.1, 1.3, 0];
v1 = [0, 0, 0];
v2 = [-1, -1, 0];
d1 = 1;
d2 = 1;
m = 0.05;
dt = 1e-4;
n = 1/dt;
t_col =0;
count =0;

% Solution arrays

x1e = zeros(2,3);
x2e = zeros(2,3);
v1e = zeros(2,3);
v2e = zeros(2,3);

%% Initialize the solution arrays

x1e(1,:) = x1;
x2e(1,:) = x2;
v1e(1,:) = v1;
v2e(1,:) = v2;

x1e(2,:) = x1;
x2e(2,:) = x2;
v1e(2,:) = v1;
v2e(2,:) = v2;

contact_point_bf = zeros(2,3);
contact_point_af = zeros(2,3);

%% Solution loop
%% 

for j = 1:2
    if j == 1
        e = e1;
    elseif j == 2
        e = e2;
    end
    k_unload = e^2 * k_load;
    % Leap frog
    [delta,normal] = overlap(x1e(j,:),x2e(j,:),d1,d2);
    v12e = v1e(j,:) - v2e(j,:);
    
    if dot(v12e,normal)<=0
        k = k_load;
    else
        k = k_unload;
    end
    F = force(delta,normal,k_load);
    
    

    
    for i = 1:n
        if dot(v12e,normal)<0
            k = k_load;
        elseif dot(v12e,normal)>0
            k = k_unload;
        end
        
        
        % 
        vpred2 = v2e(j,:) + 0.5*dt*(F/m);
        x2e(j,:) = x2e(j,:) + vpred2*dt;
    
        % particle 1
        vpred1 = v1e(j,:) + 0.5*dt*(-F/m);
        x1e(j,:) = x1e(j,:) + vpred1*dt;
    
        [delta,normal] = overlap(x1e(j,:),x2e(j,:),d1,d2);
        
        if delta>0
            t_col = t_col + dt;  %duration of collision for e2
        end

        if t_col == 0
            count = count +1; % count the steps time before the start of collision
        end
        t_bcol = count * dt; % duration before coll
        t_acol = t_bcol + t_col;
        count_acol = t_acol/dt;

        if i == 5391     % contact point
            contact_point_af = 0.5*(x1e(j,:) +x2e(j,:)); % contact point after collision
        end

        
        F = force(delta,normal,k);
    
        % velocity corrector
        v2e(j,:) = vpred2 + 0.5*dt*(F/m);
        v1e(j,:) = vpred1 + 0.5*dt*(-F/m);
        v12e = v1e(j,:) - v2e(j,:);
        
    end
end

%% Overlap/normal
function [delta,normal] = overlap(x1,x2,d1,d2)
normal = (x1-x2)/norm(x1-x2);
delta = (0.5*d1 + 0.5*d2) - norm(x1-x2);
if delta>0
    delta = delta;
else
    delta = 0;
end
end

%% Force

function [F] = force(delta,normal,k)
if delta>0
    F = -k*delta*normal;
else
    F = [0 0 0];
end
end

%% DPM
function [v1f, v2f, dissip] = performcollision(m1,m2,e,x1,x2,v1i,v2i)
kEi = 0.5*m1*dot(v1i,v1i) + 0.5*m2*dot(v2i,v2i);
mstar = (m1*m2)/(m1+m2);
x12 = x1 - x2;
v12 = v1i - v2i;
n = x12/sqrt(dot(x12,x12));
J = -(1+e)*dot(v12,n)*mstar*n;
v1f = v1i + J/m1
v2f = v2i - J/m2
kEf = 0.5*m1*dot(v1f,v1f) + 0.5*m2*dot(v2f,v2f);
dissip = kEi - kEf
end
%% collision time
function[tcol] = collisiontime(x1,x2,v1,v2,d1,d2)
x12 = x1-x2;
v12 = v1-v2;
A = dot(v12,v12);
B = 2*dot(x12,v12);
C = dot(x12,x12) - 0.25*(d1+d2)^2;
D = B^2 - 4*A*C;
if (D<=0)
    tcol = -1;
else
    tcol1 = (-B-sqrt(D))/(2*A);
    tcol2 = (-B+sqrt(D))/(2*A);
    if (tcol2<0)
        tcol = -1;
    else
        if (tcol1>0)
            tcol = tcol1;
        else
            tcol = tcol2;
        end
    end
end
end


