function Part1_Func(delV)

%Constants
m_o = 9.11e-31; %Rest masss of electron [kg]
m_ef = 0.26*m_o; %Effective mass of an electron
kB = 1.3806e-23; %Bolzmann's constant [m^2 kg s^-2 k^-1]
T_lat = 300; %Lattice Temperature [k]
tow_mn = 0.2e-12; %Mean time between collisions [sec]
q_e = 1.602e-19; %Charge of an electron [C]

v_th = sqrt((2*kB*T_lat)/(m_ef));

region_x = 200e-9; %Width of region [m]
region_y = 100e-9; %Height of region [m]

num_p = 1000; %Number of particles
num_dt = 1000; %Number of time steps
num_plot = 10; %Number of particles to plot

dt = ((1/100)*region_x)/v_th; %Time step [sec]
dx(1,:) = zeros(1,num_p); %X Position change
dx(1,:) = zeros(1,num_p); %X Position change


E_p1 = delV/region_x;
F_p1 = q_e*E_p1;
a_p1 = F_p1/m_ef;
%Position and Velocity Arrays
Px(1,:) = zeros(1,num_p); %X-coordinate
Py(1,:) = zeros(1,num_p); %Y-coordinate
Vx(1,:) = zeros(1,num_p); %X-velocity
Vy(1,:) = zeros(1,num_p); %Y-velocity

%Top and bottom spectral vs. diffusive
top_bc = 1; %0 = diff, 1 = spec
bottom_bc = 1; %0 = diff, 1 = spec
%Generate colours for number of plotted particles
colours = zeros(3,num_plot);
for d=1:num_plot
    colours(:,d) = [rand;rand;rand];
end
%colours = {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'};

V_rand_store(1,:) = zeros(1,num_p); 


Pscat = 1 - exp(-(dt)/(tow_mn)); %Scattering Probability Function
for j=1:num_p
   %Assign each particle a random position to start not inside the boxes
   keep_trying = 1;
   
  Px(j) = rand*region_x;
  Py(j) = rand*region_y;

 
   %Random velocity from dist with random direction
   %Random Velocity
   v_rand = sqrt(kB*T_lat/m_ef)*randn+v_th; %Normal Dist
   V_rand_store(j) = v_rand;
   theta = rand*360; %random direction
   
   Vx(j) = v_rand*cosd(theta);
   Vy(j) = v_rand*sind(theta);
   
end
figure('Name', 'Velocity Histogram')
histogram(V_rand_store);
xlabel('Velocity (m/s)')
ylabel('Number of Particles with Velocity')
title('Particle Velocity Histogram')
%Next position of plotted particles
Px_next = zeros(1,num_p);
Py_next = zeros(1,num_p);

%Old position of plotted particles
Px_old = zeros(1,num_p);
Py_old = zeros(1,num_p);
T = zeros(1,num_p);
T_avg = zeros(1,num_dt);
time = zeros(1,num_dt);

%Time since last scattered
mft = zeros(1,num_p);
avg_time_array = zeros(1,num_p);

%Concentration
n = 10^15 / 0.0001; %[cm^2] to [m^2]
Jx = zeros(1,num_dt);
time_curr = zeros(1,num_dt);
running_time = 0;
figure('Name', 'P1: Particle Traj')
for i=1:num_dt
    %Current
    Jx(i) = q_e.*n.*mean(Vx);
    time_curr(i) = running_time;
    running_time = running_time + dt;
    %Scattering
    rand_scat = rand(1,num_p);

    Iscat = find(Pscat > rand_scat);
    for z=1:length(Iscat)
        v_rand = sqrt(kB*T_lat/m_ef)*randn+v_th; %Normal Dist
        theta = rand*360; %random direction
        Vx(Iscat(z)) = v_rand*cosd(theta);
        Vy(Iscat(z)) = v_rand*sind(theta);

        if(avg_time_array(Iscat(z)) ~= 0) %Not first value
            avg_time_array(Iscat(z)) = mean([avg_time_array(Iscat(z)) mft(Iscat(z))]);
        else
            avg_time_array(Iscat(z)) = mft(Iscat(z));
        end
        mft(Iscat(z)) = 0;
    end

    
    
    %Acceleration
    ax = a_p1;
    ay = 0;
    %Get next position
    Vx_next = Vx + ax.*dt;
    Px_next = Px + (Vx_next).*dt;
    Py_next = Py + (Vy + ay.*dt).*dt;
    
%X Boundary Conditions
    %Left and Right
    Px(Px_next > region_x) = 0;
    Px(Px_next < 0) = region_x;
    
    %Top
    Vy(Py >= region_y) = -1*(Vy(Py >= region_y) + ay.*mft(Py >= region_y));
    
    %Bottom
    Vy(Py <=0) = -1*(Vy(Py <=0) + ay.*mft(Py <=0));
    
    Px_old = Px;
    Vx = Vx + ax.*dt;
    Px = Px + (Vx).*dt; %Update x position
    Py_old = Py;
    Py = Py + Vy.*dt; %Update y position
    
    mft = mft + dt; %Add time since last scatter or bounce
    
    %Calculate the semiconductor temperature
    v_avg = mean(sqrt(Vx.^2 + Vy.^2));
    T_avg(i) = (v_avg^2*m_ef)/(2*kB);
    
    %Calculate the temperature of each electron
    T = ((sqrt(Vx.^2 + Vy.^2)).^2.*m_ef)./(2*kB);
    
    time(i) = i*dt;
    %Plot subet of particles
    for k=1:num_plot
        hold on;
        plot([Px_old(k) Px(k)],[Py_old(k) Py(k)], 'Color', colours(:,k));
    end
    if(i == 1)
        xlabel('X Position')
        ylabel('Y Position')
        xlim([0 region_x])
        ylim([0 region_y])
    end
    title(['Particle Trajectories (Avg Temp = ',num2str(T_avg(i)), ' K)'])
    

    fprintf('Iteration = %i \n', i);
    pause(0.05);
    
end
%Average time to scatter
tau_avg = mean(avg_time_array);
fprintf('Average time between collisions = %d \n', tau_avg);

%Measure the Mean Free Path
mfp_avg = v_avg*tau_avg;
fprintf('Average mean free path = %d \n', mfp_avg);

%Plot Current
figure('Name', 'P1: Current')
plot(time_curr, Jx);
xlabel('Time (seconds)')
ylabel('Current')
title('Current Plot')

%Electron Density Map
figure('Name', 'P1: Density')
[den, cent] = hist3([Px',Py'], 'Nbins', [40 20]);
[X,Y] = meshgrid(cent{1}, cent{2});
surf(X,Y,den')
colorbar
xlabel('Position X')
ylabel('Position Y')
title('Electron Density')

%Temperature Map
%Bin the electrons
bin_length = 200e-9/40;
bin_width = 100e-9/20;
temp_binned = zeros(40,20);
for i=1:num_p
    for j=1:length(cent{1}) %X bin
         if(j == 1)
            x_in_bin = (Px > 0) & (Px <= (cent{1}(j)+bin_length/2));
         else
            x_in_bin = (Px > (cent{1}(j-1)-bin_length/2)) & (Px <= (cent{1}(j)+bin_length/2));
         end
       for k = 1:length(cent{2}) %Y bin
          if(k == 1)
              y_in_bin = (Py > 0) & (Py <= (cent{2}(k)+bin_width/2));
          else
              y_in_bin = (Py > (cent{2}(k-1)-bin_width/2)) & (Py <= (cent{2}(k)+bin_width/2));
          end
          
          curr_sum = sum(T(x_in_bin & y_in_bin));
          if(curr_sum == 0)
              temp_binned(j,k) = 0;
          else
              temp_binned(j,k) = mean(T(x_in_bin & y_in_bin));
          end
          
       end
    end
end
figure('Name', 'P1: Temperature')

surf(X,Y,temp_binned')
colorbar
xlabel('Position X')
ylabel('Position Y')
title('Temperature Map')

end

