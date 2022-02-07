%ELEC4700 Assignment 1
%Zachary St. Pierre
%101094217



%-------------------------------------------------------------------------%
clc;
clear;
close all;




%-------------------------------------------------------------------------%
%Variables%

N_e = 1000;  %number of electrons
N_plot = randi(N_e, 5, 1);  %number of electrons to plot

m_0 = 9.1093837015e-31;
m_n = 0.26 * m_0;      % mass of electrons

c_bolt = 1.381e-23;
Temp = 300; 

size_x = 200e-9; %length of region
size_y = 100e-9; %width of region


%duration and time step
nT = 5e-15;
T = 100*nT;


%velocity

Vth = sqrt(2.*c_bolt.*Temp./m_n)

Vx = [];
Vy = [];

% position
Px = [];
Py = []; 



%-------------------------------------------------------------------------%
%Create random points%
Px = rand(N_e,1)*size_x;
Py = rand(N_e,1)*size_y;

Vx = randn(N_e,1) * Vth;
Vy = randn(N_e,1) * Vth;


figure(1)
hold on
axis([0 size_x 0 size_y]);

for i = 1:1000

   Px_old = Px;
   Py_old = Py;


    

   Px = Px_old + (Vx * nT);
   Py = Py_old + (Vy * nT);

   Px_outside_right = Px > size_x;
   Px_outside_left  = Px < 0;

   Py_top = Py > size_y;
   Py_bottom = Py < 0;

   Px(Px_outside_right) = Px(Px_outside_right) - size_x;
   Px(Px_outside_left) = Px(Px_outside_left) + size_x;

   Px_old(Px_outside_right) = Px_old(Px_outside_right) - size_x;
   Px_old(Px_outside_left) = Px_old(Px_outside_left) + size_x;

   Vy(Py_top) = (Vy(Py_top) + 180) * -1;
   Vy(Py_bottom) = (Vy(Py_bottom) + 180) * -1;
  


   for j = 1: length(N_plot)
    plot([Px_old(N_plot(j)) Px(N_plot(j))]',[Py_old(N_plot(j)) Py(N_plot(j))]','SeriesIndex',j)
    hold on
   end 


   avg_velocity = mean(mean(sqrt(Vx.^2 + Vy.^2)));
   avg_temperature(i) = (avg_velocity.^2 .* m_n)./(2.*c_bolt);
   t_temp(i) = nT*i;

    title("T = " + avg_temperature(i) + "K");

    pause(0.01)
end

figure(2)
plot(t_temp,avg_temperature)
title('Average temperature of Silicon Plane')
xlabel('Time')
ylabel('Temperature (K)')
