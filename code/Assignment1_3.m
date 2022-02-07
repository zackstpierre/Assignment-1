%ELEC4700 Assignment 1
%Zachary St. Pierre
%101094217



%-------------------------------------------------------------------------%
clc;
clear;
close all;




%-------------------------------------------------------------------------%
%Variables%

N_e = 100;  %number of electrons
N_plot = randi(N_e, 5, 1);  %number of electrons to plot

m_0 = 9.1093837015e-31;
m_n = 0.26 * m_0;      % mass of electrons

c_bolt = 1.381e-23;
Temp = 300;

size_x = 200e-9; %length of region
size_y = 100e-9; %width of region

t_mn = 0.2e-12;

%define box regions

box_length = 25e-9;
box_height = 40e-9;

rectangle('Position', [(size_x/2) - (box_length/2) , 0 , box_length, box_height])
hold on
rectangle('Position', [(size_x/2) -  (box_length/2) , (size_y - box_height) , box_length, box_height])
hold on

axis([0 size_x 0 size_y])

%duration and time step
nT = 5e-15;
T = 100*nT;


%probability of scatter

Pscat = 1 - exp(-nT/t_mn);
%Pscat = 0.5;
Pscatter = [];

%velocity

Vth = sqrt(2.*c_bolt.*Temp./m_n)

Vx = [];
Vy = [];

% position
Px = [];
Py = []; 

%Maxwell-Boltzman deviation calculations
deviation = sqrt((c_bolt*Temp)/m_n);
rand_speed = deviation + Vth;


%-------------------------------------------------------------------------%
%Create random points%
Px = rand(N_e,1)*size_x;
Py = rand(N_e,1)*size_y;


Vx = randn(N_e,1) * rand_speed;
Vy = randn(N_e,1) * rand_speed;


P_inside_box_x = ((size_x/2) - (box_length/2) < Px) &  (Px < (size_x/2) + box_length/2) & (box_height > Py)  | (Py > box_height + (size_y - 2* box_height));

while sum(P_inside_box_x) > 0
    Px(P_inside_box_x) = rand(1)*size_x;
    Py(P_inside_box_x) = rand(1)*size_y;

    P_inside_box_x = ((size_x/2) - (box_length/2) < Px) &  (Px < (size_x/2) + box_length/2) & (box_height > Py)  | (Py > box_height + (size_y - 2* box_height));
end


counter = 1;
for i = 1:1000

   Pscatter = rand(N_e,1);

   Px_old = Px;
   Py_old = Py;

   P_scatter_locations = Pscatter < Pscat;
   Vy(P_scatter_locations) = randn(1) * Vth;
   Vx(P_scatter_locations) = randn(1) * Vth;
    

   Px = Px_old + (Vx * nT);
   Py = Py_old + (Vy * nT);

   Px_outside_right = Px > size_x;
   Px_outside_left  = Px < 0;

   P_inside_box_x = ((size_x/2) - (box_length/2) < Px) &  (Px < (size_x/2) + box_length/2) & ((box_height > Py)  | (Py > box_height + (size_y - 2* box_height)));
  
   
   Py_top_box = Py > box_height + (size_y - 2* box_height) ;
   Py_bottom_box = box_height > Py;

   P_xy_top = P_inside_box_x + Py_top_box;
   P_xy_bottom = P_inside_box_x + Py_bottom_box;

   P_gap_top = P_xy_top > 1;
   P_gap_bottom = P_xy_bottom > 1;
   
   
   Vy(P_gap_top) = Vy(P_gap_top) * -1;
   Vy(P_gap_bottom) =   Vy(P_gap_bottom) * -1;
   Py(P_gap_top) = Py_old(P_gap_top);
   Py(P_gap_bottom) = Py_old(P_gap_bottom);

  
   Vx(P_inside_box_x) = Vx(P_inside_box_x) * -1;
   Px(P_inside_box_x) = Px_old(P_inside_box_x);

  
   Py_top = Py > size_y;
   Py_bottom = Py < 0;

   Px(Px_outside_right) = Px(Px_outside_right) - size_x;
   Px(Px_outside_left) = Px(Px_outside_left) + size_x;

   Px_old(Px_outside_right) = Px_old(Px_outside_right) - size_x;
   Px_old(Px_outside_left) = Px_old(Px_outside_left) + size_x;

  

   Vy(Py_top) = Vy(Py_top) * -1;
   Vy(Py_bottom) = Vy(Py_bottom) * -1;
  
  


   for j = 1: length(N_plot)
    plot([Px_old(N_plot(j)) Px(N_plot(j))]',[Py_old(N_plot(j)) Py(N_plot(j))]','SeriesIndex',j)
    hold on
   end 

   avg_velocity = mean(mean(sqrt(Vx.^2 + Vy.^2)));
   avg_temperature(i) = (avg_velocity.^2 .* m_n)./(2.*c_bolt);
   t_temp(i) = nT*i;

    title(avg_temperature(i)); title("T = " + avg_temperature(i) + "K");
    


    pause(0.01)
end


figure(2)
plot(t_temp,avg_temperature)
title('Average temperature of Silicon Plane')
xlabel('Time')
ylabel('Temperature (K)')