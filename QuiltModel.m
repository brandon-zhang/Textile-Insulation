function QuiltModel

H_sample = 0.1016; %Sample height, m
L_sample = 0.1524; %Sample length, m
W_seam = 0.005; %Width of seam depression, m
L_seam = 0.001; %Thickness at seam, m
n_seam_horz = 1; %Number of horizontal seams
n_seam_vert = 1; %Number of vertical seams
k_batt = 0.05; %Thermal conductivity of batting, W/m K

k_thread = 0.05; %Thermal conductivity of thread, W/m K
A_thread = (0.05/1000)^2 * pi; %Thickness of thread, m2
tpi = 15; %Threads per inch

rho_air = 1.225; %Density of air, kg/m3
cp_air = 1000; %Heat capacity of air, J/kg K

Q_heater = 18.22*0.32; %Heat FLOW into the metal plate
rho_al = 2.7; %Aluminum density, g/cm3
cp_al = 0.9; %Heat capacity of aluminum, J/g K
V_al = 0.0154*0.4*100*100; %Volume of aluminum heater plate, cm3

A_cond = (H_sample*L_sample) - W_seam * (n_seam_horz*L_sample + n_seam_vert*H_sample);

k_air = 0.026; %Thermal conductivity of air, W/mK
L_metal_fabric = 0.0001; %Metal-Fabric air gap, m
R_metal_fabric = L_metal_fabric/(k_air*A_cond); %Thermal resistance of metal-fabric air gap, K/W

k_fabric = 0.08; %Thermal conductivity of fabric, W/mK
L_fabric = 0.0056/6; %Thickness of fabric, m
R_fabric = L_fabric/(k_fabric*A_cond); %Thermal resistance of fabric, K/W

L_fabric_batting = 0; %Fabric-batting air gap, m
R_fabric_batting = L_fabric_batting/(k_air*A_cond); %Thermal resistance of fabric-batting air gap, K/W

k_batting = 0.05; %Thermal conductivity of batting, W/mK
L_batting = 0.003; %Thickness of batting, m
R_batting = L_batting/(k_batting*A_cond); %Thermal resistance of batting, K/W

R_cond_total = 2*(R_metal_fabric + R_fabric + R_fabric_batting) + R_batting; %Total thermal resistance of non-seam area

g = 9.81; %Acceleration due to gravity, m/s2
visc = 1.7*10^-5; %Kinematic viscosity, m2/s. Replace with non-constant version!
diffus = 1.9*10^-5; %Thermal diffusivity, m2/s.
Pr = 0.71; %Prantdl number, unitless

L_total = 2*L_fabric+L_batting;
theta = atan((L_total-L_seam)/W_seam); %Seam angle
theta_conv = theta * 180/pi
V_seam_vert = (W_seam^2 * tan(theta))/4 * H_sample; %Vertical seam volume, m3
V_seam_horz = (W_seam^2 * tan(theta))/4 * L_sample; %Horizontal seam volume, m3

rho_fabric = 1.225; %Density of fabric, kg/m3
cp_fabric = 1000; %Heat capacity of fabric, J/kg K
V_fabric_vert = W_seam/cos(theta) * L_fabric * H_sample; %Volume of fabric in seam
V_fabric_horz = W_seam/cos(theta) * L_fabric * L_sample; %Volume of fabric in seam

T_cold = 27.5; %Temperature of cold plate, C

tspan = [0: 3600];
Z0 = [T_cold, T_cold, T_cold, T_cold, T_cold, T_cold, T_cold, T_cold, T_cold];

[T, zout] = ode45(@slope_func, tspan, Z0);
[T2, zout2] = ode45(@slope_func2, tspan, Z0);

plot(T, zout(:,1));
hold on
plot(T2, zout2(:,1));
legend({'Conductivity + Convection','Only Conductivity'},'Location','southeast')
hold off

cond_conv_dT = zout(end, 1) - T_cold
cond_dT = zout2(end,1) - T_cold

    function result = Nusselt(dT, L_char, T_air)
        %Based on Churchill and Chu's approx. for the Nusselt number
        Ra = (g*(1/T_air)*abs(dT)*L_char^3)/(visc*diffus);
        result = 0.68 + 0.67*Ra^(1/4) * (1 + (0.492/Pr)^(9/16))^(-4/9);
    end

    function dzdt = slope_func(T, Z)
        
        dzdt = [0;0;0;0;0;0;0;0;0];
        
        T_heater = Z(1); %Heater temperature
        
        T_air_1_vert = Z(2); %Temperature of the vertical first air gap
        T_fabric_1_vert = Z(3); %Temperature of the first fabric in the vertical seam
        T_fabric_2_vert = Z(4); %Temperature of the second fabric in the vertical seam
        T_air_2_vert = Z(5); %Temperature of the second vertical air gap
        
        T_air_1_horz = Z(6); %Temperature of the horizontal first air gap
        T_fabric_1_horz = Z(7); %Temperature of the first fabric in the horizontal seam
        T_fabric_2_horz = Z(8); %Temperature of the second fabric in the horizontal seam
        T_air_2_horz = Z(9); %Temperature of the second horizontal air gap
        
        %For vertical seams:
        
        R_metal_air_1_vert = H_sample/(k_air*H_sample*W_seam*Nusselt(T_heater - T_air_1_vert, H_sample, T_air_1_vert));
        Q_metal_air_1_vert = (T_heater - T_air_1_vert)/R_metal_air_1_vert;
        %Convective heating resistance, with h_sample as L_char
        R_air_fabric_1_vert = H_sample/(k_air*H_sample*2*(W_seam/(2*cos(theta)))*Nusselt(T_air_1_vert - T_fabric_1_vert, H_sample, T_air_1_vert));
        Q_air_fabric_1_vert = (T_air_1_vert - T_fabric_1_vert)/R_air_fabric_1_vert;
        %Convective cooling resistance into the fabric
        dT_air_1_vert = 1/(V_seam_vert*rho_air*cp_air) * (Q_metal_air_1_vert - Q_air_fabric_1_vert);
        dzdt(2) = dT_air_1_vert;
        %Change in temperature of the air pocket based on net heat flow
        
        R_batt_comp_vert = 2*tan(theta)/(k_batt*H_sample*(-log(L_seam) + log(L_seam + W_seam*tan(theta))));
        %2 * the integration of the conductivity for the compressed
        %batting, assuming a triangle
        R_stitch_vert = (H_sample*tpi*L_seam)/(k_thread * A_thread);
        %Resistance of the stitches, modeling as cylinders that pass
        %through the fabric
        R_seam_vert = 1/(1/R_stitch_vert + 1/(2*R_fabric + R_batt_comp_vert));
        %Resistance of the seam, assuming parallel conductance through
        %stitches and compressed batting
        Q_seam_vert = (T_fabric_1_vert - T_fabric_2_vert)/R_seam_vert;
        %Heat flow through seam from Fourier's Law
        dT_fabric_1_vert = 1/(V_fabric_vert*rho_fabric*cp_fabric) * (Q_air_fabric_1_vert - Q_seam_vert);
        dzdt(3) = dT_fabric_1_vert;
        %Change in temperature of the first fabric layer based on net heat
        %flow
        
        R_air_fabric_2_vert = H_sample/(k_air*H_sample*W_seam*Nusselt(T_fabric_2_vert - T_air_2_vert, H_sample, T_air_2_vert));
        Q_air_fabric_2_vert = (T_fabric_2_vert-T_air_2_vert)/R_air_fabric_2_vert;
        %Convective heating resistance, with h_sample as L_char
        dT_fabric_2_vert = 1/(V_fabric_vert*rho_fabric*cp_fabric) * (Q_seam_vert-Q_air_fabric_2_vert);
        dzdt(4) = dT_fabric_2_vert;
        R_metal_air_2_vert = H_sample/(k_air*H_sample*2*(W_seam/(2*cos(theta)))*Nusselt(T_air_2_vert - T_cold, H_sample, T_air_1_vert));
        Q_metal_air_2_vert = (T_air_2_vert - T_cold)/R_metal_air_2_vert;
        %Convective cooling resistance into the fabric
        dT_air_2_vert = 1/(V_seam_vert*rho_air*cp_air) * (Q_air_fabric_2_vert - Q_metal_air_2_vert);
        dzdt(5) = dT_air_2_vert;
        %Change in temperature of the air pocket based on net heat flow
        
        %For horizontal seams:
        
        R_metal_air_1_horz = L_sample/(k_air*L_sample*W_seam*Nusselt(T_heater - T_air_1_horz, L_sample, T_air_1_horz));
        Q_metal_air_1_horz = (T_heater - T_air_1_horz)/R_metal_air_1_horz;
        %Convective heating resistance, with L_sample as L_char
        R_air_fabric_1_horz = L_sample/(k_air*L_sample*2*(W_seam/(2*cos(theta)))*Nusselt(T_air_1_horz - T_fabric_1_horz, L_sample, T_air_1_horz));
        Q_air_fabric_1_horz = (T_air_1_horz - T_fabric_1_horz)/R_air_fabric_1_horz;
        %Convective cooling resistance into the fabric
        dT_air_1_horz = 1/(V_seam_horz*rho_air*cp_air) * (Q_metal_air_1_horz - Q_air_fabric_1_horz);
        dzdt(6) = dT_air_1_horz;
        %Change in temperature of the air pocket based on net heat flow
        
        R_batt_comp_horz = 2*tan(theta)/(k_batt*L_sample*(-log(L_seam) + log(L_seam + W_seam*tan(theta))));
        %2 * the integration of the conductivity for the compressed
        %batting, assuming a triangle
        R_stitch_horz = (L_sample*tpi*L_seam)/(k_thread * A_thread);
        %Resistance of the stitches, modeling as cylinders that pass
        %through the fabric
        R_seam_horz = 1/(1/R_stitch_horz + 1/(2*R_fabric + R_batt_comp_horz));
        %Resistance of the seam, assuming parallel conductance through
        %stitches and compressed batting
        Q_seam_horz = (T_fabric_1_horz - T_fabric_2_horz)/R_seam_horz;
        %Heat flow through seam from Fourier's Law
        dT_fabric_1_horz = 1/(V_fabric_horz*rho_fabric*cp_fabric) * (Q_air_fabric_1_horz - Q_seam_horz);
        dzdt(7) = dT_fabric_1_horz;
        %Change in temperature of the first fabric layer based on net heat
        %flow
        
        R_air_fabric_2_horz = L_sample/(k_air*L_sample*W_seam*Nusselt(T_fabric_2_horz - T_air_2_horz, L_sample, T_air_2_horz));
        Q_air_fabric_2_horz = (T_fabric_2_horz-T_air_2_horz)/R_air_fabric_2_horz;
        %Convective heating resistance, with L_sample as L_char
        dT_fabric_2_horz = 1/(V_fabric_horz*rho_fabric*cp_fabric) * (Q_seam_horz-Q_air_fabric_2_horz);
        dzdt(8) = dT_fabric_2_horz;
        R_metal_air_2_horz = L_sample/(k_air*L_sample*2*(W_seam/(2*cos(theta)))*Nusselt(T_air_2_horz - T_cold, L_sample, T_air_1_horz));
        Q_metal_air_2_horz = (T_air_2_horz - T_cold)/R_metal_air_2_horz;
        %Convective cooling resistance into the fabric
        dT_air_2_horz = 1/(V_seam_horz*rho_air*cp_air) * (Q_air_fabric_2_horz - Q_metal_air_2_horz);
        dzdt(9) = dT_air_2_horz;
        %Change in temperature of the air pocket based on net heat flow
        
        R_vert = n_seam_vert/(1/R_metal_air_1_vert + 1/R_air_fabric_1_vert + 1/R_seam_vert + 1/R_air_fabric_2_vert + 1/R_metal_air_2_vert);
        R_horz = n_seam_horz/(1/R_metal_air_1_horz + 1/R_air_fabric_1_horz + 1/R_seam_horz + 1/R_air_fabric_2_horz + 1/R_metal_air_2_horz);
        R_total = 1/(1/R_cond_total + 1/R_vert + 1/R_horz);
        Q_heater_out = (T_heater - T_cold)/R_total;
        dT_heater = (1/(V_al*rho_al*cp_al)) * (Q_heater - Q_heater_out);
        dzdt(1) = dT_heater;
        for i = 1:9
            if (dzdt(i) < 0)
                dzdt(i) = 0;
            end
        end
    end

    function dzdt = slope_func2(T, Z)
        dzdt = [0;0;0;0;0;0;0;0;0];
        R_cond_adjusted = R_cond_total * A_cond/(H_sample*L_sample);
        T_heater = Z(1); %Heater temperature
        Q_heater_out = (T_heater - T_cold)/R_cond_adjusted;
        dT_heater = (1/(V_al*rho_al*cp_al)) * (Q_heater - Q_heater_out);
        dzdt(1) = dT_heater;
    end
end

