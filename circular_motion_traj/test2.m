clc
clear all
load('data.mat')
[begin_x, begin_y] = path();
time = begin_x.x;
pre_traj = [];
pre_traj(1:3,:)= begin_x.y(1:3,:);
pre_traj(4,:) = -begin_x.y(6,:);
figure
plot(time, pre_traj(1,:), 'k-');xlabel('t')
hold on
plot(time, pre_traj(2,:)); 
hold on
plot(time, pre_traj(3,:)); 
plot(time, pre_traj(4,:));

t_f = 6;
new_time = 0:0.1:t_f;
traj_spline(1,:) = spline(time, pre_traj(1,:),new_time);
disp(traj_spline(1,:))
traj_spline(2,:) = spline(time, pre_traj(2,:),new_time);
traj_spline(3,:) = spline(time, pre_traj(3,:),new_time);
traj_spline(4,:) = spline(time, pre_traj(4,:),new_time);
final_traj(1,:) = [traj_spline(1,:) x(2:end) ];
final_traj(2,:) = [traj_spline(2,:) vx(2:end)];
final_traj(3,:) = [traj_spline(3,:) acc_x(2:end)];
final_traj(4,:) = [traj_spline(4,:) jerk_x(2:end)];

figure
plot(new_time, traj_spline(1,:));
hold on
plot(new_time, traj_spline(2,:));
hold on
plot(new_time, traj_spline(3,:));
hold on
plot(new_time, traj_spline(4,:));

append_time = t_f + 0.1:0.1:t_f + 37.9 - 0.1;
whole_time = [new_time append_time];
figure
plot(whole_time, final_traj(1,:));
hold on
plot(whole_time, final_traj(2,:));
hold on
plot(whole_time, final_traj(3,:));
hold on
plot(whole_time, final_traj(4,:));

legend('pos', 'vel', 'acc', 'jerk')

save("traj_x", 'final_traj', 'whole_time')

function [sol_x, sol_y] = path

tf = 6;
ti = linspace(0, tf, 61);
solinit = bvpinit(ti, [0 0 0 0.5 0.5 0.5]);

opts= bvpset('Stats', 'on', 'RelTol', 1e-6);

sol_x = bvp4c(@ode, @bc_x, solinit, opts);
sol_y = bvp4c(@ode, @bc_y, solinit, opts);
% t = sol.x;
% y = sol.y;


    function dydx = ode(t,y,p)
        M = 9;
        u = -y(6);
        u = max(-M, u);
        u = min(M, u);
        dydx = [y(2);
            y(3);
            u;
            0;
            -y(4);
            -y(5)];
    end

    function res = bc_x(ya, yb, p)
        res = [ya(1);
            ya(2);
            ya(3);
            yb(1);
            yb(2) + 7;
            yb(3)];
    end
    function res = bc_y(ya, yb, p)
        res = [ya(1);
            ya(2);
            ya(3);
            yb(1) + 15;
            yb(2);
            yb(3) + 3.2666]
    end

%     function y = init(t, p)
%         alpha = -0.0513;
%         beta = 0.1128;
%         gamma = 0.0501;
%         y =[alpha * (t).^5 + beta * (t).^4 + gamma *(t).^3;
%             alpha * (t).^4 + beta *(t).^3 + gamma * (t).^2;
%             alpha * (t).^3 + beta *(t).^2 + gamma * (t);
%             -2 * alpha;
%             2 & alpha * t + 2 * beta;
%             -alpha * t.^2 -2 * beta*t -2*gamma];
%     end
end

