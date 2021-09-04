clc
clear all
load('data_outdoor.mat')
[begin_x, begin_y, end_x, end_y] = path();

begin_time_x = begin_x.x;                                  %% time for x axis
begin_time_y = begin_y.x;                                  %% time for y axis

pre_traj_x = [];                                           %% guidance begin traj for x         
pre_traj_y = [];                                           %% guidance begin traj for y
pre_traj_x(1:3,:)= begin_x.y(1:3,:);
pre_traj_x(4,:) = -begin_x.y(6,:);
pre_traj_y(1:3, :) = begin_y.y(1:3, :);
pre_traj_y(4, :) = -begin_y.y(6, :);
% figure
% plot(begin_time_x, pre_traj_x(1,:), 'k-');xlabel('t')
% hold on
% plot(begin_time_x, pre_traj_x(2,:)); 
% hold on
% plot(begin_time_x, pre_traj_x(3,:)); 
% plot(begin_time_x, pre_traj_x(4,:));

t_f = 6;
new_time = 0:0.1:t_f;
traj_spline_x(1,:) = spline(begin_time_x, pre_traj_x(1,:),new_time);       %% fitting the guidance begin in x
traj_spline_x(2,:) = spline(begin_time_x, pre_traj_x(2,:),new_time);
traj_spline_x(3,:) = spline(begin_time_x, pre_traj_x(3,:),new_time);
traj_spline_x(4,:) = spline(begin_time_x, pre_traj_x(4,:),new_time);

traj_spline_y(1, :) = spline(begin_time_y, pre_traj_y(1, :), new_time);    %% fitting the guidance begin in y
traj_spline_y(2, :) = spline(begin_time_y, pre_traj_y(2, :), new_time);
traj_spline_y(3, :) = spline(begin_time_y, pre_traj_y(3, :), new_time);
traj_spline_y(4, :) = spline(begin_time_y, pre_traj_y(4, :), new_time);

final_traj_x(1,:) = [traj_spline_x(1,:) x(2:end) ];                        %% output guidance traj in x
final_traj_x(2,:) = [traj_spline_x(2,:) vx(2:end)];
final_traj_x(3,:) = [traj_spline_x(3,:) acc_x(2:end)];
final_traj_x(4,:) = [traj_spline_x(4,:) jerk_x(2:end)];

final_traj_y(1, :) = [traj_spline_y(1, :) y(2:end)];                        %% output guidance traj in y
final_traj_y(2, :) = [traj_spline_y(2, :) vy(2:end)];
final_traj_y(3, :) = [traj_spline_y(3, :) acc_y(2:end)];
final_traj_y(4, :) = [traj_spline_y(4, :) jerk_y(2:end)];

% figure
% plot(new_time, traj_spline_x(1,:));
% hold on
% plot(new_time, traj_spline_x(2,:));
% hold on
% plot(new_time, traj_spline_x(3,:));
% hold on
% plot(new_time, traj_spline_x(4,:));

circular_motion_end_time = 37.9;

append_time = t_f + 0.1:0.1:t_f + circular_motion_end_time - 0.1;    %% append time is the time in circular motion
whole_time = [new_time append_time];             %% guidance time + circular motion time + 
                                                 %       possible later tailing time
% figure
% plot(whole_time, final_traj_x(1,:));
% hold on
% plot(whole_time, final_traj_x(2,:));
% hold on
% plot(whole_time, final_traj_x(3,:));
% hold on
% plot(whole_time, final_traj_x(4,:));
% legend('pos', 'vel', 'acc', 'jerk')
% figure
% plot(final_traj_x(1,:), final_traj_y(1, :))
end_time_x = end_x.x;                                                       %% compute the tail traj
end_time_y = end_y.x;                                                       
tail_traj_x = [];
tail_traj_y = [];
tail_traj_x(1:3, :) = end_x.y(1:3, :);
tail_traj_x(4, :) = -end_x.y(6, :);
tail_traj_y(1:3, :) = end_y.y(1:3, :);
tail_traj_y(4, :) = end_y.y(6, :);

tail_traj_spline_x(1, :) = spline(end_time_x, tail_traj_x(1, :), new_time);%% fitting the tail traj for x and y
tail_traj_spline_x(2, :) = spline(end_time_x, tail_traj_x(2, :), new_time);
tail_traj_spline_x(3, :) = spline(end_time_x, tail_traj_x(3, :), new_time);
tail_traj_spline_x(4, :) = spline(end_time_x, tail_traj_x(4, :), new_time);

tail_traj_spline_y(1, :) = spline(end_time_y, tail_traj_y(1, :), new_time);
tail_traj_spline_y(2, :) = spline(end_time_y, tail_traj_y(2, :), new_time);
tail_traj_spline_y(3, :) = spline(end_time_y, tail_traj_y(3, :), new_time);
tail_traj_spline_y(4, :) = spline(end_time_y, tail_traj_y(4, :), new_time);


final_traj_x = [final_traj_x tail_traj_spline_x(:, 2:end)];                %% concantenation the whole traj points
final_traj_y = [final_traj_y tail_traj_spline_y(:, 2:end)];

tail_time = new_time + t_f + circular_motion_end_time;
whole_time = [whole_time tail_time(1:end - 1)];                                        %% concantenation the time
% figure
% plot(end_time_x, tail_traj_x(1, :));
% hold on
% plot(end_time_x, tail_traj_x(2, :));
% hold on
% plot(end_time_x, tail_traj_x(3, :));
% hold on
% plot(end_time_x, tail_traj_x(4, :));
% legend('pos_tail', 'vel_tail', 'acc_tail', 'jerk_tail');
figure
plot(whole_time, final_traj_x(1, :))
hold on
plot(whole_time, final_traj_x(2, :))
hold on
plot(whole_time, final_traj_x(3, :))
hold on
plot(whole_time, final_traj_x(4, :))
legend('pos_x', 'vel_x', 'acc_x', 'jerk_x')


figure
plot(final_traj_x(1,:), final_traj_y(1,:))


save("outdoor", 'final_traj_x', 'final_traj_y', 'whole_time')

function [sol_x, sol_y, end_sol_x, end_sol_y] = path

tf = 6;
ti = linspace(0, tf, 61);
solinit = bvpinit(ti, [0 0 0 0.5 0.5 0.5]);

opts= bvpset('Stats', 'on', 'RelTol', 1e-6);

sol_x = bvp4c(@ode, @bc_x, solinit, opts);
sol_y = bvp4c(@ode, @bc_y, solinit, opts);
end_sol_x = bvp4c(@ode, @end_bc_x, solinit, opts);
end_sol_y = bvp4c(@ode, @end_bc_y, solinit, opts);
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
            yb(1) - 15;
            yb(2);
            yb(3) + 3.2666]
    end


    function res = end_bc_x(ya, yb, p)
        res = [ya(1) - 14.031;
            ya(2) + 2.474;
            ya(3) + 3.0558;
            yb(1);
            yb(2);
            yb(3)];
    end
    function res = end_bc_y(ya, yb, p)
        res = [ya(1) - 5.3015;
            ya(2) - 6.5482;
            ya(3) + 1.1545;
            yb(1);
            yb(2);
            yb(3)]
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

