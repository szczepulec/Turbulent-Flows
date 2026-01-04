
function [u_mean, v_mean, uu_rms, vv_rms, uv_rms, y] = compute_velocities_rms(Channel, U_0, base)


y = [Channel.y]'; 
h = Channel(1).h; 
nu= Channel(1).nu; 

y = [0; y; h];   

% velocities
% u
u_mean = zeros(23,1);

for i = 1:23
    u_mean(i) = sum(Channel(i).u .* Channel(i).tt) / sum(Channel(i).tt); % eqn 10.1 
    
    % u' = u - u_mean
    u_fluc{i} = Channel(i).u - u_mean(i);
end
u_mean = [0; u_mean; U_0]; 

% v
v_mean = zeros(23,1);

for i = 1:23
    v_mean(i) = sum(Channel(i).v .* Channel(i).tt) / sum(Channel(i).tt); % eqn 10.1 
    
    v_fluc{i} = Channel(i).v - v_mean(i);
end
v_mean = [0; v_mean; U_0]; 

% velocities_rms
% u'^2
uu_rms = zeros(23, 1);
for i = 1:23
    uu_rms(i) = (sum(u_fluc{i}.^2 .* Channel(i).tt) / sum(Channel(i).tt)).^(0.5); % eqn 10.2
end

% v'^2
vv_rms = zeros(23, 1);
for i = 1:23
    vv_rms(i) = (sum(v_fluc{i}.^2 .* Channel(i).tt) / sum(Channel(i).tt)).^(0.5); % eqn 10.2
end

% -uv'^2
uv_rms = zeros(23, 1);
for i = 1:23
    uv_rms(i) = (sum(-v_fluc{i}.*u_fluc{i} .* Channel(i).tt) / sum(Channel(i).tt)).^(0.5); % eqn 10.2
end %% Ask if this should be negative

end