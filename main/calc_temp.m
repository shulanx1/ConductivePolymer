function temp = calc_temp(I_temp_down, dt_down)

if nargin < 2, dt_down = 0.02; end

Fs_down = 1/dt_down;
Wn = 10./(Fs_down /2);
[b,a] = butter(6, Wn, 'low');

V = filtfilt(b, a, I_temp_down)*10; % in mV

temp = (V-1250)/5;
end