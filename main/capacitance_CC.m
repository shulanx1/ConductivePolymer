function [Rm,Cm, Ihold, RMP, R_square] = capacitance_CC(base_path, date, cell_idx, idx_f, if_plot, if_save)
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 0; end


preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
filename = sprintf('%s.capacitance_CCstep.%d.wcp', preflix,idx_f);

try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.capacitance_CCstep.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end


n_channel = out.channel_no;
n_recording = length(out.rec_index);


dt = out.T(2)-out.T(1); % sample time in s

   
Vm = out.S{3};
Im = out.S{4};
t = 0:dt:dt*(size(Vm,1)-1);

Ihold = zeros(1, n_recording);
RMP = zeros(1, n_recording);

for i = 1:n_recording
    pulse_start = find(diff(Im(:,i))<-25);
    pulse_start = pulse_start(1);
    RMP(i) = mean(Vm(1:pulse_start,i));
    Ihold(i) = mean(Im(1:pulse_start,i));
    Vm(:,i) = Vm(:,i) - RMP(i);
end

Vm_mean = mean(Vm,2);
Istep = -50; %pA
Rin_est = mean(Vm_mean(pulse_start+floor(0.6/dt):pulse_start+floor(1/dt)))/Istep*1e3;
% ft = fittype('a1*(1-exp(-x/tau1))+a2*(1-exp(-x/tau2))', 'independent', 'x');
% V_temp = Vm(floor(0.5/dt):floor(0.65/dt), i_fit_range(i))-mean(Vm(1:floor(0.5/dt),i_fit_range(i)));
% t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [Istep_sub(i_fit_range(i))*Rin_pooled/1e3,Istep_sub(i_fit_range(i))*Rin_pooled/1e3, 5,5];
% opts.Lower = [min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), 0, 0];
% opts.Upper = [max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), 500,500];
% [curve,gof,output] = fit(t_temp', V_temp, ft, opts);
ft = fittype('a1*(1-exp(-x/tau1))+b', 'independent', 'x');
V_temp = Vm_mean(pulse_start:pulse_start+floor(0.15/dt));
t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
[curve,gof,output] = fit(t_temp', V_temp, ft, 'Start',[Istep*Rin_est/1e3,0, 5], 'Lower',[Istep*200/1e3, -5, 0],'Upper',[Istep/1e3, 0, 500]);  
R_square = gof.adjrsquare;
Rm = curve.a1/Istep*1e3;
Cm = curve.tau1/Rm*1e3;


if if_plot
figure,scatter(t_temp,V_temp)
hold on
plot(t_temp,curve(t_temp))
end

if if_save
    save(fullfile(path,sprintf('capacitance_CCstep_%d.mat', idx_f)),'Cm','Rm', 'Im', 'Vm','Ihold','RMP','Vm_mean');
end

RMP = mean(RMP);
Ihold = mean(Ihold);
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end