clear all;
close all;

%% Constants
h = 6.62607004e-34;
q = 1.60217662e-19;
c = 299792458;

%% User inputs
%This parameter indicates whether the values of the parameters in the FB+ model should be plotted.
plot_param = 0; 

%This parameter expresses the degree with which the parameter should be
%fitted (1 = linear, 2 = quadratic, etc.).
deg_fitting = 1;

%The filenames of all measured samples and there bandgap energies
name = ["perov_1_557.csv" "perov_1_615.csv" "perov_1_665.csv"];
Eg = [1.557 1.615 1.665];

%This parameter indicates how many oscillators are used. 
%Typical values are 2 or 3, which can be observed by visual inspection. 
N_osc = 3;

%The bandgap energy that is wanted
Eg_wanted = 1.6:0.01:1.7;


%% Loading files from the folder data
Nsamples = length(Eg);
wav = 300:5:1200;
n_wav = zeros(length(wav),Nsamples);
k_wav = zeros(length(wav),Nsamples);

% The data is read from the folder
cd('Data')
for i = 1:Nsamples
    data = table2array(readtable(name(i)));
    n_wav(:,i) = interp1(data(:,1),data(:,2),wav);
    k_wav(:,i) = interp1(data(:,1),data(:,3),wav);
end
cd ..


%% Fitting with the FB+ method
%S, E_c, W, n_const, and E_min represent the parameters in the FB+ model.
%For each sample, these parameters are found by fitting the curves. For S,
%E_c, and W, multiple values are needed for each sample, as there are
%multiple oscilators in one material.
S = zeros(Nsamples,N_osc);
E_c = zeros(Nsamples,N_osc);
W = zeros(Nsamples,N_osc);
n_const = zeros(Nsamples,1);
E_min = zeros(Nsamples,1);

Error_k = zeros(Nsamples,1);
Error_n = zeros(Nsamples,1);

%This are the initial guesses for all the parameters.
S_0= [0.149,0.078,0.056];       S_0 = S_0(1:N_osc);
E_c_0 = [1.597,2.418,3.392];    E_c_0 = E_c_0(1:N_osc);
W_0 = [0.08,0.387,0.448];   W_0 = W_0(1:N_osc);
E_min_0 = Eg;

%An energy axis is made, such that the refractive index can be expressed as
%a function of energy.
E_int = 1e-21;
E = h*c/(wav(end)*1e-9):E_int:h*c/(wav(1)*1e-9);

%To obtain a higher accuracy in the simulated absorption, some wavelengths
%are more important to predict accurately than others. Therefore, a weight
%is given to the wavelength region 600 nm to 900 nm.
weight_wav = ones(1,length(wav));
[~,weight_ind_min] = min(abs(wav-600));
[~,weight_ind_max] = min(abs(wav-900));
weight_wav(weight_ind_min:weight_ind_max) = 4;
weight_E = interp1(wav*1e-9,weight_wav,h*c./E);

% The options for the fminsearch function
options = optimset('fzero');
options.MaxFunEvals = 1e7;
options.MaxIter = 1e7;
error = 1e-7;
options.TolFun = error;
options.TolX = error;

%For each sample, the parameters are fitted.
for i = 1:Nsamples
    k_E = interp1(wav*1e-9,k_wav(:,i),h*c./E);
    n_E = interp1(wav*1e-9,n_wav(:,i),h*c./E);
    x0 = [S_0, E_c_0,W_0,E_min_0(i)];
    fun_k = @(x)estimate_error_Buildin_k_inc_Eg(x,E,k_E,N_osc,weight_E);

    [x,~] = fminsearch(fun_k,x0,options);

    S(i,:) = abs(x(1:N_osc));
    E_c(i,:) = abs(x(N_osc+1:2*N_osc));
    W(i,:) = abs(x(2*N_osc+1:3*N_osc));
    E_min(i,:) = abs(x(end));

    k_E_fit = Make_k_E_curve(S(i,:),E_c(i,:),W(i,:),E_min(i,:),E);
    k_fit = interp1(E,k_E_fit,h*c./(wav*1e-9),'linear','extrap')';


    funmin = @(x)error_constant_n(wav,k_fit,n_wav(:,i),x);
    n_const(i) = fminsearch(funmin,1);

end


%% Show the fitted NK data
figure
for i = 1:Nsamples
    if Nsamples == 6 || Nsamples == 5
        subplot(3,2,i)
    else
        subplot(Nsamples,1,i)
    end
    yyaxis left
    hold on
    plot(wav,k_wav(:,i))
    k_E_fitted = Make_k_E_curve(S(i,:),E_c(i,:),W(i,:),E_min(i),E);
    k_wav_fitted = interp1(E,k_E_fitted,h*c./(wav*1e-9),"linear","extrap");
    plot(wav,k_wav_fitted);
    ylabel('k [-]')
    yyaxis right
    plot(wav,n_wav(:,i))
    n_wav_fitted = KK_integration(wav,k_wav_fitted)+n_const(i);
    plot(wav,n_wav_fitted);
    Error_k(i) = rms(k_wav(:,i)-k_wav_fitted');
    Error_n(i) = rms(n_wav(:,i)-n_wav_fitted');
    disp(append('Fitting error k for sample ',num2str(i),' is ',num2str(Error_k(i))));
    disp(append('Fitting error n for sample ',num2str(i),' is ',num2str(Error_n(i))));
    xlim([300, 1000])
    xlabel('Wavelength [nm]')
    title(append('E_g = ',num2str(Eg(i)),' eV'))
    legend('Actual','Fitted')
end
disp(append('Average fitting error k is ',num2str(mean(Error_k))));
disp(append('Average fitting error n is ',num2str(mean(Error_n))));


%% Make the fit of all parameters

fit_S = zeros(N_osc,deg_fitting+1);
fit_Ec = zeros(N_osc,deg_fitting+1);
fit_W = zeros(N_osc,2);

for i = 1:N_osc
    fit_S(i,:) = polyfit(Eg(:),S(:,i),deg_fitting);
    fit_Ec(i,:) = polyfit(Eg(:),E_c(:,i),deg_fitting);
    fit_W(i,:) = polyfit(Eg(:),W(:,i),deg_fitting);
end
fit_Emin = polyfit(Eg(:),E_min(:),deg_fitting);
fit_n_const = polyfit(Eg(:),n_const(:),deg_fitting);

%Plot parameters
if plot_param
    figure
    for i = 1:N_osc
        subplot(4,N_osc,i)
        hold on
        plot(Eg,E_c(:,i))
        plot(X_range,polyval(fit_Ec(i,:),X_range))
        xlabel(xlabel_str)
        title(append('E_j (',num2str(i),')'))
    end
    for i = 1:N_osc
        subplot(4,N_osc,i+N_osc)
        hold on
        plot(Eg,gamma(:,i))
        plot(X_range,polyval(fit_W(i,:),X_range))
        xlabel(xlabel_str)
        title(append('\gamma (',num2str(i),')'))
    end
    for i = 1:N_osc
        subplot(4,N_osc,i+2*N_osc)
        hold on
        plot(Eg,f(:,i))
        plot(X_range,polyval(fit_S(i,:),X_range))
        xlabel(xlabel_str)
        title(append('f (',num2str(i),')'))
    end
    subplot(4,N_osc,3*N_osc+1)
    hold on
    plot(Eg,n_const)
    plot(X_range,polyval(fit_n_const,X_range))
    xlabel('E_g [eV]')
    title('n_{const}')
    subplot(4,N_osc,3*N_osc+2)
    hold on
    plot(Eg,Eg)
    plot(X_range,polyval(fit_Emin,X_range))
    xlabel(xlabel_str)
    title('E_g optical [eV]')
end


%% Make new samples
n_wanted = zeros(length(wav),length(Eg_wanted));
k_wanted = zeros(length(wav),length(Eg_wanted));

color_begin = [1,0,0]; %Red
color_end = [0,0,1]; %Blue
color = [linspace(color_begin(1),color_end(1),length(Eg_wanted))',linspace(color_begin(2),color_end(2),length(Eg_wanted))',linspace(color_begin(3),color_end(3),length(Eg_wanted))'];

lgd_txt = cell(1,length(Eg_wanted));
fig = figure;
fig.Position = fig.Position + [-300, 0, 600,0];
sub1 = subplot(1,2,1); hold on; box on;
sub2 = subplot(1,2,2); hold on; box on;
for i = 1:length(Eg_wanted)
    S_insert = [polyval(fit_S(1,:),Eg_wanted(i));polyval(fit_S(2,:),Eg_wanted(i));polyval(fit_S(3,:),Eg_wanted(i))];
    E_c_insert = [polyval(fit_Ec(1,:),Eg_wanted(i));polyval(fit_Ec(2,:),Eg_wanted(i));polyval(fit_Ec(3,:),Eg_wanted(i))];
    W_insert = [polyval(fit_W(1,:),Eg_wanted(i));polyval(fit_W(2,:),Eg_wanted(i));polyval(fit_W(3,:),Eg_wanted(i))];
    E_min_insert = polyval(fit_Emin(1,:),Eg_wanted(i));
    n_c = polyval(fit_n_const(1,:),Eg_wanted(i));

    k_E_fit = Make_k_E_curve(S_insert,E_c_insert,W_insert,E_min_insert,E);
    k_wanted(:,i) = interp1(E,k_E_fit,h*c./(wav*1e-9),'linear','extrap')';
    n_wanted(:,i) = KK_integration(wav,k_wanted(:,i))+n_c;


    nm = wav'*1e-9;
    n = n_wanted(:,i);
    k = max(k_wanted(:,i),0);
    text = table(nm,n,k);
    writetable(text,append('perov_',num2str(Eg_wanted(i)*100),'.csv'),'FileType','text')

    subplot(1,2,1);
    plot(wav,n,'color',color(i,:),'LineWidth',2);

    subplot(1,2,2);
    plot(wav,k,'color',color(i,:),'LineWidth',2);
    lgd_txt{i} = append('E_g = ',num2str(Eg_wanted(i),3),' eV'); 
end

subplot(1,2,1);
xlim([300, 900]);
ylim([1.5, 3]);
xlabel('Wavelength [nm]')
ylabel('n [-]')
ax = gca;
ax.FontSize = 15;
subplot(1,2,2);
xlim([300, 900]);
ylim([0, 1.5]);
xlabel('Wavelength [nm]')
ylabel('k [-]')
ax = gca;
ax.FontSize = 15;
Lgd = legend(lgd_txt,'position',[0.47,0.2,0.05,0.7],'box','off','FontSize',15);
sub1.Position = sub1.Position + [-0.06,0.04,0,-0.04];
sub2.Position = sub2.Position + [0.06,0.04,0,-0.04];


%% Functions

function [error] = estimate_error_Buildin_k_inc_Eg(x,E,k_E,N_peaks,weight_E)
% This function estimates the error between the measured data and the fitted data
% Input:
%   x:
%       This vector contains the parameter that describe the k data
%       -S: represents the height of each peak
%       -E_c: represents the location of each peak
%       -W: represents the width of each peak
%       -E_min: represents the minimum energy
%   E:
%       This is the range of energies for which the values are given
%   k_E:
%       This is the k-data which is fitted. It is expressed as a function
%       of energy.
% Output:
%   error:
%       This is the deviation between the measured data and the actual data
%   Written by Y. Blom
S = abs(x(1:N_peaks));
E_c = abs(x(N_peaks+1:2*N_peaks));
W = abs(x(2*N_peaks+1:3*N_peaks));
E_min = abs(x(end));
[k_E_fit] = Make_k_E_curve(S,E_c,W,E_min,E);
error = rms((k_E-k_E_fit).*weight_E);
end

function [k_E_fit] = Make_k_E_curve(S,E_c,W,E_min,E)
% This function makes a fitted k_E curve based on the fitting parameters
% Input:
%   S:
%       represents the height of each peak
%   E_c:
%       represents the location of each peak
%   W:
%       represents the width of each peak
%   E_min:
%       is the minimum energy that has a nonzero k-values
%   E:
%       This is the range of energies for which the values are given
% Output:
%   k_E_fit:
%       This is the fitted data
%   Written by Y. Blom
q = 1.60217662e-19;
k_E_fit = zeros(1,length(E));
for i = 1:length(E)
    if E(i) > E_min*q
        k_E_fit(i) = sum(S.*(E(i)/q-E_min).^2./((E(i)/q-E_c).^2+W.^2));
    end
end
end

function [n_KK] = KK_integration(wav,k)
%This function applies the KK based on the second method described in:
% [1] K. Ohta et al, Applied Spectroscopy, 42(6), 952-957, (1988).
% Input:
%   wav:
%       The wavelength range of the input data
%   k:
%       The k values
% Output:
%   n_KK:
%       The n data after the KK transformation
% Written by Y. Blom

% Constants
h = 6.62607004e-34;
c = 299792458;

%The k data is transfered to a function of energy
E_int = 1e-21;
E = h*c/(wav(end)*1e-9):E_int:h*c/(wav(1)*1e-9);
k_E = interp1(wav*1e-9,k,h*c./E);
n_KK_E = zeros(size(k_E));

%A numerical integration is done for each wavelength/ energy
for i = 1:length(n_KK_E)
    %The integral is split into three parts as described by [1].
    j = 1:length(n_KK_E);
    j = j(j~=i);
    int1 = 2/pi*E_int*sum(E(j).*k_E(j)./(E(j).^2-E(i)^2));

    int2 = 1/pi*E_int*k_E(i)/(2*E(i));
    if i == 1
        int3 = 1/pi*E_int*(k_E(i+1)-k_E(i))/E_int;
    elseif i == length(n_KK_E)
        int3 = 1/pi*E_int*(k_E(i)-k_E(i-1))/E_int;
    else
        int3 = 1/pi*E_int*(k_E(i+1)-k_E(i-1))/(2*E_int);
    end
    n_KK_E(i) = max(1+int1+int2+int3,0);
end
% The data is transferred to a function of wavelength
n_KK = interp1(E,n_KK_E,h*c./(wav*1e-9),'linear','extrap');
end

function [error] = error_constant_n(wav,k,n,c)
% This function calculates the error between the integrated n and measured n
% for a certain additional constant c
% Input:
%   wav:
%       The wavelength range of the input data
%   k:
%       The k values
%   n:
%       The n values
%   c:
%       The constant that will be added to the n-data
% Output:
%   error:
%       The difference between the measured and simulated curve
% Written by Y. Blom
n_test = KK_integration(wav,k)';
error = rms(n_test+c-n);
end
