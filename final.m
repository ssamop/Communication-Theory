clc
clear all
close all
a = 1;
b = 9;
t = [0.0:0.1:30];
x = a * sin(0.5*b*pi*t);

Ns = 2;

ts = t(1:Ns:end);
xs = x(1:Ns:end);

levels = [2,4,8,16,32];

meanError = zeros(size(levels));
xq = zeros(size(xs));
errors = zeros(size(xs));

pos = zeros(size(xs));
minimum= 1000;
diff=0;
eff= zeros(size(levels));
CR=zeros(size(levels));

for i = 1:length(levels)
    step = (max(x) - min(x)) / levels(i);
    
    %Axs = (min(xs):step:max(xs));
    B = (min(xs) + (step/2):step:max(xs) - (step/2));
   
    
    for k = 1:length(xs)
        minimum= 1000;
        for p = 1:length(B)
            diff= abs(xs(k)-B(p));
            if(diff<=minimum)
                minimum=diff;
                pos(k) = p-1; %encoded msg
                xq(k)= B(p);
            end
        end
    end
    
    errors = abs(xs - xq);
    meanError(i) = sum(errors) / length(xs);
    
    varpract(i) = var(abs(xs - xq));
    vartheo(i) = (step^2) / 12;
    SQNRprac(i) = max(x)^2 / varpract(i);
    SQNRthe(i) = 3 * levels(i)^2;

    prob=[];
    symbols = [];
    symbols = unique(pos);
    frequencies= histc(pos, symbols);
    prob=frequencies/length(pos);
    dict = huffmandict(symbols, prob);
    encode = huffmanenco(pos, dict);
    data=encode;
    decode = huffmandeco(data, dict); %source decoding
    
    dec=zeros(size(xs));
    for p=1:length(xs)
        dec(p)=B(decode(p)+1); %decoding (shabah el xq)
    end
    I= -log2(prob)
    H= sum(-prob.*log2(prob))
    V=ceil(I)
    L=sum(prob.*V)
    eff(i)=(H/L)
    CR(i)= (ceil(log2(length(xs))))/ levels(i);

    figure;
    plot(t, x, 'o-', 'DisplayName', 'Original Signal');
    hold on;
    plot(ts, dec, 's-', 'DisplayName', 'Decoded');
    title(['Decode - Level ', num2str(levels(i))]);
    legend('show');
    grid on;

end

% Plot the original signal
figure;
plot(t, x, 'DisplayName', 'Original');
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;

% Plot the sampled signal
figure;
plot(ts, xs, 'DisplayName', 'Sampled Signal');
title('Sampled Signal');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;

figure;
plot(levels, meanError, 'DisplayName', 'Mean');
title('Mean');
xlabel('Levels');
ylabel('Errors');
legend('show');
grid on;

figure;
plot(levels, varpract, 'o-', 'DisplayName', 'Variance Practical');
hold on;
plot(levels, vartheo, 's-', 'DisplayName', 'Variance Theo');
title('Variance Practical vs Theoretical');
xlabel('Number of Levels (L)');
ylabel('Variance');
legend('show');
grid on;

figure;
plot(levels, SQNRprac, 'o-', 'DisplayName', 'Practical SQNR');
hold on;
plot(levels, SQNRthe, 's-', 'DisplayName', 'Theoritical SQNR');
title('SQNR');
xlabel('Number of Levels (L)');
ylabel('SQNR');
legend('show');
grid on;

