%clear all
ITap3=0;
y=0;
f=0;
f0=0;
f1=0;
f2=0;
f3=0;
x0=0;
zl=0;
x=0;
s=0;
outbit=0;
outvalue3TAP=0;


alpha=0.04; %db/cm/GHz%
v0=1;  %1st TAP
tau = (1/20); %nanosec%
l=10*2.54; %cm%
k= 0.0366*alpha*l/(tau);
zu=100;
dif = 0.001;
n = (zu - zl)/dif;
for x =1:1:20
f0=0;
f1=0;
zl=0;
s=0;
f1 = (2*v0/3.14)*exp(-k*zl)*cos((2*x-2*x0-1)*zl);
y(1)=f1;
for i=2:n
zl = zl+dif;
f1 = (2*v0/3.14)*((sin(zl))/zl)*exp(-k*zl)*cos((2*x-2*x0-1)*zl);
f=f1;
y(i)=f;
s =(dif/2)*(y(i)+y(i-1))+s;
end
JTap1(x)= s;
end
JTap1;

Inbit=0;
n=10000; %no. of BIT
Inbit = randi([0,1],[1,n]); %generating random Bit
inbit = Inbit ;
if inbit(1)==0
    inbit(1)=1;
else
    inbit(1)=1;
end
inbit;


inbit;
sm=0;
sm1=0;
k= 1; %*floor(tau/0.01);
j=0;
i=0;
p=length(JTap1)/k;
for i=1:p
    y(i)= JTap1(1+(i-1).*k);
    sm(i)=JTap1(1+(i-1).*k);
end
y;
sm;
for j=2:(length(inbit)+1)
    if inbit(j-1)==1
    for i=2:p+1
   while (i+(j-2))>length(sm)
            sm(i+(j-2))=0;
   end   
        sm1(i+(j-2))= y(i-1)+sm(i+(j-2));
    end 
    for i=1:j-1
      sm1(i)= sm(i);
    end
    else 
        for i=2:p+1
   while (i+(j-2))>length(sm)
            sm(i+(j-2))=0;
   end   
        sm1(i+(j-2))=sm(i+(j-2));
   end 
    for i=1:j-1
      sm1(i)= sm(i);
    end
    end
    sm=sm1;
    outvalue(j-1)=sm1(j);
end
outvalue;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **ISI-Aware Data Reconstruction Algorithm**
oneval = [];
zeroval = [];
mo = 1;
no = 1;

outbit = zeros(1, length(outvalue)); % Preallocate Output Bits
JTap1_extended = zeros(1, length(outvalue)); % Initialize with zeros
JTap1_extended(1:20) = JTap1; % Keep the first 20 values unchanged
JTap1 = JTap1_extended; % Replace original JTap1

for j = 1:length(outvalue)
    % Compute ISI using only previous logic '1' values (JTap1)
    if j > 1
        ISI_j = sum(JTap1(2:j)); % ISI contribution from previous bits
    else
        ISI_j = 0; % No ISI for the first bit
    end

    % Decision Rule
    if (outvalue(j) - ISI_j) <= 0
        outbit(j) = 0;
        zeroval(mo) = outvalue(j);
        mo = mo + 1;
    else
        outbit(j) = 1;
        oneval(no) = outvalue(j);
        no = no + 1;
    end
end

% **Compute BER Using Probability Distribution**
miu1 = mean(oneval);
sigma1 = std(oneval);
miu0 = mean(zeroval);
sigma0 = std(zeroval);

x = -0.5:0.01:1.5;
y1 = @(x) exp(-((miu1 - x).^2) / (2 * sigma1^2)) / sqrt(2 * pi * sigma1^2);
y2 = @(x) exp(-((miu0 - x).^2) / (2 * sigma0^2)) / sqrt(2 * pi * sigma0^2);
yf = @(x) y1(x) - y2(x);

x0 = 0.1;
xint = fzero(yf, x0);

% Integration for BER
v1L = -2.9;
v1u = xint;
dif1 = 0.001;
n = (v1u - v1L) / dif1;

fn = exp(-((miu1 - v1L)^2) / (2 * sigma1^2)) / sqrt(2 * pi * sigma1^2);
y(1) = fn;
s1 = 0;

for i = 2:n
    v1L = v1L + dif1;
    fn = exp(-((miu1 - v1L)^2) / (2 * sigma1^2)) / sqrt(2 * pi * sigma1^2);
    y(i) = fn;
    s1 = (dif1 / 2) * (y(i) + y(i-1)) + s1;
end

I1 = s1;

v0L = xint;
v0u = 2.5;
dif0 = 0.001;
n = (v0u - v0L) / dif0;

fn = exp(-((miu0 - v0L)^2) / (2 * sigma0^2)) / sqrt(2 * pi * sigma0^2);
y(1) = fn;
s0 = 0;

for i = 2:n
    v0L = v0L + dif0;
    fn = exp(-((miu0 - v0L)^2) / (2 * sigma0^2)) / sqrt(2 * pi * sigma0^2);
    y(i) = fn;
    s0 = (dif0 / 2) * (y(i) + y(i-1)) + s0;
end

I0 = s0;
probBer = I1 + I0;

% **Compute SNR**
Signalval = miu1 - miu0;
Noiseval = sqrt(sigma1^2 + sigma0^2);
SNRval = Signalval / Noiseval;

% **BER Calculation by Bit Errors**
noe = sum(Inbit ~= outbit);
BER = noe / length(Inbit);

% **Final Output Values**
fnloutptxl = [miu1, sigma1, miu0, sigma0, xint, max(zeroval), min(oneval), sm(1), probBer, SNRval];




