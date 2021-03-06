clear

N = 10^6; % broj bita po simbolu
Eb_N0_dB = [0:40]; 
nTx = 2; %dvije predajne antene
nRx = 2; %dvije prijemne antene 
         %MIMO 2X2
         
d = 2;
pathLossEksp = 2;
K = 4;

%%MMSE - Rayleighev kanal
for i = 1:length(Eb_N0_dB)

    % Predajnik
    ip = rand(1,N)>0.5; % generi�u se 0 i 1 sa jednakom vjerovatno?om
    s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

    %kreiranje matrice [nRx,nTx,N/nTx]
    sMod = kron(s,ones(nRx,1));
    sMod = reshape(sMod,[nRx,nTx,N/nTx]);

   
    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; % Rayleigh-ev kanal
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % bijeli Gaussov �um, varijanca 0
   

    % Kanal i dodavanje �uma
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(i)/20)*n;

    % Prijemnik

    % Formiranje MMSE ekvalizatorske matrice W = inv(H^H*H+sigma^2*I)*H^H
    % H^H*H je dizmenzija [nTx x nTx], kod nas [2 x 2] 
    % Potrebno je za [2x2] matricu odrediti [a b; c d] = 1/(ad-bc)[d -b;-c a] 
    hCof = zeros(2,2,N/nTx)  ; 
    hCof(1,1,:) = sum(h(:,2,:).*conj(h(:,2,:)),1) + 10^(-Eb_N0_dB(i)/10);  
    hCof(2,2,:) = sum(h(:,1,:).*conj(h(:,1,:)),1) + 10^(-Eb_N0_dB(i)/10); 
    hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); 
    hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1);
    hDen = ((hCof(1,1,:).*hCof(2,2,:)) - (hCof(1,2,:).*hCof(2,1,:))); 
    hDen = reshape(kron(reshape(hDen,1,N/nTx),ones(2,2)),2,2,N/nTx);  
    hInv = hCof./hDen; % inv(H^H*H)

    hMod =  reshape(conj(h),nRx,N); % H^H 
    
    yMod = kron(y,ones(1,2)); % formatting prijemnih simbola za ekvalizaciju
    yMod = sum(hMod.*yMod,1); % H^H * y 
    yMod =  kron(reshape(yMod,2,N/nTx),ones(1,2)); % formatiranje
    yHat = sum(reshape(hInv,2,N).*yMod,1); % inv(H^H*H)*H^H*y
   
    % prijemnik - hard decision dekodiranje
    ipHat = real(yHat)>0;

    % brojanje gre�ki
    nErr(i) = size(find([ip- ipHat]),2);

end
simBer_MMSE_rejli = nErr/N; % simulacija BER za MMSE - Rayleigh-ev kanal

%%MMSE 
for i = 1:length(Eb_N0_dB)

    % Predajnik
    ip = rand(1,N)>0.5; % generi�u se 0 i 1 sa jednakom vjerovatno?om
    s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

    %kreiranje matrice [nRx,nTx,N/nTx]
    sMod = kron(s,ones(nRx,1));
    sMod = reshape(sMod,[nRx,nTx,N/nTx]);

    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx)*(1/K) + j*randn(nRx,nTx,N/nTx)]; 
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % bijeli Gaussov �um, varijanca 0
   % h = h*d^(-pathLossEksp/2)% Rice-ev kanal

    % Kanal i dodavanje �uma
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(i)/20)*n;

    % Prijemnik

    % Formiranje MMSE ekvalizatorske matrice W = inv(H^H*H+sigma^2*I)*H^H
    % H^H*H je dizmenzija [nTx x nTx], kod nas [2 x 2] 
    % Potrebno je za [2x2] matricu odrediti [a b; c d] = 1/(ad-bc)[d -b;-c a] 
    hCof = zeros(2,2,N/nTx)  ; 
    hCof(1,1,:) = sum(h(:,2,:).*conj(h(:,2,:)),1) + 10^(-Eb_N0_dB(i)/10);  
    hCof(2,2,:) = sum(h(:,1,:).*conj(h(:,1,:)),1) + 10^(-Eb_N0_dB(i)/10); 
    hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); 
    hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1);
    hDen = ((hCof(1,1,:).*hCof(2,2,:)) - (hCof(1,2,:).*hCof(2,1,:))); 
    hDen = reshape(kron(reshape(hDen,1,N/nTx),ones(2,2)),2,2,N/nTx);  
    hInv = hCof./hDen; % inv(H^H*H)

    hMod =  reshape(conj(h),nRx,N); % H^H 
    
    yMod = kron(y,ones(1,2)); % formatting prijemnih simbola za ekvalizaciju
    yMod = sum(hMod.*yMod,1); % H^H * y 
    yMod =  kron(reshape(yMod,2,N/nTx),ones(1,2)); % formatiranje
    yHat = sum(reshape(hInv,2,N).*yMod,1); % inv(H^H*H)*H^H*y
   
    % prijemnik - hard decision dekodiranje
    ipHat = real(yHat)>0;

    % brojanje gre�ki
    nErr(i) = size(find([ip- ipHat]),2);

end
simBer_MMSE_rice = nErr/N; % simulacija BER za MMSE - Riceov kanal

close all
figure
semilogy(Eb_N0_dB,simBer_MMSE_rejli,'ro-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer_MMSE_rice,'ko-','LineWidth',2);
axis([0 40 10^-4 0.5])
grid on
legend('(nTx=2,nRx=2, MMSE(Rayleigh-ev kanal))', '(nTx=2,nRx=2, MMSE(Rice-ov kanal))');
xlabel('Eb/No,dB');
ylabel('Bit Error Rate');
title('BER za BPSK modulaciju sa MMSE ekvalajzerom\newline u Rayleigh-evom i Rice-ovom kanalu');
