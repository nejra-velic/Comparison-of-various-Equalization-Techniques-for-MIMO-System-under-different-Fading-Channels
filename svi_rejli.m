clear

N = 10^6; % broj bita po simbolu
Eb_N0_dB = [0:40]; 
nTx = 2; %dvije predajne antene
nRx = 2; %dvije prijemne antene 
         %MIMO 2X2
         
%%ML
for z = 1:length(Eb_N0_dB)

     % Predajnik
    ip = rand(1,N)>0.5; % generišu se 0 i 1 sa jednakom vjerovatno?om
    s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

    %kreiranje matrice [nRx,nTx,N/nTx]
    sMod = kron(s,ones(nRx,1));
    sMod = reshape(sMod,[nRx,nTx,N/nTx]);

    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; % Rayleigh-ev kanal
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % bijeli Gaussov šum, varijanca 0

    % Kanal i dodavanje šuma
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(z)/20)*n;

    % Maximum Likelihood prijemnik
    % ----------------------------
    % [s1 s2 ] = [+1,+1 ]
    sHat1 = [1 1];	
    sHat1 = repmat(sHat1,[1 ,N/2]);
    sHat1Mod = kron(sHat1,ones(nRx,1));	
    sHat1Mod = reshape(sHat1Mod,[nRx,nTx,N/nTx]);	
    zHat1 = squeeze(sum(h.*sHat1Mod,2)) ;
    J11 = sum(abs(y - zHat1),1);
    
    % [s1 s2 ] = [+1,-1 ]
    sHat2 = [1 -1];	
    sHat2 = repmat(sHat2,[1 ,N/2]);
    sHat2Mod = kron(sHat2,ones(nRx,1));	
    sHat2Mod = reshape(sHat2Mod,[nRx,nTx,N/nTx]);	
    zHat2 = squeeze(sum(h.*sHat2Mod,2)) ;
    J10 = sum(abs(y - zHat2),1);
    
    % [s1 s2 ] = [-1,+1 ]
    sHat3 = [-1 1];	
    sHat3 = repmat(sHat3,[1 ,N/2]);
    sHat3Mod = kron(sHat3,ones(nRx,1));	
    sHat3Mod = reshape(sHat3Mod,[nRx,nTx,N/nTx]);	
    zHat3 = squeeze(sum(h.*sHat3Mod,2)) ;
    J01 = sum(abs(y - zHat3),1);
    
    % [s1 s2 ] = [-1,-1 ]
    sHat4 = [-1 -1];	
    sHat4 = repmat(sHat4,[1 ,N/2]);
    sHat4Mod = kron(sHat4,ones(nRx,1));	
    sHat4Mod = reshape(sHat4Mod,[nRx,nTx,N/nTx]);	
    zHat4 = squeeze(sum(h.*sHat4Mod,2)) ;
    J00 = sum(abs(y - zHat4),1);
    
    % pronalaženje minimuma od svih kombinacija 
    rVec = [J11;J10;J01;J00];
    [jj dd] = min(rVec,[],1);

    % mapiranje u bite
    ref = [1 1; 1 0; 0 1; 0 0 ];
    ipHat = zeros(1,N);
    ipHat(1:2:end) = ref(dd,1);
    ipHat(2:2:end) = ref(dd,2);

    % brojanje greški
    nErr(z) = size(find([ip- ipHat]),2);

end

simBer_ML = nErr/N; % simulacija BER za ML 
         
%%MRC - nelinearni ekvalajzer, treba dati najbolje rezultate

    for m = 1:length(Eb_N0_dB)
        
        % Predajnik
        ip = rand(1,N)>0.5; % generišu se 0 i 1 sa jednakom vjerovatno?om
        s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

        n = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % bijeli Gaussov šum, varijanca 0
        h = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % Rayleigh-ev kanal
        
        % Kanal i dodavanje šuma
        sD = kron(ones(nRx,1),s);
        y = h.*sD + 10^(-Eb_N0_dB(m)/20)*n;

        % ekvalizacija MRC-om 
        yHat =  sum(conj(h).*y,1)./sum(h.*conj(h),1); 

        % prijemnik - hard decision dekodiranje
        ipHat = real(yHat)>0;

        % brojanje greški
        nErr(m) = size(find([ip- ipHat]),2);

    end

simBer_MRC = nErr/N; % simulacija BER za MRC

%%MMSE 
for i = 1:length(Eb_N0_dB)

    % Predajnik
    ip = rand(1,N)>0.5; % generišu se 0 i 1 sa jednakom vjerovatno?om
    s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

    %kreiranje matrice [nRx,nTx,N/nTx]
    sMod = kron(s,ones(nRx,1));
    sMod = reshape(sMod,[nRx,nTx,N/nTx]);

    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; % Rayleigh-ev kanal
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % bijeli Gaussov šum, varijanca 0
   

    % Kanal i dodavanje šuma
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

    % brojanje greški
    nErr(i) = size(find([ip- ipHat]),2);

end
simBer_MMSE = nErr/N; % simulacija BER za MMSE

%%ZF
for j = 1:length(Eb_N0_dB)

    % Predajnik
    ip = rand(1,N)>0.5; % generišu se 0 i 1 sa jednakom vjerovatno?om
    s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

    %kreiranje matrice [nRx,nTx,N/nTx]
    sMod = kron(s,ones(nRx,1));
    sMod = reshape(sMod,[nRx,nTx,N/nTx]);

    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; % Rayleigh-ev kanal
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % bijeli Gaussov šum, varijanca 0

    % Kanal i dodavanje šuma
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(j)/20)*n;

    % Prijemnik

    % Formiranje ZF ekvalizacijske matrice W = inv(H^H*H)*H^H
    % H^H*H je dimenzija [nTx x nTx], kod nas [2 x 2] 
    % inverzna matrica [a b; c d] = 1/(ad-bc)[d -b;-c a]
    hCof = zeros(2,2,N/nTx)  ; 
    hCof(1,1,:) = sum(h(:,2,:).*conj(h(:,2,:)),1) ;  
    hCof(2,2,:) = sum(h(:,1,:).*conj(h(:,1,:)),1) ;  
    hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); 
    hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1); 
    hDen = ((hCof(1,1,:).*hCof(2,2,:)) - (hCof(1,2,:).*hCof(2,1,:))); 
    hDen = reshape(kron(reshape(hDen,1,N/nTx),ones(2,2)),2,2,N/nTx);  
    hInv = hCof./hDen; % inv(H^H*H)

    hMod =  reshape(conj(h),nRx,N); % H^H 
    
    yMod = kron(y,ones(1,2)); % primljeni simboli za ekvalizaciju
    yMod = sum(hMod.*yMod,1); % H^H * y 
    yMod =  kron(reshape(yMod,2,N/nTx),ones(1,2)); 
    yHat = sum(reshape(hInv,2,N).*yMod,1); % inv(H^H*H)*H^H*y
    
    % prijemnik - hard decision dekodiranje
    ipHat = real(yHat)>0;

    % brojanje greški
    nErr(j) = size(find([ip- ipHat]),2);

end

simBer_ZF = nErr/N; % simulacija BER za ZF


close all
figure
semilogy(Eb_N0_dB,simBer_MMSE,'bo-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer_ZF,'mo-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer_ML,'ro-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer_MRC,'ko-','LineWidth',2);
axis([0 40 10^-4 0.5])
grid on
legend('(nTx=2,nRx=2, MMSE)', '(nTx=2,nRx=2, ZF)', '(nTx=2, nRx=2, ML)', '(nRx=2, MRC)');
xlabel('Eb/No,dB');
ylabel('Bit Error Rate');
title('BER za BPSK modulaciju sa MMSE, ZF, ML i MRC ekvalajzerima\newline u Rayleigh-evom kanalu');