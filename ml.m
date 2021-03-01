clear

N = 10^6; % broj bita po simbolu
Eb_N0_dB = [0:40]; 
nTx = 2; %dvije predajne antene
nRx = 2; %dvije prijemne antene 
         %MIMO 2X2
         
K = 4; 

%%ML - rejli
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

simBer_ML_rejli = nErr/N; % simulacija BER za ML u Rayleigh-evom kanalu

%%ML-rajs
for z = 1:length(Eb_N0_dB)

     % Predajnik
    ip = rand(1,N)>0.5; % generišu se 0 i 1 sa jednakom vjerovatno?om
    s = 2*ip-1; % BPSK modulacija 0 -> -1; 1 -> 0

    %kreiranje matrice [nRx,nTx,N/nTx]
    sMod = kron(s,ones(nRx,1));
    sMod = reshape(sMod,[nRx,nTx,N/nTx]);

    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx)*(1/K) + j*randn(nRx,nTx,N/nTx)]; 
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

simBer_ML_rice = nErr/N; % simulacija BER za ML u Rice-ovom kanalu

close all
figure
semilogy(Eb_N0_dB,simBer_ML_rejli,'ro-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer_ML_rice,'ko-','LineWidth',2);
axis([0 30 10^-4 0.5])
grid on
legend('(nTx=2,nRx=2, ML(Rayleigh-ev kanal))', '(nTx=2,nRx=2, ML(Rice-ov kanal))');
xlabel('Eb/No,dB');
ylabel('Bit Error Rate');
title('BER za BPSK modulaciju sa ML ekvalajzerom\newline u Rayleigh-evom i Rice-ovom kanalu');

