function states = db_to_dstates(signal_freq_bins, ...
                                noise_freq_bins, ...
                                resolution, ...
                                cont_vol, ...
                                cont_int)

    res = checkForHarmonics(signal_freq_bins, noise_freq_bins);
    if res~= false
        error(['Signal and Noise frequencies contain harmonics in noise_freq=',...
            num2str(res(1))])
    end

    % SET THIS PATH:
    path = 'C:\PPGEE\SBEB_CBA_24\ASSR - Coleta OFFLINE';
    
    %vetor dos voluntários 
    Vvoluntario = {'Abdon';'Ana';'BBB';'Colatina';'Erick';'Luciana';...
        'Sombra';'Quenaz';'Vinicius';'Sacola';'Wreikson'}; %vetor dos voluntário 
    
    %vetor da intensidade 
    Vintensidade = {'70';'60';'50';'40';'30'}; 
    load('eletrodos.mat')
    ganho  = 200;
    
    remoc = [.1]/ganho; 

    voluntario = cell2mat(Vvoluntario(cont_vol,:));

    if cont_int<0
         load([voluntario 'ESP'], 'x','Fs','binsM','freqEstim') 
    else
        intensidade = cell2mat(Vintensidade(cont_int,:));
        load([voluntario '_'  intensidade 'dB'], 'x','Fs','binsM','freqEstim') 
    end

    nfft = Fs;%1segundo de sinal 
     
    %retirar componente DC por janela (fiz isso pq no processamento em
    %tempo real é por janela)
    x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção
    


    fcInferior = 70;
    fcSuperior = 110;
    [b,a] = butter(2,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);

    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
    % freqz(d,[],Fs)
    x = filtfilt(d,x);

    d2 = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',fcInferior,'HalfPowerFrequency2',fcSuperior, ...
        'DesignMethod','butter','SampleRate',Fs);
    % freqz(d2,[],Fs)
    x = filtfilt(d2,x);
    
    % x = filter(b,a,x); 
    % excluir os dois primeiros segundos do inicio da coleta 
    x(:,1:2,:) =[]; 
    
    
    
    %encontrar o valor máximo por canal 
    Vmax = squeeze(max(max(abs(x)),[],3));
    ind = Vmax>remoc;
    
    xmedia = squeeze(mean(x(:,~ind,:),2));

      %   eletrodos =
      % 
      % 16×2 char array
      % 
      %   'FC'
      %   'F4'
      %   'T6'
      %   'P4'
      %   'T4'
      %   'Oz'
      %   'C4'
      %   'T5'
      %   'P3'
      %   'F7'
      %   'F3'
      %   'T3'
      %   'C3'
      %   'Fz'
      %   'Pz'
      %   'Cz'
    pos_eletrodo= 1;
    xmedia = x(:,:,pos_eletrodo);

    freq=0:(Fs-1); 

    xmedia = xmedia./std(xmedia);
    
    SIGNALS = fft(xmedia(:,:));%*2/nfft*1e9;
    FS = Fs;
    NFFT = nfft;
    SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
    f = FS/2*linspace(0,1,NFFT/2+1)'; % only half the FFT spectrum is valid
    max_length = size(SIGNALS,2);

    MSC = nan(numel(f),max_length);
    CSM=nan(numel(f),max_length);
    GFT=nan(numel(f),max_length);

    all_freqs = [signal_freq_bins noise_freq_bins];

    for idx_episodio = 1:max_length 
        
        % para cada frequencia em f
        for idx_f = 1:1:numel(f)
            
            
            M = idx_episodio;
            if M>40
                M=40;

                X_atual = SIGNALS(:,idx_episodio-M+1:idx_episodio);

            else

                X_atual = SIGNALS(:,1:idx_episodio);             
                
            end
            
            GFT(idx_f,idx_episodio) = sum(abs(X_atual(idx_f,1:M,:)).^2)./...
                        (sum(abs(X_atual(idx_f,1:M,:)).^2)+sum(abs(X_atual(numel(f),1:M,:)).^2));
    
            c1_csm = (sum(cos(angle(X_atual(idx_f,:))))./M).^2;
            c2_csm = (sum(sin(angle(X_atual(idx_f,:))))./M).^2;
            CSM(idx_f,idx_episodio) = c1_csm+c2_csm;

            num_msc = abs(sum(X_atual(idx_f,:)))^2;
            den_msc= M*sum(abs(X_atual(idx_f,:)).^2);
            MSC(idx_f,idx_episodio) = num_msc/den_msc;
        end
    end
    
    % cria matriz dos estados continuos:
    c_states = nan(numel(signal_freq_bins)+numel(noise_freq_bins), ... % freqs
                    4,... % states
                    max_length); % windows

    all_freqs = [signal_freq_bins,noise_freq_bins];
    
    c_states(:,1,:) = log10(CSM(all_freqs,:));
    c_states(:,2,:) = log10(GFT(all_freqs,:)); 
    c_states(:,3,:) = log10(MSC(all_freqs,:));
    c_states(:,4,:) = repmat(100.*[1:max_length]./max_length,numel(all_freqs),1);

    % discretiza a matriz de estados
    resolucao = resolution; % 10 niveis discretos para cada estado
    states =  nan(size(c_states));

    for idx_freq = 1:numel(all_freqs)
        for idx_state = 1:size(c_states,2)
            values = c_states(idx_freq,idx_state,:);

            min_val = min(c_states(:,idx_state,:),[],'all');
            max_val = max(c_states(:,idx_state,:),[],'all');
    
            states(idx_freq,idx_state,:) = discretize_val( ...
                                            values, min_val, max_val, resolucao);
        end
    end

  
end