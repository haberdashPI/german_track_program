function E_stimuli(run_main,run_train)

    % ---------------------------------------------------------------------
    % running config ------------------------------------------------------
    % ---------------------------------------------------------------------
    
    if nargin==0
        run_main = 1;
        run_train = 1;
        runexample = 0;
    else
        runexample = 0;
    end
    
    % ---------------------------------------------------------------------
    % experiment/token setup ----------------------------------------------
    % ---------------------------------------------------------------------

    experiment_cfg = [];
    if 1
        book_name = 'passages1';
        readers = [10,21];
        [all_sentences,fs] = read_sentences(book_name,readers);
    end
    path = 'Experiments/E1/';
    rep_num = 1;
    normval = 30;
    tmr = 7;
    
    experiment_cfg.rep_num = rep_num;
    experiment_cfg.fs = fs;
    experiment_cfg.all_sentences = all_sentences;
    experiment_cfg.normval = normval;
    experiment_cfg.tmr = tmr;
    experiment_cfg.dev_len = 0.5; 
    
    test_block_cfg = [];
    test_block_cfg.dev_cases = [1 1; 2 1; 1 2; 2 2; -1 -1];
    test_block_cfg.dev_probs = [3; 2; 2; 1; 2]/10;
    test_block_cfg.num_trials = 50;  
    
    train_block_cfg = [];
    train_block_cfg.dev_cases = [1 1; -1 -1; 2 2; -1 -1; 1 2; 2 1; -1 -1; 1 1; 2 2; -1 -1];
    train_block_cfg.dev_probs = ones(size(train_block_cfg.dev_cases,1),1)/size(train_block_cfg.dev_cases,1);
    train_block_cfg.num_trials = size(train_block_cfg.dev_cases,1);  
    
    % ---------------------------------------------------------------------
    % running of the actual sections for stims ----------------------------
    % ---------------------------------------------------------------------
    
    if run_main
        experiment_cfg.block_cfg = test_block_cfg;
        get_all_exp_stuff(0);
        experiment_cfg.block_cfg.trial_sentences = trial_sentences;
        experiment_cfg.block_cfg.trial_dev_speakers = trial_dev_speakers;
        experiment_cfg.block_cfg.trial_dev_direction = trial_dev_direction;
        experiment_cfg.block_cfg.target_times = target_times;
        experiment_cfg.test_block_cfg = experiment_cfg.block_cfg;
    end
    
    if run_train
        experiment_cfg.block_cfg = train_block_cfg;
        get_all_exp_stuff(1);
        experiment_cfg.block_cfg.trial_sentences = trial_sentences;
        experiment_cfg.block_cfg.trial_dev_speakers = trial_dev_speakers;
        experiment_cfg.block_cfg.trial_dev_direction = trial_dev_direction;
        experiment_cfg.block_cfg.target_times = target_times;
        experiment_cfg.train_block_cfg = experiment_cfg.block_cfg;
    end
    
    if runexample
        experiment_cfg = make_all_stim(experiment_cfg,path);
    end
    
    if run_main==1 && run_train==1
        experiment_cfg = rmfield(experiment_cfg,'block_cfg');
        save([path 'experiment_record.mat']);
        save_for_experiment(experiment_cfg,path);
    end
    
    % ---------------------------------------------------------------------
    % nest functions ------------------------------------------------------
    % ---------------------------------------------------------------------
    
    function get_all_exp_stuff(is_training)
        sentence_perms = get_sentences_per_trial(experiment_cfg.all_sentences);
        [trial_sentences,trial_dev_speakers,trial_dev_direction] = get_trial_info(sentence_perms,experiment_cfg);
        target_times = zeros(size(trial_dev_direction));
        for trial_idx=1:size(trial_sentences,1)
            thispath = sprintf('trial_%d',trial_idx); 
            [stim,tt,h] = make_stim(experiment_cfg,trial_sentences(trial_idx,:),trial_dev_speakers(trial_idx),trial_dev_direction(trial_idx));
            target_times(trial_idx) = tt;
%             save_trial(stim,h,experiment_cfg.fs,experiment_cfg.normval,path,thispath,is_training);
        end
    end
    
end

function [all_sentences,fs] = read_sentences(book_name,readers) 

    all_sentences = cell(0,2);
    all_sentence_files = dir(['Experiments/E1/' book_name '/']);
    for file_idx=1:length(all_sentence_files)
        file_name = all_sentence_files(file_idx).name;
        if length(file_name)<3 || ~strcmp(file_name(end-2:end),'wav')
            continue;
        end
        [passage,fs] = audioread(['Experiments/E1/' book_name '/' file_name]);
        reader_id = str2double(file_name(strfind(file_name,'_dickens')-2:strfind(file_name,'_dickens')-1));
        if ~ismember(reader_id,readers)
            continue;
        end
        all_sentences_idx = find(cell2mat(all_sentences(:,2))==reader_id);
        if isempty(all_sentences_idx)
            all_sentences = [all_sentences; {[{passage} {length(passage)/fs}] reader_id}];
        else
            all_sentences{all_sentences_idx} = [all_sentences{all_sentences_idx}; {passage} {length(passage)/fs}];
        end
    end
    
end

function save_for_experiment(experiment_cfg,path)
%%
    block_cfg = experiment_cfg.train_block_cfg;
    bpath = 'stim_training/';
    save_this;
    
    block_cfg = experiment_cfg.test_block_cfg;
    bpath = '';
    save_this;    
    
    function save_this
        trial_dev_speakers = block_cfg.trial_dev_speakers;
        [~,~,ac] = unique(trial_dev_speakers,'stable');
        trial_dev_speakers(trial_dev_speakers>0) = ac(trial_dev_speakers>0);

        trial_dev_direction = block_cfg.trial_dev_direction;
        [~,~,ac] = unique(trial_dev_direction,'stable');
        trial_dev_direction = ac;
        trial_dev_direction(trial_dev_direction==3) = trial_dev_speakers(trial_dev_direction==3);

        sal = trial_dev_speakers>-1;
        this_info = [block_cfg.target_times sal trial_dev_speakers trial_dev_direction];
        dlmwrite([path bpath 'target_info_all.txt'],this_info,'delimiter',' ');

        sal = trial_dev_speakers==1;
        this_info = [block_cfg.target_times sal trial_dev_speakers trial_dev_direction];
        dlmwrite([path bpath 'target_info_obj.txt'],this_info,'delimiter',' ');

        sal = trial_dev_direction==1;
        this_info = [block_cfg.target_times sal trial_dev_speakers trial_dev_direction];
        dlmwrite([path bpath 'target_info_dir.txt'],this_info,'delimiter',' ');
    end

end

function sentence_perms = get_sentences_per_trial(all_sentences)
%%
    num_sentences_per_speaker = cellfun(@length,all_sentences(:,1));
    sentence_perms = zeros(prod(num_sentences_per_speaker),length(num_sentences_per_speaker));
    for speaker_idx=1:length(num_sentences_per_speaker)
        n = num_sentences_per_speaker(speaker_idx);
        if speaker_idx==length(num_sentences_per_speaker)
            num_rep_within = 1;
        else
            num_rep_within = prod(num_sentences_per_speaker(speaker_idx+1:end));
        end
        this_col = ones(num_rep_within,n);
        this_col = this_col .* repmat(1:n,num_rep_within,1);
        this_col = this_col(:);
        sentence_perms(:,speaker_idx) = repmat(this_col,size(sentence_perms,1)/length(this_col),1);
    end    

end

function [trial_sentences,trial_dev_speakers,trial_dev_direction] = get_trial_info(sentence_perms,experiment_cfg)
%%
    all_speakers = cell2mat(experiment_cfg.all_sentences(:,2))';
    dev_cases = experiment_cfg.block_cfg.dev_cases;
    dev_probs = experiment_cfg.block_cfg.dev_probs;
    num_trials = experiment_cfg.block_cfg.num_trials;    
    
    trial_sentences = [];
    perm_pointer = 1;
    for i=1:size(dev_cases,1)
        num_this = num_trials * dev_probs(i);
        trial_sentences = [trial_sentences; sentence_perms(mod(perm_pointer+(0:num_this-1)-1,size(sentence_perms,1))+1,:) repmat(dev_cases(i,:),num_this,1)];
        perm_pointer = mod(perm_pointer+num_this-1,size(sentence_perms,1))+1;
    end
    trial_dev_speakers = trial_sentences(:,3);
    trial_dev_speakers(trial_dev_speakers>0) = all_speakers(trial_dev_speakers(trial_dev_speakers>0));
    trial_dev_direction = categorical(trial_sentences(:,4),[-1 1 2],{'none','left','right'});
    trial_sentences = trial_sentences(:,1:2);

end

function [stim,target_time,h] = make_stim(experiment_cfg,sentence_idxs,dev_speaker,dev_direction)
%%
    % ---------------------------------------------------------------------
    % all config params ---------------------------------------------------
    % ---------------------------------------------------------------------
    
    all_sentences = experiment_cfg.all_sentences;
    fs = experiment_cfg.fs;
    tmr = experiment_cfg.tmr;
    dev_len = experiment_cfg.dev_len;
    
    % ---------------------------------------------------------------------
    % making of the stim (function calls mostly) --------------------------
    % ---------------------------------------------------------------------
    
    s1 = all_sentences{1}{sentence_idxs(1)};
    s2 = all_sentences{2}{sentence_idxs(2)};
    normalize_volume;
    len_stim = equalize_sentence_lengths;
    [t1,t2] = make_directionality;
    add_deviant;
    ear1 = s1.*direction_wave1 + s2.*direction_wave2;
    ear2 = s1.*(1-direction_wave1) + s2.*(1-direction_wave2);
    stim = [ear1, ear2];
    h=[];
%     h = show_stim(experiment_cfg,sentence_idxs,direction_wave1,direction_wave2,target_time,dev_speaker);
    
    % ---------------------------------------------------------------------
    % nest functions ------------------------------------------------------
    % ---------------------------------------------------------------------
    
    function normalize_volume
        s1 = s1/rms(s1);
        s2 = s2/rms(s2);
    end

    function len_stim = equalize_sentence_lengths
        len_stim = max(length(s1),length(s2));
        s1 = [s1; zeros(len_stim-length(s1),1)];
        s2 = [s2; zeros(len_stim-length(s2),1)];
    end

    function [t1,t2] = make_directionality
        % the switching section
        wave_time = len_stim/fs/6;
        wavelet = (1+sin(2*pi/(wave_time*2)*[0:1/fs:wave_time]-pi/2)')/2;
        % times of switches
        min_stay_len = 1; %seconds
        t1 = round(randsample(min_stay_len:0.1:len_stim/fs-min_stay_len*3-wave_time*2,1)*fs);
        t2 = round(randsample(t1/fs+wave_time+min_stay_len*2:0.1:len_stim/fs-min_stay_len-wave_time,1)*fs);
        % construction (1 = left ear)
        direction_wave1 = [zeros(t1,1); wavelet; ones(t2-t1-length(wavelet),1); flipud(wavelet); zeros(len_stim-t2-length(wavelet),1)];
        % second wave
        if 0
            t1 = round(randsample(min_stay_len:0.1:len_stim/fs-min_stay_len*3-wave_time*2,1)*fs);
            t2 = round(randsample(t1/fs+wave_time+min_stay_len:0.1:len_stim/fs-min_stay_len-wave_time,1)*fs);
            direction_wave2 = [ones(t1,1); flipud(wavelet); zeros(t2-t1-length(wavelet),1); wavelet; ones(len_stim-t2-length(wavelet),1)];
        else
            t1 = t1 + round(wave_time*fs/2);
            t2 = t2 - round(wave_time*fs/2);
            direction_wave2 = [ones(t1,1); flipud(wavelet); zeros(t2-t1-length(wavelet),1); wavelet; ones(len_stim-t2-length(wavelet),1)];
        end
    end

    function add_deviant
        if dev_direction=='right'
            target_t = randsample(find(direction_wave1(t1:end-fs)>0.8),1);
        elseif dev_direction=='left'
            target_t = randsample(find(direction_wave1(t1:end-fs)<0.2),1);
        else
            target_t = 0;
        end
        target_time = target_t/fs;
        if dev_speaker>0
            idxs = target_t:target_t+dev_len*fs;
            if dev_speaker==all_sentences{1,2}
                s1(idxs) = s1(idxs)*10^(tmr/20);
            elseif dev_speaker==all_sentences{2,2}
                s2(idxs) = s2(idxs)*10^(tmr/20);
            end
        end
    end

end

function save_trial(stim,handl,fs,normval,path,thispath,istraining)
%%
    stim = stim/normval;
    if ~istraining
        audio_path = [path 'stim_wav/' thispath '.wav'];
        img_path = [path 'stim_records/' thispath '.jpg'];
    else
        audio_path = [path 'stim_training/' thispath '.wav'];
        img_path = [path 'stim_training/' thispath '.jpg'];
    end
    audiowrite(audio_path,stim,fs);
    r = 150; % pixels per inch
    print(handl,'-dpng',sprintf('-r%d',r),img_path); 
    close;
end

function h = show_stim(experiment_cfg,sentence_idxs,direction_wave1,direction_wave2,target_time,dev_speaker)
%%
    all_sentences = experiment_cfg.all_sentences;
    fs = experiment_cfg.fs;
    colors = ['b','r','g'];
    
    h = figure; hold on
    
    t = (1:length(direction_wave1))/fs;
    l1 = plot(t,direction_wave1,'color',colors(1),'displayname',[num2str(all_sentences{1,2}) '-' num2str(sentence_idxs(1))]); 
    l2 = plot(t,direction_wave2,'color',colors(2),'displayname',[num2str(all_sentences{2,2}) '-' num2str(sentence_idxs(2))]);
    legend([l1 l2],'location','northwest');
    ylim([-0.5 1.5]);
    xlim([0 length(direction_wave1)/fs]);
    set(gca,'ytick',[0 1],'yticklabel',{'right','left'});
    xlabel('Time (s)');
    
    t = target_time:1/fs:target_time+experiment_cfg.dev_len;
    if dev_speaker==all_sentences{1,2}
        plot(t,direction_wave1(round(t(1)*fs):round(t(end)*fs)),'linewidth',10,'color',colors(1));
    elseif dev_speaker==all_sentences{2,2}
        plot(t,direction_wave2(round(t(1)*fs):round(t(end)*fs)),'linewidth',10,'color',colors(2));
    end

end














