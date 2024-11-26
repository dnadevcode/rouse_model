function [hFig,f] = rouse_gui(visible)

    % GUI for rouse model

    if nargin < 1
        visible = 'on';
    end
        
    %% Add source files to path and initialize settings
    mFilePath = mfilename('fullpath');
    mfolders = split(mFilePath, {'\', '/'});
    [fd, fe] = fileparts(mFilePath);
    % read settings txt
    setsTable  = readtable(fullfile(fd,'data','rouse_model_settings.txt'),'Format','%s%s%s');

    outputRes = []; % for selecting good/bad
    savePath = [];

    processFolders = 1; % whether to process single files or folders
    
    warning(''); % empty warning
    addpath(genpath(fullfile(mfolders{1:end - 1})));
    
    [~, lwid] = lastwarn;
    
    if strcmp(lwid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
        error('Unexpected error when asserting source folder path.')
    end

    %% Generate UI
    versionRM = importdata(fullfile(fd,'src','VERSION'));
        
    % create tabbed figure
    hFig = figure('Name', ['Rouse model GUI v',num2str(versionRM)], ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 0.3 0.8], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'figure', ...
        'Visible',visible ... 
    );

    hPanel = uipanel('Parent', hFig);
    h = uitabgroup('Parent',hPanel);
    t1 = uitab(h, 'title', 'Rouse model');
    tsHCC = uitabgroup('Parent',t1);
    hPanelImport = uitab(tsHCC, 'title', 'Settings tab');

    % Checklist as a loop
%     checkItems = [];%setsTable.Var2(14:20);
% 
%    % checkbox for things to plot and threshold
%     for i = 1:length(checkItems)
%         itemsList{i} = uicontrol('Parent', hPanelImport, 'Style', 'checkbox','Value', str2double(setsTable.Var1{13+i}),'String',{checkItems{i}},'Units', 'normal', 'Position', [0.45 .83-0.05*i 0.3 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
%     end
% 
%     itmsLst2 = [23 33 34];
%     checkItems2 = setsTable.Var2(itmsLst2 );
%     for i = 1:length(checkItems2)
%         itemsList2{i} = uicontrol('Parent', hPanelImport, 'Style', 'checkbox','Value', str2double(setsTable.Var1{itmsLst2(i)}),'String',{checkItems2{i}},'Units', 'normal', 'Position', [0.6 .83-0.05*i 0.3 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
%     end
%     
    % parameters with initial values
    textItems =  setsTable.Var2([1:8]);    
    values = setsTable.Var1([1:8]);

    for i=[1:8] % these will be in two columns
        positionsText{i} =   [0.5-0.5*mod(i,2) .93-0.1*ceil(i/2) 0.3 0.03];
        positionsBox{i} =   [0.5-0.5*mod(i,2) .88-0.1*ceil(i/2) 0.3 0.05];
    end
    

    for i=1:length(textItems)
        textListT{i} = uicontrol('Parent', hPanelImport, 'Style', 'text','String',{textItems{i}},'Units', 'normal', 'Position', positionsText{i},'HorizontalAlignment','Left');%, 'Max', Inf, 'Min', 0);  [left bottom width height]
        textList{i} = uicontrol('Parent', hPanelImport, 'Style', 'edit','String',{strip(values{i})},'Units', 'normal', 'Position', positionsBox{i});%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    end
   
     runButton = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Run'},'Callback',@run,'Units', 'normal', 'Position', [0.2 0.2 0.4 0.15]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    
    


    function run(src, event)
        display(['Started Rouse model simulation v',num2str(versionRM)]) %todo: make settings read-off automatic in order to make updating easier

        
        % Nb number of beads
        Nb =  str2double(textList{1}.String);
        
        zeta = eval(textList{2}.String{1}); % Friction coefficient of a bead, Ns/m
        ks =  eval(textList{3}.String{1}); % Spring constant N/m
        Xs =  eval(textList{4}.String{1}); % equilibrium rest length of the spring m
        kT = eval(textList{5}.String{1}); % J
        h=  eval(textList{6}.String{1}); % timestep
        
        numFrames =  str2double(textList{7}.String);
        frameLen =  str2double(textList{8}.String);       % time difference (s) between two frames
        fTime = numFrames*frameLen; % simulation stop time

        mex rouse.c;

        [sol] = rouse(Xs, Nb, zeta, ks, kT, h, fTime ,frameLen,numFrames)';

        f=figure;plot(sol)
        xlabel('time')
        ylabel('$x_i(t)$','Interpreter','latex')
        set(gca,'ytick',[])

        
        % Get the current time as a string
        currentTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
        
        % Create a filename with the current time
        filename = sprintf('rouse_trajectories%s.txt', currentTime);
        
        % Save the matrix to the file
        writematrix(sol, filename);
        
        % Display confirmation
        fprintf('Rouse bead trajectories saved to %s\n', filename);



    end
end


