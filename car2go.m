function car2go( varargin)
% car2go visualize and analyse
clear all
global  CarsData haxMap hTblCarData haxGraph Year0 car1Day1 ...
        hpupTable hpupGraph hpupGrType LastGraph LastTable LastGrType ...
        histData histTitle xTitle

[LastGraph, LastTable, LastGrType] = deal(NaN);    
Year0 = 2014;
figName = 'car2go';
singleton = true;
% Note: if non-singleton then take care about 1) Global vars 2) timer 
%


% ========================================================================
% ====    Keep Single     ================================================
% ========================================================================
if singleton
    allGUIs = findall(0,'name',figName);
    isGUIopen = allGUIs > 0;
    if numel(isGUIopen) > 0
        lastGUI = allGUIs(end); 
        figure(lastGUI)
        choice = questdlg('Only one instance','', ...
            'Stay in Existing','Close Existing and Open New','Stay in Existing');
        % Handle response
        switch choice
            case 'Stay in Existing'
                disp([mfilename,'=>',choice ])
                return;
            case 'Close Existing and Open New'
                disp([mfilename,'=>',choice])
                delete(allGUIs)
        end

    end
end

% ========================================================================
% ====    Size and Colours     ===========================================
% ========================================================================
colors = bone(20);
bgc = colors(12,:);
colors = repmat(colors(8,:),10,1);
% colors = colors(8:end,:);
% bgc = colors(4,:);

set(0, 'unit','pixels')
% Default units
tmp = get(0,'screensize');
if tmp(3) > 1400
    defaultFontsize = 8;
    figW = 0.36;
else
    defaultFontsize = 7;
    figW = 0.60;
end

% ========================================================================
% ====    MAIN FIGURE     ================================================
% ========================================================================

set(0, 'unit','normalized')
hF = figure(...    'numbertitle', 'off',...
    'WindowStyle','normal',...
    'name', figName,...
    'units', 'normalized',...
    'color', bgc,...
    'position',[0.02 0.2 figW 0.75],... %ceil(get(0,'screensize') .* [1 1 0.975 0.65]),...
    'menubar', 'none','toolbar','none','visible','on');

set(hF,'DefaultUicontrolUnits','normalized',...
    'DefaultUicontrolFontSize',defaultFontsize);

% ========================================================================
% ====    TAB PANEL     ==================================================
% ========================================================================
MainPnlNames = {'Table','Map','Graph'};

set(0, 'unit','pixels')
[~,hMainPnl,~] = tabPanel(hF,MainPnlNames,...
    'panelpos',[0.02 0.05 0.95 0.80],...
    'tabpos','Top',...
    'tabHeight',60,...
    'colors',colors,...
    'highlightColor','w',...'c'
    'tabCardPVs',{'bordertype','etchedin','fontsize',defaultFontsize},...
    'tabLabelPVs',{'fontsize',11,'Rotation',0});


% ========================================================================
% ====  Get Config Params  ===============================================
% ========================================================================

if nargin == 0
    [InputData] = ioCarCfgF('r');
else
    CarCfgFname = varargin{1};
    [InputData] = ioCarCfgF('r',CarCfgFname);
end

currCarInpFileNo = 1;            
currCarInpF = InputData{currCarInpFileNo,2};            

% Sampling time
% P.D. 21/04/2014: Think to get it from the data file
SmplTimeStr = InputData{currCarInpFileNo,3};
SmplTime = (1/60)*str2double(SmplTimeStr);

ScaleLat = InputData{currCarInpFileNo,6}{1};
ScaleLon = InputData{currCarInpFileNo,6}{2};
TankPerKm = 1/550;

% ========================================================================
% ====    Info Table  ===================================================
% ========================================================================
htableInpF = uitable(hF,'Tag','tableInpF',...
    'units','normalized','pos',[0.02 0.90 0.95 0.085],'enable','on');

isVarEdit = false(1,4);

CityStr = InputData{currCarInpFileNo,1};  
PeriodStr = 'Not defined yet';

inpFLegend = {'City','input_Fname','period','dT'};
InfoStr = {CityStr, currCarInpF, PeriodStr, SmplTimeStr};
ColW_inpData = colWidth(htableInpF, InfoStr);

set(htableInpF,'ColumnEditable',isVarEdit,...
    'ColumnName',inpFLegend,'ColumnWidth', ColW_inpData, ...
    'Data', InfoStr,'FontSize',11);
% ========================================================================
% ====    Tab Panels Content     =========================================
% ========================================================================

for iTab = 1:length(MainPnlNames)
    setupPanel(MainPnlNames{iTab},iTab);
end

% ========================================================================
% ====    Menu     =======================================================
% ========================================================================
MenuFile =  uimenu(hF,'Label','Files');
            uimenu(MenuFile,'Label','CfgFile ...','Callback',{@getCfgFile});
            uimenu(MenuFile,'Label','Select Input File','Callback',{@CBinpFile});


% ========================================================================
% ====  Get Car Data    ==================================================
% ========================================================================
isHr = true;
CarsData = getCarData(currCarInpF,isHr);


% ========================================================================
% ====    Update Table  ===================================================
% ========================================================================
PeriodStr = CarsData.Period;
InfoStr = get(htableInpF, 'data');
InfoStr{3} = PeriodStr;
ColW_inpData = colWidth(htableInpF, InfoStr);

set(htableInpF,'ColumnWidth', ColW_inpData,'Data', InfoStr);

% ========================================================================
% ====    Display     ====================================================
% ========================================================================
currGraphNo = 1;
currTableNo = 1;
UpdatePanel(currGraphNo, currTableNo);

% ========================================================================
% ========================================================================
% ====   NESTED  FUNCTIONS     ===========================================
% ========================================================================
% ========================================================================

    function CarsData = getCarData(carInpF,isHr)
    %   
    %  formT = {'hr' 'full'}
        
        fid = fopen(carInpF,'r');
        nums = fread(fid);
        strCfg = char(nums');
        fclose(fid);
        
        CarsData = struct; %struct([]);
        
    % Get column names from the first line 
        tmpStr = textscan(strCfg,'%s',1, 'delimiter', '\n');
        Legends = regexp(tmpStr{1},',','split');
        CarsData.FileLegend = strtrim(Legends{1});
    
        
    % --- Collect data by car's license No  -----------------------------
    % ---- A) Reg Expr path -------------------
        tic
        if isHr
            expr_data = ...
                ['"(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):\d\d"\,'...   % capture (Day+hr+min) skip seconds
                ' *(\w{3,6})\,'...     % license
                ' *(-?\d{1,3}\.\d{0,6})\, *(-?\d{1,3}\.\d{0,6})\,'... %lat, lon
                ' *(\d{1,3})\, *([\w]{1,6})\, *([\w]{1,6})\s+']; % fuel, ext, int
            
            data = regexpi(strCfg,expr_data,'tokens');
            data = data';
            %CarsData.RawData = ;
            CarsLic = cellfun(@(x) x{6}, data,'uni',false); % List of all logged Licenses
            CarsY = cellfun(@(x) x{1}, data,'uni',false);
            CarsM = cellfun(@(x) x{2}, data,'uni',false);
            CarsD = cellfun(@(x) x{3}, data,'uni',false);
            CarsHr = str2double(cellfun(@(x) x{4}, data,'uni',false));
            CarsMin = str2double(cellfun(@(x) x{5}, data,'uni',false));
            Date1str   = [CarsY{1},'/',CarsM{1},'/',CarsD{1}];
            DateEndstr = [CarsY{end},'/',CarsM{end},'/',CarsD{end}];
            PeriodStr = [Date1str,'--',DateEndstr];
            CarsData.Period = PeriodStr;
            
            CarsD = str2double(CarsD);
            CarsDay1 = CarsD(1);
            TimeLine = (CarsD-CarsDay1)*24 + CarsHr + CarsMin/60;
        else
            %CarsTime = datenum(CarsTime,'yyyy-mm-dd HH:MM');
        end
        
        toc

        car1Day1 = CarsD(1);        
        DayLine = unique(CarsD);
        numelD = numel(DayLine);
        
        CarsLicUniq =  unique(CarsLic);
        CarsData.CarsLic = CarsLicUniq;

        TimeLine = unique(TimeLine);
        CarsData.TimeLine = TimeLine;

        numelT = numel(TimeLine);
        numelCars = length(CarsLicUniq);
        
        MatTOC = {'Lat', 'Lon', 'Fuel'}; 
        CarsData.MatTOC = MatTOC;

        [CarsLocInd, CarsGoDist, CarsGoTime, DataStat,...
            CarsGoDistCum, CarsGoTimeCum, CarsFuelCum] = deal(NaN(numelT,numelCars));
        
        [CarsGoDistDaily, CarsGoTimeDaily] = deal(NaN(numelCars, numelD));
        
        iLic = 0;
        while ~isempty(CarsLic)
            iLic = iLic +1;
            Lic1 = CarsLicUniq{iLic};  % First license in the list
            indxLic = strcmp(Lic1,CarsLic); 
            car1Data = data(indxLic);

            % Stat<int8> characterizes data quality: 
            % 1 - High, 2 - Average, 3 - Poor
            [carTime, carDayLine, GoDist, GoTime, LocInd, stat,FuelCum, ...
                GoDistCum, GoTimeCum, GoDistDaily, GoTimeDaily] = CarState(car1Data);
            
            [~, ind1_inAll] = ismember(carTime, TimeLine);
            
            for i = 1:numel(carTime)
                CarsLocInd(ind1_inAll(i),iLic) = LocInd(i);
                CarsGoDist(ind1_inAll(i),iLic) = GoDist(i);
                CarsGoTime(ind1_inAll(i),iLic) = GoTime(i);
                DataStat(ind1_inAll(i),iLic) = stat(i);
                CarsGoDistCum(ind1_inAll(i),iLic) = GoDistCum(i);
                CarsGoTimeCum(ind1_inAll(i),iLic) = GoTimeCum(i);
                CarsFuelCum(ind1_inAll(i),iLic) = FuelCum(i);
            end
            
            [~, indDays] = ismember(carDayLine, DayLine);
            for i = 1:numel(carDayLine)
                CarsGoDistDaily(iLic,indDays(i)) =  GoDistDaily(i);
                CarsGoTimeDaily(iLic,indDays(i)) =  GoTimeDaily(i);
            end
            
            % Removes all lines with data for the Marked sensor
            data(indxLic) = [];
            CarsLic(indxLic) = [];
            
        end % While
        
        toc
        CarsData.CarsLocInd = CarsLocInd;
        CarsData.CarsGoDist = CarsGoDist;
        CarsData.CarsGoTime = CarsGoTime;
        CarsData.CarsGoDistCum = CarsGoDistCum;
        CarsData.CarsGoTimeCum = CarsGoTimeCum;
        CarsData.CarsFuelCum = CarsFuelCum;
        CarsData.DataStat = DataStat;
        CarsData.CarsGoDistDaily = CarsGoDistDaily;
        CarsData.CarsGoTimeDaily = CarsGoTimeDaily;
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [carTime, carDayLine, GoDist, GoTime, LocInd, stat, carFuel, ...
            GoDistCumT, GoTimeCumT, GoDistDaily, GoTimeDaily] = CarState(carData)
    %
    % stat digit value:     1   2   3   4   5   6   7   8   9
    %  1-pos [dT,min]      <2, <4, <8, <16,                >16 
    %  2-pos [Vel=Dist/dT] <2, <4, <8, <16,<32             >32
        carD = str2double(cellfun(@(x) x{3}, carData,'uni',false));
        carHr = str2double(cellfun(@(x) x{4}, carData,'uni',false));
        carMin = str2double(cellfun(@(x) x{5}, carData,'uni',false));
        carTime = (carD-car1Day1)*24 + carHr + carMin/60;
        carLat = str2double(cellfun(@(x) x{7}, carData,'uni',false));
        carLon = -str2double(cellfun(@(x) x{8}, carData,'uni',false));
        carFuel = 0.01*str2double(cellfun(@(x) x{9}, carData,'uni',false));
        %car1Time = datenum(car1Time,'yyyy-mm-dd HH:MM');
        dT = SmplTime*ones(size(carTime)); %deltaF(carTime);

        dLat = deltaF(carLat);
        dLon = deltaF(carLon);
        dX = ScaleLat*dLat; % km
        dY = ScaleLon*dLon; % km
        GoDist = sqrt(dX.^2 + dY.^2);

        dFuel = deltaF(carFuel);
        %[usedFuel addedFuel] = deal(zeros(numel(carFuel),1));
        isGasAdded = dFuel > 0;
        usedFuel = (dFuel - abs(dFuel))/2; % Negative matrix decomposition
        addedFuel = (dFuel + abs(dFuel))/2;% Positive matrix decomposition
        usedWhenAdded = GoDist(isGasAdded)*TankPerKm;
        addedFuelCor = addedFuel;
        addedFuelCor(isGasAdded) = addedFuel(isGasAdded)+ usedWhenAdded;
        Fuel = usedFuel + addedFuelCor;

        isGoDist = GoDist > 0.0;
        isFuelChanged = Fuel ~= 0;
        isGoTime = isGoDist | isFuelChanged;
        
        GoTime = dT .* double(isGoTime);
        
        stat = 1111 * ones(numel(carFuel),1);
        LocInd = NaN(numel(carFuel),1);
        
        carDayLine = unique(carD);
        carDay = floor(carD);
        [~, markInd] = ismember(carDayLine, carDay);
        Nmark = numel(markInd);
        dInd = zeros(Nmark,2);
        for iM = 1:Nmark-1
            dInd(iM,:) = [markInd(iM), markInd(iM+1)-1];
        end
        dInd(Nmark,:) = [markInd(Nmark), numel(carDay)];
        
        [GoDistCumT, GoDistDaily]  = AccumF(GoDist, dInd); %, markInd);
        [GoTimeCumT, GoTimeDaily] = AccumF(GoTime, dInd);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dF = deltaF(F)
           dF = zeros(numel(F),1);
           dF(2:end) = F(2:end)-F(1:end-1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Fcum, Fmark] = AccumF(F, dInd) %, markInd)
           n = numel(F); 
           Fcum = zeros(n,1);
           nMark = size(dInd,1);
           Fmark = zeros(nMark,1);
           Fbuf = 0;
           for j = 1:nMark
               for i = dInd(j,1):dInd(j,2)
                   Fcum(i) = Fbuf + F(i);
                   Fbuf = Fcum(i);
               end
               Fmark(j) = Fbuf;
               if j > 1
                   Fmark(j) = Fmark(j) - Fmark(j-1);
               end
           end
           
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function UpdatePanel(varargin)
        % Handles to be updated
        % haxMap hTblCarData haxGraph

        [isGraphNew, isGrTypeNew] = deal(false);
        currGraph = get(hpupGraph,'val');
        if currGraph ~= LastGraph
            switch currGraph
                case 1
                    histData = CarsData.CarsGoTimeDaily;
                    histTitle = 'Car GoTime daily (hr)';
                    xTitle = 'Time, hr';
                    LastGraph = 1;
                case 2
                    histData = CarsData.CarsGoDistDaily;
                    histTitle = 'Car GoDist daily (km)';
                    xTitle = 'Distance, km';
                    LastGraph = 2;
            end
            isGraphNew = true;
        end
        
        currGrType = get(hpupGrType,'val');
        if currGrType ~= LastGrType
            isGrTypeNew = true;
        end        
        
        if isGraphNew || isGrTypeNew
            [nelements, xcenters] = hist(haxGraph, histData);
            nPlots = size(nelements,2);
            LgndStr = {[]};
            for i = 1:nPlots
                LgndStr{i} = ['day',num2str(i)];
            end

            switch currGrType
                case 1
                    bar(xcenters,nelements)
                    LastGrType = 1;
                case 2
                    plot(xcenters,nelements,'--','LineWidth',2,...
                        'MarkerEdgeColor','k','MarkerFaceColor','w',...
                        'Marker','o','MarkerSize',7)
                    LastGrType = 2;
            end
            title(histTitle,'FontSize',12) %,'FontWeight','bold')
            legend(LgndStr)
            xlabel(xTitle,'FontSize',10,'FontWeight','bold')
            grid

        end
        
        
        currTable = get(hpupTable,'val');
        if currTable ~= LastTable
            switch currTable
                case 1
                    TableData = CarsData.CarsGoTime;
                    LastTable = 1;
                case 2
                    TableData = CarsData.CarsGoDist;
                    LastTable = 2;
                case 3
                    TableData = CarsData.CarsFuelCum;
                    LastTable = 3;
            
%         CarsData.CarsLocInd = CarsLocInd;
%         CarsData.CarsGoDist = CarsGoDist;
%         CarsData.CarsGoTime = CarsGoTime;
%         CarsData.CarsGoDistCum = CarsGoDistCum;
%         CarsData.CarsGoTimeCum = CarsGoTimeCum;
%         CarsData.CarsFuelCum = CarsFuelCum;

            end
            isVarEdit = false(1,length(CarsData.CarsLic));
            set(hTblCarData,'ColumnEditable',isVarEdit,...
                'ColumnName',['time',CarsData.CarsLic'],...
                'Data', num2cell([CarsData.TimeLine, TableData]));
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupPanel(MainPnlNames,iTab)
        % Handles to be updated
        % haxMap hTblCarData haxGraph
        parent = hMainPnl{1}(iTab);
        bgc = get(parent,'backgroundcolor');
        switch MainPnlNames
            case 'Table'
                % Initially, displays raw data from input file
                hTblCarData = uitable(parent,'Tag','CarData','enable','inactive',...
                    'units','normalized','pos',[0.05 0.05 0.9 0.77]);

                hpupTable = uicontrol(parent,'style','popup',...
                    'Units','normalized','Tag','pupTable',...
                    'position',[0.05 0.85 0.2 0.12],...
                    'string',{'GoTime', 'GoDist', 'Fuel'},...
                    'value', 1,'visible','on',...
                    'Callback', @UpdatePanel);

                
            case 'Map'
                % Initially, displays static map (png OR jpg) with known
                % GEO corners such that each point can be properly
                % coordinated
                haxMap = axes('parent',parent,'units','normalized',...
                    'pos',[0.05 0.2 0.9 0.7]);

            case 'Graph'
                haxGraph = axes('parent',parent,'units','normalized',...
                    'pos',[0.05 0.2 0.9 0.63]);
                
                hpupGraph = uicontrol(parent,'style','popup',...
                    'Units','normalized','Tag','pupGraph',...
                    'position',[0.05 0.85 0.2 0.12],...
                    'string',{'GoTime', 'GoDist'},...
                    'value', 1,'visible','on',...
                    'Callback', @UpdatePanel);
                
                hpupGrType = uicontrol(parent,'style','popup',...
                    'Units','normalized','Tag','pupGraph',...
                    'position',[0.3 0.85 0.2 0.12],...
                    'string',{'bars', 'line'},...
                    'value', 1,'visible','on',...
                    'Callback', @UpdatePanel);

        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% -----------------------------------------------------------------
% ----   colWidth    ----------------------------------------------
% -----------------------------------------------------------------
function colW = colWidth(hTable, strCell)
    % Calculates column width for the given cell Array of strings

    posUnits = get(hTable,'units');
    set(hTable,'units','pixels')
    tmppos = get(hTable,'pos');
    maxLpxl = tmppos(3);
    set(hTable,'units',posUnits)

    [nr,nc] = size(strCell);
    maxLpxl = maxLpxl - (29+1.25*nc);
    maxW = 15*ones(1,nc);
    colW = num2cell(maxW);
    for j = 1:nc
        for i = 1:nr
            colW_ = length(strCell{i,j})*5;
            if colW_ > colW{j}
                colW{j} = colW_;
            end
        end
    end

    colWmat = cell2mat(colW);
    Lpxl = sum(colWmat);
    scale = maxLpxl/Lpxl;
    colWmat = round(scale*colWmat);
    colW = num2cell(colWmat); % mat2cell more complicated

end

% ========================================================================
% ========================================================================
% ====   Callback  FUNCTIONS     ===========================================
% ========================================================================
% ========================================================================
    function CBinpFile(varargin)
        
        LstStr = sprintf('(%s) file: %s ( %s min)',InputData);
        [selection, ok] = listdlg('ListString', LstStr,'SelectionMode','single' );
        
        if ok == 0
            selection = 1;
        end
        CarInpFnew = InputData{selection,2};
        
        if ~strcmpi(CarInpFnew, currCarInpF)
            currCarInpF = CarInpFnew;
            isHr = true;
            CarsData = getCarData(currCarInpF,isHr);
        end
    end
        
% hTopPnl = uipanel('parent',hF,'units','normalized',...
%             'pos',[0.02 0.72 0.95 0.25],'backgroundcolor',bgc);
%         inpFLegend = {'city','input_file_name','frq'};
%         htableInpF = uitable(hTopPnl,'Tag','tableInpF',...
%             'units','normalized','pos',[0.02 0.35 0.68 0.6],'enable','on');
%         ColW_inpData = colWidth(htableInpF, inpFLegend);
%         isVarEdit = false(1,3);
%         set(htableInpF,'ColumnEditable',isVarEdit,...
%             'ColumnName',inpFLegend,'ColumnWidth', ColW_inpData,...
%             'Data', InputData);
%         isNew_inpF = false; %#ok<NASGU>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function CB_graph_sel(varargin)
%         h = varargin{1};
%         currGraph = get(h,'val');
%         UpdatePanel(currGraph, LastTable)
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function CB_table_sel(varargin)
%         h = varargin{1};
%         currTable = get(h,'val');
%         UpdatePanel(LastGraph, currTable)
%     end

% ========================================================================
% ====     End of Functions     ==========================================
% ========================================================================

    
end    

% ========================================================================
% ====     THE END      =================================================
% ========================================================================
