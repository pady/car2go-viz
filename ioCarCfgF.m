function [Data] = ioCarCfgF(oper, varargin )
% Reads Config file and parses it into structure
%   Input: 1) starting path     
%

% == 1) Gets config file
switch oper
    case 'r'
        if ~isempty(varargin)
            cfgFile = varargin{1};
        else
            defPath = pwd;
            cfgFile = 'car2goV01.cfg';
            % cfgFile = [pwd,filesep,'car2goV01.cfg'];
            answerID = exist(cfgFile,'file');
            if answerID ~= 2
                [cfgFname, cfgPath] = uigetfile({'*.cfg'},'Building config file',defPath);
                cfgFile = [cfgPath,cfgFname];
            end
            
        end
        
        fid = fopen(cfgFile,'r');
        nums = fread(fid);
        strCfg = char(nums');
        fclose(fid);
        
       [Data] = readCfgData(strCfg);
    case 'w'
       Data = varargin{1}; 
       cfgFname = varargin{2}; 
       [Data] = writeCfgData(Data,cfgFname);
        
end
end

% ========================================================================
% ========================================================================
% ====     FUNCTIONS     =================================================
% ========================================================================
% ========================================================================

function [Data] = readCfgData(Str)
%
exprData = ['< *([\w]*) *>'...   % Openning tag is city name 
    '\s*inpfile[= ]*"([^\f\n\r\t\v]*)"',...
    '\s*freq[= ]*"([\d]*)"',...
    '\s*Center[= ]*"([\d\.]*),[ ]*([\d\.]*)"',...
    '\s*Map[= ]*"([^\f\n\r\t\v]*)"',...
    '\s*Scale[= ]*"([\d\.E-]*),[ ]*([\d\.E-]*)"',...
    '\s*< */[\w ]*>'];            % closing tag
StrData = regexpi(Str,exprData,'tokens');
nData = length(StrData);


Data = cell(1,1);
for iD = 1:nData
    Data{iD,1} = StrData{1,iD}{1,1};  % City name
    Data{iD,2} = StrData{1,iD}{1,2};  % inpfile
    Data{iD,3} = StrData{1,iD}{1,3};  % freq
    Data{iD,4}{1} = str2double(StrData{1,iD}{1,4});  % Center Lat
    Data{iD,4}{2} = str2double(StrData{1,iD}{1,5});  % Center Lon
    Data{iD,5} = StrData{1,iD}{1,6};  % Map file
    Data{iD,6}{1} = str2double(StrData{1,iD}{1,7});  % Scale Lat
    Data{iD,6}{2} = str2double(StrData{1,iD}{1,8});  % Scale Lon
end

end

function status = writeCfgData(Data,cfgFname)
cfgFname %#ok<NOPRT>
Data{1}; %#ok<VUNUS>
status = true;
end



