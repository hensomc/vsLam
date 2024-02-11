function eigB = getBaseEig(bcType)
%getBaseEig()
%    eigB  = returned Baseline eigenvalues

switch bcType
    case {'cfff','fffc'}
        eigB = [207.7; 507.1; 1274.8; 1625.7; 1850.6; 3238.4];  %C-F-F-F
    case 'sfsf'
        eigB = [572.89; 960.09; 2184.88; 2316.94; 2780.33; 4207.97];  %S-F-S-F
    case 'cfcf'
        eigB = [1324.73; 1578.13; 2597.11; 3656.54; 4018.21; 4752.85];  %C-F-C-F
    case 'ssss'
        eigB = [1174.23; 2935.58; 2935.58; 4696.93; 5871.16; 5871.16];  %S-S-S-S
    case 'cccc'
        eigB = [2140.86; 4366.99; 4366.79; 6442.22; 7828.22; 7863.91];  %C-C-C-C
    otherwise
        %warning('Unexpected bcType')
end
end