function startmatlabpool(varargin)
if  nargin ==1
    defsize = varargin{1};
else
    defsize = feature('numCores');
end

pools = gcp('nocreate');
if isempty(pools)
    poolsize = 0;
else
    poolsize = pools.NumWorkers;
end
    
if poolsize == 0  
    if nargin == 0  
        %parpool('local');
        c = parcluster;
        parpool(c);
    else  
        try  
            parpool('local',defsize);  
        catch ce  
            c = parcluster;
            parpool(c); 
            fail_p = gcp('nocreate');  
            fail_size = fail_p.NumWorkers;  
            display(ce.message);  
            display(strcat('输入的核数不正确，采用的默认配置defsize=',num2str(fail_size)));  
        end  
    end  
else  
    disp('parpool start ...');  
    if poolsize ~= defsize  
        closematlabpool();  
        startmatlabpool(defsize);  
    end  
end  