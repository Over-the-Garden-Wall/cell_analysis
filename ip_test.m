function ip_test(script, varargin)
    p = inputParser;    
    p.addRequired('script', @ischar);
    p.addOptional('format', 'html', @ischar);
    
    p.parse(script, varargin{:});
    
    s = p.Results;
end