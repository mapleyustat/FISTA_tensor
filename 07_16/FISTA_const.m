function opt = FISTA_const(x, f , g ,lambda, Prox_g, Grad_f, L , options)
% ************* FISTA Algorithm with constant step size ****************
% Input:
%   x         :  inital guess
%   f         :  function: 0.5*||Ax-b||^2
%   g         :  complexity function
%   lambda    :  trade-off parameter
%   Prox_g    :  proximal function of g(x)
%   Grad_f    :  gradient of f(x)
%   L         :  Lipschitz constant of Grad_f
%   options   :  control structure:
%                options.maxite : maximum iteration
%                options.quiet  : display result on each iteration
%                options.err    : error tolerance
% Output:
%   opt       :  output structure
%                opt.x          : optimal x
%                opt.value      : optimal value of objective function
%                opt.err        : error of each iteration
%                opt.t          : t of each iteration
%
% **********************************************************************
% by Jamie Zemin Zhang
% 07/16/2014

% 
% if isempty(options.maxite)
%     options.maxite = 1000;
% end
% 
% if isempty(options.quiet)
%     options.quiet = 0;
% end
% 
% if isempty(options.err)
%     options.err = 1e-2;
% end
maxite = options.maxite;

if ~options.quiet
        fprintf('%3s\t%16s\t%7s\n',...
            'iter','complexity','err')         ;
end

F = @(x)f(x)+g(x) ;
t = 1;  
y = x;
opt.x          = zeros(size(x,1),maxite);
opt.value      = zeros(maxite,1);
opt.complexity = zeros(maxite,1);

for iter = 1:maxite
    x_old     =  x;
    t_old     = t;
    [x,value] = Prox_g( y - 1/L*Grad_f(y), 1/L );
    t         = (1+sqrt(1+4*t_old^2))/2;
    y         = x + (t_old-1)/(t)*(x-x_old);
    
    opt.t           = t                  ;
    opt.x(:,iter)   = x                  ;
    opt.value(iter) = F(x)               ;
    opt.err(iter)   = 2*f(x)             ;
    opt.complexity(iter) = value/lambda   ;
    if ~options.quiet
        if mod(iter-1,10) ==0
            fprintf('%3d\t%10.4f\t%10.4f\n',...
                iter,opt.complexity(iter),opt.err(iter))  ;
        end
    end

    if f(x) < options.err
        break;
    end
end


