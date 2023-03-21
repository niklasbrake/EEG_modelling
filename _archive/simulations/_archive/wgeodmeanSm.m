function vnew = wgeodmeanSm(data,w,Error)
    % geodmeanS2    weighted geodesic mean of data on S^m using Exponential map and its %               inverse map. Iteratively finds the weighted mean.
    %
    % Inputs:   vdata   - a matrix (m+1)-by-n : a column vector is a point on S^m
    %           w :     - a vector of weights (default = equal weights)
    %           Error   - a error allowance for optimization  (default = 1e-10)
    %
    % Output:
    %           wgeodmean    - weighted geodesic mean
    %
    %
    % 2008.07.08, 2008.09.28, Aug 10, 2009
    % 10/12/2011
    % Sungkyu Jung

    [m, n] = size(data);

    if nargin == 1
        w = 1/n*ones(1,n);
        Error = 1e-10;
    end

    if nargin == 2
        Error = 1e-10;
    end

    vini = wmean(w,data,2);
    vini = vini/norm(vini);
    % Initial candidate for geodesic mean

    diff = 1;
    cnt = 0;
    while diff > Error
        rot= rotMat(vini);
        logdata = LogNPd(rot*data) ;
        u = wmean(w,logdata,2);
        vnew = rot\ExpNPd(u);

        diff = norm(vnew - vini);
        vini = vnew;
        cnt = cnt+1;
        % if cnt > 50
        %     disp(['cnt > 50 ' num2str(diff)]);
        %     break;
        % end;
    end
    function rot = rotMat(b,a,alpha)
        % ROTMAT returns a rotation matrix that rotates unit vector b to a
        %
        %   rot = rotMat(b) returns a d x d rotation matrix that rotate
        %   unit vector b to the north pole (0,0,...,0,1)
        %
        %   rot = rotMat(b,a ) returns a d x d rotation matrix that rotate
        %   unit vector b to a
        %
        %   rot = rotMat(b,a,alpha) returns a d x d rotation matrix that rotate
        %   unit vector b towards a by alpha (in radian)
        %
        %    See also .

        % Last updated Nov 7, 2009
        % Sungkyu Jung
        [s1 s2]=size(b);
        d = max(s1,s2);
        b= b/norm(b);
        if (min(s1,s2) ~= 1 || nargin==0)
            help rotMat
            return
        end

        if s1<=s2;
            b = b';
        end

        if nargin == 1;
            a = [zeros(d-1,1); 1];
            alpha = acos(a'*b);
        end
        if nargin == 2;
            alpha = acos(a'*b);
        end
        if abs(a'*b - 1) < 1e-15;
            rot = eye(d);
            return,
        end
        if abs(a'*b + 1) < 1e-15;
            rot = -eye(d);
            return,
        end

        c = b - a * (a'*b); c = c / norm(c);
        A = a*c' - c*a' ;

        rot = eye(d) + sin(alpha)*A + (cos(alpha) - 1)*(a*a' +c*c');
    end
    function Exppx = ExpNPd(x)
        % EXPNP Riemannian exponential map at North pole of S^k
        %       ExpNP(v) returns (k+1) x n matrix where each column is a point on a
        %                sphere and the input v is k x n matrix where each column
        %                is a point on tangent  space at north pole.
        %
        %
        %   See also LogNPd.

        % Last updated Oct 20, 2009
        % Sungkyu Jung

        [d n] = size(x);
        nv = sqrt(sum(x.^2));
        Exppx = [repmat(sin(nv)./nv,d,1).*x ;cos(nv)];
        Exppx(:,nv < 1e-16) = repmat([zeros(d,1);1],1,sum(nv < 1e-16));
    end
    function Logpx = LogNPd(x)
        % LOGNP Riemannian log map at North pole of S^k
        %       LogNP(x) returns k x n matrix where each column is a point on tangent
        %       space at north pole and the input x is (k+1) x n matrix where each column
        %       is a point on a sphere.
        %
        %
        %   See also ExpNPd.
        % Last updated Oct 20, 2009
        % Sungkyu Jung
        [d n] = size(x);
        scale = acos(x(end,:))./sqrt(1-x(end,:).^2);
        scale(isnan(scale)) =1;
        Logpx = repmat(scale,d-1,1).*x(1:(end-1),:);
    end
    function y = wmean(w,data,dim)
        if(nargin<3)
            dim = 1;
        end
        y = sum(w.*data,dim)./sum(w,dim);
    end
end