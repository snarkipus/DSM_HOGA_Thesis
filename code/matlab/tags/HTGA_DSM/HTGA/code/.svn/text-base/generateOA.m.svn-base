function OA = generateOA(N)

%//////////////////////////////////////////////////////////////////////////
%///
%/// OA = generateOA(N)
%///
%/// Description:
%/// Create a 2-level orthogonal matrix of width N [Leung,Wang,2001]
%///
%/// Arguments:     [int]   N:  Required Number of Columns        
%///
%/// Returns:       [mat]   OA: Orthogonal Matrix
%///
%//////////////////////////////////////////////////////////////////////////


% Process Input
if( mod(log2(N),1) == 0 ) % Power of 2 Check
    J = log2(N)+1;
else
    J = ceil(log2(N));
end

Q = 2; % Number of Levels

% Create the basic columns
for k = 1:J    
    j = (Q^(k-1)-1)/(Q-1)+1;    
    for i = 1:Q^J 
        OA(i,j) = mod(floor( (i-1)/Q^(J-k) ),Q);        
    end    
end

% Create the non-basic columns
for k = 2:J    
    j = (Q^(k-1)-1)/(Q-1)+1;    
    for s = 1:j-1        
        for t = 1:Q-1            
            OA(:,j+(s-1)*(Q-1)+t) = mod(OA(:,s).*t+OA(:,j),Q);            
        end        
    end    
end

% Truncate columns for non-symmetric matrices
for k = (2^J-1):-1:N+1
    OA(:,k) = [];
end

OA = OA + 1;    