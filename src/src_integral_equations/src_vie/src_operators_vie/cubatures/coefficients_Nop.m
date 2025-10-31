function f = coefficients_Nop(dx,n,np,p,q,l,lp,ker_type)

%% inputs and outputs 

% dx       : resolution
% n        : (3x1) normal vector on the surface of the observation voxel
% np       : (3x1) normal vector on the surface of the source voxel
% p        : \in {1,2,3} indicating the component of the unit vector associated with the testing function
% q        : \in {1,2,3} indicating the component of the unit vector associated with the basis function
% l        : \in {1,2,3,4} indicating the scalar part of the linear/pulse testing function for each component of the vector field
% lp       : \in {1,2,3,4} indicating the scalar part of the linear/pulse basis function for each component of the vector field
% ker_type : \in {1,2,3,4} indicating the reduced surface-surface kernel to be calculated
% f        : double value of the coefficient


%% calculate one coefficient

% testing function's unit vector
pp = zeros(3,1);
pp(p) = 1;

% basis function's unit vector
qq = zeros(3,1);
qq(q) = 1;


switch ker_type 
    
    case 1
        %% surface-surface kernel 1
        
        f = dot(cross(n,pp),cross(np,qq));
        
        
    case 2
        %% surface-surface kernel 2
       
        switch l
            case 1
                h = [0 0 0];
            case 2
                h = [1/dx 0 0];
            case 3
                h = [0 1/dx 0];
            case 4
                h = [0 0 1/dx];
        end
        
        f = dot ( cross(cross(h,pp),qq) , n );
        
        
    case 3
        %% surface-surface kernel 3
       
        switch lp
            case 1
                hp = [0 0 0];
            case 2
                hp = [1/dx 0 0];
            case 3
                hp = [0 1/dx 0];
            case 4
                hp = [0 0 1/dx];
        end

        f = dot ( cross(pp,cross(hp,qq)) , n );

        
    case 4
        %% surface-surface kernel 4
        
        switch l
            case 1
                h = [0 0 0];
            case 2
                h = [1/dx 0 0];
            case 3
                h = [0 1/dx 0];
            case 4
                h = [0 0 1/dx];
        end

        switch lp
            case 1
                hp = [0 0 0];
            case 2
                hp = [1/dx 0 0];
            case 3
                hp = [0 1/dx 0];
            case 4
                hp = [0 0 1/dx];
        end

        a = cross(h,pp);

        q1 = cross(a,[1 0 0]);
        q2 = cross(a,[0 1 0]);
        q3 = cross(a,[0 0 1]);

        w1 = hp * dot(qq,[1 0 0]);
        w2 = hp * dot(qq,[0 1 0]);
        w3 = hp * dot(qq,[0 0 1]);

        f = ( dot(q1,n)*dot(w1,np) + dot(q2,n)*dot(w2,np) + dot(q3,n)*dot(w3,np) ); 
        
end

end