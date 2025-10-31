function [U,V] = aca(A, tol)

	[M,N] = size(A);

	I = 1;

	U = [];
	V = [];
	
	SF2 = 0;

	for k=1:min(M,N)
		% find value and col J of largest element
		% of R in row I
		Ri = A(I, :);
		if k > 1
			Ri = Ri - U(I, :)*V';
		end
		[~,J] = max(abs(Ri));

		% set Vk to Linf normalized copy
		% of Ith row of Residual
		Vk = (Ri / Ri(J))';
		% set Uk to copy of
		% Jth column of Residual
		Uk = A(:,J);
		if k > 1
			Uk = Uk - U*(V(J,:))';
		end
		
		% update norm of approximation
		nVk = norm(Vk);
		nUk = norm(Uk);

		SF2(k+1) = SF2(k) + (nVk*nUk)^2;
		if k > 1
			SF2(k+1) = SF2(k+1) + 2*sum(real((U'*Uk).*(V'*Vk))); 
		end
		
		V = [V Vk];
		U = [U Uk];
		
		% terminate if added norm is small relative
		% to norm of entire approximation
		if nUk*nVk <= tol*sqrt(SF2(k+1))
			break;
		end
		
		% calculate index of next row to zero out
		Uk = abs(Uk);
		Uk(I) = 0;
		[~,I] = max(Uk);
	end

end
