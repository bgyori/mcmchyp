function [gamma,gammas,ming,eta] = getGammaIter(f)
	[n,nf] = size(f);
	i = 1;
	eta(1) = 1;
	for j=1:nf
		Vf(j) = getVf(f(:,j));
	end
	while true
		for j=1:nf
			g(i,j) = getGammaC(f(:,j),eta(i),eta(i));
			gammas(i,j) = 1-(g(i,j)/Vf(j)).^(1/eta(i));
		end
		ming(i) = min(gammas(i,:));
		if i>1
			if ming(i)>=ming(i-1)
				break;
			end
		end
		eta(i+1) = ceil(log(n*ming(i))/(4*log(1/(1-ming(i)))));
		i = i + 1;
	end
	gamma = ming(end);
end