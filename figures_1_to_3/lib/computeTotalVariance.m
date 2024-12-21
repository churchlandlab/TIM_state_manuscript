function cMap = computeTotalVariance(Vc, Vm, U)
Vc = reshape(Vc,size(Vc,1),[]);
Vm = reshape(Vm,size(Vm,1),[]);
if length(size(U)) == 3
    U = arrayShrink(U, squeeze(isnan(U(:,:,1))));
end
covVc = cov(Vc');  % S x S
covVm = cov(Vm');  % S x S
cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
covP = sum((U * cCovV) .* U, 2)';  % 1 x P
varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
cMap = gather((covP ./ stdPxPy)') .^ 2;
end