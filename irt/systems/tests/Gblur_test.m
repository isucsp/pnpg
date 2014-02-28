% Gblur_test.m
% Test the Gblur object

if ~isvar('A'),	printm 'setup Gblur_test'
	psf = [0 1 2 1 0; 1 2 4 3 1; 0 2 3 1 0];
	psf = psf / sum(psf(:));
	idim = [64 70];

	mask = true(idim);
	A = Gblur(mask, 'psf', psf, 'type', 'conv,same');
	Af = Gblur(mask, 'psf', psf, 'type', 'fft,same');

	test_adjoint(A, 'big', 1e-11, 'complex', 1)
	test_adjoint(Af, 'big', 1e-11, 'complex', 1)

	im clf; im pl 3 3
	im(1, A.arg.psf, 'psf'), cbar
	im(2, A.arg.mask, 'mask'), cbar
prompt
end

% test A and A'
if 1
	x = shepplogan(idim(1), idim(2), 1);
	y1 = A * x;
	y2 = Af * x;

	x1 = A' * y1;
	x2 = Af' * y1;

	im(3, x, 'x')
	im(4, y1, 'A * x')
	im(5, x1, 'A'' * y1')

	equivs(x1, x2)
	equivs(y1, y2)
prompt
end

if 1
	Fatrix_test_basic(A, mask) % paces
end

% check adjoint
if 1, printm 'test adjoint'
	As = Gblur(true(21,22), 'psf', psf, 'type', 'fft,same');
%	test_adjoint(As);
	test_adjoint(As, 'big', 1e-12, 'complex', 1)
end
