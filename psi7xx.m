function z = psi7xx(x,y)
z = -640*y.^4-240*log(x).*y.^4+2160*x.^2.*log(x).*y.^2+2160*x.^2.*y.^2-450*x.^4.*log(x)-165*x.^4;