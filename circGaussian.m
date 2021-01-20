function f = circGaussian(phistim,phipref,sigma,ampl)

f = ampl * exp(-0.5*(phistim-phipref).^2/sigma^2) / (sigma*sqrt(2*pi));


