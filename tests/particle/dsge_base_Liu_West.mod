@#include "dsge_base.inc"

options_.TeX=true;

estimation(order=3,posterior_sampling_method='online',filter_algorithm=nlkf,proposal_approximation=montecarlo,posterior_sampler_options=('particles',100));

collect_latex_files;
[status, cmdout]=system(['pdflatex -halt-on-error -interaction=nonstopmode ' M_.fname '_TeX_binder.tex']);
if status
    cmdout
    error('TeX-File did not compile.')
end
