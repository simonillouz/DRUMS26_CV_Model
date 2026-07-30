%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------
function [rout,J] = model_wrap(pars,data)
%global ALLPARS INDMAP

tpars =  data.ALLPARS;
tpars(data.INDMAP') = pars;

[rout, J] = SolveModel(tpars,data);
