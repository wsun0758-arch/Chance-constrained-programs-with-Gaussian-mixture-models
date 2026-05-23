function y=myfun(x)
global GMModel;

y=-sum(GMModel.ComponentProportion*GMModel.mu*x);