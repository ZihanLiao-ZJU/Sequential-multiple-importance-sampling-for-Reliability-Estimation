classdef SeMIS < SeMCS
    methods
        function obj = SeMIS(IntBay,sampler)
            sampler.BayStc = IntBay;
            obj@SeMCS(sampler);
        end
    end
end
