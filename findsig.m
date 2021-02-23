function contour_level = findsig (input_density, siglevel)
%Inputs 
%input_density:     meshgrid of probability density 
%siglevel:          significance level for contour
%
%Outputs
%contour_level:     the z-value that creates a contour that encompasses the 
%                   percent of the data specified by the significance level 
%                   (so you can contour at that z-value).
contour_level = fzero(@fun, [0 1]);
    function y = fun(x)
        base_density = input_density-(x);
        base_density(base_density<0) = 0;
        base_density = nonzeros(base_density)+x;
        y = siglevel - sum(base_density);
    end
end