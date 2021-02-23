%% Calculate Kuiper and KS from densities
function [V, V_1_to_1, mean_V_sink, range_V_sink, D, D_1_to_1, mean_D_sink, range_D_sink]=...
    KS_from_density(input_density1, input_density2, dimensions, hObject,Y1, input_variables)
% Takes as input two sets of density distributions of either univariate or
% bivariate data. Calculates KS D value for 1 and 2 dimensional samples.
% Calculates Kuiper V Value for 1 dimensional samples.  Returns an empty
% matrix for V values for 2 dimensional samples.  


%Get/set variables
if nargin == 4
    data=guidata(hObject);
    input_variables = data.input_variables;
    X1=data.X1;
    Y1=data.Y1;
elseif nargin == 6
    %use input variables for X1, Y1, input_variables
    X1 = hObject;
    Y1 = Y1;
    input_variables = input_variables;
else 
    err_dlg=errordlg('Check input variables.','Error!');
	waitfor(err_dlg);
end
if length(size(input_density1)) ~= length(size(input_density2))
    err_dlg=errordlg('Input distributions must have the same dimensions');
	waitfor(err_dlg);
end
 

if length(size(input_density1))==1;
    mean_V_sink = {'N/A'};
    range_V_sink = {'N/A'};
end

if length(size(input_density1))==2 || dimensions == 1
        pdps_active1=input_density1;
        pdps_active2=input_density2;
        nsamples1 = size(pdps_active1, 2);
        nsamples2 = size(pdps_active2, 2);
    
        for i=1:nsamples1
        for j=1:nsamples2
            V(i,j)=max(cumsum(pdps_active1(:,i))-...
                cumsum(pdps_active2(:,j)))+max(cumsum(pdps_active2(:,j))...
                -cumsum(pdps_active1(:,i)));
            D(i,j) = max(abs(cumsum(pdps_active1(:,i))-cumsum(pdps_active2(:,j))));
        end
        end
        mean_V_sink_temp=[];
        mean_D_sink_temp=[];
        V_1_to_1 = diag(V);
        for i=1:nsamples1
             mean_V_sink_temp = vertcat(mean_V_sink_temp, diag(V,i));
             mean_D_sink_temp = vertcat(mean_D_sink_temp, diag(D,i));
        end
        mean_V_sink=round(mean(mean_V_sink_temp),3);
        range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
        mean_D_sink=round(mean(mean_D_sink_temp),3);
        range_D_sink=round(max(mean_D_sink_temp)-min(mean_D_sink_temp),3);
end

if length(size(input_density1))==3 || dimensions == 0
        density_active1 = input_density1;
        density_active2 = input_density2;
        nsamples1 = size(input_density1, 3);
        nsamples2 = size(input_density2,3);
        
        %enable for testing% plot_2D_distribution(density_active1, X1, Y1, figure, input_variables);
        %enable for testing% plot_2D_distribution(density_active2, X1, Y1, figure, input_variables);
        
        for i=1:nsamples1
            density_active_CDF1(:,:,i) = cumsum(cumsum(density_active1(:,:,i),1),2);    % take the CDF of x at y = max (y)
            CDF1_Q1(:,:,i) = density_active_CDF1(:,:,i);                                            %Local CDF for Quadrant 1
            CDF1_Q2(:,:,i) = cumsum(cumsum(density_active1(:,:,i), 1,'reverse'), 2);            %Local CDF for Quadrant 2
            CDF1_Q3(:,:,i) = cumsum(cumsum(density_active1(:,:,i), 1,'reverse'), 2, 'reverse'); %Local CDF for Quadrant 3
            CDF1_Q4(:,:,i) = cumsum(cumsum(density_active1(:,:,i), 1), 2, 'reverse');           %Loval CDF for Quadrant 4
            
        end
        CDF1 = cat(3,CDF1_Q1, CDF1_Q2, CDF1_Q3, CDF1_Q4);
        
        for i=1:nsamples2
            density_active_CDF2(:,:,i) = cumsum(cumsum(density_active2(:,:,i),1),2);    % take the CDF of x at y = max (y)
            CDF2_Q1(:,:,i) = density_active_CDF2(:,:,i);                                            %Local CDF for Quadrant 1
            CDF2_Q2(:,:,i) = cumsum(cumsum(density_active2(:,:,i), 1,'reverse'), 2);            %Local CDF for Quadrant 2
            CDF2_Q3(:,:,i) = cumsum(cumsum(density_active2(:,:,i), 1,'reverse'), 2, 'reverse'); %Local CDF for Quadrant 3
            CDF2_Q4(:,:,i) = cumsum(cumsum(density_active2(:,:,i), 1), 2, 'reverse');           %Loval CDF for Quadrant 4
        end
        CDF2 = cat(3,CDF2_Q1, CDF2_Q2, CDF2_Q3, CDF2_Q4);
        
               
        %Kolmogorov-Smirnov distance
        for i=1:nsamples1
            for j=1:nsamples2
            D_Q1(i,j) = max(max(abs(CDF1_Q1(:,:,i) - CDF2_Q1(:,:,j)))); %Maximum absolute difference for Quadrant 1
            D_Q2(i,j) = max(max(abs(CDF1_Q2(:,:,i) - CDF2_Q2(:,:,j)))); %Maximum absolute difference for Quadrant 2
            D_Q3(i,j) = max(max(abs(CDF1_Q3(:,:,i) - CDF2_Q3(:,:,j)))); %Maximum absolute difference for Quadrant 3
            D_Q4(i,j) = max(max(abs(CDF1_Q4(:,:,i) - CDF2_Q4(:,:,j)))); %Maximum absolute difference for Quadrant 4
            
            end
        end
        D_temp = cat(3, D_Q1, D_Q2, D_Q3, D_Q4);
        D = max(D_temp, [], 3);
        
        %Kuiper distance
        for i=1:nsamples1
            for j=1:nsamples2
            V_Q1(i,j) = max(max((CDF1_Q1(:,:,i) - CDF2_Q1(:,:,j)))) + max(max((CDF2_Q1(:,:,j) - CDF1_Q1(:,:,i)))); %Maximum sum of differences for Quadrant 1
            V_Q2(i,j) = max(max((CDF1_Q2(:,:,i) - CDF2_Q2(:,:,j)))) + max(max((CDF2_Q2(:,:,j) - CDF1_Q2(:,:,i)))); %Maximum sum of differences for Quadrant 2
            V_Q3(i,j) = max(max((CDF1_Q3(:,:,i) - CDF2_Q3(:,:,j)))) + max(max((CDF2_Q3(:,:,j) - CDF1_Q3(:,:,i)))); %Maximum sum of differences for Quadrant 3
            V_Q4(i,j) = max(max((CDF1_Q4(:,:,i) - CDF2_Q4(:,:,j)))) + max(max((CDF2_Q4(:,:,j) - CDF1_Q4(:,:,i)))); %Maximum sum of differences for Quadrant 4
            end
        end
        V_temp = cat(3, V_Q1, V_Q2, V_Q3, V_Q4);
        V = max(V_temp, [], 3);
     
    mean_V_sink_temp=[];
    mean_D_sink_temp=[];
    
    if nsamples1==nsamples2
        for i=1:nsamples1
            mean_V_sink_temp = vertcat(mean_V_sink_temp,diag(V,i));
            mean_D_sink_temp = vertcat(mean_D_sink_temp,diag(D,i));
        end
        
    mean_V_sink  = round(mean(mean_V_sink_temp),3);
    range_V_sink = round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
    mean_D_sink  = round(mean(mean_D_sink_temp),3);
    range_D_sink = round(max(mean_D_sink_temp)-min(mean_D_sink_temp),3);
    else
        mean_V_sink  = [];
        range_V_sink = [];
        mean_D_sink  = [];
        range_D_sink = [];
    end
    
end
V_1_to_1 = diag(V);
D_1_to_1 = diag(D);

if length(size(input_density1))>3; 
    err_dlg=errordlg('Error, data dimensions >3');
	waitfor(err_dlg);    
end