%%    PLOT 2D Distributions %%%
function [] = plot_2D_distribution ...
    (input_densities, X1, Y1, input_uipanel, input_variables, color, contour, contour_percentile);
%clear uipanels if they exist
try child = allchild(input_uipanel);
    %clf; 
    clear map
catch; 
end;
tic
density_active = input_densities;
nsamples = size(density_active,3);
%parent = strcat('H.',input_uipanel);
x_int = ((max(X1(1,:))-min(X1(1,:)))/10);
xtics = min(X1(1,:)):x_int:max(X1(1,:));
y_int = ((max(Y1(:,1))-min(Y1(:,1)))/5);
ytics = min(Y1(:,1)):y_int:max(Y1(:,1));
%set(input_uipanel,'Visible', 'off');
f = waitbar(0, 'Plotting 2D densities', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0)


map=jet(100);


for i=1:nsamples
            left = 0.075;
            bottom = 0.075 + (((nsamples-i)/nsamples)*0.90);
            width = 0.9;
            height = 0.85/nsamples;
            pos = [left bottom width height];
            
            ax(i)=axes(input_uipanel, 'Position', pos);
               
            z_at_sig_level = findsig(density_active(:,:,i),contour_percentile/100);
            contour_level = z_at_sig_level;
            
            %%Plot color ramp, contours, or both
            if color == 1 && contour == 0
                %colormap(map);
                s = surf(ax(i),X1,Y1,density_active(:,:,i), 'EdgeColor','none');
                colormap(ax(i),map);
            end
            if color == 0 && contour == 1
                
                map = [1 1 1; 0 0 0];
                %colormap(ax(i),map);
                contourf(ax(i),X1,Y1,density_active(:,:,i),[0 contour_level 0.99],...
                   'black','LineWidth',1);
               colormap(ax(i),map);
                                
            end
            if color == 1 && contour ==1
                
                s = surf(ax(i),X1,Y1,density_active(:,:,i), 'EdgeColor','none');
                colormap(ax(i),map);
                hold(ax(i),'on')
                z = get(s,'ZData');
                set(s,'ZData',z-1e-10*z)
                c = contour3(ax(i),X1,Y1,density_active(:,:,i),[contour_level contour_level],'white','LineWidth',1);
                hold off
            end
            %%End Plot
            
            view(ax(i),[0, 90]);
            
            ax(i).TickDir = 'out';
            ax(i).TickLength = [0.01 0.02];
            ax(i).XTickLabel = {};
            xlim(ax(i),[min(X1(1,:)) max(X1(1,:))]);
            ylim(ax(i),[min(Y1(:,1)) max(Y1(:,1))]);
            xticks(ax(i), xtics)
            yticks(ax(i),ytics)
            %axis(ax(i),[lower_Xlim  upper_Xlim  lower_Ylim  upper_Ylim]);
           
            if getappdata(f,'canceling')
                break
            end
            waitbar(i/nsamples, f);
end
        delete(f);
        f=msgbox('Scaling and labeling','Please wait');
            xticklabels(ax(nsamples),'auto');
            xlabel(ax(nsamples), string(input_variables(1,1)),'Visible','on');
            yticklabels('auto')
            plot_top = ax(1);
            plot_bottom = ax(nsamples);
            p1=get(plot_top,'position');
            p2=get(plot_bottom,'position');
            height = -p2(2)+p1(2)+p1(4);
toc
tic
            h3=axes(input_uipanel, 'position',[p2(1) p2(2) p2(3) height],'visible','off');
            
            ylab = ylabel(h3, string(input_variables(1,2)), 'Visible','on');
            uistack(h3,'bottom');
toc
        delete(f);
        set(input_uipanel,'Visible', 'on');