%% Gaps info
function salida=gaps_stat_sample( ind,time_step, time_slots, gapMinTime )

Col2=[42.7/100 69.8/100 29.8/100; 9.8/100 29.8/100 0 ]   
LnWBar=2.3;
BarW=0.8;

CT2=[0.4, 0.44, 0];
timeTot=time_step*time_slots;


 rep=0;
     
    ctAux=0; 
    InterfaceRep=false;
    VertexRep=false;
				
	simple_max=114*40;
	cluster_max=simple_max*8;

    
  	ngapMinTime=floor(gapMinTime/time_step);
	ngapMinTime=10;
      
      clearvars gapVertex; 
      clearvars gapInt;
      
      ct_real=0;
      ct=0;
      ct_adh=0;
      if exist(sprintf('GapStat_%d.txt',ind), 'file')
      
        InterfaceRep=true;
        VertexRep=true;

        loop_continue=true;
      
        Vertex=true;
        Interface=true;
	ct_adh=ct_adh+1;
	


        while (loop_continue)
          
          if exist(sprintf('GapVertex_%d_%d.txt', ind,ct), 'file')
        %
            data = load (sprintf('GapVertex_%d_%d.txt', ind,ct));
            [n,m]=size(data);
            if (n > ngapMinTime )
      
              ct_real=ct_real+1;
              ct=ct+1;
               
              gapVertex(ct_real,1)=(data(n,1)-data(1,1))/60;   %time (min)
              gapVertex(ct_real,2)=mean(data(:,2))/1000000;   %size
              gapVertex(ct_real,3)=max(data(:,2))/1000000;    %max size
              gapVertex(ct_real,4)=mean(data(:,3));    %	Stress
         
              clearvars data;
            else  
              ct=ct+1;
            end
  
          else
        
            loop_continue=false 
            if ct==0
      		    Vertex=false;	
      	    end   
              
          end
      
        end  
        
        if ct_real==0
      		    Vertex=false;
        end
        
        %% Interface Gaps
        loop_continue=true;
        ct=0;
        ct_real=0;
        while (loop_continue)
          
          if exist(sprintf('R_%d/sample_%d/GapInterface_%d_%d.txt', ind,ct), 'file')

        %   

            data = load (sprintf('R_%d/sample_%d/GapInterface_%d_%d.txt', ind,ct));
            [n,m]=size(data);
            if (n > ngapMinTime)
      
              ct_real=ct_real+1;
              ct=ct+1;
        
              gapInt(ct_real,1)=(data(n,1)-data(1,1))/60;   %time (min)
              gapInt(ct_real,2)=mean(data(:,2))/1000000;   %size  (um2)
              gapInt(ct_real,3)=max(data(:,2))/1000000;    %max size (um2)
              gapInt(ct_real,4)=mean(data(:,3));    %	Stress
         
              clearvars data;
            else
                 ct=ct+1;
            end	
       
          else
        	  if (ct==0)
      		    Interface=false;	
      	    end   
              loop_continue=false    
              
          end
          
        end
        
         
        if ct_real==0
      		    Interface=false;
        end
       %% Plot
    
    
    
     
  
  
         
        if Vertex==true
          gapVertexStat(rep+1-ctAux,1)=mean(gapVertex(:,1));
       
        
        else
          gapVertexStat(rep+1-ctAux,1)=0;
        end
        if Interface==true
          gapIntStat(rep+1-ctAux,1)= mean(gapInt(:,1));
         
        else
          gapIntStat(rep+1-ctAux,1)=0;
        end
        
      
        if Vertex==true
          gapVertexStat(rep+1-ctAux,2)=mean(gapVertex(:,2));
         
        else
          gapVertexStat(rep+1-ctAux,2)=0;
        end
        if Interface==true
          gapIntStat(rep+1-ctAux,2)= mean(gapInt(:,2));
        else
          gapIntStat(rep+1-ctAux,2)=0;
        end
        
        
     
        gapVertexStat(rep+1-ctAux,3)=0;
        gapIntStat(rep+1-ctAux,3)=0;
       
        if Vertex==true

         [n,m]=size(gapVertex);
         gapVertexStat(rep+1-ctAux,3)=n/timeTot;
        
        end
        if Interface==true
         [n2,m]=size(gapInt);
         gapIntStat(rep+1-ctAux,3)=n2/timeTot;
       
        end
        grid off
        

    
        
       
        if Vertex==true
          gapVertexStat(rep+1-ctAux,4)=mean(gapVertex(:,4));
         
        else
          gapVertexStat(rep+1-ctAux,4)=0;
        end
        if Interface==true
          gapIntStat(rep+1-ctAux,4)= mean(gapInt(:,4));
      
        else
          gapIntStat(rep+1-ctAux,4)=0;
        end
        grid off
        
        
        %axis([0 6 0 Vmax])
       

        
      else % exists
        ctAux=ctAux+1;
  
     end  % end if exists
  
      

   
   
    
    h=figure()

    
    if VertexRep==true
 	[aux1,aux3]=size(gapVertexStat);
   	 SE_Vertex=std(gapVertexStat(:,1))/sqrt(aux1);
     bar(2, mean(gapVertexStat(:,1)),'FaceColor', [0, 0.44, 0.9],'linewidth', LnWBar, 'barwidth', BarW);
     hold on;
     errorbar(2, mean(gapVertexStat(:,1)),std(gapVertexStat(:,1)),'k','linewidth', LnWBar);
    end
    if InterfaceRep==true
   	 [aux2,aux3]=size(gapIntStat);
    	SE_Int=std(gapIntStat(:,1))/sqrt(aux2);
     bar(3.5, mean(gapIntStat(:,1)),'FaceColor', CT2 ,'linewidth', LnWBar, 'barwidth', BarW);
     hold on;
     errorbar(3.5, mean(gapIntStat(:,1)),std(gapIntStat(:,1)),'k','linewidth', LnWBar);
    end
    grid off
    
    
    
    set(gca, 'xscale', 'lin', 'yscale', 'lin','YGrid', 'off' ,'xminorgrid', 'off', 'yminorgrid', 'of','box', 'off', ...
         'fontsize',18 ,'xtick', [2 3.5], 'xticklabel', {'Vertex','Interface'},  ...
         'ticklength', [0.015 0.015], 'linewidth', 1.1, 'Position', [0.22 0.21 0.55 0.73]);
     
    
    
    ylabel('Time (min)','fontsize',21);
    a = get(gca,'XTickLabel');
    set(gca,'fontsize',19)
    % 
    % 
    % lgd=legend('Vertex','Interface', 'Best');
    % set(lgd,'FontSize',17);
    %   legend('boxoff')
    % 
    if (VertexRep==true ||  InterfaceRep==true)
        saveas(h,sprintf('Time_%d.jpg',ind));
    end    
    hold off;
    
    
    h=figure()
    


    if VertexRep==true
     SE_Vertex=std(gapVertexStat(:,2))/sqrt(aux1);
     bar(2, mean(gapVertexStat(:,2)),'FaceColor',  [0, 0.44, 0.9],'linewidth',LnWBar,  'barwidth', BarW);
     hold on;
     errorbar(2, mean(gapVertexStat(:,2)),std(gapVertexStat(:,2)),'k','linewidth', LnWBar);
    end
    if InterfaceRep==true
     SE_Int=std(gapIntStat(:,2))/sqrt(aux2);
     bar(3.5, mean(gapIntStat(:,2)),'FaceColor', CT2 ,'linewidth',LnWBar, 'barwidth', BarW);
     hold on;
     errorbar(3.5, mean(gapIntStat(:,2)),std(gapIntStat(:,2)),'k','linewidth', LnWBar);
    end
    grid off
    
    
    %axis([0 6 0 Vmax])
    set(gca, 'xscale', 'lin', 'yscale', 'lin','YGrid', 'off' ,'xminorgrid', 'off', 'yminorgrid', 'of','box', 'off', ...
         'fontsize',15  ,'xtick', [2 3.5], 'xticklabel', {'Vertex','Interface'},  ...
         'ticklength', [0.015 0.015], 'linewidth', 1.1, 'Position', [0.22 0.21 0.55 0.73]);
     
    
    
    ylabel('Size (um2)','fontsize',17);
    a = get(gca,'XTickLabel');
    set(gca,'fontsize',17)
    
    % 
    % lgd=legend('Vertex','Interface', 'Best');
    % set(lgd,'FontSize',17);
    %   legend('boxoff')
    
     if (VertexRep==true ||  InterfaceRep==true)
         saveas(h,sprintf('Size_%d.jpg',ind));
     end
    
    hold off;
    
    h=figure()
    

    if VertexRep==true
     SE_Vertex=std(gapVertexStat(:,3))/sqrt(aux1);
     bar(2, mean(gapVertexStat(:,3)),'FaceColor', [0, 0.44, 0.9],'linewidth',LnWBar, 'barwidth', BarW);
     hold on;
     errorbar(2, mean(gapVertexStat(:,3)),std(gapVertexStat(:,3)),'k','linewidth', LnWBar);
    end
    if InterfaceRep==true
     SE_Int=std(gapIntStat(:,3))/sqrt(aux2);
     bar(3.5, mean(gapIntStat(:,3)),'FaceColor', CT2 ,'linewidth',LnWBar, 'barwidth', BarW);
     hold on;
     errorbar(3.5, mean(gapIntStat(:,3)),std(gapIntStat(:,3)),'k','linewidth', LnWBar);
    end
    grid off
    
    
    
    %axis([0 6 0 Vmax])
    set(gca, 'xscale', 'lin', 'yscale', 'lin','YGrid', 'off' ,'xminorgrid', 'off', 'yminorgrid', 'of','box', 'off', ...
         'fontsize',18  ,'xtick', [2 3.5], 'xticklabel', {'Vertex','Interface'},  ...
         'ticklength', [0.015 0.015], 'linewidth', 1.1, 'Position', [0.22 0.21 0.55 0.73]);
    
    ylabel('Gap Number/hour','fontsize',21);
    a = get(gca,'XTickLabel');
    set(gca,'fontsize',19)
    
if VertexRep==true
    [nn,mm]=size(gapVertexStat);
        
    [nn2,mm2]=size(gapIntStat);
    str = sprintf('Sample %d',nn);
    t=text(2.5,1,str,'FontSize',16)
end	 

  
  %  annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
%      lgd=legend(sprintf('Vertex (%d)',nn),sprintf('Interface (%d)',nn2));
%      set(lgd,'FontSize',17);
%        legend('boxoff')
    
     if (VertexRep==true ||  InterfaceRep==true)
         saveas(h,sprintf('Gap_Number_%d.jpg',ind));
     end
    hold off;



%     h=figure()
%     if VertexRep==true
%      bar(2, mean(gapVertexStat(:,4)),'FaceColor', Col2(1,:),'linewidth',1.9);
%      hold on;
%      errorbar(2, mean(gapVertexStat(:,4)),std(gapVertexStat(:,4)),'k','linewidth', LnWBar);
%     end
%     if InterfaceRep==true
%      bar(3.5, mean(gapIntStat(:,4)),'FaceColor', CT2 ,'linewidth',1.9);
%      hold on;
%      errorbar(3.5, mean(gapIntStat(:,4)),std(gapIntStat(:,4)),'k','linewidth', LnWBar);
%     end
%     grid off
%     
%     
%     
%     %axis([0 6 0 Vmax])
%     set(gca, 'xscale', 'lin', 'yscale', 'lin','YGrid', 'on' ,'xminorgrid', 'off', 'yminorgrid', 'of','box', 'off', ...
%          'fontsize',18  ,'xtick', [2 3.5], 'xticklabel', {'Vertex','Interface'},  ...
%          'ticklength', [0.015 0.015], 'linewidth', 1.1, 'Position', [0.22 0.21 0.55 0.73]);
%     
%     ylabel('Stress','fontsize',21);
%     a = get(gca,'XTickLabel');
%     set(gca,'fontsize',19)
%     
%     % 
%     % lgd=legend('Vertex','Interface', 'Best');
%     % set(lgd,'FontSize',17);
%     %   legend('boxoff')
%     
%      if (VertexRep==true ||  InterfaceRep==true)
%          saveas(h,sprintf('R_%d/Stress_%d.jpg',glob,ind-1));
%      end
    hold off;
clearvars gapVertexStat;
clearvars gapIntStat;
   
 quit;
end
    
    
    
