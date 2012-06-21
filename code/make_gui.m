uictrls = {};
uictrls.fig     = figure('Position',[10,400,450,250],...
                        'NumberTitle','off','DoubleBuffer','on',...
                        'BackingStore','on','Renderer','OpenGL',...
                        'Name','Granule Matcher','MenuBar','none');
    
clf(uictrls.fig);    


uictrls.quitbutton         = uicontrol(uictrls.fig,'Style','pushbutton','Position',[280 10 150 20],...
            'String','Quit','Value',0,'Callback',{'quitbutton'});
global quitsim;
quitsim=0;


uictrls.GC_to_use_label = uicontrol('style','text',...
                        'Position',[280 120 100 20],'String','GC to model: ','FontSize',14);  
uictrls.GC_to_use = uicontrol('style','edit',...
                        'Position',[400 120 30 20],'String',num2str(1));    
                    
uictrls.runs_to_avg_label = uicontrol('style','text',...
                        'Position',[280 80 100 20],'String','Runs to avg: ','FontSize',14);  
uictrls.runs_to_avg = uicontrol('style','edit',...
                        'Position',[400 80 30 20],'String',num2str(10));    

spikes          = uicontrol(uictrls.fig,'Style','checkbox','Position',[50 10 150 20],...
            'String','Spiking On','Value',0);
        
uictrls.spikes          = spikes;

synapse_sliders = {};
synapse_numbers = {};
for i=1:3

    synapse_sliders{i}  = uicontrol('Style', 'slider',...
                            'Min',0,'Max',10,'Value',0.5,...
                            'Position', [130 80+40*(i-1) 100 40]);
                        
    synapse_numbers{i}  = uicontrol('Style','popup','position',[20 80+40*(i-1) 100 40],...
    'string',{'no input', mftypes{:}},'Value',1);
        

end

uictrls.synapse_sliders = synapse_sliders;
uictrls.synapse_numbers = synapse_numbers;


synapse_title  = uicontrol('Style','text',...
                    'Position',[40 300 165 20],'FontSize',14,...
                    'String','Synapse number/weight');
                
uictrls.synapse_title = synapse_title;

uictrls.tau_syn_slider  = uicontrol('Style', 'slider',...
                        'Min',0,'Max',5,'Value',0.2,...
                        'Position', [280 150 100 20]);
                    
uictrls.tau_syn_label  = uicontrol('style','edit',...
                        'Position',[400 150 30 20],'String',num2str(0));       


uictrls.tau_mem_slider  = uicontrol('Style', 'slider',...
                        'Min',0,'Max',30,'Value',8.7,...
                        'Position', [280 180 100 20]);                    

uictrls.tau_mem_label  = uicontrol('style','edit',...
                        'Position',[400 180 30 20],'String',num2str(0));       

