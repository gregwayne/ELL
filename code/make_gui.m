uictrls = {};
uictrls.fig     = figure('Position',[10,400,400,250],...
                        'NumberTitle','off','DoubleBuffer','on',...
                        'BackingStore','on','Renderer','OpenGL',...
                        'Name','Granule Matcher','MenuBar','none');
    
clf(uictrls.fig);    

resample        = uicontrol(uictrls.fig,'Style','pushbutton','Position',[220 120 150 20],...
            'String','Resample from this synapse','Value',0,'Callback',{'resample_synapses'});
compute         = uicontrol(uictrls.fig,'Style','pushbutton','Position',[220 40 150 20],...
            'String','Compute','Value',0,'Callback',{'compute_voltage'});
        
uictrls.resample        = resample;
uictrls.compute         = compute;

granule         = uicontrol(uictrls.fig,'Style','pushbutton','Position',[220 80 150 20],...
            'String','Get Granule','Value',0,'Callback',{'choose_granule'});
        
uictrls.granule         = granule;
        

spikes          = uicontrol(uictrls.fig,'Style','checkbox','Position',[50 10 150 20],...
            'String','Spiking On','Value',0);
        
uictrls.spikes          = spikes;

synapse_sliders = {};
synapse_numbers = {};
for i=1:3

    synapse_sliders{i}  = uicontrol('Style', 'slider',...
                            'Min',0,'Max',1,'Value',0.5,...
                            'Position', [30+50*(i-1) 50 40 100],...
                            'Callback', {@set_weight,i});
            
    synapse_numbers{i}  = uicontrol('style','edit',...
                            'Position',[50+50*(i-1) 160 30 20],...
                            'Callback',{@Temp_call,i});            

end

uictrls.synapse_sliders = synapse_sliders;
uictrls.synapse_numbers = synapse_numbers;


synapse_title  = uicontrol('Style','text',...
                    'Position',[40 200 165 20],'FontSize',14,...
                    'String','Synapse number/weight');
                
uictrls.synapse_title = synapse_title;




