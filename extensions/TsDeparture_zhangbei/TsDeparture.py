# -*- coding: utf-8 -*-
# pass in Python3.6

def run_departure(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_TsDeparture')
    gi.log.info('TsDeparture START')    
    ####--------------End PART ONE--------------####    
    # part 2 ----------读取前节点界面参数---------- #
    Methods = data['pars']['Methods']
    print('debug: Threshold//'+str(Methods))            
    #data_path = data['OutputPath']
    orig_file = data['tsdata'] #pathlib.Path(data_path,'tsdata.txt')
    na_values = None
    date_format = None
    print('debug: orig_file//'+str(orig_file))
    ####--------------End PART TWO--------------####
    # part 3 ----------设置程序运行参数---------- #
    #Res_flag = data['GetResult'] #Print,JSON
    Mod_flag = data['DataMode']  #FileList,Url,DataTable
    # parameters for saving data
    tmppath = pathlib.Path(__file__).parent
    path = pathlib.Path(tmppath,'tmpData')
    if not path.exists():
        path.mkdir()
    res_file = pathlib.Path(tmppath,'tmpData','tsdata.csv')
    png_file = pathlib.Path(tmppath,'tmpData','tsdata.png')
    html_file = pathlib.Path(tmppath,'tmpData','tsdata.html')
    ####--------------End PART THREE--------------####
    # part 4 ----------专业逻辑实现部分---------- #
    # data processing
    import geoist.snoopy.tsa as tsa
    d=pd.read_csv(pathlib.Path(orig_file),delimiter=';',parse_dates=True,index_col=[0],na_values=na_values)
    res= tsa.departure(d['origin_data'].interpolate(),freq=Methods) 
    d['departure'] = res['departure_mean']
    d['residual'] = res['departure']
    d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
    gi.log.info('method = ' + str(Methods))
    gi.log.info('data length = ' + str(len(d['origin_data'])))
    # output html
    if Mod_flag == 'Url':
        import bokeh.plotting as bp
        from bokeh.palettes import Spectral4
        from bokeh.models import ColumnDataSource    
        p = bp.figure(title="The preliminary result by method = "+str(Methods)
                     , plot_width=535, plot_height=350, x_axis_type="datetime")
        bp.output_file(html_file, mode = 'inline')
        source = ColumnDataSource(d)
        for data, name, color in zip([source.column_names[1], source.column_names[2],source.column_names[3]], ["ORI", "Departure","Residual"], Spectral4):
           p.line(source.column_names[0], data, source=source,color=color, alpha=0.8, legend=name)
        p.xaxis.axis_label = 'Date'
        p.yaxis.axis_label = 'Value'
        p.legend.location = "top_left"
        p.legend.click_policy="hide"	
        bp.save(p)
        print(str(html_file)) #输出网络地址
        gi.log.info('output html finished')
     #output static pic
    if Mod_flag == 'FileList':
        SMALL_SIZE = 12
        MEDIUM_SIZE = 15
        BIGGER_SIZE = 18
        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=MEDIUM_SIZE,loc='upper left')    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title        
        ax = d.plot(figsize=(15,12),y=[d.columns[0],'departure','residual'])
        ax.set_xlabel('Date')
        ax.set_ylabel('Value')
        plt.grid()
        plt.legend()
        plt.title("The preliminary result by method={}".format(Methods),loc='left')
        plt.savefig(str(png_file),format='png') 
        gi.log.info('output png file finished')
        print(png_file) #输出一个图片
    ####--------------End PART FOUR--------------####
    # part 5 ------------输出数据------------------ #
    print(res_file) #输出数据表格文件zai 
    ####--------------End PART FIVE--------------####
        
if __name__ == '__main__':
    print('debug: starting//TsDeparture.py by zhangbei')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_departure(data)  #业务实现
    print('debug: finished//TsDeparture is finished and output now...')
