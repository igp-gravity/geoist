# -*- coding: utf-8 -*-
# pass in Python3.6
def run_TsSeason(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_TsSeason')
    gi.log.info('TsSeason START')    
    ####--------------End PART ONE--------------####    
    # part 2 ----------读取前节点界面参数---------- #
    period_str = data['pars']['Period']
    print('debug: Period//'+str(period_str))
    period=int(period_str)    
    na_values_output = 0.0            
    orig_file = data['tsdata'] 
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
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)
    d = d.interpolate() 
    decomposition = tsa.seasonal_decompose(d,freq=period,extrapolate_trend='freq')      
    d['trend'] = decomposition.trend.fillna(na_values_output)
    d['seasonal'] = decomposition.seasonal.fillna(na_values_output)
    d['residual'] = decomposition.resid.fillna(na_values_output)
    d.to_csv(pathlib.Path(res_file),sep=';')
    gi.log.info('period = ' + str(period))
    gi.log.info('data length = ' + str(len(d['trend'])))
    # output html
    if Mod_flag == 'Url':
        import bokeh.plotting as bp
        from bokeh.palettes import Spectral4
        from bokeh.models import ColumnDataSource    
        p = bp.figure(title="The preliminary result by period = "+str(period)
                     , plot_width=535, plot_height=350, x_axis_type="datetime")
        bp.output_file(html_file, mode = 'inline')
        source = ColumnDataSource(d)
        
        for data, name, color in zip([source.column_names[1], source.column_names[2], source.column_names[3], source.column_names[4]], ["ORI", "Trend", "Seasonal", "Residual"], Spectral4):
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
        #ax = d.plot(figsize=(15,12),y=[d.columns[0],'despiked'])
        fig = decomposition.plot()
        fig.set_size_inches(15,12)
        plt.grid()
        plt.legend()
        plt.title("The preliminary result by period={}".format(period),loc='left')
        plt.savefig(str(png_file),format='png') 
        gi.log.info('output png file finished')
        print(png_file) #输出一个图片
    ####--------------End PART FOUR--------------####
    # part 5 ------------输出数据------------------ #
    print(res_file) #输出数据表格文件zai 
    ####--------------End PART FIVE--------------####
        
if __name__ == '__main__':
    print('debug: starting//TsSeason.py by zhangbei')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_TsSeason(data)  #业务实现
    print('debug: finished//TsSeason is finished and output now...')



