# -*- coding: utf-8 -*-
# pass in Python3.6

def run_tsTools(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_tsTools')
    gi.log.info('tsTools START')    
    ####--------------End PART ONE--------------####    
    # part 2 ----------读取前节点界面参数---------- #
    algorithm_name = data['pars']['Tools'] 
    print('debug: algorithm_name//'+str(algorithm_name))
    
    WindowSize = data['pars']['WindowSize'] 
    Periods = data['pars']['Periods'] 
    Factor = data['pars']['Factor'] 
    FitNum = data['pars']['FitNum'] 
    Methods = data['pars']['Methods'] 
    WindowSize=int(WindowSize)    
    na_values = None
    nan_value = 0
    date_format = None
    center = False
    data_path = data['OutputPath']
    orig_file = data['tsdata'] #pathlib.Path(data_path,'gradata.txt')
    print('debug: orig_file//'+str(orig_file))
    #####--------------End PART TWO--------------####   
    ## part 3 ----------设置程序运行参数---------- #
    Res_flag = data['GetResult'] #Print,JSON
    Mod_flag = data['DataMode']  #FileList,Url,DataTable
    # parameters for saving data
    tmppath = pathlib.Path(__file__).parent
    path = pathlib.Path(tmppath,'tmpData')
    if not path.exists():
        path.mkdir()
    res_file = pathlib.Path(tmppath,'tmpData','tsdata.csv')
    png_file = pathlib.Path(tmppath,'tmpData','tsdata.png')
    html_file = pathlib.Path(tmppath,'tmpData','tsdata.html')
    #####--------------End PART THREE--------------####
    ## part 4 ----------专业逻辑实现部分---------- #
    # data processing
    import geoist.snoopy.tsa as tsa
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)    
    if algorithm_name == '去线性':
        res = tsa.detrend(d[d.columns[0]].interpolate())
        d['trans_data'] = d[d.columns[0]].interpolate() - res
        d['residual'] = res
        d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
        params = 'linear'
    elif algorithm_name == '指数平滑':
        d['trans_data'] = d.interpolate().ewm(alpha=Factor).mean()
        d['residual'] = d['origin_data'].interpolate()-d['trans_data']
        d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
        params = Factor
    elif algorithm_name == '滑动平均':
        d['trans_data'] = d.interpolate().rolling(window=WindowSize,center=center,min_periods=1).mean()
        d['residual'] = d['origin_data'].interpolate()-d['trans_data']
        d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';',header=False)
        params = WindowSize
    elif algorithm_name == '一阶差分':
        res = d['origin_data'].interpolate()
        for i in range(1):
            res = res.diff(periods=int(Periods))
        d['trans_data'] = res.fillna(nan_value)
        d['residual'] = d['origin_data'].interpolate()-d['trans_data']
        d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
        params = Periods
    elif algorithm_name == '二阶差分':
        res = d['origin_data'].interpolate()
        for i in range(2):
            res = res.diff(periods=int(Periods))
        d['trans_data'] = res.fillna(nan_value)
        d['residual'] = d['origin_data'].interpolate()-d['trans_data']
        d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
        params = Periods
    else: #多项式拟合
        import numpy as np
        res = d['origin_data'].interpolate()
        x = np.arange(len(res))
        z1 = np.polyfit(x, res, int(FitNum))
        print('debug: fit params//'+str(z1))
        d['trans_data'] = np.polyval(z1,x)
        d['residual'] = d['origin_data'].interpolate()-d['trans_data']
        d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
        params = FitNum
        
    gi.log.info('algorithm_name = ' + str(algorithm_name))
    gi.log.info('data length = ' + str(len(d['origin_data'])))
    # output html
    if Mod_flag == 'Url':
        import bokeh.plotting as bp
        from bokeh.palettes import Spectral4
        from bokeh.models import ColumnDataSource    
        p = bp.figure(title="The preliminary result by algorithm = "+str(algorithm_name)
                     , plot_width=535, plot_height=350, x_axis_type="datetime")
        bp.output_file(html_file, mode = 'inline')
        source = ColumnDataSource(d)
        
        for data, name, color in zip([source.column_names[1], source.column_names[2], source.column_names[3]], ["ORI", "Trans", "Residual"], Spectral4):
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
        ax = d.plot(figsize=(15,12),y=[d.columns[0],'trans_data','residual'])
        ax.set_xlabel('Date')
        ax.set_ylabel('Value')
        plt.grid()
        plt.legend()
        plt.title("The preliminary result by params={}".format(params),loc='left')
        plt.savefig(str(png_file),format='png') 
        gi.log.info('output png file finished')
        print(png_file) #输出一个图片
    #####--------------End PART FOUR--------------####
    #
    ## part 5 ------------输出数据------------------ #
    print(res_file) #输出数据表格文件 
    #####--------------End PART FIVE--------------####

if __name__ == '__main__':
    print('debug: starting//TsTools.py by zhangbei')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_tsTools(data)  #业务实现
    print('debug: finished//TsTools is finished and output now...')
