# -*- coding: utf-8 -*-
# pass in Python3.6
def run_EQCatlog(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_RandGen')
    gi.log.info('RandGen START')    
    ####--------------End PART ONE--------------####    
    #import shutil  
    #shutil.copy(argfile,  'debug.JSON') #仅在调试时候用    
    # part 2 ----------读取前节点界面参数---------- #
    windowsize = int(data['pars']['windowsize']) 
    Chan1 = data['pars']['Chan1'] 
    Chan2 = data['pars']['Chan2'] 
    Chan3 = data['pars']['Chan3'] 
    print('debug: 生成信号长度//'+str(windowsize))
    #####--------------End PART TWO--------------####
    ## part 3 ----------设置程序运行参数---------- #
    #Res_flag = data['GetResult'] #Print,JSON
    Mod_flag = data['DataMode']  #FileList,Url,DataTable
    # parameters for saving data
    tmppath = pathlib.Path(__file__).parent
    path = pathlib.Path(tmppath,'tmpData')
    if not path.exists():
        path.mkdir()    
    csv_file = pathlib.Path(tmppath,'tmpData','tsdata.csv')
    png_file = pathlib.Path(tmppath,'tmpData','tsdata.png')
    html_file = pathlib.Path(tmppath,'tmpData','catalog.html')
    #####--------------End PART THREE--------------####

    ## part 4 ----------专业逻辑实现部分---------- #
    # data processing "正态分布;0-1分布;卡方分布;泊松分布;指数分布"
    import numpy as np
    #np.random.seed(20190801) #随机数可以预测
    chandata = pd.DataFrame(columns=['id', 'chan1', 'chan2', 'chan3'])
    chandata['id'] = np.arange(windowsize)
    if Chan1 == '正态分布':
        chandata['chan1'] = np.random.randn(windowsize)  #dataframe
    elif Chan1 == '0-1分布':
        chandata['chan1'] = np.random.rand(windowsize)  #dataframe
    elif Chan1 == '卡方分布':
        chandata['chan1'] = np.random.chisquare(1, windowsize)  #dataframe
    elif Chan1 == '指数分布':
        chandata['chan1'] = np.random.standard_exponential(windowsize)  #dataframe
    elif Chan1 == '泊松分布':
        chandata['chan1'] = np.random.poisson(5, windowsize)                 
    else:
        chandata['chan1'] = np.random.randint(1, windowsize)  #dataframe

    if Chan2 == '正态分布':
        chandata['chan2'] = np.random.randn(windowsize)  #dataframe
    elif Chan2 == '0-1分布':
        chandata['chan2'] = np.random.rand(windowsize)  #dataframe
    elif Chan2 == '卡方分布':
        chandata['chan2'] = np.random.chisquare(1, windowsize)  #dataframe
    elif Chan2 == '指数分布':
        chandata['chan2'] = np.random.standard_exponential(windowsize)  #dataframe
    elif Chan2 == '泊松分布':
        chandata['chan2'] = np.random.poisson(5, windowsize)                 
    else:
        chandata['chan2'] = np.random.randint(1, windowsize)   #dataframe

    if Chan3 == '正态分布':
        chandata['chan3'] = np.random.randn(windowsize)  #dataframe
    elif Chan3 == '0-1分布':
        chandata['chan3'] = np.random.rand(windowsize)  #dataframe
    elif Chan3 == '卡方分布':
        chandata['chan3'] = np.random.chisquare(1, windowsize)  #dataframe
    elif Chan3 == '指数分布':
        chandata['chan3'] = np.random.standard_exponential(windowsize)  #dataframe
    elif Chan3 == '泊松分布':
        chandata['chan3'] = np.random.poisson(5, windowsize)                 
    else:
        chandata['chan3'] = np.random.randint(1, windowsize)   #dataframe     
        
    chandata.to_csv(csv_file,sep=';',index=False,encoding="utf-8")
    
    gi.log.info('data length = ' + str(windowsize))

    # output html
    if Mod_flag == 'Url':
        import bokeh.plotting as bp
        from bokeh.palettes import Spectral4
        from bokeh.models import ColumnDataSource    
        p = bp.figure(title="The preliminary result by windowsize = "+str(windowsize)
                     , plot_width=535, plot_height=350)
        bp.output_file(html_file, mode = 'inline')
        source = ColumnDataSource(chandata)
        #print(source.column_names)
        for data, name, color in zip([source.column_names[2], source.column_names[3]
            , source.column_names[4]], [Chan1, Chan2, Chan3], Spectral4):
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
        ax = d.plot(figsize=(15,12),y=[d.columns[0],'chan1','chan2','chan3'])
        ax.set_xlabel('Date')
        ax.set_ylabel('Value')
        plt.grid()
        plt.legend()
        plt.title("The preliminary result by windowsize={}".format(windowsize),loc='left')
        plt.savefig(str(png_file),format='png') 
        gi.log.info('output png file finished')
        print(png_file) #输出一个图片
#    #####--------------End PART FOUR--------------####
#    #
#    ## part 5 ------------输出数据------------------ #
    print(csv_file) #输出数据表格文件zai 
#    #####--------------End PART FIVE--------------####

if __name__ == '__main__':
    print('debug: starting//RandGen.py by chenshi')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_EQCatlog(data)  #业务实现
    print('debug: RandGen is finished.')


