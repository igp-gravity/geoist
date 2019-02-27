# -*- coding: utf-8 -*-
# pass in Python3.6
def run_TsEDA(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_TsEDA')
    gi.log.info('TsEDA START')    
    ####--------------End PART ONE--------------####    
    #import shutil  
    #shutil.copy(argfile,  'debug.JSON') #仅在调试时候用    
    # part 2 ----------读取前节点界面参数---------- #
    Method = data['pars']['method'] 
    print('debug: EDA分析方案//'+str(Method))
    orig_file = data['inputdata']
    print('debug: 分析数据文件//'+str(orig_file))    
    #####--------------End PART TWO--------------####
    ## part 3 ----------设置程序运行参数---------- #
    #Res_flag = data['GetResult'] #Print,JSON
    Mod_flag = data['DataMode']  #FileList,Url,DataTable
    # parameters for saving data
    tmppath = pathlib.Path(__file__).parent
    path = pathlib.Path(tmppath,'tmpData')
    if not path.exists():
        path.mkdir()    
    #csv_file = pathlib.Path(tmppath,'tmpData','tsdata.csv')
    png_file1 = pathlib.Path(tmppath,'tmpData','tsdata1.png')
    png_file2 = pathlib.Path(tmppath,'tmpData','tsdata2.png')
    png_file3 = pathlib.Path(tmppath,'tmpData','tsdata3.png')
    #html_file = pathlib.Path(tmppath,'tmpData','catalog.html')
    #####--------------End PART THREE--------------####

    ## part 4 ----------专业逻辑实现部分---------- #
    # data processing "正态分布;0-1分布;卡方分布;泊松分布;指数分布"
    import seaborn as sns
    #np.random.seed(20190801) #随机数可以预测
    d=pd.read_csv(pathlib.Path(orig_file),delimiter=';')
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
    Method = 'AN3'
    if Mod_flag == 'FileList': 
        if Method == 'AN1':
            sns.pairplot(d)
            #plt.title("The EDA result by method={}".format(Method),loc='left')
            plt.savefig(str(png_file1),format='png')
            plt.figure(figsize=(15,12))
            sns.heatmap(d.corr(), annot=True, fmt=".2f")
            plt.title("The Heatmap of correlation by method={}".format(Method),loc='left')           
            plt.savefig(str(png_file1),format='png')  
            print(png_file1)
        elif Method == 'AN2':
            #g = sns.JointGrid(x='chan1', y='chan2',data = d)
            #g = g.plot(sns.regplot, sns.distplot)             
            sns.jointplot(data = d, x='chan1', y='chan2', kind='reg', color='r' )        
            plt.savefig(str(png_file1),format='png')
            sns.jointplot(data = d, x='chan1', y='chan3', kind='reg', color='g' )
            plt.savefig(str(png_file2),format='png')
            sns.jointplot(data = d, x='chan1', y='chan4', kind='reg', color='b' )
            #sns.distplot(d['chan2'].dropna())
            #sns.distplot(d['chan3'].dropna())
            #sns.distplot(d['chan4'].dropna())
            plt.savefig(str(png_file3),format='png') 
            print(png_file1)
            print(png_file2)
            print(png_file3)            
        else:
            fig, axes = plt.subplots(2, 2,figsize=(15, 12)) #, sharex=True)  #ax = axes[0]
            sns.distplot(d['chan1'].dropna(), ax = axes[0, 0])
            sns.distplot(d['chan2'].dropna(), ax = axes[0, 1])
            sns.distplot(d['chan3'].dropna(), ax = axes[1, 0])
            sns.distplot(d['chan4'].dropna(), ax = axes[1, 1])
            sns.plt.savefig(str(png_file1),format='png')  
            print(png_file1)
        gi.log.info('output png file finished={}'.format(Method))
         #输出一个图片
    #gi.log.info('data length = ' + str(len(d)))

if __name__ == '__main__':
    print('debug: starting//TsEDA.py by chenshi')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_TsEDA(data)  #业务实现
    print('debug: TsEDA is finished.')
