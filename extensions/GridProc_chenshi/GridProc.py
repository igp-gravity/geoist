# -*- coding: utf-8 -*-
# pass in Python3.6
def run_GridProc(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_GridProc')
    gi.log.info('GridProc START')    
    ####--------------End PART ONE--------------####    
    import shutil  
    shutil.copy(argfile,  'debug.JSON') #仅在调试时候用    
    # part 2 ----------读取前节点界面参数---------- #
    dataName = data['pars']['data'] 
    formatName = data['pars']['format'] 
    print('debug: dataName:{1}and{2}'.format(1=dataName,2=formatName))
#    #####--------------End PART TWO--------------####
#    ## part 3 ----------设置程序运行参数---------- #
#    #Res_flag = data['GetResult'] #Print,JSON
#    Mod_flag = data['DataMode']  #FileList,Url,DataTable
    # parameters for saving data
    tmppath = pathlib.Path(__file__).parent
    path = pathlib.Path(tmppath,'tmpData')
    if not path.exists():
        path.mkdir()    
    file1 = pathlib.Path(tmppath,'tmpData','catalog.csv')
    fileinfo1 = pathlib.Path(tmppath,'tmpData','catsample.csv')
    filedesc1 = pathlib.Path(tmppath,'tmpData','catdesc.csv')
    downfile = pathlib.Path(tmppath,'tmpData','Hash.txt')
    html_file = pathlib.Path(tmppath,'tmpData','catalog.html')
    #####--------------End PART THREE--------------####

#    ## part 4 ----------专业逻辑实现部分---------- #
#    # data processing
#    from geoist.others import fetch_data
#    import requests, os
#    #先获得HASH信息从Github到本地
#    url = "https://raw.githubusercontent.com/gravity-igpcea/dataset/master/"
#    filename = url + 'regHash.txt'
#    f = open(downfile,'wb+')    
#    with f:
#        r = requests.get(filename, stream=True)
#        r.raise_for_status()
#        for chunk in r.iter_content(chunk_size=1024):
#            if chunk:
#                f.write(chunk)
#
#    fetch_data.Drepo1.load_registry(downfile)
#    
#    if CatalogName == 'CN':
#        cncat = fetch_data.fetch_catalogCN()  #dataframe
#    elif CatalogName == 'GEM':
#        cncat = fetch_data.fetch_catalogGEM()  #dataframe
#    else:
#        cncat = fetch_data.fetch_catalog('catalog_xc.TXT')  #dataframe
#
#    cncat.to_csv(file1,sep=';',index=False,encoding="utf-8")
#    catinfo = pd.concat([cncat.head(),cncat.tail()])    
#    catinfo.to_csv(fileinfo1,sep=';',index=True,encoding="utf-8")
#    catdesc = cncat.describe()
#    catdesc.to_csv(filedesc1,sep=';',index=True,encoding="utf-8", float_format = '%.2f')
#    
#    gi.log.info('Download = ' + str(CatalogName))
#    gi.log.info('data length = ' + str(len(cncat['id'])))
#
#    # output html
#    if Mod_flag == 'Url':
#        import bokeh.plotting as bp
#        from bokeh.palettes import Spectral4
#        from bokeh.models import ColumnDataSource    
#        p = bp.figure(title="The catalog with index = "+str(CatalogName)
#                     , plot_width=535, plot_height=350)
#        bp.output_file(html_file, mode = 'inline')
#        p.vbar(x = cncat['id'].values, width=1, bottom=0, top = cncat['mag'].values, color='firebrick', alpha=0.8)
#        p.xaxis.axis_label = 'Index'
#        p.yaxis.axis_label = 'Magnitude'
#        bp.save(p)
#        print(str(html_file)) #输出网络地址
#        gi.log.info('output html finished')
##     #output static pic
#    if Mod_flag == 'FileList':
#        print(fileinfo1) #输出数据表格文件zai 
#        print(filedesc1) #输出数据表格文件zai 
#    #####--------------End PART FOUR--------------####
#    #
#    ## part 5 ------------输出数据------------------ #
    print(file1) #输出数据表格文件zai 
#    #####--------------End PART FIVE--------------####

if __name__ == '__main__':
    print('debug: starting//GridProc.py by chenshi')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_GridProc(data)  #业务实现
    print('debug: GridProc is finished.')



#
#import pandas as pd
#from pandas import DataFrame
#import numpy as np    
#
#
##接入节点界面中所设参数，并转换成固定代码。
#Intpol_str = '[$Interpolation$]'              #插值方法：线性;样条，分别赋值0;1 字符型
#Xaxis_str = '[$Xaxis$]'                       #新坐标X：手动输入数值,正则校验，字符型
#Yaxis_str = '[$Yaxis$]'                       #新坐标Y：手动输入数值,正则校验，字符型
#RowNum_str = '[$RowNum$]'                     #新行数：手动输入整数,正则校验，字符型
#ConNum_str = '[$ConNum$]'                     #新列数：手动输入整数,正则校验，字符型
#RowGap_str = '[$RowGap$]'                     #新行距：手动输入数值,正则校验，字符型
#ConGap_str = '[$ConGap$]'                     #新列距：手动输入数值,正则校验，字符型
#
##转换数据类型
#Intpol = int(Intpol_str)             
#Xaxis = float(Xaxis_str)                
#Yaxis = float(Yaxis_str) 
#RowNum = int(RowNum_str)        
#ConNum = int(ConNum_str)              
#RowGap = float(RowGap_str)             
#ConGap = float(ConGap_str)                     
#
#basepath = GetNodeBasePath()                      #获取自定义节点路径
#df = CsvFile2DataFrame(basepath + 'grddata.txt')  #获取grd文件全路径，包含文件名，grdfile = df['DocName']
#
#
#
#OutputText(Intpol_str + ' | ' + Xaxis_str + ' | ' + Yaxis_str + ' | ' + RowNum_str +  ' | ' +ConNum_str + ' | ' + RowGap_str + ' | ' + ConGap_str)
#OutputTable(df, basepath + 'output.txt')
#
#
