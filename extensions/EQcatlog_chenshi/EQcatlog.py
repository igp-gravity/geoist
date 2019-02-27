# -*- coding: utf-8 -*-
# pass in Python3.6
def run_EQCatlog(data):
    # part 1 ----------import and set loginfo---------- #  
    import pandas as pd
    import matplotlib.pyplot as plt
    import pathlib
    # output logging 
    import geoist as gi      
    gi.log.setname('Datist_EQCatlog')
    gi.log.info('EQCatlog START')    
    ####--------------End PART ONE--------------####    
    #import shutil  
    #shutil.copy(argfile,  'debug.JSON') #仅在调试时候用    
    # part 2 ----------读取前节点界面参数---------- #
    CatalogName = data['pars']['Dataset'] 
    print('debug: CatalogName//'+str(CatalogName))
    #####--------------End PART TWO--------------####
    ## part 3 ----------设置程序运行参数---------- #
    #Res_flag = data['GetResult'] #Print,JSON
    Mod_flag = data['DataMode']  #FileList,Url,DataTable
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

    ## part 4 ----------专业逻辑实现部分---------- #
    # data processing
    from geoist.others import fetch_data
    import requests, os
    #先获得HASH信息从Github到本地
    url = "https://raw.githubusercontent.com/gravity-igpcea/dataset/master/"
    filename = url + 'regHash.txt'
    f = open(downfile,'wb+')    
    with f:
        r = requests.get(filename, stream=True)
        r.raise_for_status()
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)

    fetch_data.Drepo1.load_registry(downfile)
    
    if CatalogName == 'CN':
        cncat = fetch_data.fetch_catalogCN()  #dataframe
    elif CatalogName == 'GEM':
        cncat = fetch_data.fetch_catalogGEM()  #dataframe
    else:
        cncat = fetch_data.fetch_catalog('catalog_xc.TXT')  #dataframe

    cncat.to_csv(file1,sep=';',index=False,encoding="utf-8")
    catinfo = pd.concat([cncat.head(),cncat.tail()])    
    catinfo.to_csv(fileinfo1,sep=';',index=True,encoding="utf-8")
    catdesc = cncat.describe()
    catdesc.to_csv(filedesc1,sep=';',index=True,encoding="utf-8", float_format = '%.2f')
    
    gi.log.info('Download = ' + str(CatalogName))
    gi.log.info('data length = ' + str(len(cncat['id'])))

    # output html
    if Mod_flag == 'Url':
        import bokeh.plotting as bp
        from bokeh.palettes import Spectral4
        from bokeh.models import ColumnDataSource    
        p = bp.figure(title="The catalog with index = "+str(CatalogName)
                     , plot_width=535, plot_height=350)
        bp.output_file(html_file, mode = 'inline')
        p.vbar(x = cncat['id'].values, width=1, bottom=0, top = cncat['mag'].values, color='firebrick', alpha=0.8)
        p.xaxis.axis_label = 'Index'
        p.yaxis.axis_label = 'Magnitude'
        bp.save(p)
        print(str(html_file)) #输出网络地址
        gi.log.info('output html finished')
#     #output static pic
    if Mod_flag == 'FileList':
        print(fileinfo1) #输出数据表格文件zai 
        print(filedesc1) #输出数据表格文件zai 
#    #####--------------End PART FOUR--------------####
#    #
#    ## part 5 ------------输出数据------------------ #
    print(file1) #输出数据表格文件zai 
#    #####--------------End PART FIVE--------------####

if __name__ == '__main__':
    print('debug: starting//EQCatlog.py by chenshi')    
    import sys, json
    argfile=sys.argv[1] #json参数  
    print('debug: open//'+argfile) 
    with open(argfile,'rb') as f:
        data = json.load(f)    
    run_EQCatlog(data)  #业务实现
    print('debug: EQCatlog is finished.')


