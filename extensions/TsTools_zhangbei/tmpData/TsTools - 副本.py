"""
    此代码分为四个部分，第一部分为参数替换区，本部分代码的目的为接入节点页面中设置的参数，
并将之替换为代码中使用的参数变量。用户在使用时，可以只修改此处的参数，不用去修改代码中的参数。
    第二部分为数据转换区，功能为接入前节点传入的数据，并转化为代码可以处理的格式。
	第三部分为代码执行区，用户可以在此编辑将要执行的代码。在实际使用时可以重点关注此部分。
	第四部分为结果导出区，此部分的功能为导出符合之后节点使用的结果文件。

"""

import pandas as pd
import matplotlib.pyplot as plt
import geoist.snoopy.tsa as tsa
import sys
import pathlib
import matplotlib.pyplot as plt
import shutil
import io
import base64

#接入节点界面中所设参数，并转换成固定代码。
algorithm_name  = '[$Tools$]'    #算法类型：去线性;指数平滑;滑动平均;一阶差分;二阶差分;距平分析;   传入的值为选项文字本身。
period_str = '[$Periods$]'       #差分步长：≥1的整数;
window_str = '[$WindowSize$]'    #窗口大小：≥1的整数;
factor_str = '[$Factor$]'        #系数：0.00~1小数;

Methods = '[$Methods$]'          #分析方法：日平距;月平距;   传入的值为选项文字本身。
RunModes = '[$RunModes$]'        #运行方式：静态图;交互图;反计算;   传入的值为选项文字本身。



         
#获取前节点数据
data_path = GetNodeBasePath()
orig_file = pathlib.Path(data_path,'gradata.txt')
nan_value = 0
na_values = None
date_format = None

# parameters for saving data
html_file = AllocateFileName('html')

# parameters for saving data
res_file = AllocateTableFileName()

na_values_output = 0.0


# parameters for processing data
factor = float(factor_str)

if (__name__ == '__main__') and ( algorithm_name == '指数平滑') :
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)
    d['ewm'] = d.interpolate().ewm(alpha=factor).mean()
    d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
    
    image = io.BytesIO()
    d.plot(figsize=(15,12),y=[d.columns[0],'ewm'])
    plt.savefig(image,format='png')
    data=base64.b64encode(image.getvalue()).decode('utf8')
    picdata = "<img src=data:image/png;base64,{}>".format(data)
    htmldata = "<html>"+picdata+"</html>"
    with open(pathlib.Path(html_file),'w') as f:
        f.write(htmldata) 

# parameters for processing data
window_size = int(window_str)
center = False

if (__name__ == '__main__') and ( algorithm_name == '滑动平均') :
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=';',index_col=[0],na_values=na_values)
    d['ma'] = d.interpolate().rolling(window=window_size,center=center,min_periods=1).mean()
    d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';',header=False)
    
    image = io.BytesIO()
    d.plot(figsize=(15,12),y=[d.columns[0],'ma'])
    plt.savefig(image,format='png')
    data=base64.b64encode(image.getvalue()).decode('utf8')
    picdata = "<img src=data:image/png;base64,{}>".format(data)
    htmldata = "<html>"+picdata+"</html>"
    with open(pathlib.Path(html_file),'w') as f:
        f.write(htmldata)  


if (__name__ == '__main__') and ( algorithm_name == '去线性') :
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)
    d['detrend'] = tsa.detrend(d[d.columns[0]].interpolate())
    d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')

    image = io.BytesIO()
    d.plot(figsize=(15,12),y=[d.columns[0],'detrend'])
    plt.savefig(image,format='png')
    data=base64.b64encode(image.getvalue()).decode('utf8')
    picdata = "<img src=data:image/png;base64,{}>".format(data)
    htmldata = "<html>"+picdata+"</html>"
    with open(pathlib.Path(html_file),'w') as f:
        f.write(htmldata)

# parameters for processing data
periods = int(period_str)

if (__name__ == '__main__') and ( algorithm_name == '一阶差分') :
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)
    res = d['origin_data'].interpolate()
    for i in range(1):
        res = res.diff(periods=int(periods))
    d['diff'] = res.fillna(nan_value)
    d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
    
    image = io.BytesIO()
    d.plot(figsize=(15,12),y=[d.columns[0],'diff'])
    plt.savefig(image,format='png')
    data=base64.b64encode(image.getvalue()).decode('utf8')
    picdata = "<img src=data:image/png;base64,{}>".format(data)
    htmldata = "<html>"+picdata+"</html>"
    with open(pathlib.Path(html_file),'w') as f:
        f.write(htmldata)

if (__name__ == '__main__') and ( algorithm_name == '二阶差分') :
    d=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)
    res = d['origin_data'].interpolate()
    for i in range(2):
        res = res.diff(periods=int(periods))
    d['diff'] = res.fillna(nan_value)
    d.to_csv(pathlib.Path(res_file),date_format=date_format,sep=';')
    
    image = io.BytesIO()
    d.plot(figsize=(15,12),y=[d.columns[0],'diff'])
    plt.savefig(image,format='png')
    data=base64.b64encode(image.getvalue()).decode('utf8')
    picdata = "<img src=data:image/png;base64,{}>".format(data)
    htmldata = "<html>"+picdata+"</html>"
    with open(pathlib.Path(html_file),'w') as f:
        f.write(htmldata)

OpenHttpUrl('file:///' + html_file)
