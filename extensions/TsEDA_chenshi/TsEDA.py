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
import shutil
import io
import base64

#接入节点界面中所设参数，并转换成固定代码。
windowsize_str = '[$windowsize$]'     

#转换数据类型
window_size = int(windowsize_str)                 


#获取前节点数据
data_path = GetNodeBasePath()
orig_file = pathlib.Path(data_path,'gradata.txt')
na_values = None
date_format = None


# parameters for saving data
res_file = AllocateTableFileName()
html_file = AllocateFileName('html')

na_values = None
date_format = None



# parameters for saving data
res_data = AllocateFileName('txt')

if __name__ == '__main__':
    data=pd.read_csv(pathlib.Path(orig_file),parse_dates=True,delimiter=";",index_col=[0],na_values=na_values)
    data['origin_data'] = data[data.columns[0]].interpolate()  
    res = tsa.adfuller(data['origin_data'].values)
    with open(pathlib.Path(res_data),'w') as f:
        tsa.print_adf(res,'original data',file=f)
    
    image = io.BytesIO()
    data['mean'] = data['origin_data'].rolling(window=window_size).mean()
    data['std'] = data['origin_data'].rolling(window=window_size).std()
    data.plot(figsize=(15,12),y=['origin_data','mean','std'])
    plt.savefig(image,format='png')
    data=base64.b64encode(image.getvalue()).decode('utf8')
    picdata = "<img src=data:image/png;base64,{}>".format(data)
    htmldata = "<html>"+picdata+"</html>"
    with open(pathlib.Path(html_file),'w') as f:
        f.write(htmldata)

OpenHttpUrl('file:///' + html_file)

