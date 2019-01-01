# coding: utf-8
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import re
import pandas as pd

#url_res = ['http://projects.biotec.tu-dresden.de/metapocket/result/1532488901_86/3tms-out.html',
#'http://projects.biotec.tu-dresden.de/metapocket/result/1532496555_47/8adh-out.html']
url_res = []
q1 = []
q2 = []

def get_url(url):
    pattern = re.compile(r'http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+')   
    url = re.findall(pattern,url)
    return url
    
def get_pockets(str1):
    pattern = re.compile(r"\'(.*?)\'")    
    str1 = re.findall(pattern,str1)
    return str1
    
def main():
    count = 0
    
    for line in  open('/Users/congwanwang/Desktop/48-bound-unbound/48_pdb_ids_copy.txt'):
        res = line.strip().split()
        a = res[0]
        b = res[1]
        q1.append(a)
        q2.append(b)
        
    for i in q1[:3]:
        count += 1
        print(i)
        params = OrderedDict([("pdbid", (None,i)),
                      ("chain", (None, '')),
                      ('filename', (None, '')),
                      ('number', (None, '3')),
                      ('email', (None, ''))])
        
        res = requests.post('http://projects.biotec.tu-dresden.de/mpcgi/metapocket.cgi', data=params)
         # 打印的结果：
        soup = BeautifulSoup(res.text, "lxml")
        url = str(soup.select("meta")[1])
        url = get_url(url)
        url_res.append(str(url[0]))
        print(url)
        print(count)
        
#to csv 文件
def get_res():
    
    res = []
    for i in url_res:
        r1 = []
        r2 = []
        print(i)
        res_1 = requests.get(i)
        soup = BeautifulSoup(res_1.text, "lxml")
        #print(soup)
        a = str(soup.select("pre")[1]).split("\n")[-5:-2]
        print("+++++++"*10)
        for  cc in a:
            r1.append(get_pockets(cc))
        b = str(soup.select("pre")[2]).split("\n")[-4:-1]
        for j in b:
            r2.append(j.split())
        for k in range(len(r1)):
            r2[k].append(",".join(r1[k]))
            r2[k].append(i[-13:-9])
        res.extend(r2)
        print(res)
    a = pd.DataFrame(res)
    print(a)
    a.to_csv("res.csv",index=False)
    
if __name__ == '__main__':
    main()
    #get_res()
