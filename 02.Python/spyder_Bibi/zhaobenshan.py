#!/usr/bin/env python3

#%%
import requests
from fake_useragent import UserAgent
import json
import time
import pandas as pd


# 下载网页评论数据
def get_page_json(url):
    try:
        ua = UserAgent(verify_ssl=False)
        headers = {"User-Agent": ua.random}
        json_comment = requests.get(url, headers=headers).text
        return json_comment
    except:
        return None


# 解析网页评论数据
def parse_page_json(json_comment):
    try:
        comments = json.loads(json_comment)
    except:
        return "error"

    comments_list = []
    #获取当页数据有多少条评论（一般情况下为20条）
    num = len(comments['data']['replies'])

    for i in range(num):
        comment = comments['data']['replies'][i]
        comment_list = []
        floor = comment['floor']
        ctime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(comment['ctime'])) # 时间转换
        likes = comment['like']
        author = comment['member']['uname']
        sex = comment['member']['sex']
        level = comment['member']['level_info']['current_level']
        content = comment['content']['message'].replace('n', '') # 将评论内容中的换行符去掉
        #print(content)
        rcount = comment['rcount']
        comment_list.append(floor)
        comment_list.append(ctime)
        comment_list.append(likes)
        comment_list.append(author)
        comment_list.append(sex)
        comment_list.append(level)
        comment_list.append(content)
        comment_list.append(rcount)

        comments_list.append(comment_list)

    save_to_csv(comments_list)


def save_to_csv(comments_list):
    data = pd.DataFrame(comments_list)
    # 注意存储文件的编码为utf_8_sig，不然会乱码，后期会单独深入讲讲为何为这样（如果为utf-8）
    data.to_csv('春晚鬼畜_1.csv', mode='a', index=False, sep=',',\
            header=False, encoding='utf_8_sig')


def main():
    base_url = "https://api.bilibili.com/x/v2/reply?&type=1&oid=19390801&pn=1"
    json_comment = get_page_json(base_url)
    parse_page_json(json_comment)


    # 通过首页获取评论总页数
#    pages = int(json.loads(get_page_json(base_url))['data']['page']['count'])//20
#    for page in range(pages):
#        url = "https://api.bilibili.com/x/v2/reply?&type=1&oid=19390801&pn=" + str(page)
#        json_comment = get_page_json(url)
#        parse_page_json(json_comment)
#        print("正在保存第%d页" % int(page+1))
#
#        if page % 20 == 0:
#            time.sleep(5)


main()


#%%
import requests
from lxml import etree
import time
import re

# 获取某市区域的所有链接
def get_areas(url):
    print('start grabing areas')
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.108 Safari/537.36'}
    resposne = requests.get(url, headers=headers)
    content = etree.HTML(resposne.text)
    areas = content.xpath("//dd[@data-index = '0']//div[@class='option-list']/a/text()")
    areas_link = content.xpath("//dd[@data-index = '0']//div[@class='option-list']/a/@href")
    for i in range(1, len(areas)):
        area = areas[i]
        area_link = areas_link[i]
        link = 'https://bj.lianjia.com' + area_link
        print("开始抓取页面")
        get_pages(area, link)


# 通过获取某一区域的页数，来拼接某一页的链接
def get_pages(area, area_link):
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.108 Safari/537.36'}
    resposne = requests.get(area_link, headers=headers)
    pages = int(re.findall("page-data=\'{\"totalPage\":(\d+),\"curPage\"", resposne.text)[0])
    print("这个区域有" + str(pages) + "页")
    for page in range(1, pages+1):
        url = 'https://bj.lianjia.com/zufang/dongcheng/pg' + str(page)
        print("开始抓取" + str(page) + "的信息")
        get_house_info(area, url)


# 获取某一区域某一页的详细房租信息
def get_house_info(area, url):
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.108 Safari/537.36'}
    time.sleep(2)
    try:
        resposne = requests.get(url, headers=headers)
        content = etree.HTML(resposne.text)
        info = []
        for i in range(30):
            title = content.xpath("//div[@class='where']/a/span/text()")[i]
            room_type = content.xpath("//div[@class='where']/span[1]/span/text()")[i]
            square = re.findall("(\d+)", content.xpath("//div[@class='where']/span[2]/text()")[i])[0]
            position = content.xpath("//div[@class='where']/span[3]/text()")[i].replace(" ", "")
            try:
                detail_place = re.findall("([\u4E00-\u9FA5]+)租房", content.xpath("//div[@class='other']/div/a/text()")[i])[0]
            except Exception as e:
                detail_place = ""
            floor =re.findall("([\u4E00-\u9FA5]+)\(", content.xpath("//div[@class='other']/div/text()[1]")[i])[0]
            total_floor = re.findall("(\d+)",content.xpath("//div[@class='other']/div/text()[1]")[i])[0]
            try:
                house_year = re.findall("(\d+)",content.xpath("//div[@class='other']/div/text()[2]")[i])[0]
            except Exception as e:
                house_year = ""
            price = content.xpath("//div[@class='col-3']/div/span/text()")[i]
            with open('链家北京租房.txt', 'a', encoding='utf-8') as f:
                f.write(area + ',' + title + ',' + room_type + ',' + square + ',' + position + ',' + detail_place + ',' + floor + ',' + total_floor + ',' + price + ',' + house_year + '\n')
            print('writing work has done!continue the next page')
    except Exception as e:
        print('ooops! connecting error, retrying.....')
        time.sleep(20)
        return get_house_info(area, url)


def main():
    print('start!')
    url = 'https://bj.lianjia.com/zufang'
    get_areas(url)


if __name__ == '__main__':
    main()
