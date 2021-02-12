#! /Users/tsuno/.pyenv/shims/python3
# -*- coding:utf-8 -*-
import os, sys
#import Image
#import urllib2
#from cStringIO import StringIO


#zipアーカイブからファイルを読み込むため。通常は必要ないはず。
#sys.path.insert(0, 'reportlab.zip')

import reportlab
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.units import cm

#
import linecache
#

class Report():

    def __init__(self):
        #self.FONT_NAME = "Helvetica"
        self.FONT_NAME = "GenShinGothic"
        GEN_SHIN_GOTHIC_MEDIUM_TTF = "./fonts/GenShinGothic-Monospace-Medium.ttf"
        # フォント登録
        pdfmetrics.registerFont(TTFont('GenShinGothic', GEN_SHIN_GOTHIC_MEDIUM_TTF))
        #font_size = 20
        #c.setFont('GenShinGothic', font_size)

    ########################################################################
    # 文字と画像を配置
    def create_row(self,c, index, data):
        #y_shift = -240 * index
        y_shift = -180 * index
        c.setFont(self.FONT_NAME, 9)
        """
        for i in range(0,len(data)):
            # txt
            c.drawString(300, 720-(i-1)*10 + y_shift, data[i])
        """
        c.drawString(55, self.ypos(0,y_shift), data[0].encode('utf-8'))
        c.drawString(55, self.ypos(1,y_shift), data[1].encode('utf-8'))

        # Slab Condition
        lx = "{:.2f}".format(float(data[2]))
        ly = "{:.2f}".format(float(data[3]))
        t  = data[4]
        dt = data[5]
        w  = data[6]
        creep  = data[7]
        fc     = data[8]
        gamma  = data[9]
        # Rebar Arrangement
        lx1Pitch  = data[10]
        lx2Pitch  = data[11]
        ly1Pitch  = data[12]
        ly2Pitch  = data[13]
        #
        ind_bound = int( data[14] )
        # Combo
        lx1bar  = data[19]
        lx2bar  = data[20]
        ly1bar  = data[21]
        ly2bar  = data[22]
        # 出力
        mx1 = data[23]
        mx2 = data[24]
        my1 = data[25]
        my2 = data[26]
        #
        atx1 = data[27]
        atx2 = data[28]
        aty1 = data[29]
        aty2 = data[30]
        #
        reqt = data[31]
        tbyl = data[32]
        #
        sfx1   = data[33]
        sfx2   = data[34]
        sfy1   = data[35]
        sfy2   = data[36]
        defl   = data[37]
        dBySan = data[38]

        # add
        ft     = data[39]

        # Design Condition
        c.drawString(260, self.ypos(0,y_shift),
                     "Lx = "
                     )
        c.drawString(290, self.ypos(0,y_shift),
                     lx + "m, ")
        c.drawString(350, self.ypos(0,y_shift),
                     "Ly = "
                     )
        c.drawString(380, self.ypos(0,y_shift),
                     ly + " m"\
                     )
        #
        c.drawString(260, self.ypos(1,y_shift),
                     "t = "
                     )
        c.drawString(290, self.ypos(1,y_shift),
                     t + " mm, "\
                     )
        c.drawString(350, self.ypos(1,y_shift),
                     "dt = "
                     )
        c.drawString(380, self.ypos(1,y_shift),
                     dt + " mm"\
                     )
        #
        c.drawString(260, self.ypos(2,y_shift),
                     "fc = "
                     )
        c.drawString(290, self.ypos(2,y_shift),
                     fc + " N/mm2, "\
                     )
        c.drawString(350, self.ypos(2,y_shift),
                     "( γ = " + gamma + " kN/m3)"\
                     )
        c.drawString(430, self.ypos(2,y_shift),
                     "ft = " + ft + " N/mm2"\
                     )
        #
        c.drawString(260, self.ypos(3,y_shift),
                     "w = " 
                     )
        c.drawString(290, self.ypos(3,y_shift),
                     w + " kN/m2, "\
                     )
        c.drawString(350, self.ypos(3,y_shift),
                     "K = "
                     )
        c.drawString(380, self.ypos(3,y_shift),
                     creep \
                     )

        # Result for the stress
        """
        c.drawString(330, self.ypos(5,y_shift), "S,Ext. " )
        c.drawString(380, self.ypos(5,y_shift), "S,Cen. " )
        c.drawString(430, self.ypos(5,y_shift), "L,Ext. " )
        c.drawString(480, self.ypos(5,y_shift), "L,Cen. " )
        """
        c.drawString(330, self.ypos(5,y_shift), "短辺端部 " )
        c.drawString(380, self.ypos(5,y_shift), "短辺中央 " )
        c.drawString(430, self.ypos(5,y_shift), "長辺端部" )
        c.drawString(480, self.ypos(5,y_shift), "長辺中央" )

        c.drawString(260, self.ypos(6,y_shift), "M " )
        c.drawString(260, self.ypos(7,y_shift), "at " )
        c.drawString(260, self.ypos(8,y_shift), "上端" )
        c.drawString(260, self.ypos(9,y_shift), "下端" )
        c.drawString(260, self.ypos(10,y_shift), "検定比" )

        c.drawString(290, self.ypos(6,y_shift), "kN.m/m " )
        c.drawString(290, self.ypos(7,y_shift), "mm2 " )

        c.drawString(330, self.ypos(6,y_shift), mx1)
        c.drawString(380, self.ypos(6,y_shift), mx2 )
        c.drawString(430, self.ypos(6,y_shift), my1 )
        c.drawString(480, self.ypos(6,y_shift), my2 )

        c.drawString(330, self.ypos(7,y_shift), atx1 )
        c.drawString(380, self.ypos(7,y_shift), atx2 )
        c.drawString(430, self.ypos(7,y_shift), aty1 )
        c.drawString(480, self.ypos(7,y_shift), aty2 )

        c.drawString(330, self.ypos(8,y_shift), lx1bar + "@" + lx1Pitch )
        c.drawString(380, self.ypos(9,y_shift), lx2bar + "@" + lx2Pitch )
        c.drawString(430, self.ypos(8,y_shift), ly1bar + "@" + ly1Pitch )
        c.drawString(480, self.ypos(9,y_shift), ly2bar + "@" + ly2Pitch )

        c.drawString(330, self.ypos(10,y_shift), sfx1 )
        c.drawString(380, self.ypos(10,y_shift), sfx2 )
        c.drawString(430, self.ypos(10,y_shift), sfy1 )
        c.drawString(480, self.ypos(10,y_shift), sfy2 )

        # Deflection
        c.drawString(260, self.ypos(12,y_shift), "変形")
        c.drawString(260, self.ypos(13,y_shift), "Req.t = ")
        c.drawString(300, self.ypos(13,y_shift), reqt + " mm,")
        c.drawString(350, self.ypos(13,y_shift), "( t/Lx = " )
        c.drawString(430, self.ypos(13,y_shift), tbyl + ")")

        c.drawString(260, self.ypos(14,y_shift), "δ = " )
        c.drawString(300, self.ypos(14,y_shift), defl + " mm, ")
        c.drawString(350, self.ypos(14,y_shift), "δ/Lx = " )
        c.drawString(430, self.ypos(14,y_shift), dBySan)
        """
        for i in range(2,len(data)):
            # txt
            #c.drawString(300, 720-(i-1)*10 + y_shift, data[i])
            c.drawString(500, 720-(i-1)*10 + y_shift, data[i])
        """
        # png
        imagefile=self.boundimage(ind_bound)
        c.drawImage(imagefile, 70,  y_shift + 540, width=5*cm , preserveAspectRatio=True)

    def boundimage(self,index):
        if index == 0 or index == 1: # ４辺固定
            image_data = "./images/4sideFix.jpg"
        elif index == 2: # ３辺固定
            image_data = "./images/m3_1.jpg"
        elif index == 3:# ３辺固定
            image_data = "./images/m3_2.jpg"
        elif index == 4:# ２辺固定
            image_data = "./images/m2.jpg"
        elif index ==5: # 4辺支持
            image_data = "./images/m4pin.jpg"
        elif index ==6: # 3辺支持長辺支持
            image_data = "./images/m3_1pin.jpg"
        elif index ==7: # ３辺固定短辺支持
            image_data = "./images/m3-1pin2.jpg"
        elif index ==8: # 2辺固定2辺支持
            image_data = "./images/m2_2pin.jpg"
        elif index ==9: # 2辺固定2辺支持
            image_data = "./images/m2_2pin2.jpg"
        elif index ==10: # 2辺固定2辺支持
            image_data = "./images/m2_2pin3.jpg"
        elif IdBound ==11: # 短辺1辺固定3辺支持
            image_data = "./images/m1-3pin1.jpg"
        elif IdBound ==12: # 長辺1辺固定3辺支持
            image_data = "./images/m1-3pin2.jpg"
        else:
            print("Error report/def")

        return image_data

    def ypos(self,ipos,y_shift):
        return 730-(ipos-1)*10 + y_shift

    ########################################################################
    # pdfの作成
    def print_page(self, c, index, nCase):


        #タイトル描画
        c.setFont(self.FONT_NAME, 20)
        #c.drawString(50, 795, u"Design of the twoway slab")
        c.drawString(50, 795, u"スラブの設計")

        #グリッドヘッダー設定
        xlist = [40, 240, 560]
        ylist = [760, 780]
        c.grid(xlist, ylist)

        #sub title
        c.setFont(self.FONT_NAME, 12)
        c.drawString(55, 765, u"符号")
        c.drawString(260, 765, u"設計")

        #データを描画
        ########################################################################
        #for i, data in range(0,int(nCase)):

        for i in range(0,nCase):
#            f = open(inputf[index+i],'r')

            """
            tmpData = []
            while True:
                line = f.readline()
                if line:
                    if line != '\n':
                        tmpData.append(line.replace('\n',''))
                    else:
                        tmpData.append('')
                else:
                    break
            """
            line = linecache.getline('./db/rcslab.txt', index+i+1 )
            data = line.split(', ')
            linecache.clearcache()
            #f.close()
            #data = tmpData
            self.create_row( c, i, data )

        #最後にグリッドを更新
        #ylist = [40,  280,  520,  760]
        #ylist = [40,  160, 280, 400, 520, 640, 760]
        ylist = [40,  220,  400, 580, 760]
        c.grid(xlist, ylist[4 - nCase:])
        #ページを確定
        c.showPage()

    ########################################################################
    # pdfの作成
    def print_head(self, c , title):

        #title = 'Sample Project'

        #タイトル描画
        c.setFont(self.FONT_NAME, 20)
        c.drawString(50, 795, title)

        #sub title
        c.setFont(self.FONT_NAME, 12)

        #データを描画
        ########################################################################
        inputf = './db/input.txt'
        f = open(inputf,'r')
        tmpData = []
        while True:
            line = f.readline()
            if line:
                if line != '\n':
                    tmpData.append(line.replace('\n',''))
                else:
                    tmpData.append('')
            else:
                break
        f.close()
        data = tmpData
        #c.setFont(self.FONT_NAME, 9)
        for i in range(0,len(data)):
            # txt
            c.drawString(55, 720-(i-1)*14, data[i])
        """
        # Model Diagram
        imagefile = './db/model.png'
        c.drawImage(imagefile, 60,  -300, width=18*cm , preserveAspectRatio=True)
        """
        #ページを確定
        c.showPage()
    ########################################################################
    # whole control
    def create_pdf(self, dataNum, pdfFile, title):

        # Parameter -------
        # inputf   : path to text file
        # imagefile: path to png file
        # pdfFile  : name of making pdf file

        #フォントファイルを指定して、フォントを登録
        #folder = os.path.dirname(reportlab.__file__) + os.sep + 'fonts'
        #pdfmetrics.registerFont(TTFont(FONT_NAME, os.path.join(folder, 'ipag.ttf')))
        #出力するPDFファイル
        c = canvas.Canvas(pdfFile)

        # ページ数
        ########################################################################
        #dataNum = len(inputf)
        numPage = dataNum // 4
        numMod = dataNum % 4
        #print(numPage,numMod)
        if numMod >= 1:
            numPage = numPage + 1

        # pdfの作成
        ########################################################################
        #self.print_head( c , title)

        for i in range(0,numPage):
            index = 4*i # index: 参照データのインデックス
            if numPage == 1:
                self.print_page( c, index, dataNum)
            elif i != numPage-1 and numPage != 1:
                self.print_page( c, index, 4)
            else:
                if numMod != 0:
                    self.print_page( c, index, numMod)
                else:
                    self.print_page( c, index, 4 )

        #pdfファイル生成
        ########################################################################
        c.save()
        print ("repot.py is Okay!!.")

########################################################################
# END CLASS


"""
########################################################################
# test script

pathname = "./test.pdf"
obj = Report()
# テキストの読み込み
########################################################################
inputf = []
inputf.append("./db/rcslab.txt")
inputf.append("./db/rcslab.txt")
inputf.append("./db/rcslab.txt")
inputf.append("./db/rcslab.txt")
inputf.append("./db/rcslab.txt")
inputf.append("./db/rcslab.txt")

title = "sample"

obj.create_pdf(3,pathname,title)
"""
