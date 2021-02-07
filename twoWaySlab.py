#! /Users/tsuno/.pyenv/shims/python3
# -*- coding: utf-8 -*-
#
# generated by wxGlade 0.9.6 on Thu Oct 29 22:41:10 2020
#

import numpy, matplotlib
if matplotlib.__version__ < '2.2':
    raise ValueError("Minimum Matplotlib version required: 2.2")

#
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

#
import aijRc
#
import higashi

import wx
import os
import csv
import pandas as pd
import linecache
#from shutil import make_archive

# read from glade
import gui

# excel format
import openpyxl
from openpyxl.utils import get_column_letter # 列幅の指定 2020/05/27

# pdf export
import report

#import plate
# begin wxGlade: extracode
# end wxGlade


class MyFrame2(gui.MyFrame2):


    def OnChangeBound(self,event):
        IdBound = self.combo_box_bound.GetSelection()
        if IdBound == 0 or IdBound == 1: # ４辺固定
            image_data = "./images/4sideFix.jpg"
        elif IdBound == 2: # ３辺固定
            image_data = "./images/m3_1.jpg"
        elif IdBound == 3:# ３辺固定
            image_data = "./images/m3_2.jpg"
        elif IdBound == 4:# ２辺固定
            image_data = "./images/m2.jpg"
        elif IdBound ==5: # 4辺支持
            image_data = "./images/m4pin.jpg"
        elif IdBound ==6: # 3辺支持長辺支持
            image_data = "./images/m3_1pin.jpg"
        elif IdBound ==7: # ３辺固定短辺支持
            image_data = "./images/m3-1pin2.jpg"
        elif IdBound ==8: # 2辺固定2辺支持
            image_data = "./images/m2_2pin.jpg"
        elif IdBound ==9: # 2辺固定2辺支持
            image_data = "./images/m2_2pin2.jpg"
        elif IdBound ==10: # 2辺固定2辺支持
            image_data = "./images/m2_2pin3.jpg"
        else:
            dlg = wx.MessageDialog(self, 'Bound, input Erro',
                                   'Bound Error',
                                   wx.OK | wx.ICON_INFORMATION
                                   )
            dlg.ShowModal()
            dlg.Destroy()

        image = wx.Image(image_data)
        image = image.Scale(150,150,wx.IMAGE_QUALITY_BICUBIC)
        bitmap = image.ConvertToBitmap()
        wx.StaticBitmap(self.panel_bound, -1, bitmap, pos=(0,0) )

        #self.Clear_R()

    def OnClear(self,event):
        self.Clear_R()

    def Clear_R(self):
        #
        self.text_ctrl_mx1.SetValue('')
        self.text_ctrl_mx2.SetValue('')
        self.text_ctrl_my1.SetValue('')
        self.text_ctrl_my2.SetValue('')
        #
        self.text_ctrl_atx1.SetValue('')
        self.text_ctrl_atx2.SetValue('')
        self.text_ctrl_aty1.SetValue('')
        self.text_ctrl_aty2.SetValue('')
        #
        self.text_ctrl_rebarx1.SetValue('')
        self.text_ctrl_rebarx2.SetValue('')
        self.text_ctrl_rebary1.SetValue('')
        self.text_ctrl_rebary2.SetValue('')
        #
        self.text_ctrl_reqt.SetValue('')
        self.text_ctrl_tbyl.SetValue('')
        #
        self.text_ctrl_sfx1.SetValue('')
        self.text_ctrl_sfx2.SetValue('')
        self.text_ctrl_sfy1.SetValue('')
        self.text_ctrl_sfy2.SetValue('')
        self.text_ctrl_def.SetValue('')
        self.text_ctrl_dBySpan.SetValue('')
        #
        self.text_ctrl_lx1PitchOut.SetValue('')
        self.text_ctrl_lx2PitchOut.SetValue('')
        self.text_ctrl_ly1PitchOut.SetValue('')
        self.text_ctrl_ly2PitchOut.SetValue('')

    def ListShow(self):
        self.list_ctrl_output.DeleteAllItems();
        idTotal = self.text_ctrl_idView.GetValue()
        outFile = './db/rcslab.txt'
        index = 0
        if idTotal != '0':
            with open(outFile) as f:
                for row in csv.reader(f):
                    self.list_ctrl_output.InsertItem(index, index)
                    self.list_ctrl_output.SetItem(index, 0, str(int(index+1)))
                    self.list_ctrl_output.SetItem(index, 1, row[0])
                    self.list_ctrl_output.SetItem(index, 2, row[1])
                    index += 1
            f.close()
    # Remove Button

    def OnMove(self,event):

        # read
        idTotal = int(self.text_ctrl_idTotal.GetValue())
        id_move1 = int(self.text_ctrl_move1.GetValue())
        id_move2 = int(self.text_ctrl_move2.GetValue())
        data = []

        f = open('./db/rcslab.txt')

        for i in range(0,idTotal):
            line = f.readline()
            data.append(line)
        f.close()

        outFile = './db/rcslab.txt'
        fout = open(outFile, "w")

        if id_move2 < id_move1:
            for i in range(0,idTotal):
                if id_move2 == i+1:
                    fout.writelines(str(data[id_move1-1]))
                if id_move1 != i+1:
                    fout.writelines(str(data[i]))
        else:
            for i in range(0,idTotal):
                if id_move1 != i+1:
                    fout.writelines(str(data[i]))
                if id_move2 == i+1:
                    fout.writelines(str(data[id_move1-1]))

        fout.close()

        # output
        self.ListShow()

    def OnRemove(self,event):


        if self.text_ctrl_remove.GetValue() == '':
            print(self.text_ctrl_remove.GetValue())
            dlg = wx.MessageDialog(self, 'Pls, input Remove No.',
                                   'Error',
                                   wx.OK | wx.ICON_INFORMATION
                                   )
            dlg.ShowModal()
            dlg.Destroy()

        else:
            # read
            idTotal = int(self.text_ctrl_idTotal.GetValue())
            id_remove = int(self.text_ctrl_remove.GetValue())

            data = []

            f = open('./db/rcslab.txt')

            for i in range(0,idTotal):
                line = f.readline()
                data.append(line)
            f.close()

            outFile = './db/rcslab.txt'
            fout = open(outFile, "w")

            for i in range(0,idTotal):
                if id_remove != i+1:
                    fout.writelines(str(data[i]))

            fout.close()

            # output
            self.text_ctrl_idTotal.SetValue( str(idTotal-1))
            self.ListShow()

    def BarShow(self):
        # Read
        ####################
        # Rebar Arrangement show to display
        lx1Pitch  = self.text_ctrl_lx1Pitch.GetValue()
        lx2Pitch  = self.text_ctrl_lx2Pitch.GetValue()
        ly1Pitch  = self.text_ctrl_ly1Pitch.GetValue()
        ly2Pitch  = self.text_ctrl_ly2Pitch.GetValue()
        #
        lx1bar  = self.combo_box_lx1bar.GetValue()
        lx2bar  = self.combo_box_lx2bar.GetValue()
        ly1bar  = self.combo_box_ly1bar.GetValue()
        ly2bar  = self.combo_box_ly2bar.GetValue()
        # Display
        ####################
        self.text_ctrl_rebarx1.SetValue(lx1bar)
        self.text_ctrl_rebarx2.SetValue(lx2bar)
        self.text_ctrl_rebary1.SetValue(ly1bar)
        self.text_ctrl_rebary2.SetValue(ly2bar)
        #
        self.text_ctrl_lx1PitchOut.SetValue(lx1Pitch)
        self.text_ctrl_lx2PitchOut.SetValue(lx2Pitch)
        self.text_ctrl_ly1PitchOut.SetValue(ly1Pitch)
        self.text_ctrl_ly2PitchOut.SetValue(ly2Pitch)

    def SubShow(self,data):

        # substitute
        title = data[0]
        subTitle = data[1]
        # Slab Condition
        lx = data[2]
        ly = data[3]
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
        # Combo
        ind_bound = int( data[14] )
        ind_lx1 = int( data[15] )
        ind_lx2 = int( data[16] )
        ind_ly1 = int( data[17] )
        ind_ly2 = int( data[18] )
    #
        lx1bar  = data[19]
        lx2bar  = data[20]
        ly1bar  = data[21]
        ly2bar  = data[22]
        #
        """
        lx1bar = data[23]
        lx2bar = data[24]
        ly1bar = data[25]
        ly2bar = data[26]
        """
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
        ft     = data[39]

        # Title
        self.text_ctrl_title.SetValue(title)
        self.text_ctrl_subtitle.SetValue( subTitle)
        # Slab Condition
        self.text_ctrl_lx.SetValue( lx)
        self.text_ctrl_ly.SetValue( ly)
        self.text_ctrl_t.SetValue( t)
        self.text_ctrl_dt.SetValue( dt)
        self.text_ctrl_w.SetValue( w)
        self.text_ctrl_creep.SetValue( creep)
        self.text_ctrl_fc.SetValue( fc)
        #add
        self.text_ctrl_ft.SetValue( ft)
        #
        self.text_ctrl_gamma.SetValue( gamma)
        # Rebar Arrangement
        self.text_ctrl_lx1Pitch.SetValue( lx1Pitch)
        self.text_ctrl_lx2Pitch.SetValue( lx2Pitch)
        self.text_ctrl_ly1Pitch.SetValue( ly1Pitch)
        self.text_ctrl_ly2Pitch.SetValue( ly2Pitch)
        # Combo
        self.combo_box_bound.SetSelection( ind_bound)
        self.combo_box_lx1bar.SetSelection( ind_lx1)
        self.combo_box_lx2bar.SetSelection( ind_lx2)
        self.combo_box_ly1bar.SetSelection( ind_ly1)
        self.combo_box_ly2bar.SetSelection( ind_ly2)
        #
        """
        self.combo_box_lx1bar.SetValue( lx1bar)
        self.combo_box_lx2bar.SetValue( lx2bar)
        self.combo_box_ly1bar.SetValue( ly1bar)
        self.combo_box_ly2bar.SetValue( ly2bar)
        """
        #
        # output side
        self.text_ctrl_mx1.SetValue( mx1)
        self.text_ctrl_mx2.SetValue( mx2)
        self.text_ctrl_my1.SetValue( my1)
        self.text_ctrl_my2.SetValue( my2)
        #
        self.text_ctrl_atx1.SetValue( atx1)
        self.text_ctrl_atx2.SetValue( atx2)
        self.text_ctrl_aty1.SetValue( aty1)
        self.text_ctrl_aty2.SetValue( aty2)
        #
        self.text_ctrl_reqt.SetValue( reqt)
        self.text_ctrl_tbyl.SetValue( tbyl)
        #
        self.text_ctrl_sfx1.SetValue( sfx1)
        self.text_ctrl_sfx2.SetValue( sfx2)
        self.text_ctrl_sfy1.SetValue( sfy1)
        self.text_ctrl_sfy2.SetValue( sfy2)
        self.text_ctrl_def.SetValue( defl)
        self.text_ctrl_dBySpan.SetValue( dBySan)


    # Show Button
    ####################
    def OnShow(self,event):

        # read
        id_show = self.list_ctrl_output.GetFirstSelected() + 1
        line = linecache.getline('./db/rcslab.txt', id_show )
        data = line.split(', ')
        # test
        """
        print('Id = ' + str(id_show))
        print('line' + line)
        print(data)
        """
        # output
        idView = id_show
        self.text_ctrl_idView.SetValue(str(idView))
        self.SubShow(data)
        self.BarShow()
        linecache.clearcache()
        #
        # output remove and move
        self.text_ctrl_remove.SetValue(str(idView))
        self.text_ctrl_move1.SetValue(str(idView))
        #
        self.OnChangeBound(event)

    # Header Button
    def OnPre(self,event):

        # read
        id_show = int(self.text_ctrl_idView.GetValue())
        idTotal = int(self.text_ctrl_idTotal.GetValue())
        #
        if id_show == 1 :
            id_show_next = idTotal
        else:
            id_show_next = id_show - 1

        self.text_ctrl_idView.SetValue(str(id_show_next))
        line = linecache.getline('./db/rcslab.txt', id_show_next )
        data = line.split(', ')
        self.SubShow(data)
        self.BarShow()
        linecache.clearcache()

    def OnNext(self,event):

        # read
        id_show = int(self.text_ctrl_idView.GetValue())
        idTotal = int(self.text_ctrl_idTotal.GetValue())
        #
        if id_show >= idTotal :
            id_show_next = 1
        else:
            id_show_next = id_show + 1

        self.text_ctrl_idView.SetValue(str(id_show_next))
        line = linecache.getline('./db/rcslab.txt', id_show_next )
        data = line.split(', ')
        self.SubShow(data)
        self.BarShow()
        linecache.clearcache()

    def OnPlus(self,event):
        idTotal = int(self.text_ctrl_idTotal.GetValue())
        id_View = idTotal + 1
        self.text_ctrl_idView.SetValue(str(id_View))
        self.Clear_R()
        self.text_ctrl_title.SetValue('No.'+str(id_View))

    def OnQuit(self, event):
        self.Close()

    def rebarAt(self,index):
        return at

    # 断面算定: 引張鉄筋比以下としてRC基準で算定
    # md:    Design Bending moment = []
    # ast**: Area of reinforcement, mm2
    # ft:    permissible bar stress, N/mm2
    # t :    Slab thickness
    # dt:    distance tension bar center from the compressibe fiber
    def sf(self,md,astx1,astx2,asty1,asty2,ft,t,dt):
        #
        j = ( t - dt ) * 7.0 / 8.0
        # demand area
        reqatx1  = md[0] * 10**6 / ( ft * j )
        reqatx2  = md[1] * 10**6 / ( ft * j )
        reqaty1  = md[2] * 10**6 / ( ft * j )
        reqaty2  = md[3] * 10**6 / ( ft * j )
        #
        sfx1 = reqatx1/astx1
        sfx2 = reqatx2/astx2
        sfy1 = reqaty1/asty1
        sfy2 = reqaty2/asty2
        #
        return reqatx1,reqatx2,reqaty1,reqaty2,sfx1,sfx2,sfy1,sfy2

    # ４辺固定版、略算の設計応力と変形
    def side4(self,lx,ly,t,dt,w,creep,ec,wp):
        #
        # 応力
        wx = ly**4 / ( lx**4+ly**4 ) * w
        mx1 = 1.0/12.0 * wx * lx**2
        mx2 = 2.0/3.0  * mx1
        my1 = 1.0/24.0 * w  * lx**2
        my2 = 2.0/3.0  * my1

        # 必要スラブ
        lamda = ly/lx
        reqt = 0.02 * ( lamda - 0.7 ) / ( lamda - 0.6 )
        reqt = reqt * ( 1.0 + wp/10.0 + lx/10.0 ) * lx * 1000.0

        return mx1,my2,my1,my2,reqt


    def OnStore(self,event):

        #
        idView  = self.text_ctrl_idView.GetValue()
        idTotal = self.text_ctrl_idTotal.GetValue()
        #

        data = []
        # Title
        title = self.text_ctrl_title.GetValue()
        data.append(title)
        subTitle = self.text_ctrl_subtitle.GetValue()
        data.append(subTitle)
        # Slab Condition
        lx = self.text_ctrl_lx.GetValue()
        data.append(lx)
        ly = self.text_ctrl_ly.GetValue()
        data.append(ly)
        t  = self.text_ctrl_t.GetValue()
        data.append(t)
        dt = self.text_ctrl_dt.GetValue()
        data.append(dt)
        w  = self.text_ctrl_w.GetValue()
        data.append(w)
        creep  = self.text_ctrl_creep.GetValue()
        data.append(creep)
        fc     = self.text_ctrl_fc.GetValue()
        data.append(fc)
        gamma  = self.text_ctrl_gamma.GetValue()
        data.append(gamma)
        # Rebar Arrangement
        lx1Pitch  = self.text_ctrl_lx1Pitch.GetValue()
        data.append(lx1Pitch )
        lx2Pitch  = self.text_ctrl_lx2Pitch.GetValue()
        data.append(lx2Pitch )
        ly1Pitch  = self.text_ctrl_ly1Pitch.GetValue()
        data.append(ly1Pitch )
        ly2Pitch  = self.text_ctrl_ly2Pitch.GetValue()
        data.append(ly2Pitch )
        # Combo
        ind_bound = self.combo_box_bound.GetSelection()
        data.append(ind_bound)
        ind_lx1 = self.combo_box_lx1bar.GetSelection()
        data.append(ind_lx1)
        ind_lx2 = self.combo_box_lx2bar.GetSelection()
        data.append(ind_lx2)
        ind_ly1 = self.combo_box_ly1bar.GetSelection()
        data.append(ind_ly1)
        ind_ly2 = self.combo_box_ly2bar.GetSelection()
        data.append(ind_ly2)
        #
        lx1bar  = self.combo_box_lx1bar.GetValue()
        data.append(lx1bar )
        lx2bar  = self.combo_box_lx2bar.GetValue()
        data.append(lx2bar )
        ly1bar  = self.combo_box_ly1bar.GetValue()
        data.append(ly1bar )
        ly2bar  = self.combo_box_ly2bar.GetValue()
        data.append(ly2bar )


        #
        """
        lx1bar = lx1bar +'\n' + str(lx1Pitch)
        data.append(lx1bar)
        lx2bar = lx2bar +'\n' + str(lx2Pitch)
        data.append(lx2bar)
        ly1bar = ly1bar +'\n' + str(ly1Pitch)
        data.append(ly1bar)
        ly2bar = ly2bar +'\n' + str(ly2Pitch)
        """
        # 出力
        mx1 = self.text_ctrl_mx1.GetValue()
        data.append(mx1)
        mx2 = self.text_ctrl_mx2.GetValue()
        data.append(mx2)
        my1 = self.text_ctrl_my1.GetValue()
        data.append(my1)
        my2 = self.text_ctrl_my2.GetValue()
        data.append(my2)
        #
        atx1 = self.text_ctrl_atx1.GetValue()
        data.append(atx1)
        atx2 = self.text_ctrl_atx2.GetValue()
        data.append(atx2)
        aty1 = self.text_ctrl_aty1.GetValue()
        data.append(aty1)
        aty2 = self.text_ctrl_aty2.GetValue()
        data.append(aty2)
        #
        reqt = self.text_ctrl_reqt.GetValue()
        data.append(reqt)
        tbyl = self.text_ctrl_tbyl.GetValue()
        data.append(tbyl)
        #
        sfx1   = self.text_ctrl_sfx1.GetValue()
        data.append(sfx1  )
        sfx2   = self.text_ctrl_sfx2.GetValue()
        data.append(sfx2  )
        sfy1   = self.text_ctrl_sfy1.GetValue()
        data.append(sfy1  )
        sfy2   = self.text_ctrl_sfy2.GetValue()
        data.append(sfy2  )
        defl   = self.text_ctrl_def.GetValue()
        data.append(defl  )
        dBySan = self.text_ctrl_dBySpan.GetValue()
        data.append(dBySan)

        # add
        ft     = self.text_ctrl_ft.GetValue()
        data.append(ft)


        ####################
        outFile = './db/rcslab.txt'
        if idTotal == '0':
            fout = open(outFile, "w")
        else:
            fout = open(outFile, "a")

        for i in range(len(data)):
            fout.writelines(str(data[i]))
            fout.writelines(', ')
        fout.writelines('\n')
        fout.close()

        #idView_next = int(idView) + 1
        idView_next = int(idTotal) + 2
        idTotal_next = int(idTotal) + 1
        self.text_ctrl_idView.SetValue(str(idView_next))
        self.text_ctrl_idTotal.SetValue(str(idTotal_next))
        #
        self.text_ctrl_title.SetValue('No.'+str(idView_next))
        #

        self.ListShow()
        self.Clear_R()

        # Pandas Check
        #df = pd.read_csv("./db/rcslab.txt")
        #print(df)


    # Make report
    def OnReport(self,event):

        with wx.FileDialog(self, "Save Pdf File", wildcard="Input File (*.pdf)|*.pdf",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
            # save the current contents in the file
            pathname = fileDialog.GetPath() # pdf file

            try:
                num = int( self.text_ctrl_idTotal.GetValue() )
                title = "none"
                #for i in range(0,num):
                obj = report.Report()
                obj.create_pdf(num,pathname,title)

            except IOError:
                wx.LogError("Cannot save current data in file '%s'." % pathname)


    # Export Csv Sheet
    ########################################################################
    def OnExport(self,event):
        with wx.FileDialog(self, "Save Csv File", wildcard="Output File (*.csv)|*.csv",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
            # save the current contents in the file
            pathname = fileDialog.GetPath()
            try:
                input_path = './db/rcslab.txt'
                with open(input_path) as f:
                    s = f.read()
                with open(pathname, 'w') as file:
                    file.write(s)
            except IOError:
                wx.LogError("Cannot save current data in file '%s'." % pathname)


    # Import Csv Sheet
    ########################################################################
    def OnImport(self,event):
        with wx.FileDialog(self, "Open Input csv file", wildcard="Input files (*.csv)|*.csv",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
            # save the current contents in the file
            pathname = fileDialog.GetPath()
            try:
                with open(pathname) as f:
                    s = f.read()

                input_path = './db/rcslab.txt'
                with open(input_path, 'w') as file:
                    file.write(s)

                # Listの読み込み
                line_count = 0
                with open(input_path) as f:
                    for line in f:
                        line_count += 1
                print("Success Import!!!")

                self.text_ctrl_idTotal.SetValue(str(line_count))
                #self.text_ctrl_idView.SetValue(str(1))
                self.ListShow()
                self.list_ctrl_output.Select(0)
                self.OnShow(event)

            except IOError:
                wx.LogError("Cannot import current data in file '%s'." % pathname)



    def OnCal(self, event):  # wxGlade: MyFrame.<event_handler>

        # Read Data
        ########################################################################
        # Title
        title = self.text_ctrl_title.GetValue()
        subTitle = self.text_ctrl_subtitle.GetValue()
        # Slab Condition
        lx = float( self.text_ctrl_lx.GetValue() )
        ly = float( self.text_ctrl_ly.GetValue() )
        t  = float( self.text_ctrl_t.GetValue() )
        dt = float( self.text_ctrl_dt.GetValue() )
        w  = float( self.text_ctrl_w.GetValue() )
        creep  = float( self.text_ctrl_creep.GetValue() )
        fc     = float( self.text_ctrl_fc.GetValue() )
        gamma  = float( self.text_ctrl_gamma.GetValue() )
        # Rebar Arrangement
        lx1Pitch  = float( self.text_ctrl_lx1Pitch.GetValue() )
        lx2Pitch  = float( self.text_ctrl_lx2Pitch.GetValue() )
        ly1Pitch  = float( self.text_ctrl_ly1Pitch.GetValue() )
        ly2Pitch  = float( self.text_ctrl_ly2Pitch.GetValue() )
        # Combo
        ind_bound = self.combo_box_bound.GetSelection()
        ind_lx1 = self.combo_box_lx1bar.GetSelection()
        ind_lx2 = self.combo_box_lx2bar.GetSelection()
        ind_ly1 = self.combo_box_ly1bar.GetSelection()
        ind_ly2 = self.combo_box_ly2bar.GetSelection()
        #
        lx1bar  = self.combo_box_lx1bar.GetValue()
        lx2bar  = self.combo_box_lx2bar.GetValue()
        ly1bar  = self.combo_box_ly1bar.GetValue()
        ly2bar  = self.combo_box_ly2bar.GetValue()
        #
        """
        lx1bar = lx1bar +'\n' + str(lx1Pitch)
        lx2bar = lx2bar +'\n' + str(lx2Pitch)
        ly1bar = ly1bar +'\n' + str(ly1Pitch)
        ly2bar = ly2bar +'\n' + str(ly2Pitch)
        """
        # Test
#        print(astx1,astx2,asty1,asty2)
        """
        print(title)
        print(subTitle)
        print(lx,ly)
        print(t,dt)
        print(w)
        print(creep)
        print(fc,gamma)
        print(lx1Pitch)
        print(lx2Pitch)
        print(ly1Pitch)
        print(ly2Pitch)
        print(ind_bound)
        print(ind_lx1)
        print(ind_lx2)
        print(ind_ly1)
        print(ind_ly2)
        """
        # add
        ft     = float( self.text_ctrl_ft.GetValue() )
        #
        # Preparation Cal
        ########################################################################
        wp = w - gamma * t / 1000.0
        # 配筋断面積
        astx1 = aijRc.Aij_rc_set().Ra_p(lx1bar,lx1Pitch)
        astx2 = aijRc.Aij_rc_set().Ra_p(lx2bar,lx2Pitch)
        asty1 = aijRc.Aij_rc_set().Ra_p(ly1bar,ly1Pitch)
        asty2 = aijRc.Aij_rc_set().Ra_p(ly2bar,ly2Pitch)
        # Solve
        ########################################################################
        ec = aijRc.Aij_rc_set().Ec(fc,gamma-1.0) # Young Modulus
        #ft = 195.0
        nu = 0.0
        # 応力
        IdBound = ind_bound
        #self.combo_box_bound.GetSelection()
        if IdBound == 0: # ４辺固定
            md = self.side4(lx,ly,t,dt,w,creep,ec,wp)
        else:
            md = higashi.Higashi().solve(IdBound,lx,ly,t,w,creep,ec,nu,5,5)
        print(md)
        # 断面算定
        sfRes = self.sf(md,astx1,astx2,asty1,asty2,ft,t,dt)

        # output to gui
        ########################################################################
        #
        # 出力
        self.text_ctrl_mx1.SetValue("{:.1f}".format(md[0]))
        self.text_ctrl_mx2.SetValue("{:.1f}".format(md[1]))
        self.text_ctrl_my1.SetValue("{:.1f}".format(md[2]))
        self.text_ctrl_my2.SetValue("{:.1f}".format(md[3]))

        self.text_ctrl_atx1.SetValue("{:.0f}".format(sfRes[0]))
        self.text_ctrl_atx2.SetValue("{:.0f}".format(sfRes[1]))
        self.text_ctrl_aty1.SetValue("{:.0f}".format(sfRes[2]))
        self.text_ctrl_aty2.SetValue("{:.0f}".format(sfRes[3]))

        self.text_ctrl_sfx1.SetValue("{:.2f}".format(sfRes[4]))
        self.text_ctrl_sfx2.SetValue("{:.2f}".format(sfRes[5]))
        self.text_ctrl_sfy1.SetValue("{:.2f}".format(sfRes[6]))
        self.text_ctrl_sfy2.SetValue("{:.2f}".format(sfRes[7]))

        if IdBound == 0: # 4辺固定
            self.text_ctrl_def.SetValue("-")
            self.text_ctrl_dBySpan.SetValue("-")
            self.text_ctrl_reqt.SetValue("{:.0f}".format(md[4]))
            self.text_ctrl_tbyl.SetValue('1/'+"{:.0f}".format(lx*1000/t))
        else:
            self.text_ctrl_def.SetValue("{:.0f}".format(md[4]))
            self.text_ctrl_dBySpan.SetValue('1/'+"{:.0f}".format(lx*1000/md[4]))
            self.text_ctrl_reqt.SetValue("-")
            self.text_ctrl_tbyl.SetValue('1/'+"{:.0f}".format(lx*1000/t))

        self.BarShow()
        #

        """
        obj = plate.slab()
        lamda = ly/lx
        print(obj.mm3(lamda))
        """

class MyFrame(gui.MyFrame):

    def OnQuit(self, event):
        self.Close()

    def OnRcslab(self,event):
        MyFrame2(None, wx.ID_ANY, "").Show()
        return True

    def OnTest(self,event):
        MyFrame_test(None, wx.ID_ANY, "").Show()
        return True
"""
class MyFrame_test(gui.MyFrame_test):
    print('h')
"""

# end of class MyFrame

# end of class MyFrame

class MyApp(wx.App):

    def OnInit(self):
        new_dir_path = './db'
        try:
            os.mkdir(new_dir_path)
        except FileExistsError:
            pass
        #self.frame = MyFrame(None, wx.ID_ANY, "")
        self.frame = MyFrame2(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        #MyFrame2(None, wx.ID_ANY, "").Show()
        return True
# 上のふたつの#をとって、frameをコメントアウトすれば、ホーム画面を立ち上げる。
# end of class MyApp

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
