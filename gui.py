#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# generated by wxGlade 0.9.6 on Fri Feb 12 14:16:43 2021
#

import wx

# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode
# end wxGlade


class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        
        # Menu Bar
        self.frame_menubar = wx.MenuBar()
        wxglade_tmp_menu = wx.Menu()
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "&NEW\tCtrl+N", "")
        self.Bind(wx.EVT_MENU, self.OnOpen, id=item.GetId())
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "E&xit", "")
        self.Bind(wx.EVT_MENU, self.OnQuit, id=item.GetId())
        self.frame_menubar.Append(wxglade_tmp_menu, "&File\tCtrl+F")
        self.SetMenuBar(self.frame_menubar)
        # Menu Bar end
        self.panel_1 = wx.Panel(self, wx.ID_ANY)
        self.panel_2 = wx.Panel(self, wx.ID_ANY)
        self.button_1 = wx.Button(self.panel_2, wx.ID_ANY, "Design of two way slab")
        self.button_2 = wx.Button(self.panel_2, wx.ID_ANY, "Design of one way slab")
        self.panel_3 = wx.Panel(self, wx.ID_ANY)
        self.panel_4 = wx.Panel(self, wx.ID_ANY)

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.OnRcslab, self.button_1)
        self.Bind(wx.EVT_BUTTON, self.OnTest, self.button_2)
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: MyFrame.__set_properties
        self.SetTitle("RcAji")
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: MyFrame.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_1.Add(self.panel_1, 1, wx.EXPAND, 0)
        sizer_2.Add(self.button_1, 0, wx.EXPAND, 0)
        sizer_2.Add(self.button_2, 0, 0, 0)
        label_1 = wx.StaticText(self.panel_2, wx.ID_ANY, "label_1")
        sizer_2.Add(label_1, 0, 0, 0)
        label_2 = wx.StaticText(self.panel_2, wx.ID_ANY, "label_2")
        sizer_2.Add(label_2, 0, 0, 0)
        self.panel_2.SetSizer(sizer_2)
        sizer_1.Add(self.panel_2, 1, wx.EXPAND, 0)
        sizer_1.Add(self.panel_3, 1, wx.EXPAND, 0)
        sizer_1.Add(self.panel_4, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade

    def OnOpen(self, event):  # wxGlade: MyFrame.<event_handler>
        print("Event handler 'OnOpen' not implemented!")
        event.Skip()

    def OnQuit(self, event):  # wxGlade: MyFrame.<event_handler>
        print("Event handler 'OnQuit' not implemented!")
        event.Skip()

    def OnRcslab(self, event):  # wxGlade: MyFrame.<event_handler>
        print("Event handler 'OnRcslab' not implemented!")
        event.Skip()

    def OnTest(self, event):  # wxGlade: MyFrame.<event_handler>
        print("Event handler 'OnTest' not implemented!")
        event.Skip()

# end of class MyFrame

class MyFrame2(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame2.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.panel_6 = wx.Panel(self, wx.ID_ANY)
        self.panel_7 = wx.Panel(self.panel_6, wx.ID_ANY)
        self.button_3 = wx.Button(self.panel_7, wx.ID_ANY, "Pre")
        self.button_8 = wx.Button(self.panel_7, wx.ID_ANY, "Next")
        self.text_ctrl_idView = wx.TextCtrl(self.panel_7, wx.ID_ANY, "1", style=wx.TE_READONLY)
        self.text_ctrl_idTotal = wx.TextCtrl(self.panel_7, wx.ID_ANY, "0", style=wx.TE_READONLY)
        self.button_9 = wx.Button(self.panel_7, wx.ID_ANY, "+")
        self.panel_8 = wx.Panel(self.panel_6, wx.ID_ANY)
        self.panel_10 = wx.Panel(self.panel_8, wx.ID_ANY)
        self.text_ctrl_title = wx.TextCtrl(self.panel_10, wx.ID_ANY, "No.1")
        self.text_ctrl_subtitle = wx.TextCtrl(self.panel_10, wx.ID_ANY, "S20")
        self.text_ctrl_lx = wx.TextCtrl(self.panel_10, wx.ID_ANY, "3.0", style=wx.TE_RIGHT)
        self.text_ctrl_ly = wx.TextCtrl(self.panel_10, wx.ID_ANY, "4.0", style=wx.TE_RIGHT)
        self.text_ctrl_t = wx.TextCtrl(self.panel_10, wx.ID_ANY, "200", style=wx.TE_RIGHT)
        self.text_ctrl_dt = wx.TextCtrl(self.panel_10, wx.ID_ANY, "40", style=wx.TE_RIGHT)
        self.text_ctrl_w = wx.TextCtrl(self.panel_10, wx.ID_ANY, "10.0", style=wx.TE_RIGHT)
        self.combo_box_bound = wx.ComboBox(self.panel_10, wx.ID_ANY, choices=[u"４辺固定（略算）", u"４辺固定（精算）", u"３辺固定長辺自由", u"３辺固定短辺自由", u"２隣辺固定２辺自由", u"４辺支持", u"３辺固定長辺支持", u"３辺固定短辺支持", u"２隣辺固定２辺支持", u"２対辺固定長辺２辺支持", u"２対辺固定短辺２辺支持"], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.text_ctrl_creep = wx.TextCtrl(self.panel_10, wx.ID_ANY, "16", style=wx.TE_RIGHT)
        self.text_ctrl_fc = wx.TextCtrl(self.panel_10, wx.ID_ANY, "21", style=wx.TE_RIGHT)
        self.text_ctrl_ft = wx.TextCtrl(self.panel_10, wx.ID_ANY, "195")
        self.text_ctrl_gamma = wx.TextCtrl(self.panel_10, wx.ID_ANY, "24", style=wx.TE_RIGHT)
        self.combo_box_lx1bar = wx.ComboBox(self.panel_10, wx.ID_ANY, choices=["D10", "D10+D13", "D13", "D13+D16", "D16", "D16+D19", "D19", "D19+D22", "D22"], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.text_ctrl_lx1Pitch = wx.TextCtrl(self.panel_10, wx.ID_ANY, "200", style=wx.TE_RIGHT)
        self.combo_box_lx2bar = wx.ComboBox(self.panel_10, wx.ID_ANY, choices=["D10", "D10+D13", "D13", "D13+D16", "D16", "D16+D19", "D19", "D19+D22", "D22"], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.text_ctrl_lx2Pitch = wx.TextCtrl(self.panel_10, wx.ID_ANY, "200", style=wx.TE_RIGHT)
        self.combo_box_ly1bar = wx.ComboBox(self.panel_10, wx.ID_ANY, choices=["D10", "D10+D13", "D13", "D13+D16", "D16", "D16+D19", "D19", "D19+D22", "D22"], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.text_ctrl_ly1Pitch = wx.TextCtrl(self.panel_10, wx.ID_ANY, "200", style=wx.TE_RIGHT)
        self.combo_box_ly2bar = wx.ComboBox(self.panel_10, wx.ID_ANY, choices=["D10", "D10+D13", "D13", "D13+D16", "D16", "D16+D19", "D19", "D19+D22", "D22"], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.text_ctrl_ly2Pitch = wx.TextCtrl(self.panel_10, wx.ID_ANY, "200", style=wx.TE_RIGHT)
        self.panel_11 = wx.Panel(self.panel_8, wx.ID_ANY)
        self.panel_boundtmp = wx.Panel(self.panel_11, wx.ID_ANY)
        self.panel_bound = wx.Panel(self.panel_boundtmp, wx.ID_ANY)
        self.panel_5 = wx.Panel(self.panel_boundtmp, wx.ID_ANY)
        self.text_ctrl_mx1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_mx2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_my1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_my2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_atx1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_atx2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_aty1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_aty2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_rebarx1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_rebarx2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_rebary1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_rebary2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_lx1PitchOut = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_lx2PitchOut = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_ly1PitchOut = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_ly2PitchOut = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_sfx1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_sfx2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_sfy1 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_sfy2 = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_reqt = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_tbyl = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_def = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.text_ctrl_dBySpan = wx.TextCtrl(self.panel_11, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_RIGHT)
        self.panel_9 = wx.Panel(self.panel_6, wx.ID_ANY)
        self.button_6 = wx.Button(self.panel_9, wx.ID_ANY, "Store")
        self.button_7 = wx.Button(self.panel_9, wx.ID_ANY, "Cal")
        self.panel_14 = wx.Panel(self.panel_6, wx.ID_ANY)
        self.list_ctrl_output = wx.ListCtrl(self.panel_14, wx.ID_ANY, style=wx.LC_HRULES | wx.LC_REPORT | wx.LC_VRULES)
        self.panel_15 = wx.Panel(self.panel_14, wx.ID_ANY)
        self.panel_17 = wx.Panel(self.panel_15, wx.ID_ANY)
        self.text_ctrl_remove = wx.TextCtrl(self.panel_15, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.button_11 = wx.Button(self.panel_15, wx.ID_ANY, "Remove")
        self.text_ctrl_move1 = wx.TextCtrl(self.panel_15, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.text_ctrl_move2 = wx.TextCtrl(self.panel_15, wx.ID_ANY, "")
        self.button_15 = wx.Button(self.panel_15, wx.ID_ANY, "Move")
        self.panel_16 = wx.Panel(self.panel_15, wx.ID_ANY)
        self.button_14 = wx.Button(self.panel_15, wx.ID_ANY, "Report")
        self.button_12 = wx.Button(self.panel_15, wx.ID_ANY, "Import")
        self.button_13 = wx.Button(self.panel_15, wx.ID_ANY, "Export")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.OnPre, self.button_3)
        self.Bind(wx.EVT_BUTTON, self.OnNext, self.button_8)
        self.Bind(wx.EVT_BUTTON, self.OnPlus, self.button_9)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_title)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_subtitle)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_lx)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_ly)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_t)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_dt)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_w)
        self.Bind(wx.EVT_COMBOBOX, self.OnChangeBound, self.combo_box_bound)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_creep)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_fc)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_gamma)
        self.Bind(wx.EVT_COMBOBOX, self.OnClear, self.combo_box_lx1bar)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_lx1Pitch)
        self.Bind(wx.EVT_COMBOBOX, self.OnClear, self.combo_box_lx2bar)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_lx2Pitch)
        self.Bind(wx.EVT_COMBOBOX, self.OnClear, self.combo_box_ly1bar)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_ly1Pitch)
        self.Bind(wx.EVT_COMBOBOX, self.OnClear, self.combo_box_ly2bar)
        self.Bind(wx.EVT_TEXT, self.OnClear, self.text_ctrl_ly2Pitch)
        self.Bind(wx.EVT_BUTTON, self.OnStore, self.button_6)
        self.Bind(wx.EVT_BUTTON, self.OnCal, self.button_7)
        self.Bind(wx.EVT_LIST_DELETE_ITEM, self.OnRemove, self.list_ctrl_output)
        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnShow, self.list_ctrl_output)
        self.Bind(wx.EVT_BUTTON, self.OnRemove, self.button_11)
        self.Bind(wx.EVT_BUTTON, self.OnMove, self.button_15)
        self.Bind(wx.EVT_BUTTON, self.OnReport, self.button_14)
        self.Bind(wx.EVT_BUTTON, self.OnImport, self.button_12)
        self.Bind(wx.EVT_BUTTON, self.OnExport, self.button_13)
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: MyFrame2.__set_properties
        self.SetTitle("Desighn of RC Slab")
        _icon = wx.NullIcon
        _icon.CopyFromBitmap(wx.Bitmap("./icons/twoWay_Icon.jpg", wx.BITMAP_TYPE_ANY))
        self.SetIcon(_icon)
        self.text_ctrl_idView.SetMinSize((40, 22))
        self.text_ctrl_idTotal.SetMinSize((40, 22))
        self.text_ctrl_title.SetMinSize((150, 22))
        self.text_ctrl_subtitle.SetMinSize((150, 22))
        self.text_ctrl_lx.SetMinSize((60, 22))
        self.text_ctrl_ly.SetMinSize((60, 22))
        self.text_ctrl_t.SetMinSize((60, 22))
        self.text_ctrl_dt.SetMinSize((60, 22))
        self.text_ctrl_w.SetMinSize((60, 22))
        self.combo_box_bound.SetMinSize((170, 25))
        self.combo_box_bound.SetSelection(0)
        self.text_ctrl_creep.SetMinSize((40, 22))
        self.text_ctrl_fc.SetMinSize((30, 22))
        self.text_ctrl_ft.SetMinSize((30, 22))
        self.text_ctrl_gamma.SetMinSize((30, 22))
        self.combo_box_lx1bar.SetMinSize((100, 25))
        self.combo_box_lx1bar.SetSelection(0)
        self.text_ctrl_lx1Pitch.SetMinSize((50, 22))
        self.combo_box_lx2bar.SetMinSize((100, 25))
        self.combo_box_lx2bar.SetSelection(0)
        self.text_ctrl_lx2Pitch.SetMinSize((50, 22))
        self.combo_box_ly1bar.SetMinSize((100, 25))
        self.combo_box_ly1bar.SetSelection(0)
        self.text_ctrl_ly1Pitch.SetMinSize((50, 22))
        self.combo_box_ly2bar.SetMinSize((100, 25))
        self.combo_box_ly2bar.SetSelection(0)
        self.text_ctrl_ly2Pitch.SetMinSize((50, 22))
        self.panel_bound.SetMinSize((150, 150))
        self.text_ctrl_mx1.SetMinSize((65, 22))
        self.text_ctrl_mx2.SetMinSize((65, 22))
        self.text_ctrl_my1.SetMinSize((65, 22))
        self.text_ctrl_my2.SetMinSize((65, 22))
        self.text_ctrl_atx1.SetMinSize((65, 22))
        self.text_ctrl_atx2.SetMinSize((65, 22))
        self.text_ctrl_aty1.SetMinSize((65, 22))
        self.text_ctrl_aty2.SetMinSize((65, 22))
        self.text_ctrl_rebarx1.SetMinSize((65, 22))
        self.text_ctrl_rebarx2.SetMinSize((65, 22))
        self.text_ctrl_rebary1.SetMinSize((65, 22))
        self.text_ctrl_rebary2.SetMinSize((65, 22))
        self.text_ctrl_lx1PitchOut.SetMinSize((65, 22))
        self.text_ctrl_lx2PitchOut.SetMinSize((65, 22))
        self.text_ctrl_ly1PitchOut.SetMinSize((65, 22))
        self.text_ctrl_ly2PitchOut.SetMinSize((65, 22))
        self.text_ctrl_sfx1.SetMinSize((65, 22))
        self.text_ctrl_sfx2.SetMinSize((65, 22))
        self.text_ctrl_sfy1.SetMinSize((65, 22))
        self.text_ctrl_sfy2.SetMinSize((65, 22))
        self.text_ctrl_reqt.SetMinSize((65, 22))
        self.text_ctrl_tbyl.SetMinSize((65, 22))
        self.text_ctrl_def.SetMinSize((65, 22))
        self.text_ctrl_dBySpan.SetMinSize((65, 22))
        self.list_ctrl_output.AppendColumn("ID", format=wx.LIST_FORMAT_LEFT, width=40)
        self.list_ctrl_output.AppendColumn("Title", format=wx.LIST_FORMAT_LEFT, width=150)
        self.list_ctrl_output.AppendColumn("Sub Titile", format=wx.LIST_FORMAT_LEFT, width=150)
        self.text_ctrl_remove.SetMinSize((40, 22))
        self.text_ctrl_move1.SetMinSize((40, 22))
        self.text_ctrl_move2.SetMinSize((40, 22))
        self.panel_14.SetMinSize((625, 150))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: MyFrame2.__do_layout
        sizer_5 = wx.BoxSizer(wx.VERTICAL)
        sizer_6 = wx.BoxSizer(wx.VERTICAL)
        sizer_32 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_36 = wx.BoxSizer(wx.VERTICAL)
        sizer_41 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_38 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_37 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_24 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_8 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_21 = wx.BoxSizer(wx.VERTICAL)
        sizer_28 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_27 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_29 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_33 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_26 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_25 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_23 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_22 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_34 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_9 = wx.BoxSizer(wx.VERTICAL)
        sizer_20 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_19 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_18 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_17 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_16 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_35 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_14 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_15 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_13 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_12 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_4 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_3 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_11 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_10 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7.Add(self.button_3, 0, wx.ALL, 2)
        sizer_7.Add(self.button_8, 0, wx.ALL, 2)
        sizer_7.Add((20, 20), 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_7.Add(self.text_ctrl_idView, 0, wx.ALL, 1)
        label_46 = wx.StaticText(self.panel_7, wx.ID_ANY, "/")
        sizer_7.Add(label_46, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_7.Add(self.text_ctrl_idTotal, 0, wx.ALL, 1)
        sizer_7.Add(self.button_9, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        self.panel_7.SetSizer(sizer_7)
        sizer_6.Add(self.panel_7, 0, wx.ALL | wx.EXPAND, 3)
        static_line_1 = wx.StaticLine(self.panel_6, wx.ID_ANY)
        sizer_6.Add(static_line_1, 0, wx.EXPAND, 0)
        label_3 = wx.StaticText(self.panel_10, wx.ID_ANY, "Title")
        label_3.SetMinSize((70, 16))
        sizer_10.Add(label_3, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_10.Add(self.text_ctrl_title, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_10, 0, wx.ALL | wx.EXPAND, 1)
        label_4 = wx.StaticText(self.panel_10, wx.ID_ANY, "Sub Title")
        label_4.SetMinSize((70, 16))
        sizer_11.Add(label_4, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_11.Add(self.text_ctrl_subtitle, 0, 0, 0)
        sizer_9.Add(sizer_11, 0, wx.ALL | wx.EXPAND, 1)
        label_5 = wx.StaticText(self.panel_10, wx.ID_ANY, "Lx")
        label_5.SetMinSize((20, 16))
        sizer_3.Add(label_5, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_3.Add(self.text_ctrl_lx, 0, 0, 0)
        label_6 = wx.StaticText(self.panel_10, wx.ID_ANY, "m")
        label_6.SetMinSize((50, 16))
        sizer_3.Add(label_6, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_7 = wx.StaticText(self.panel_10, wx.ID_ANY, "Ly")
        label_7.SetMinSize((20, 16))
        sizer_3.Add(label_7, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_3.Add(self.text_ctrl_ly, 0, 0, 0)
        label_8 = wx.StaticText(self.panel_10, wx.ID_ANY, "m")
        label_8.SetMinSize((30, 16))
        sizer_3.Add(label_8, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_3, 0, wx.ALL | wx.EXPAND, 1)
        label_9 = wx.StaticText(self.panel_10, wx.ID_ANY, "t")
        label_9.SetMinSize((20, 16))
        sizer_4.Add(label_9, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_4.Add(self.text_ctrl_t, 0, 0, 0)
        label_10 = wx.StaticText(self.panel_10, wx.ID_ANY, "mm")
        label_10.SetMinSize((50, 16))
        sizer_4.Add(label_10, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_11 = wx.StaticText(self.panel_10, wx.ID_ANY, "dt")
        label_11.SetMinSize((20, 16))
        sizer_4.Add(label_11, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_4.Add(self.text_ctrl_dt, 0, 0, 0)
        label_12 = wx.StaticText(self.panel_10, wx.ID_ANY, "mm")
        label_12.SetMinSize((30, 16))
        sizer_4.Add(label_12, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_4, 0, wx.ALL | wx.EXPAND, 1)
        label_13 = wx.StaticText(self.panel_10, wx.ID_ANY, "w")
        label_13.SetMinSize((20, 16))
        sizer_12.Add(label_13, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_12.Add(self.text_ctrl_w, 0, 0, 0)
        label_14 = wx.StaticText(self.panel_10, wx.ID_ANY, "kN/m2")
        label_14.SetMinSize((50, 16))
        sizer_12.Add(label_14, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_12, 0, wx.ALL | wx.EXPAND, 1)
        label_15 = wx.StaticText(self.panel_10, wx.ID_ANY, u"支持条件")
        label_15.SetMinSize((80, 16))
        sizer_13.Add(label_15, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_13.Add(self.combo_box_bound, 0, 0, 0)
        sizer_9.Add(sizer_13, 0, wx.ALL | wx.EXPAND, 1)
        label_16 = wx.StaticText(self.panel_10, wx.ID_ANY, u"変形増大率")
        label_16.SetMinSize((80, 16))
        sizer_15.Add(label_16, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_15.Add(self.text_ctrl_creep, 0, 0, 0)
        sizer_9.Add(sizer_15, 0, wx.ALL | wx.EXPAND, 1)
        label_17 = wx.StaticText(self.panel_10, wx.ID_ANY, "Fc")
        label_17.SetMinSize((150, 16))
        sizer_14.Add(label_17, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_14.Add(self.text_ctrl_fc, 0, 0, 0)
        label_18 = wx.StaticText(self.panel_10, wx.ID_ANY, "N/mm2")
        label_18.SetMinSize((50, 16))
        sizer_14.Add(label_18, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_14, 0, wx.ALL | wx.EXPAND, 1)
        label_47 = wx.StaticText(self.panel_10, wx.ID_ANY, "ft")
        label_47.SetMinSize((150, 16))
        sizer_35.Add(label_47, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_35.Add(self.text_ctrl_ft, 0, 0, 0)
        label_48 = wx.StaticText(self.panel_10, wx.ID_ANY, "N/mm2")
        label_48.SetMinSize((50, 16))
        sizer_35.Add(label_48, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_35, 0, wx.ALL | wx.EXPAND, 1)
        label_19 = wx.StaticText(self.panel_10, wx.ID_ANY, u"単位重量")
        label_19.SetMinSize((150, 16))
        sizer_16.Add(label_19, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_16.Add(self.text_ctrl_gamma, 0, 0, 0)
        label_20 = wx.StaticText(self.panel_10, wx.ID_ANY, "kN/m3")
        label_20.SetMinSize((50, 16))
        sizer_16.Add(label_20, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_9.Add(sizer_16, 0, wx.ALL | wx.EXPAND, 1)
        label_21 = wx.StaticText(self.panel_10, wx.ID_ANY, u"配筋　間隔@mm")
        sizer_9.Add(label_21, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        label_22 = wx.StaticText(self.panel_10, wx.ID_ANY, u"短辺端部")
        label_22.SetMinSize((70, 16))
        sizer_17.Add(label_22, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_17.Add(self.combo_box_lx1bar, 0, 0, 0)
        label_23 = wx.StaticText(self.panel_10, wx.ID_ANY, "@")
        sizer_17.Add(label_23, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_17.Add(self.text_ctrl_lx1Pitch, 0, 0, 0)
        sizer_9.Add(sizer_17, 0, wx.ALL | wx.EXPAND, 1)
        label_24 = wx.StaticText(self.panel_10, wx.ID_ANY, u"短辺中央")
        label_24.SetMinSize((70, 16))
        sizer_18.Add(label_24, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_18.Add(self.combo_box_lx2bar, 0, 0, 0)
        label_25 = wx.StaticText(self.panel_10, wx.ID_ANY, "@")
        sizer_18.Add(label_25, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_18.Add(self.text_ctrl_lx2Pitch, 0, 0, 0)
        sizer_9.Add(sizer_18, 0, wx.ALL | wx.EXPAND, 1)
        label_26 = wx.StaticText(self.panel_10, wx.ID_ANY, u"長辺端部")
        label_26.SetMinSize((70, 16))
        sizer_19.Add(label_26, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_19.Add(self.combo_box_ly1bar, 0, 0, 0)
        label_27 = wx.StaticText(self.panel_10, wx.ID_ANY, "@")
        sizer_19.Add(label_27, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_19.Add(self.text_ctrl_ly1Pitch, 0, 0, 0)
        sizer_9.Add(sizer_19, 0, wx.ALL | wx.EXPAND, 1)
        label_28 = wx.StaticText(self.panel_10, wx.ID_ANY, u"長辺中央")
        label_28.SetMinSize((70, 16))
        sizer_20.Add(label_28, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_20.Add(self.combo_box_ly2bar, 0, 0, 0)
        label_29 = wx.StaticText(self.panel_10, wx.ID_ANY, "@")
        sizer_20.Add(label_29, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_20.Add(self.text_ctrl_ly2Pitch, 0, 0, 0)
        sizer_9.Add(sizer_20, 0, wx.ALL | wx.EXPAND, 0)
        self.panel_10.SetSizer(sizer_9)
        sizer_8.Add(self.panel_10, 0, wx.ALL | wx.EXPAND, 2)
        sizer_34.Add((80, 20), 0, 0, 0)
        sizer_34.Add(self.panel_bound, 1, wx.EXPAND, 0)
        sizer_34.Add(self.panel_5, 1, wx.EXPAND, 0)
        self.panel_boundtmp.SetSizer(sizer_34)
        sizer_21.Add(self.panel_boundtmp, 4, wx.EXPAND, 0)
        sizer_22.Add((80, 20), 0, 0, 0)
        label_32 = wx.StaticText(self.panel_11, wx.ID_ANY, u"短辺端部")
        label_32.SetMinSize((65, 16))
        label_32.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_22.Add(label_32, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_33 = wx.StaticText(self.panel_11, wx.ID_ANY, u"短辺中央")
        label_33.SetMinSize((65, 16))
        label_33.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_22.Add(label_33, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_34 = wx.StaticText(self.panel_11, wx.ID_ANY, u"長辺端部")
        label_34.SetMinSize((65, 16))
        label_34.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_22.Add(label_34, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_35 = wx.StaticText(self.panel_11, wx.ID_ANY, u"長辺中央")
        label_35.SetMinSize((65, 16))
        label_35.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_22.Add(label_35, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_21.Add(sizer_22, 0, wx.ALL | wx.EXPAND, 2)
        label_30 = wx.StaticText(self.panel_11, wx.ID_ANY, "M")
        label_30.SetMinSize((30, 16))
        label_30.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_23.Add(label_30, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_31 = wx.StaticText(self.panel_11, wx.ID_ANY, "kN.m")
        label_31.SetMinSize((50, 16))
        label_31.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_23.Add(label_31, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_23.Add(self.text_ctrl_mx1, 0, 0, 0)
        sizer_23.Add(self.text_ctrl_mx2, 0, 0, 0)
        sizer_23.Add(self.text_ctrl_my1, 0, 0, 0)
        sizer_23.Add(self.text_ctrl_my2, 0, 0, 0)
        sizer_21.Add(sizer_23, 0, wx.ALL | wx.EXPAND, 1)
        label_36 = wx.StaticText(self.panel_11, wx.ID_ANY, "at")
        label_36.SetMinSize((30, 16))
        label_36.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_25.Add(label_36, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_37 = wx.StaticText(self.panel_11, wx.ID_ANY, "mm2")
        label_37.SetMinSize((50, 16))
        label_37.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_25.Add(label_37, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_25.Add(self.text_ctrl_atx1, 0, 0, 0)
        sizer_25.Add(self.text_ctrl_atx2, 0, 0, 0)
        sizer_25.Add(self.text_ctrl_aty1, 0, 0, 0)
        sizer_25.Add(self.text_ctrl_aty2, 0, 0, 0)
        sizer_21.Add(sizer_25, 0, wx.ALL | wx.EXPAND, 1)
        label_38 = wx.StaticText(self.panel_11, wx.ID_ANY, u"配筋")
        label_38.SetMinSize((80, 16))
        label_38.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_26.Add(label_38, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_26.Add(self.text_ctrl_rebarx1, 0, 0, 0)
        sizer_26.Add(self.text_ctrl_rebarx2, 0, 0, 0)
        sizer_26.Add(self.text_ctrl_rebary1, 0, 0, 0)
        sizer_26.Add(self.text_ctrl_rebary2, 0, 0, 0)
        sizer_21.Add(sizer_26, 0, wx.ALL | wx.EXPAND, 1)
        sizer_33.Add((80, 20), 0, 0, 0)
        sizer_33.Add(self.text_ctrl_lx1PitchOut, 0, 0, 0)
        sizer_33.Add(self.text_ctrl_lx2PitchOut, 0, 0, 0)
        sizer_33.Add(self.text_ctrl_ly1PitchOut, 0, 0, 0)
        sizer_33.Add(self.text_ctrl_ly2PitchOut, 0, 0, 0)
        sizer_21.Add(sizer_33, 0, wx.ALL | wx.EXPAND, 1)
        label_44 = wx.StaticText(self.panel_11, wx.ID_ANY, u"検定比")
        label_44.SetMinSize((80, 16))
        label_44.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_29.Add(label_44, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_29.Add(self.text_ctrl_sfx1, 0, 0, 0)
        sizer_29.Add(self.text_ctrl_sfx2, 0, 0, 0)
        sizer_29.Add(self.text_ctrl_sfy1, 0, 0, 0)
        sizer_29.Add(self.text_ctrl_sfy2, 0, 0, 0)
        sizer_21.Add(sizer_29, 0, wx.ALL | wx.EXPAND, 1)
        label_39 = wx.StaticText(self.panel_11, wx.ID_ANY, u"必要スラブ厚")
        label_39.SetMinSize((80, 16))
        label_39.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_27.Add(label_39, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_40 = wx.StaticText(self.panel_11, wx.ID_ANY, "mm")
        label_40.SetMinSize((40, 16))
        label_40.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_27.Add(label_40, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_27.Add(self.text_ctrl_reqt, 0, 0, 0)
        label_41 = wx.StaticText(self.panel_11, wx.ID_ANY, "t/Lx")
        label_41.SetMinSize((30, 16))
        label_41.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_27.Add(label_41, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_27.Add(self.text_ctrl_tbyl, 0, 0, 0)
        sizer_21.Add(sizer_27, 0, wx.ALL | wx.EXPAND, 1)
        label_42 = wx.StaticText(self.panel_11, wx.ID_ANY, u"最大変位")
        label_42.SetMinSize((80, 16))
        label_42.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_28.Add(label_42, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        label_43 = wx.StaticText(self.panel_11, wx.ID_ANY, "mm")
        label_43.SetMinSize((40, 16))
        label_43.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_28.Add(label_43, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_28.Add(self.text_ctrl_def, 0, 0, 0)
        sizer_28.Add(self.text_ctrl_dBySpan, 0, 0, 0)
        sizer_21.Add(sizer_28, 0, wx.ALL | wx.EXPAND, 1)
        self.panel_11.SetSizer(sizer_21)
        sizer_8.Add(self.panel_11, 0, wx.ALL | wx.EXPAND, 3)
        self.panel_8.SetSizer(sizer_8)
        sizer_6.Add(self.panel_8, 1, wx.ALL | wx.EXPAND, 3)
        static_line_2 = wx.StaticLine(self.panel_6, wx.ID_ANY)
        sizer_6.Add(static_line_2, 0, wx.EXPAND, 0)
        sizer_24.Add((20, 20), 1, 0, 0)
        sizer_24.Add(self.button_6, 0, wx.ALL, 2)
        sizer_24.Add(self.button_7, 0, wx.ALL, 2)
        self.panel_9.SetSizer(sizer_24)
        sizer_6.Add(self.panel_9, 0, wx.ALL | wx.EXPAND, 3)
        static_line_3 = wx.StaticLine(self.panel_6, wx.ID_ANY)
        sizer_6.Add(static_line_3, 0, wx.EXPAND, 0)
        sizer_32.Add(self.list_ctrl_output, 0, wx.ALL | wx.EXPAND, 3)
        sizer_36.Add(self.panel_17, 1, wx.EXPAND, 0)
        label_53 = wx.StaticText(self.panel_15, wx.ID_ANY, "Edit")
        label_53.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_36.Add(label_53, 0, wx.ALL, 2)
        label_49 = wx.StaticText(self.panel_15, wx.ID_ANY, "No.")
        sizer_37.Add(label_49, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_37.Add(self.text_ctrl_remove, 0, 0, 0)
        sizer_37.Add((20, 20), 1, 0, 0)
        sizer_37.Add(self.button_11, 0, 0, 0)
        sizer_36.Add(sizer_37, 0, wx.ALL | wx.EXPAND, 2)
        label_50 = wx.StaticText(self.panel_15, wx.ID_ANY, "No.")
        sizer_38.Add(label_50, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_38.Add(self.text_ctrl_move1, 0, 0, 0)
        label_51 = wx.StaticText(self.panel_15, wx.ID_ANY, "To")
        sizer_38.Add(label_51, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_38.Add(self.text_ctrl_move2, 0, 0, 0)
        sizer_38.Add((20, 20), 1, 0, 0)
        sizer_38.Add(self.button_15, 0, 0, 0)
        sizer_36.Add(sizer_38, 0, wx.ALL | wx.EXPAND, 2)
        sizer_36.Add(self.panel_16, 1, wx.EXPAND, 0)
        label_54 = wx.StaticText(self.panel_15, wx.ID_ANY, "FILE")
        label_54.SetForegroundColour(wx.Colour(0, 0, 255))
        sizer_36.Add(label_54, 0, wx.ALL, 2)
        sizer_41.Add((20, 20), 1, 0, 0)
        sizer_41.Add(self.button_14, 0, 0, 0)
        sizer_41.Add(self.button_12, 0, 0, 0)
        sizer_41.Add(self.button_13, 0, 0, 0)
        sizer_36.Add(sizer_41, 1, wx.EXPAND, 0)
        self.panel_15.SetSizer(sizer_36)
        sizer_32.Add(self.panel_15, 1, wx.EXPAND, 0)
        self.panel_14.SetSizer(sizer_32)
        sizer_6.Add(self.panel_14, 0, wx.EXPAND, 0)
        self.panel_6.SetSizer(sizer_6)
        sizer_5.Add(self.panel_6, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_5)
        sizer_5.Fit(self)
        self.Layout()
        # end wxGlade

    def OnPre(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnPre' not implemented!")
        event.Skip()

    def OnNext(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnNext' not implemented!")
        event.Skip()

    def OnPlus(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnPlus' not implemented!")
        event.Skip()

    def OnClear(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnClear' not implemented!")
        event.Skip()

    def OnChangeBound(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnChangeBound' not implemented!")
        event.Skip()

    def OnStore(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnStore' not implemented!")
        event.Skip()

    def OnCal(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnCal' not implemented!")
        event.Skip()

    def OnRemove(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnRemove' not implemented!")
        event.Skip()

    def OnShow(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnShow' not implemented!")
        event.Skip()

    def OnMove(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnMove' not implemented!")
        event.Skip()

    def OnReport(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnReport' not implemented!")
        event.Skip()

    def OnImport(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnImport' not implemented!")
        event.Skip()

    def OnExport(self, event):  # wxGlade: MyFrame2.<event_handler>
        print("Event handler 'OnExport' not implemented!")
        event.Skip()

# end of class MyFrame2

class MyFrame_test(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame_test.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((400, 300))
        self.panel_13 = wx.Panel(self, wx.ID_ANY)

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: MyFrame_test.__set_properties
        self.SetTitle("frame_1")
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: MyFrame_test.__do_layout
        sizer_30 = wx.BoxSizer(wx.VERTICAL)
        sizer_31 = wx.BoxSizer(wx.HORIZONTAL)
        label_45 = wx.StaticText(self.panel_13, wx.ID_ANY, "test")
        sizer_31.Add(label_45, 0, 0, 0)
        self.panel_13.SetSizer(sizer_31)
        sizer_30.Add(self.panel_13, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_30)
        self.Layout()
        # end wxGlade

# end of class MyFrame_test

class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyFrame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True

# end of class MyApp

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
