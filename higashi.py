#! /Users/tsuno/.pyenv/shims/python3
# -*- coding: utf-8 -*-

# Calculation for the two way slab stress by Plate Structure theory
# Coded by tsunoppy on Sunday

import math

import numpy as np

# File Control
import os

# sympy
import sympy as sym
sym.init_printing()

# matplort
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 5, 5
from sympy import var
from sympy.plotting import plot3d
from sympy import solve, simplify

# 東博士による２重級数による平板解法
# 建築構造学体系11、平板構造、東洋一、小森清司著
# named Higashi as class name, respectiong Dr. Higashi

# 定義した関数
# m_4fix: ４辺固定     --- 3 form
# > w_4fix: たわみ関数
# m_3fix: 3辺固定     --- 3 form
# > w_3fix: たわみ関数
# m_2fix: 2辺固定     --- 3 form
# > w_2fix: たわみ関数
# m_3fix_1pin: ３辺固定、１辺支持     --- 5 form
# > w_3fix_1pin: たわみ関数


class Higashi:

    ########################################################################
    # Init
    def __init__(self):
        self.tmpdata = [] # test data


    def solve(self,Id_bound,lx,ly,t,w,creep,ec,nu,nmax,mmax):

        # ４辺固定
        #   ---------
        #   |       | 2b = ly
        #   ---------
        #    2a = lx
        #
        if Id_bound == 1:
            lamda = ly/lx
            md = self.m_4fix(lamda,nu,nmax,mmax)

        # 3辺固定 - 1
        #   ---------
        #           | 2b = ly
        #   ---------
        #    a = lx
        #  lamda = b/a
        elif Id_bound == 2:
            lamda = ly/lx/2.0
            md = self.m_3fix(lamda,nu,nmax,mmax)

        # 3辺固定 - 2
        #
        #   |       | a = ly
        #   ---------
        #    2b = lx
        #  lamda = b/a
        # need replace lx & ly
        elif Id_bound == 3:
            lamda = lx/ly/2.0
            md = self.m_3fix(lamda,nu,nmax,mmax)

        # 2辺固定
        #
        #           | b = ly
        #   ---------
        #     a = lx
        #  lamda = b/a
        elif Id_bound == 4:
            lamda = ly/lx
            index2 = 1
            md = self.m_2fix(lamda,nu,nmax,mmax,index2)


        # ４辺支持
        #   ---------
        #   |       | 2b = ly
        #   ---------
        #    2a = lx
        #
        elif Id_bound == 5:
            lamda = ly/lx
            md = self.m_4pin(lamda,nu,nmax,mmax)

        # 3辺固定 - 1辺支持
        #   ---------
        #   *       | 2b = ly
        #   ---------
        #    a = lx
        #  lamda = b/a
        elif Id_bound == 6:
            lamda = ly/lx/2.0
            md = self.m_3fix_1pin(lamda,nu,nmax,mmax)

        # 3辺固定 - 1辺支持
        #   ---------
        #   *       | 2b = ly
        #   ---------
        #    a = lx
        #  lamda = b/a
        # need replace lx & ly
        elif Id_bound == 7:
            lamda = lx/ly/2.0
            md = self.m_3fix_1pin(lamda,nu,nmax,mmax)

        # 2隣辺固定 - 2辺支持
        #   *********
        #   *       | b = ly
        #   ---------
        #    a = lx
        #  lamda = b/a
        elif Id_bound == 8:
            lamda = ly/lx
            md = self.m_2nfix_2pin(lamda,nu,nmax,mmax)

        # 2対辺固定 - 2辺支持
        #   *********
        #   |       | 2b = ly
        #   *********
        #    2a = lx
        #  lamda = b/a
        elif Id_bound == 9:
            lamda = ly/lx
            md = self.m_2fix_2pin(lamda,nu,nmax,mmax)

        # 2対辺固定 - 2辺支持
        #   ---------
        #   *       * 2a = ly
        #   ---------
        #    2b = lx
        #  lamda = b/a
        # need replace lx & ly
        elif Id_bound == 10:
            lamda = lx/ly
            md = self.m_2fix_2pin(lamda,nu,nmax,mmax)

        # Show Error Dialog
        else:
            dlg = wx.MessageDialog(self, 'Bound, input Erro/ Higashi.solve',
                                   'Bound Error',
                                   wx.OK | wx.ICON_INFORMATION
                                   )
            dlg.ShowModal()
            dlg.Destroy()

        # 撓みの計算
        if Id_bound == 3 or Id_bound == 7 or Id_bound ==10:
            dv = creep*md[4]*(w/1000.0)*(ly*1000.0)**4/(ec*t**3)
            return abs(md[2]*w*ly**2), md[3]*w*ly**2, \
                abs(md[0]*w*ly**2), md[1]*w*ly**2, dv
        else:
            dv = creep*md[4]*(w/1000.0)*(lx*1000.0)**4/(ec*t**3)
            return abs(md[0]*w*lx**2), md[1]*w*lx**2, \
                abs(md[2]*w*lx**2), md[3]*w*lx**2, dv
    ########################################################################
    # ４辺固定版
    # p81, ３型による解
    def m_4fix(self,lamda,nu,nmax,mmax):

        # 初期計算
        pi = math.pi
        a = 1.0
        b = lamda * a
        y = []

        # nに関する連立方程式
        ####################


        for n in range(1,nmax+1):

            vec = []

            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a
            a1 = math.cosh(beta_a) * math.sinh(beta_a) + beta_a
            a1 = a1 / ( 2.0 * (math.cosh(beta_a))**2 ) / beta_a

            # making vector
            for j in range(1,nmax+1):
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            for m in range(1,mmax+1):

                alpha = (2.0*m-1.0)/(2.0*a) *pi
                b1 = 2.0 * ( -1.0 )**m *( -1.0 )**n * alpha * beta
                b1 = b1 / ( a*b * ( alpha**2 + beta**2 ) **2 )
                vec.append(b1)

            # make y vector
            xtmp = 2.0 * (-1.0) ** n / (beta**4 * a**3 * b )
            xtmp = xtmp*  ( math.cosh(beta_a) * math.sinh(beta_a) - beta_a ) / ( 2.0 * (math.cosh(beta_a))**2 )
            y.append(xtmp)

            # make matrix
            if n == 1:
                aaa = np.array(vec)
            else:
                aaa = np.vstack((aaa,np.array(vec)))

        # mに関する連立方程式
        ####################

        for m in range(1,mmax+1):

            vec = []

            alpha = (2.0*m-1.0)/(2.0*a) *pi
            alpha_b = alpha * b

            b1 = math.cosh(alpha_b) * math.sinh(alpha_b) + alpha_b
            b1 = b1/ ( 2.0 * ( math.cosh( alpha_b ) )**2 )
            b1 = b1/ ( alpha * a )

            for n in range(1,nmax+1):

                beta = (2.0*n-1)/(2.0*b) * pi
                a1 = 2.0 * ( -1.0 )**n *( -1.0 )**m * alpha * beta
                a1 = a1/ ( a**2 * ( alpha**2 + beta**2 )**2 )
                vec.append(a1)

            # make y vector
            xtmp = 2.0 * (-1.0)**m / ( alpha**4 * a**4 )
            xtmp = xtmp * ( math.cosh(alpha_b) * math.sinh(alpha_b) - alpha_b )/ ( 2.0* ( math.cosh(alpha_b) )**2 )
            y.append(xtmp)


            # making vector
            for j in range(1,mmax+1):
                if j == m:
                    vec.append(b1)
                else:
                    vec.append(0)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))

        # 連立方程式を解く
        aaainv = np.linalg.inv(aaa)
        y = np.array(y)
        mx = aaainv @ y

        mx_end = 0
        my_end = 0
        for n in range(1,nmax+1):
            mx_end = mx_end + mx[n-1]
            #print('n=',n,mx[n-1])

        for m in range(1,mmax+1):
            my_end = my_end + mx[nmax+m-1]
            #print('m=',m,mx[nmax+m-1])

        # print calculation log
        print('# ', 'Solve, four side fix plate')
        print()
        print('ly/lx =', lamda, 'nu=', nu)
        print('nmax =', nmax, 'mmax', mmax)
        print()
        print('A =', aaa)
        print()
        print('y =', y)
        print()
        print('mx = inv(A)@y = ', mx)
        print()

        print('mx1 = ', mx_end/4)
        print('my1 = ', my_end/4)
        print()

        # たわみ関数の呼び出
        tmpData = self.w_4fix(lamda,nmax,mmax,mx,nu)

        return mx_end/4,tmpData[0],my_end/4,tmpData[1],tmpData[2]

    # w(x,y) difinition for four side fix model
    ########################################################################
    def w_4fix(self,lamda,nmax,mmax,mx,nu):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a
        #print("lamda",lamda)

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            coef = beta*a * sym.tanh(beta*a) * sym.cosh(beta*x) - beta*x * sym.sinh(beta*x)
            coef = coef /(2.0*sym.cosh(beta*a))
            coef = coef / beta**2 * sym.cos(beta*y)
            coef = coef/ a**2 * mx[n-1]

            ww = ww + coef
            #print(n-1,mx[n-1])

        for m in range(1,mmax+1):
            alpha = ( 2.0*m -1.0 )/ (2.0*a) * pi
            coef = alpha*b * sym.tanh(alpha*b) * sym.cosh(alpha*y) - alpha*y * sym.sinh(alpha*y)
            coef = coef /(2.0*sym.cosh(alpha*b) )
            coef = coef/ alpha**2 * sym.cos(alpha*x)
            coef = coef/ a**2 * mx[nmax+m-1]

            ww = ww + coef
            #print(nmax+m-1,mx[nmax+m-1])

        # 荷重項

        """
        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            coef =  -(beta*a * sym.tanh(beta*a) + 2.0 )* sym.cosh(beta*x) + beta*x*sym.sinh(beta*x)
            coef = coef/( 2.0*sym.cosh(beta*a) )
            coef = 1.0 + coef
            coef = (-1.0)**(n-1)/(beta**5) * coef * sym.cos(beta*y)
            coef = 2.0/(b*a**4) * coef

            ww = ww + coef
        """
        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = ( 2.0*m - 1.0 )/ (2.0*a) * pi
                beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi
                coef = 4.0/(a**5*b)
                coef = coef * (-1.0)**m * (-1.0)**n
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.cos(alpha*x) * sym.cos(beta*y)
                ww = ww + coef

        """
        from sympy.plotting import plot3d
        plot3d(ww,(x,0,a),(y,0,b))
        """

        # 応力の計算
        ## Bending Moment m/(p*a**2)
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )

        mx2_x0 = mx2.subs(y,0)/4.0
        my2_0y = my2.subs(x,0)/4.0

        mx2max = self.fxy_max(mx2_x0,50,0,a,"x")
        my2max = self.fxy_max(my2_0y,50,0,b,"y")

        #self.func_max(my2)
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, x, 2, y )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.0),(y,0.0)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])

        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('w(0,0) =', w00, w00*3/4)
        print('mx2=', mx2max)
        #print('my2=', my2_00/4)
        print('my2=', my2max)
        #print('vx =', vx_a0/2)
        #print('vy =', vy_0b/2)


        #p = plot3d(-ww,(x,-a,a),(y,-b,b))

        return mx2max,my2max,w00*3/4

    ########################################################################
    # 3辺固定版
    # p91, ３型による解
    def m_3fix(self,lamda,nu,nmax,mmax):

        # 初期計算
        pi = math.pi
        a = 1.0
        b = lamda * a
        y = []

        # nに関する連立方程式
        ####################


        for n in range(1,nmax+1):

            vec = []

            # for m_an
            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a
            a1 = math.cosh(beta_a) * math.sinh(beta_a) + beta_a
            a1 = a1 / ( 2.0 * (math.cosh(beta_a))**2 ) / beta_a

            # for theta_n
            c1 = (1.0-nu) * beta_a * math.sinh(beta_a) - 2.0* math.cosh(beta_a)
            c1 = c1/(2.0*(math.cosh(beta_a))**2)

            # making vector
            for j in range(1,nmax+1):
                #print('ma_n')
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            # making vector
            for j in range(1,nmax+1):
                #print('theta_n')
                if j == n:
                    vec.append(c1)
                else:
                    vec.append(0)

            for m in range(1,mmax+1):

                #print('ma_m')

                alpha = (2.0*m-1.0)/(2.0*a) *pi
                b1 = 2.0 * ( -1.0 )**m *( -1.0 )**n * alpha * beta
                b1 = b1 / ( a*b * ( alpha**2 + beta**2 ) **2 )
                vec.append(b1)

            # make y vector
            xtmp = 2.0 * (-1.0) ** n / (beta**4 * a**3 * b )
            xtmp = xtmp*  ( math.cosh(beta_a) * math.sinh(beta_a) - beta_a ) / ( 2.0 * (math.cosh(beta_a))**2 )
            y.append(xtmp)

            # make matrix
            if n == 1:
                aaa = np.array(vec)
            else:
                aaa = np.vstack((aaa,np.array(vec)))
            #print(vec)

        # nに関する連立方程式２
        for n in range(1,nmax+1):

            vec = []

            # for m_an
            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a

            a1 =  ( 1.0 - nu ) * beta_a * math.sinh(beta_a) - 2.0 * math.cosh(beta_a)
            a1 = -a1/( 2.0 * ( math.cosh(beta_a) )**2 )

            # for theta_n
            # 本文では(1.0-nu)*beta*alapha?
            c1 = (3.0+nu)* math.cosh(beta_a) * math.sinh(beta_a) + (1.0-nu)*beta_a
            c1 = c1/( 2.0 * (math.cosh(beta_a) )**2 )
            c1 = c1 * (1.0-nu) * beta_a


            # making vector
            for j in range(1,nmax+1):
                #print('ma_n')
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            # making vector
            for j in range(1,nmax+1):
                #print('theta_n')
                if j == n:
                    vec.append(c1)
                else:
                    vec.append(0)

            for m in range(1,mmax+1):
                #print('ma_m')
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                b1 = 2.0 * ( -1.0 )**n * alpha**2 * beta
                b1 = b1 / ( b* (alpha**2 + beta**2 )**2 )
                b1 = -b1 * ( 1.0 + nu*(beta/alpha)**2 )
                vec.append(b1)

            # make y vector
            xtmp = (1.0-nu)*beta_a*math.sinh(beta_a)/(2.0*(math.cosh(beta_a))**2)
            xtmp = xtmp + nu - nu/math.cosh(beta_a)
            xtmp = ( 2.0 * (-1.0) ** n / (beta**3 * a**2 * b ) ) * xtmp
            y.append(xtmp)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))
            #print(vec)


        # mに関する連立方程式
        ####################

        for m in range(1,mmax+1):

            vec = []

            alpha = (2.0*m-1.0)/(2.0*a) *pi
            alpha_b = alpha * b

            b1 = math.cosh(alpha_b) * math.sinh(alpha_b) + alpha_b
            b1 = b1/ ( 2.0 * ( math.cosh( alpha_b ) )**2 )
            b1 = b1/ ( alpha * a )

            for n in range(1,nmax+1):
                #print('ma_n')
                beta = (2.0*n-1)/(2.0*b) * pi
                a1 = 2.0 * ( -1.0 )**n *( -1.0 )**m * alpha * beta
                a1 = a1/ ( a**2 * ( alpha**2 + beta**2 )**2 )

                vec.append(a1)

            for n in range(1,nmax+1):
                #print('theta_n')
                beta = (2.0*n-1)/(2.0*b) * pi
                #文献では
                #c1 = (-1.0)**n * 2.0 * alpha**2 * beta / alpha / ( alpha**2 + beta**2 )**2
                c1 = (-1.0)**n * 2.0 * alpha**2 * beta / a / ( alpha**2 + beta**2 )**2
                c1 = c1 *( 1.0 + nu * (beta/alpha)**2 )

                vec.append(c1)

            # make y vector
            xtmp = 2.0 * (-1.0)**m / ( alpha**4 * a**4 )
            xtmp = xtmp * ( math.cosh(alpha_b) * math.sinh(alpha_b) - alpha_b )/ ( 2.0* ( math.cosh(alpha_b) )**2 )
            y.append(xtmp)


            # making vector
            for j in range(1,mmax+1):
                #print('ma_m')
                if j == m:
                    vec.append(b1)
                else:
                    vec.append(0)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))
            #print(vec)

        # 連立方程式を解く
        aaainv = np.linalg.inv(aaa)
        y = np.array(y)
        mx = aaainv @ y

        mx_end = 0
        my_end = 0
        for n in range(1,nmax+1):
            mx_end = mx_end + mx[n-1]
            #print('n=',n,mx[n-1])

        for m in range(1,mmax+1):
            my_end = my_end + mx[2*nmax+m-1]
            #print('m=',m,mx[2*nmax+m-1],my_end)

        # print calculation log
        print('# ', 'Solve, three side fix plate')
        print()
        print('ly/lx =', lamda, 'nu=', nu)
        print('2b/a =', lamda*2)
        print('nmax =', nmax, 'mmax', mmax)
        print()
        print('A =', aaa)
        print()
        print('y =', y)
        print()
        print('mx = inv(A)@y = ', mx)
        print()

        print('mx1 = ', mx_end)
        print('my1 = ', my_end)
        print()

        # たわみ関数の呼び出し
        tmpData = self.w_3fix(lamda,nmax,mmax,mx,nu)

        return mx_end,tmpData[0],my_end,tmpData[1],tmpData[2]

    # w(x,y) difinition for three side fix model
    ########################################################################
    def w_3fix(self,lamda,nmax,mmax,mx,nu):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),theta(n=1,2,...,nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a
        #print("lamda",lamda)

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        print('-- Cal w(x,y)')
        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            beta_a = beta*a

            # for man
            coef = beta_a * sym.tanh(beta_a) * sym.cosh(beta*x) - beta*x * sym.sinh(beta*x)

            coef = coef /(2.0*sym.cosh(beta_a))
            coef = coef / beta**2 * sym.cos(beta*y)
            coef = coef/ a**2 * mx[n-1]

            ww = ww + coef

        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            beta_a = beta*a

            # for thetan
            coef = (1.0-nu) * beta_a * sym.tanh(beta_a) - 1.0 -nu
            coef = coef * sym.sinh(beta*(a-x))
            coef = coef - (1.0-nu)*beta*(a-x)*sym.cosh(beta*(a-x))

            coef = coef / ( 2.0*sym.cosh(beta_a) )
            coef = coef / beta * sym.cos(beta*y)
            coef = coef / a * mx[nmax+n-1]

            ww = ww + coef

            #print(n-1,nmax+n-1)

        for m in range(1,mmax+1):
            alpha = ( 2.0*m -1.0 )/ (2.0*a) * pi
            coef = alpha*b * sym.tanh(alpha*b) * sym.cosh(alpha*y) - alpha*y * sym.sinh(alpha*y)

            coef = coef /(2.0*sym.cosh(alpha*b) )
            coef = coef/ alpha**2 * sym.cos(alpha*x)
            coef = coef/ a**2 * mx[2*nmax+m-1]

            ww = ww + coef
            #print(2*nmax+m-1)
            #print(nmax+m-1,mx[nmax+m-1])

        # 荷重項

        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = ( 2.0*m - 1.0 )/ (2.0*a) * pi
                beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi
                coef = 4.0/(a**5*b)
                coef = coef * (-1.0)**m * (-1.0)**n
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.cos(alpha*x) * sym.cos(beta*y)
                ww = ww + coef
            """
            beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi

            coef = 2.0*(-1.0)**(n-1)/ beta**5 / a**4/ b
            coef = coef* ( 2.0*sym.cosh(beta*a) - ( beta*a*sym.tanh(beta*a) + 2.0 )*sym.cosh(beta*x) + beta*x*sym.sinh(beta*x) )
            coef = coef / 2.0/ sym.cosh(beta*a)
            coef = coef * sym.cos(beta*y)
            ww = ww + coef
            """
        # 応力の計算
        ## Bending Moment m/(p*a**2)
        print('-- Cal mx(x,y), my(x,y)')
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, y, x, 2 )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.0),(y,0.0)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])

        mx2_x0 = mx2.subs(y,0)
        my2_0y = my2.subs(x,0)

        mx2max = self.fxy_max(mx2_x0,50,0,a,"x")
        my2max = self.fxy_max(my2_0y,50,0,b,"y")

        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('w(0,0) =', w00, w00*12)
        print('mx2=', mx2_00,mx2max)
        print('my2=', my2_00,my2max)
        #print('vx_a0 =', vx_a0)
        #print('vy_0b =', vy_0b)
        #print( 'mxa0', mx2.subs([(x,a),(y,0.0)]) )
        #print( 'my0b', my2.subs([(x,0),(y,b)]) )

        #return mx2_00,my2_00,w00*12
        return mx2max,my2max,w00*12

    ########################################################################
    # 2辺固定版
    # p99, ３型による解
    # index = 0: 先端固定
    # index = 1: 先端自由
    def m_2fix(self,lamda,nu,nmax,mmax,index):

        # 初期計算
        pi = math.pi
        a = 1.0
        b = lamda * a
        y = []

        # nに関する連立方程式 (3.55)
        ####################

        for n in range(1,nmax+1):

            vec = []

            # for m_an
            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a
            a1 = math.cosh(beta_a) * math.sinh(beta_a) + beta_a
            a1 = a1 / ( 2.0 * (math.cosh(beta_a))**2 ) / beta_a

            # for theta_n
            c1 = (1.0-nu) * beta_a * math.sinh(beta_a) - 2.0* math.cosh(beta_a)
            c1 = c1/(2.0*(math.cosh(beta_a))**2)

            # making vector
            for j in range(1,nmax+1):
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            # making vector
            for j in range(1,nmax+1):
                if j == n:
                    vec.append(c1)
                else:
                    vec.append(0)

            # for m_bm
            for m in range(1,mmax+1):
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                b1 = 2.0 * ( -1.0 )**m *( -1.0 )**n * alpha * beta
                b1 = b1 / ( a*b * ( alpha**2 + beta**2 ) **2 )
                vec.append(b1)

            # for theta_bm
            for m in range(1,mmax+1):
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                d1 = 2.0 * ( -1.0 )**m * alpha * beta**2
                d1 = d1/( b* ( alpha**2 + beta**2 ) **2 )
                d1 = d1* (1.0 + nu*(alpha/beta)**2 )
                vec.append(d1)

            # for gamma
            if index == 0:
                vec.append(2.0/(beta**2*a*b))

            # make y vector
            xtmp = 2.0 * (-1.0) ** n / (beta**4 * a**3 * b )
            xtmp = xtmp*  ( math.cosh(beta_a) * math.sinh(beta_a) - beta_a ) / ( 2.0 * (math.cosh(beta_a))**2 )
            y.append(xtmp)

            # make matrix
            if n == 1:
                aaa = np.array(vec)
            else:
                aaa = np.vstack((aaa,np.array(vec)))
            #print(vec)

        # nに関する連立方程式２ (3.56)
        for n in range(1,nmax+1):

            vec = []

            # for m_an
            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a

            a1 =  ( 1.0 - nu ) * beta_a * math.sinh(beta_a) - 2.0 * math.cosh(beta_a)
            a1 = -a1/( 2.0 * ( math.cosh(beta_a) )**2 )

            # for theta_n
            # 本文では(1.0-nu)*beta*alapha?
            c1 = (3.0+nu)* math.cosh(beta_a) * math.sinh(beta_a) + (1.0-nu)*beta_a
            c1 = c1/( 2.0 * (math.cosh(beta_a) )**2 )
            c1 = c1 * (1.0-nu) * beta_a


            # making vector
            for j in range(1,nmax+1):
                #print('ma_n')
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            # making vector
            for j in range(1,nmax+1):
                #print('theta_n')
                if j == n:
                    vec.append(c1)
                else:
                    vec.append(0)

            for m in range(1,mmax+1):
                #print('ma_m')
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                b1 = 2.0 * ( -1.0 )**n * alpha**2 * beta
                b1 = b1 / ( b* (alpha**2 + beta**2 )**2 )
                b1 = -b1 * ( 1.0 + nu*(beta/alpha)**2 )
                vec.append(b1)

            for m in range(1,mmax+1):
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                d1 = 2.0 * alpha**2 * beta
                d1 = d1 / ( b* (alpha**2 + beta**2 )**2 )
                d1 = -d1* ( 1.0 - nu )**2*beta*a
                vec.append(d1)

            # for gamma
            if index == 0:
                vec.append(0.0)

            # make y vector
            xtmp = (1.0-nu)*beta_a*math.sinh(beta_a)/(2.0*(math.cosh(beta_a))**2)
            xtmp = xtmp + nu - nu/math.cosh(beta_a)
            xtmp = ( 2.0 * (-1.0) ** n / (beta**3 * a**2 * b ) ) * xtmp
            y.append(xtmp)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))
            #print(vec)


        # mに関する連立方程式 (3.57)
        ####################

        for m in range(1,mmax+1):

            vec = []

            alpha = (2.0*m-1.0)/(2.0*a) *pi
            alpha_b = alpha * b

            b1 = math.cosh(alpha_b) * math.sinh(alpha_b) + alpha_b
            b1 = b1/ ( 2.0 * ( math.cosh( alpha_b ) )**2 )
            b1 = b1/ ( alpha * a )

            d1 = (1.0-nu)*alpha_b*math.sinh(alpha_b)-2.0*math.cosh(alpha_b)
            d1 = d1/ ( 2.0 * ( math.cosh( alpha_b ) )**2 )

            for n in range(1,nmax+1):
                beta = (2.0*n-1)/(2.0*b) * pi
                a1 = 2.0 * ( -1.0 )**n *( -1.0 )**m * alpha * beta
                a1 = a1/ ( a**2 * ( alpha**2 + beta**2 )**2 )

                vec.append(a1)

            for n in range(1,nmax+1):
                beta = (2.0*n-1)/(2.0*b) * pi
                #文献では??
                #c1 = (-1.0)**n * 2.0 * alpha**2 * beta / alpha / ( alpha**2 + beta**2 )**2
                c1 = (-1.0)**n * 2.0 * alpha**2 * beta / a / ( alpha**2 + beta**2 )**2
                c1 = c1 *( 1.0 + nu * (beta/alpha)**2 )

                vec.append(c1)

            # make y vector
            xtmp = 2.0 * (-1.0)**m / ( alpha**4 * a**4 )
            xtmp = xtmp * ( math.cosh(alpha_b) * math.sinh(alpha_b) - alpha_b )/ ( 2.0* ( math.cosh(alpha_b) )**2 )
            y.append(xtmp)


            # making vector
            for j in range(1,mmax+1):
                #print('ma_m')
                if j == m:
                    vec.append(b1)
                else:
                    vec.append(0)

            # making vector
            for j in range(1,mmax+1):
                #print('ma_m')
                if j == m:
                    vec.append(d1)
                else:
                    vec.append(0)

            # for gamma
            if index == 0:
                vec.append(2.0/alpha**2/a**2)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))

        # mに関する連立方程式2 (3.58)
        ####################

        for m in range(1,mmax+1):

            vec = []

            alpha = (2.0*m-1.0)/(2.0*a) *pi
            alpha_b = alpha * b

            b1 = (1.0-nu)*alpha_b*math.sinh(alpha_b) - 2.0*math.cosh(alpha_b)
            b1 = -b1/( 2.0* (math.cosh(alpha_b))**2 )

            d1 = (3.0+nu)*math.cosh(alpha_b)*math.sinh(alpha_b) + (1.0-nu)*alpha_b
            d1 = d1*(1.0-nu)/( 2.0* (math.cosh(alpha_b))**2 ) * alpha*a

            for n in range(1,nmax+1):
                beta = (2.0*n-1)/(2.0*b) * pi
                a1 = 2.0*alpha*beta**2/( a*(alpha**2+beta**2)**2 )
                a1 = -a1*(-1.0)**m*(1.0+nu*(alpha/beta)**2 )

                vec.append(a1)

            for n in range(1,nmax+1):
                beta = (2.0*n-1)/(2.0*b) * pi
                c1 = 2.0*alpha*beta**2/( a*(alpha**2+beta**2)**2 )
                c1 = -c1*(1.0-nu)**2*alpha*a

                vec.append(c1)

            # make y vector
            xtmp = (1.0-nu)*alpha_b*math.sinh(alpha_b)/2.0/(math.cosh(alpha_b))**2
            xtmp = xtmp + nu - nu/math.cosh(alpha_b)
            xtmp = 2.0 * (-1.0)**m / ( alpha**3 * a**3 ) * xtmp
            y.append(xtmp)

            # making vector
            for j in range(1,mmax+1):
                #print('ma_m')
                if j == m:
                    vec.append(b1)
                else:
                    vec.append(0)

            # making vector
            for j in range(1,mmax+1):
                #print('ma_m')
                if j == m:
                    vec.append(d1)
                else:
                    vec.append(0)

            # for gamma
            if index == 0:
                vec.append(0.0)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))

        # gammaに関する連立方程式 (3.59)
        ####################
        if index == 0:
            vec=[]
            for n in range(1,nmax+1):
                beta = (2.0*n-1)/(2.0*b) * pi
                beta_a = beta*a
                a1 = beta_a*math.sinh(beta_a)/2.0/(math.cosh(beta_a))**2 / (beta**2*a**2)
                vec.append(a1)
            for n in range(1,nmax+1):
                beta = (2.0*n-1)/(2.0*b) * pi
                beta_a = beta*a
                c1 = (1.0+nu)*math.cosh(beta_a)*math.sinh(beta_a) + (1.0-nu) * beta_a
                c1 = -c1/2.0/(math.cosh(beta_a))**2 / (beta*a)
                vec.append(c1)
            for m in range(1,mmax+1):
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                alpha_b = alpha * b
                b1 = alpha_b*math.sinh(alpha_b)/2.0/(math.cosh(alpha_b))**2/ (alpha**2*a**2)
                vec.append(b1)
            for m in range(1,mmax+1):
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                alpha_b = alpha * b
                d1 = (1.0+nu)*math.cosh(alpha_b)*math.sinh(alpha_b) + (1.0-nu)*alpha_b
                d1 = -d1/2.0/(math.cosh(alpha_b))**2 / (alpha*a)
                vec.append(d1)
            vec.append(b/a)

            xtmp = 0.0
            for m in range(1,mmax+1):
                alpha = (2.0*m-1.0)/(2.0*a) *pi
                alpha_b = alpha * b
                xtmp = xtmp + \
                    2.0*(-1.0)**m/(alpha**5*a**5)*\
                    (1.0- (alpha_b*math.sinh(alpha_b)+2.0*math.cosh(alpha_b))/2.0/(math.cosh(alpha_b))**2 )

            y.append(xtmp)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))

        # 連立方程式を解く
        ####################
        aaainv = np.linalg.inv(aaa)
        y = np.array(y)
        mx = aaainv @ y

        mx_end = 0
        my_end = 0
        for n in range(1,nmax+1):
            mx_end = mx_end + mx[n-1]
            #print('n=',n,mx[n-1])

        for m in range(1,mmax+1):
            my_end = my_end + mx[2*nmax+m-1]
            #print('m=',m,mx[2*nmax+m-1],my_end)

        # print calculation log
        print('# ', 'Solve, two side fix plate')
        print()
        print('ly/lx =', lamda, 'nu=', nu)
        print('nmax =', nmax, 'mmax', mmax)
        print()
        print('A =', aaa)
        print()
        print('y =', y)
        print()
        print('mx = inv(A)@y = ', mx)
        print()

        print('mx1 = ', mx_end)
        print('my1 = ', my_end)
        print()

        # たわみ関数の呼び出し
        tmpData = self.w_2fix(lamda,nmax,mmax,mx,nu,index)

        #return mx_end/4,tmpData[0],my_end/4,tmpData[1],tmpData[2]
        return abs(mx_end), tmpData[0], abs(my_end), tmpData[1], tmpData[2]

    # w(x,y) difinition for two side fix model
    ########################################################################
    def w_2fix(self,lamda,nmax,mmax,mx,nu,index):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),theta(n=1,2,...,nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a
        #print("lamda",lamda)

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        print('-- Cal w(x,y)')
        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            beta_a = beta*a

            # for man
            coef = beta_a * sym.tanh(beta_a) * sym.cosh(beta*x) - beta*x * sym.sinh(beta*x)

            coef = coef /(2.0*sym.cosh(beta_a))
            coef = coef / beta**2 * sym.cos(beta*y)
            coef = coef/ a**2 * mx[n-1]

            ww = ww + coef

        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            beta_a = beta*a

            # for thetan
            coef = (1.0-nu) * beta_a * sym.tanh(beta_a) - 1.0 -nu
            coef = coef * sym.sinh(beta*(a-x))
            coef = coef - (1.0-nu)*beta*(a-x)*sym.cosh(beta*(a-x))

            coef = coef / ( 2.0*sym.cosh(beta_a) )
            coef = coef / beta * sym.cos(beta*y)
            coef = coef / a * mx[nmax+n-1]

            ww = ww + coef

            #print(n-1,nmax+n-1)

        for m in range(1,mmax+1):
            alpha = ( 2.0*m -1.0 )/ (2.0*a) * pi
            coef = alpha*b * sym.tanh(alpha*b) * sym.cosh(alpha*y) - alpha*y * sym.sinh(alpha*y)

            coef = coef /(2.0*sym.cosh(alpha*b) )
            coef = coef/ alpha**2 * sym.cos(alpha*x)
            coef = coef/ a**2 * mx[2*nmax+m-1]

            ww = ww + coef
            #print(2*nmax+m-1)
            #print(nmax+m-1,mx[nmax+m-1])


        # ここがw_3fixと異なる。
        for m in range(1,mmax+1):
            alpha = ( 2.0*m -1.0 )/ (2.0*a) * pi
            coef = (1.0-nu)*alpha*b *sym.tanh(alpha*b) - 1.0 - nu
            coef = coef * sym.sinh(alpha*(b-y))
            coef = coef - ( (1.-nu) * alpha * (b-y) ) * sym.cosh( alpha*(b-y) )

            coef = coef/ ( 2.0*sym.cosh(alpha*b) )
            coef = coef / alpha * sym.cos(alpha*x)
            coef = coef / a * mx[2*nmax + mmax + m - 1]

            ww = ww + coef

        # gamma 項
        # ここがw_3fixと異なる。
        if index == 0:
            print(mx[2*nmax+2*mmax]) # check
            ww = ww + (a-x)*(b-y) / a**2 * mx[2*nmax+2*mmax]

        # 荷重項

        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = ( 2.0*m - 1.0 )/ (2.0*a) * pi
                beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi
                coef = 4.0/(a**5*b)
                coef = coef * (-1.0)**m * (-1.0)**n
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.cos(alpha*x) * sym.cos(beta*y)
                ww = ww + coef

        # 応力の計算
        ## Bending Moment m/(p*a**2)
        print('-- Cal mx(x,y), my(x,y)')
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, y, x, 2 )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.0),(y,0.0)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])

        mx2_x0 = mx2.subs(y,0)/4.0
        my2_0y = my2.subs(x,0)/4.0

        mx2max = self.fxy_max(mx2_x0,50,0,a,"x")
        my2max = self.fxy_max(my2_0y,50,0,b,"y")


        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('w(0,0) =', w00, w00*12)
        print('mx2=', mx2_00,mx2max)
        print('my2=', my2_00,my2max)
        #print('vx_a0 =', vx_a0)
        #print('vy_0b =', vy_0b)
        #print( 'mxa0', mx2.subs([(x,a),(y,0.0)]) )
        #print( 'my0b', my2.subs([(x,0),(y,b)]) )

        #return mx2_00,my2_00,w00*12
        return mx2max,my2max,w00*12

    # w(x,y) difinition for four side pin model
    ########################################################################
    def m_4pin(self,lamda,nu,nmax,mmax):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        # 荷重項

        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = ( 2.0*m - 1.0 )/ (2.0*a) * pi
                beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi
                coef = 4.0/(a**5*b)
                coef = coef * (-1.0)**m * (-1.0)**n
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.cos(alpha*x) * sym.cos(beta*y)
                ww = ww + coef

        # 応力の計算
        ## Bending Moment m/(p*a**2)
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, x, 2, y )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.0),(y,0.0)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])
        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('# ', 'Solve, four side pin plate')
        print()
        print('ly/lx =', lamda, 'nu=', nu)
        print('nmax =', nmax, 'mmax', mmax)
        print('w(0,0) =', w00, w00*3/4)
        print('mx2=', mx2_00/4)
        print('my2=', my2_00/4)
        #print('vx =', vx_a0/2)
        #print('vy =', vy_0b/2)

        return 0.0, mx2_00/4, 0.0, my2_00/4, w00*3/4

    # 3辺固定-1辺支持版
    # p115, 5型による解
    ########################################################################
    def m_3fix_1pin(self,lamda,nu,nmax,mmax):

        # 初期計算
        pi = math.pi
        a = 1.0
        b = lamda * a
        y = []

        # nに関する連立方程式 (3.100)
        ####################

        for n in range(1,nmax+1):

            vec = []

            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a
            a1 = math.cosh(beta_a) * math.sinh(beta_a) - beta_a
            a1 = - a1 / 2.0 / (math.sinh(beta_a))**2  * beta**2 * b**2

            # making vector
            for j in range(1,nmax+1):
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            for m in range(1,mmax+1):

                alpha = m/a *pi
                b1 = 2.0 * ( -1.0 )**m *( -1.0 )**(n-1) * alpha * b
                b1 = b1 / ( ( alpha/beta) **2 + 1.0 ) **2
                vec.append(b1)


            ########################################################################

            # make y vector
            xtmp = 2.0 * (-1.0) ** (n-1) / ( beta*b )
            xtmp = xtmp *(b/a)**2
            xtmp = xtmp * \
                ( math.cosh(beta_a/2.0) * math.sinh(beta_a/2) - beta_a/2 )\
                / (2.0*(math.cosh(beta_a/2.0))**2)
            y.append(xtmp)

            # make matrix
            if n == 1:
                aaa = np.array(vec)
            else:
                aaa = np.vstack((aaa,np.array(vec)))

        # mに関する連立方程式 (3.99)
        ####################

        for m in range(1,mmax+1):

            vec = []

            alpha = m/a *pi
            alpha_b = alpha * b

            b1 = math.cosh(alpha_b) * math.sinh(alpha_b) + alpha_b
            b1 = b1/2.0/(math.cosh(alpha_b))**2
            b1 = -b1 * alpha**2 * a**2

            for n in range(1,nmax+1):

                beta = (2.0*n-1)/(2.0*b) * pi
                a1 = 2.0 * ( -1.0 )**n *( -1.0 )**(m-1) * beta * a
                a1 = a1/ ( 1.0 + (beta/alpha)**2 )**2
                vec.append(a1)

            # make y vector
            xtmp = 2.0 * ( 1.0 - (-1.0)**m )
            xtmp = xtmp/(alpha*a)
            xtmp = xtmp *\
                ( math.cosh(alpha_b)*math.sinh(alpha_b) - alpha_b )\
                / 2.0 / (math.cosh(alpha_b))**2

            y.append(xtmp)


            # making vector
            for j in range(1,mmax+1):
                if j == m:
                    vec.append(b1)
                else:
                    vec.append(0)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))

        # 連立方程式を解く
        aaainv = np.linalg.inv(aaa)
        y = np.array(y)
        mx = aaainv @ y

        mx_end = 0
        my_end = 0
        for n in range(1,nmax+1):
            mx_end = mx_end + mx[n-1]
            #print('n=',n,mx[n-1])

        for m in range(1,mmax+1):
            alpha = m*pi/a
            my_end = my_end + mx[nmax+m-1] * math.sin(alpha*0.5*a)
            #print('m=',m,mx[nmax+m-1])

        # print calculation log
        print('# ', 'Solve, three side fix + one side pin supported plate')
        print()
        print('ly/lx = b/a = ', lamda, 'nu=', nu)
        print('nmax =', nmax, 'mmax', mmax)
        print()
        print('A =', aaa)
        print()
        print('y =', y)
        print()
        print('mx = inv(A)@y = ', mx)
        print()

        print('mx1 = ', mx_end)
        print('my1 = ', my_end)
        print()

        # たわみ関数の呼び出
        tmpData = self.w_3fix_1pin(lamda,nmax,mmax,mx,nu)

        return mx_end,tmpData[0],my_end,tmpData[1],tmpData[2]

    # w(x,y) difinition for three side fix with one side pin model
    # p115
    ########################################################################
    def w_3fix_1pin(self,lamda,nmax,mmax,mx,nu):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a
        #print("lamda",lamda)

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            beta_a = beta*a

            coef = ( beta_a * sym.coth(beta_a) * sym.sinh(beta*x) - beta*x * sym.cosh(beta*x) ) \
                /2/sym.sinh(beta_a) * mx[n-1] / beta**2 / a**2 * sym.cos(beta*y)

            ww = ww + coef

        for m in range(1,mmax+1):
            alpha = m/a * pi
            alpha_b = alpha*b

            coef = ( alpha_b * sym.tanh(alpha_b) * sym.cosh(alpha*y) - alpha*y * sym.sinh(alpha*y) )\
                /2/sym.cosh(alpha_b) * mx[nmax+m-1] / alpha**2 / a**2 * sym.sin(alpha*x)

            ww = ww + coef

        # 荷重項
        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = m/a * pi
                beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi

                coef = 4.0/(a**5*b)
                coef = coef * (1.0-(-1.0)**m) * (-1.0)**(n-1)
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.sin(alpha*x) * sym.cos(beta*y)
                ww = ww + coef

        # 応力の計算
        ## Bending Moment m/(p*a**2)
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )

        mx2_x0 = mx2.subs(y,0)
        my2_05y = my2.subs(x,float(0.5*a))

        print(mx2_x0)
        mx2max = self.fxy_max(mx2_x0,50,0.0000001,a,"x")
        my2max = self.fxy_max(my2_05y,50,0.000001,b,"y")
        #mx2max = 0.0
        #my2max = 0.0

        #self.func_max(my2)
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, x, 2, y )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.5*a),(y,0.0)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])

        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('w(0,0) =', w00, w00*12)
        print('mx2=', mx2max)
        #print('my2=', my2_00/4)
        print('my2=', my2max)
        #print('vx =', vx_a0/2)
        #print('vy =', vy_0b/2)


        #p = plot3d(-ww,(x,-a,a),(y,-b,b))

        return mx2max,my2max,w00*12

    def func_max(self,fxy):

        #
        print("Cal, f_x")
        f_x = sym.diff( fxy, x, 1 )
        print("Cal, f_y")
        f_y = sym.diff( fxy, y, 1 )

        print("Cal, f_x = f_y = 0")
        points = solve([f_x,f_y],[x,y])
        #points = list(map(lambda x:simplify(x),points))
        #print(points)
        for point in points:
            #print(point[0])
            #print(point[1])
            xc = float(point[0])
            yc = float(point[1])
            #print(fxy.subs([(x,point[0]),(y,point[1])]))
            print("xc=",xc,"yc=",yc)
            print("f(xc,yc)=",fxy.subs([(x,xc),(y,yc)]))

    def fxy_max(self,fxy,ndim,x1,x2,index):

        # index = "x" or "y"

        # Preparation
        dx = (x2-x1)/ndim
        print( " dx = ",dx)
        print(fxy)

        # Calculation
        fmax = 0.0

        if index == "y":
            for i in range(0,ndim+1):
                yc = x1 + i * dx
                f = fxy.subs(y,yc)
                print("{:.2f}".format(yc),"{:.5f}".format(f))

                if f >= fmax:
                    fmax = f
                    yp = yc
            print("{:.2f}".format(yp/x2),"{:.5f}".format(fmax))

        elif index == "x":
            for i in range(0,ndim+1):
                xc = x1 + i * dx
                f = fxy.subs(x,xc)
                print("{:.2f}".format(xc),"{:.5f}".format(f))
                if f >= fmax:
                    fmax = f
                    xp = xc
            print("{:.2f}".format(xp/x2),"{:.5f}".format(fmax))
        else:
            print("Error fxy_max, please select index as 'x' or 'y' ")

        return fmax

    ########################################################################
    # 2辺対辺固定2辺支持版
    # p120, ３型による解
    def m_2fix_2pin(self,lamda,nu,nmax,mmax):
        # mmax はダミー
        # 初期計算
        pi = math.pi
        a = 1.0
        b = lamda * a
        y = []
        mx = []

        for n in range(1,nmax+1):


            beta = (2.0*n-1)/(2.0*b) * pi
            beta_a = beta*a

            a1 = 2.0*(-1.0)**n / a**2 / beta**3 / b \
                *( math.cosh(beta_a)*math.sinh(beta_a) - beta_a )\
                /( math.cosh(beta_a)*math.sinh(beta_a) + beta_a )

            mx.append(a1)

        mx_end = 0
        my_end = 0
        for n in range(1,nmax+1):
            mx_end = mx_end + mx[n-1]

        # print calculation log
        print('# ', 'Solve, two side fix w/ two side pin plate')
        print()
        print('ly/lx =', lamda, 'nu=', nu)
        print('nmax =', nmax, 'mmax = 0(No Need)')
        print()
        print('mx = ', mx)
        print('mx1 = ', mx_end/4)
        print('my1 = ', 0.0)
        print()

        # たわみ関数の呼び出
        tmpData = self.w_2fix_2pin(lamda,nmax,mmax,mx,nu)

        return mx_end/4,tmpData[0],my_end/4,tmpData[1],tmpData[2]

    # w(x,y) difinition for two side fix with two side pin model
    # p119
    ########################################################################
    def w_2fix_2pin(self,lamda,nmax,mmax,mx,nu):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a
        #print("lamda",lamda)

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        for n in range(1,nmax+1):

            beta = ( 2.0*n - 1.0 )/ (2.0*b) * pi
            beta_a = beta*a

            coef = ( beta_a * sym.tanh(beta_a) * sym.cosh(beta*x) - beta*x * sym.sinh(beta*x) ) \
                /2/sym.cosh(beta_a) * mx[n-1] / beta**2 / a**2 * sym.cos(beta*y)

            ww = ww + coef

        # 荷重項
        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = ( 2.0*m - 1.0 )/ (2.0*a) * pi
                beta  = ( 2.0*n - 1.0 )/ (2.0*b) * pi

                coef = 4.0/(a**5*b)
                coef = coef * (-1.0)**m * (-1.0)**n
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.cos(alpha*x) * sym.cos(beta*y)
                ww = ww + coef

        # 応力の計算
        ## Bending Moment m/(p*a**2)
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )

        mx2_x0 = mx2.subs(y,0)/4
        my2_0y = my2.subs(x,0)/4

        #print(mx2_x0)
        mx2max = self.fxy_max(mx2_x0,50,0,a,"x")
        my2max = self.fxy_max(my2_0y,50,0,b,"y")
        #mx2max = 0.0
        #my2max = 0.0

        #self.func_max(my2)
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, x, 2, y )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.5*a),(y,0.0)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])
        #mx2max = mx2_00/4

        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('w(0,0) =', w00, w00*4/3)
        print('mx2=', mx2max)
        #print('my2=', my2_00/4)
        print('my2=', my2max)
        #print('vx =', vx_a0/2)
        #print('vy =', vy_0b/2)


        #p = plot3d(-ww,(x,-a,a),(y,-b,b))

        return mx2max,my2max,w00*4/3

    # 2隣辺固定-2辺支持版
    # p118, 1型による解
    ########################################################################
    def m_2nfix_2pin(self,lamda,nu,nmax,mmax):

        # 初期計算
        pi = math.pi
        a = 1.0
        b = lamda * a
        y = []

        # nに関する連立方程式 (3.108)
        ####################

        for n in range(1,nmax+1):

            vec = []

            beta = n/b * pi
            beta_a = beta*a
            a1 = math.cosh(beta_a) * math.sinh(beta_a) - beta_a
            a1 = - a1 / 2.0 / (math.sinh(beta_a))**2  * beta**2 * a**2

            # making vector
            for j in range(1,nmax+1):
                if j == n:
                    vec.append(a1)
                else:
                    vec.append(0)

            for m in range(1,mmax+1):

                alpha = m/a *pi
                b1 = 2.0 * ( -1.0 )**m *( -1.0 )**(n-1) * alpha * a
                b1 = b1 / ( ( alpha/beta) **2 + 1.0 ) **2 * (a/b)
                vec.append(b1)


            ########################################################################

            # make y vector
            xtmp = 2.0 * ( 1.0 - (-1.0) ** n )/ ( beta*b )
            xtmp = xtmp * \
                ( math.cosh(beta_a/2.0) * math.sinh(beta_a/2) - beta_a/2 )\
                / (2.0*(math.cosh(beta_a/2.0))**2)
            y.append(xtmp)

            # make matrix
            if n == 1:
                aaa = np.array(vec)
            else:
                aaa = np.vstack((aaa,np.array(vec)))

        # mに関する連立方程式 (3.109)
        ####################

        for m in range(1,mmax+1):

            vec = []

            alpha = m/a *pi
            alpha_b = alpha * b

            b1 = math.cosh(alpha_b) * math.sinh(alpha_b) - alpha_b
            b1 = b1/2.0/(math.sinh(alpha_b))**2
            b1 = -b1 * alpha**2 * a**2

            for n in range(1,nmax+1):

                beta = n/b * pi
                a1 = 2.0 * ( -1.0 )**n *( -1.0 )**(m-1) * beta * a
                a1 = a1/ ( 1.0 + (beta/alpha)**2 )**2
                vec.append(a1)

            # make y vector
            xtmp = 2.0 * ( 1.0 - (-1.0)**m )
            xtmp = xtmp/(alpha*a)
            xtmp = xtmp *\
                ( math.cosh(alpha_b/2.0)*math.sinh(alpha_b/2.0) - alpha_b/2.0 )\
                / 2.0 / (math.cosh(alpha_b/2.0))**2

            y.append(xtmp)


            # making vector
            for j in range(1,mmax+1):
                if j == m:
                    vec.append(b1)
                else:
                    vec.append(0)

            # make matrix
            aaa = np.vstack((aaa,np.array(vec)))

        # 連立方程式を解く
        aaainv = np.linalg.inv(aaa)
        y = np.array(y)
        mx = aaainv @ y

        mx_end = 0
        my_end = 0
        for n in range(1,nmax+1):
            mx_end = mx_end + mx[n-1]
            #print('n=',n,mx[n-1])

        for m in range(1,mmax+1):
            my_end = my_end + mx[nmax+m-1]
            #print('m=',m,mx[nmax+m-1])

        # print calculation log
        print('# ', 'Solve, three side fix + one side pin supported plate')
        print()
        print('ly/lx = b/a = ', lamda, 'nu=', nu)
        print('nmax =', nmax, 'mmax', mmax)
        print()
        print('A =', aaa)
        print()
        print('y =', y)
        print()
        print('mx = inv(A)@y = ', mx)
        print()

        print('mx1 = ', mx_end)
        print('my1 = ', my_end)
        print()

        # たわみ関数の呼び出
        tmpData = self.w_2nfix_2pin(lamda,nmax,mmax,mx,nu)
        return mx_end,tmpData[0],my_end,tmpData[1],tmpData[2]

    # w(x,y) difinition for three side fix with one side pin model
    # p115
    ########################################################################
    def w_2nfix_2pin(self,lamda,nmax,mmax,mx,nu):

        # lamda:         ly/lx
        # nmax, mmax:    フーリエ級数の打ち切り
        # mx:            [man(n=1,2,...nmax),mam(m=1,2,....mmax)]

        # ２変数の定義
        sym.var('x y', real = True)

        # 初期条件
        pi = math.pi
        a = 1.0
        b = lamda * a
        #print("lamda",lamda)

        # たわみ関数の計算
        # 計算は w * (D/pa**4) を計算
        ww = 0.0

        for n in range(1,nmax+1):

            beta = n/b * pi
            beta_a = beta*a

            coef = ( beta_a * sym.coth(beta_a) * sym.sinh(beta*x) - beta*x * sym.cosh(beta*x) ) \
                /2/sym.sinh(beta_a) * mx[n-1] / beta**2 / a**2 * sym.sin(beta*y)

            ww = ww + coef

        for m in range(1,mmax+1):
            alpha = m/a * pi
            alpha_b = alpha*b

            coef = ( alpha_b * sym.coth(alpha_b) * sym.sinh(alpha*y) - alpha*y * sym.cosh(alpha*y) )\
                /2/sym.sinh(alpha_b) * mx[nmax+m-1] / alpha**2 / a**2 * sym.sin(alpha*x)

            ww = ww + coef

        # 荷重項
        for n in range(1,nmax+1):
            for m in range(1,mmax+1):
                alpha = m/a * pi
                beta  = n/b * pi

                coef = 4.0/(a**5*b)
                coef = coef * (1.0-(-1.0)**m) * (1.0-(-1.0)**n)
                coef = coef/( alpha*beta* (alpha**2 + beta**2)**2 )
                coef = coef * sym.sin(alpha*x) * sym.sin(beta*y)
                ww = ww + coef

        # 応力の計算
        ## Bending Moment m/(p*a**2)
        dwdx2 = sym.diff( ww, x, 2 )
        dwdy2 = sym.diff( ww, y, 2 )
        mx2 = -a**2* ( dwdx2 + nu * dwdy2 )
        my2 = -a**2* ( dwdy2 + nu * dwdx2 )

        mx2_x05 = mx2.subs(y,float(0.5*b))
        my2_05y = my2.subs(x,float(0.5*a))

        #print(mx2_x0)
        mx2max = self.fxy_max(mx2_x05,50,0.0000001,a,"x")
        my2max = self.fxy_max(my2_05y,50,0.000001,b,"y")
        #mx2max = 0.0
        #my2max = 0.0

        #self.func_max(my2)
        """
        ## Reaction force v/(p*a)
        dwdx3   = sym.diff( ww, x, 3 )
        dwdxdy2 = sym.diff( ww, x, y, 2)
        dwdy3   = sym.diff( ww, y, 3 )
        dwdydx2 = sym.diff( ww, x, 2, y )
        vx      = -a**3 * ( dwdx3 + (2.0-nu) * dwdxdy2 )
        vy      = -a**3 * ( dwdy3 + (2.0-nu) * dwdydx2 )
        """
        # 結果の整理
        w00 = ww.subs([(x,0.5*a),(y,0.5*b)])
        mx2_00 = mx2.subs([(x,0.0),(y,0.0)])
        my2_00 = my2.subs([(x,0.0),(y,0.0)])

        #vx_a0  = vx.subs([(x,a),(y,0.0)])
        #vy_0b  = vy.subs([(x,0.0),(y,b)])
        # ログの出力
        print('w(0,0) =', w00, w00*12)
        print('mx2=', mx2max)
        #print('my2=', my2_00/4)
        print('my2=', my2max)
        #print('vx =', vx_a0/2)
        #print('vy =', vy_0b/2)


        #p = plot3d(-ww,(x,-a,a),(y,-b,b))

        return mx2max,my2max,w00*12

########################################################################
# End Class

obj = Higashi()
lamda = 1.0
nu = 0.0
nmax = 5
mmax = 5
#obj.m_4fix(lamda,nu,nmax,mmax)
#obj.m_3fix(lamda,nu,nmax,mmax)
#obj.m_2fix(lamda,nu,nmax,mmax,1)
#obj.m_3fix_1pin(lamda,nu,nmax,mmax)
#obj.m_2fix_2pin(lamda,nu,nmax,mmax)
#obj.m_2nfix_2pin(lamda,nu,nmax,mmax)
"""
# ２変数の定義
sym.var('x y', real = True)
# 初期条件
fxy = 2*x**3 + 4*x*y**2 - 10*x*y + y**2
obj.func_max(fxy)
"""
