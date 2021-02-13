
# twoWaySlab/２方向スラブの設計

## GUI by wxpython
![Image](./images/twoWaySlab.png)

## Features
- Follow Japanese Code
- Calculate deflection and stress of the two way RC slab
- Solve plate eq. directly by Fourier Method.
- Import/Export by csv format.
- Export pdf file report.

![Image](./images/pdf_image.png)

## Develop Memo

 長方形版の応力は、設計実務においては主要な境界条件に対して、無次元化した数値表が用意されており、その値を用いて応力、変形を推定する事が通例となっています。既存の数値計算ソフトウェアはこの数値表を利用し、内外挿する事により求めているのが多いです。

 本プログラムは、Poisson-Kirchhoffの薄板理論で支配される偏微分方程式を東博士によって考案された境地値で表された長方形平板の基本解を利用し、板内の撓み、及び応力を直接求めているところに、他のソフトウェアとは異なる特徴があります。

 コードはPythonによって記述しました。算定にあたり、以下の手法をとっています。

- ２方向長方形板の応力を有限フーリエ変換を用いた解を利用し、直接求めたものです。
- 解の誘導は、”建築構造学体系11、平板構造、東洋一ほか”を参照してください。
- 文献で求められている数式には一部誤記があるため、適宜修正しています。
- 未定乗数を求める連立方程式はnumpyモジュールを利用して解きました。
- 撓み関数から板内の応力は、開発スピードを優先し、sympyモジュールを利用して導関数を直接求めました。
- この時、級数の打ち切りは5としました。
- 計算プログラムはクラスとして定義し、クラスの名前はHigashiとなっています。

![Image](./images/Equation.png)

## Souce
```
├── README.md
├── aijRc.py
├── report.py
├── twoWaySlab.py
├── gui.py
├── higashi.py
├── gui.wxg
├── db
│   ├── rcslab.txt
├── fonts
│   ├── GenShinGothic-Monospace-Medium.ttf
├── images
│   ├── 4sideFix.jpg
│   ├── 4sideFix.png
│   ├── m2.jpg
│   ├── m2_2pin.jpg
│   ├── m2_2pin2.jpg
│   ├── m2_2pin3.jpg
│   ├── m3-1pin.jpg
│   ├── m3-1pin2.jpg
│   ├── m3_1.jpg
│   ├── m3_1pin.jpg
│   ├── m3_2.jpg
│   ├── m4pin.jpg
```
## twoWaySlab.py
 Main Program

## Higashi.py
 Calculation of the two way plate

### how to
``` python
obj = Higashi()
obj.solve(Id_bound,lx,ly,t,w,creep,ec,nu,nmax,mmax)
```
### input parameter
- Id_bound: Boundary condition (1-5)
![Image](./images/bound.png)
- lx,ly: Shorter & Longer span (m)
- t: Slab thickness (mm)
- w: distributed load (kN/m2)
- creep: Creep Factore (16)
- nu : Poisson's ratio (-)
- nmax, mmax : Num. of the Fourier series (-)
### return
- by list = [Mx1, Mx2, My1, My2, dv]
- Mx1: Negative Momend at Ext. End for the shorter span
- Mx2: Positive Moment at Cent. for the shorter span
- My1: Negative Momend at Ext. End for the longer span
- My2: Positive Moment at Cent. for the longer span

### example console output
``` SHLL
#  Solve, four side fix plate

ly/lx = 1.3333333333333333 nu= 0.0
nmax = 5 mmax 5

A = [[ 0.5090853   0.          0.          0.          0.          0.18675521
  -0.01495861  0.0034888  -0.00129928  0.00061684]
 [ 0.          0.14293034  0.          0.          0.         -0.03721607
   0.02075058 -0.00756752  0.00327597 -0.00166207]
 [ 0.          0.          0.08489664  0.          0.          0.01004821
  -0.01285853  0.00747021 -0.00401272  0.00227042]
 [ 0.          0.          0.          0.06063058  0.         -0.00391218
   0.00716243 -0.00577602  0.00381133 -0.00243721]
 [ 0.          0.          0.          0.          0.04715702  0.0018927
  -0.00413512  0.00412076 -0.0032123   0.00230562]
 [ 0.24900694 -0.04962143  0.01339761 -0.00521623  0.00252361  0.34803985
   0.          0.          0.          0.        ]
 [-0.01994481  0.02766744 -0.0171447   0.0095499  -0.00551349  0.
   0.10611185  0.          0.          0.        ]
 [ 0.00465173 -0.01009003  0.00996028 -0.00770136  0.00549435  0.
   0.          0.06366198  0.          0.        ]
 [-0.00173238  0.00436796 -0.00535029  0.00508177 -0.00428306  0.
   0.          0.          0.04547284  0.        ]
 [ 0.00082246 -0.00221609  0.00302723 -0.00324962  0.00307416  0.
   0.          0.          0.          0.03536777]]

y = [-1.76840332e-01  4.74081318e-03 -6.22833216e-04  1.62159911e-04
 -5.93426112e-05 -1.39099658e-01  2.02765647e-03 -2.62809136e-04
  6.84113770e-05 -2.50351648e-05]

mx = inv(A)@y =  [-0.271811   -0.01676099  0.01236399 -0.00646475  0.00359169 -0.20818584
 -0.02484464  0.01004986 -0.00472525  0.00259829]

mx1 =  -0.06977026632967204
my1 =  -0.056276894239432954

 Search maximum value by x/a = dx =  0.02
 x = 0.00 mmax = 0.02969
 Search maximum value by x/a = dx =  0.026666666666666665
 y = 0.00 mmax = 0.01310
w(0,0) = 0.0314773344818457 -> 0.0236080008613843
mx2= 0.0296938294167846
my2= 0.0131015854262382
(6.279323969670483, 2.67244464751061, 5.064920481548966, 1.17914268836144, 1.76389833295548)
```

## aijRc.py
### how to
``` python
obj = Aij_rc_set()
obj.Ra(index)
obj.Ra_p(index,pitch)
obj.Ec(fc,gamma)
```
### input
- index: "D10", "D13", ......, "D41"
- pitch: "200", "100"
- fc: compressive concrete strength, N/mm2
- gamma: concrete dry density, kN/m3
### return
- obj.Ra: Area of the bar
- obj.Ra_p: Area of the bar within 1.0 m
- Ec: young modulus by AIJ Standard

## Module
- numpy
- sympy
- pandas
- reportlab

## report.py
use "./fonts/GenshinGothic-Monospace-Medium.ttf" for japanese
https://gammasoft.jp/blog/pdf-japanese-font-by-python/
``` python
self.FONT_NAME = "GenShinGothic"
GEN_SHIN_GOTHIC_MEDIUM_TTF = "./fonts/GenShinGothic-Monospace-Medium.ttf"
# フォント登録
pdfmetrics.registerFont(TTFont('GenShinGothic', GEN_SHIN_GOTHIC_MEDIUM_T
```

# Log
## 2020/02/01
Csv入出力の実装
- OnImport
- OnEport
## 2020/02/03
pdfレポートの作成実装

## Run
### For Mac & Linux
Download souse code, then,
``` shell
python3 ./twoWaySlab.py
```

### For Windows

By power shell
``` DOS
> pyinstaller main.py --onefile --noconsole --icon=icons/twoWay_Icon.ico
> ./dist/twoWaySlab/twoWaySlab.exe
> mv ./images ./dist/twoWaySlab/twoWaySlab.exe
> mv ./fonts ./dist/twoWaySlab/twoWaySlab.exe
```

Check Release, and click main/rigidWink.exe!
