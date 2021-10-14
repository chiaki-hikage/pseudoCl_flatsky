# 内容
本コードは2次元平面として近似できる天域領域において、弱い重力レンズ効果のゆがみ(シアー)場のパワースペクトルを測定するコードです。

弱い重力レンズ効果は宇宙大規模構造による重力の効果によって遠方銀河の形がゆがむ現象で、宇宙のダークマター分布を調べることができるユニークな観測量です。
銀河はもともと真円ではなく、弱い重力レンズ効果は楕円率をほんの数％程度変える程度の小さい効果であるため、
広い天域にわたって遠方にある多くの銀河の形を精密に観測し、統計的に重力レンズゆがみの大きさを引き出す処理が必要になります。
しかし、天域の方角によって精度良く観測できる銀河の数が少なかったり、明るい星が手前にあるとその背後にある銀河の形を測れないことがあり、
重力レンズゆがみの場は一様ではありません。

擬似スペクトル法は、非一様な場の影響を補正し、本来の重力レンズゆがみの角度パワースペクトルを測定する方法です。
以下は、レイトレーシングシミュレーションから実際の観測を想定した非一様な重力レンズゆがみ場を作成し、
シミュレーションに入力した重力レンズパワースペクトルが再現できることを確かめた図です。

<img width="623" alt="スクリーンショット 2021-10-14 14 02 26" src="https://user-images.githubusercontent.com/86592645/137254998-644ab80a-8409-45a9-b4b3-9b6556aa7e7a.png">

また、すばるハイパーシュプリームカムの観測データに応用し、宇宙論情報を引き出すのにも使われました。

References:  
- Shear Power Spectrum Reconstruction using Pseudo-Spectrum Method  
Chiaki Hikage, Masahiro Takada, Takashi Hamana, David Spergel  
Mon. Not. Roy. Astron. Soc., Vol.412 Issue 1 (2011), 65-74

- Cosmology from cosmic shear power spectra with Subaru Hyper Suprime-Cam first-year data  
Chiaki Hikage, Masamune Oguri, Takashi Hamana, et al.  
Publ. Astron. Soc. Japan, Vol.71, Issue 2, 43 (2019)  
